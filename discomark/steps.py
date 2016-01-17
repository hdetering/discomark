from __future__ import print_function
from discomark.models import *
from discomark import utils
import datetime
import io
import os
import re
import shutil
import subprocess
import sys
from glob import glob
from Bio import SeqIO, AlignIO
from Bio.Align import AlignInfo
from Bio.Align import MultipleSeqAlignment
#from Bio.Align.Applications import MafftCommandline # can only be used with python >=2.7
from Bio.Blast import NCBIWWW, NCBIXML
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import prifipy

################################
# 1. parse predicted orthologs #
################################
def merge_species(input_dir, ortho_dir, orthologs, log_fh=sys.stderr):
    print("Parsing input files...\n", file=log_fh)
    # combine species
    print("\nMerging orthologs for all species in folder %s" % ortho_dir, file=log_fh)
    for ortho in orthologs:
        if len(ortho.sequences) > 0:
            with open(os.path.join(ortho_dir, "%s.fasta" % ortho.id), 'wt') as f:
                for db_seq in ortho.sequences:
                    seq = Seq(db_seq.residues)
                    rec = SeqRecord(seq)
                    rec.id = db_seq.fasta_id
                    rec.description = db_seq.description
                    SeqIO.write(rec, f, 'fasta')


###########################
# 2. align ortholog files #
###########################
def align_orthologs(ortho_dir, aligned_dir, orthologs, settings, log_fh=sys.stderr):
    print("\nAligning ortholog sequences...", file=log_fh)
    # align each ortholog
    for o in orthologs:
        ortho_fn = os.path.join(ortho_dir, "%s.fasta" % o.id)
        align_fn = os.path.join(aligned_dir, '%s.fasta' % o.id)
        # alignment makes sense only if file contains >1 sequences
        if len(o.sequences) > 1:
            #cline = ['mafft','--localpair','--maxiterate','16','--inputorder','--preservecase', ortho_fn]
            cline = ['mafft'] + [x for x in sum(settings, ()) if len(x.strip())>0] + [ortho_fn]
            print("\t%s " % ' '.join(cline), file=log_fh)
            #stdout = subprocess.check_output(cline) # won't work with python <2.7!
            stdout = subprocess.Popen(cline, stdout=subprocess.PIPE, stderr=None).communicate()[0] # this works with python <2.7
            #mafft_cline = MafftCommandline(input=ortho_fn, localpair=True, maxiterate=16)
            #print("\t%s " % mafft_cline)
            # run MAFFT
            #stdout, stderr = mafft_cline()
            with open(align_fn, "wb") as handle:
                handle.write(stdout)
        # otherwise just copy the ortholog file
        elif len(o.sequences) > 0:
            shutil.copyfile(ortho_fn, align_fn)


######################
# 3. trim alignments #
######################
def trim_alignments(aligned_dir, trimmed_dir, settings, log_fh=sys.stderr):
    print("\nTrimming alignments...", file=log_fh)
    aligned_files = next(os.walk(aligned_dir))[2]
    aligned_files = [os.path.join(aligned_dir, f) for f in os.listdir(aligned_dir) if os.path.isfile(os.path.join(aligned_dir,f))]

    for f in aligned_files:
        o_id = os.path.split(f)[1].split('.')[0]
        out = os.path.join(trimmed_dir, "%s.fasta" % o_id)
        trimal_params = ['trimal', '-in', f, '-out', out, '-htmlout', "%s.html" % out, '-keepheader']
        trimal_params += [x for x in sum(settings, ()) if len(x.strip())>0]
        subprocess.call(trimal_params, stdout=log_fh, stderr=log_fh)


# create BLAST database for reference alignment
def makeblastdb(genome, log_fh=sys.stderr):
    path, filename = os.path.split(genome)
    cline = ['makeblastdb', '-in', genome, '-dbtype', 'nucl']
    print("\t%s\n" % ' '.join(cline), file=log_fh)
    subprocess.call(cline, stdout=log_fh, stderr=log_fh)

###################################
# 4. map against reference genome #
###################################
def map_to_reference(query_dir, mapped_dir, genome, settings, log_fh=sys.stderr):
    # create the BLAST db to map against
    print("\nCreating BLAST database from reference...", file=log_fh)
    makeblastdb(genome, log_fh)

    # combine query sequences into a single FASTA file
    query_files = glob(os.path.join(query_dir, '*.fasta'))
    query_fn = os.path.join(mapped_dir, 'query.fasta')
    query_file = open(query_fn, 'wt')
    for f in query_files:
        o_id = os.path.split(f)[1].split('.')[0]
        SeqIO.write(SeqIO.parse(open(f, 'rt'), 'fasta'), query_file, 'fasta')
    query_file.close()

    # run BLAST
    print("\nRunning BLAST...", file=log_fh)
    out_fn = os.path.join(mapped_dir, 'blast.out')
    # blast_options = ['-query', query_fn, '-db', genome, '-out', out_fn]
    # cline = ['blastn'] + blast_options
    cline = NcbiblastnCommandline(query=query_fn, db=genome, out=out_fn, outfmt='"6 std sstrand"', **dict(settings))
    print("\t%s\n" % cline, file=log_fh)
    stdout, stderr = cline()

    return out_fn

def add_reference(source_dir, target_dir, genome, hits, mafft_settings, log_fh):
    # copy all source alignments to target dir (so alignments without ref mapping don't get lost)
    for f in glob(os.path.join(source_dir, '*.fasta')):
        shutil.copy(f, target_dir)

    # combine ortholog and reference sequences
    for rec in SeqIO.parse(open(genome, 'rt'), 'fasta'):
        if rec.id in hits:
            for rec_hits in hits[rec.id]:
                in_fn  = os.path.join(target_dir, "%s.fasta" % rec_hits['ortholog'])
                out_fn = os.path.join(target_dir, "%s.ref.fa" % rec_hits['ortholog'])
                directions = set() # store set of mapping directions
                with open(out_fn, 'wt') as out_f:
                    # write out ortholog sequences
                    for seq in SeqIO.parse(in_fn, 'fasta'):
                        if seq.id in rec_hits['seqs']:
                            directions.add(rec_hits['seqs'][seq.id])
                        else:
                            print("[WARNING] ortholog sequence '%s' not found in Blast hits." % seq.id, file=log_fh)
                        SeqIO.write(seq, out_f, 'fasta')
                    if len(directions) > 1:
                        print("[WARNING] reference seq '%s' has ortholog seqs mapped in both directions, thus it will not be included in the alignment." % rec.id)
                    else:
                        # retrieve relevant slice of reference
                        start = max(0, rec_hits['range'][0]-100)
                        end = min(len(rec), rec_hits['range'][1]+100)
                        rec_slice = rec[start:end]

                        orientation = directions.pop()
                        # reverse complement reference if necessary
                        if orientation == 'minus':
                            rec_slice.seq = rec_slice.seq.reverse_complement()
                            rec_slice.id = rec_slice.id + '_rv'

                        SeqIO.write(rec_slice.upper(), out_f, 'fasta')

    # align combined files using MAFFT
    print("Realigning Orthologs (including reference)...", file=log_fh)
    for f in glob(os.path.join(target_dir, '*.ref.fa')):
        o_id = os.path.split(f)[1].split('.')[0]
        # run MAFFT (preserve input order, so ref seq is last)
        cline = ['mafft'] + [x for x in sum(mafft_settings, ()) if len(x.strip())>0] + [f]
        print("\t%s " % ' '.join(cline), file=log_fh)
        stdout = subprocess.Popen(cline, stdout=subprocess.PIPE, stderr=log_fh).communicate()[0] # this works with python <2.7

        with open(os.path.join(target_dir, '%s.fasta' % o_id), 'wb') as handle:
            handle.write(stdout)

# alternative handling if no reference was provided
def convertFastaToClustal(in_dir, out_dir):
    for in_path in glob(os.path.join(in_dir, '*.fasta')):
        out_fn = os.path.basename(in_path).replace('.fasta', '.aln')
        out_path = os.path.join(out_dir, out_fn)
        AlignIO.convert(in_path, 'fasta', out_path, 'clustal')


#####################
# 5. design primers #
#####################

def design_primers(source_dir, target_dir, settings, logfile):
    print("\nDesigning primers using PriFi...\n", file=logfile)
    # get rid of previous files
    utils.purge_dir(target_dir)
    aln_files = glob(os.path.join(source_dir, '*.fasta'))
    print("\tChecking for empty alignments...", file=logfile)
    for f in aln_files:
        try:
            align = AlignIO.read(f, 'fasta')
            filename = os.path.basename(f)
            shutil.copyfile(f, os.path.join(target_dir, filename))
        except Exception:
            print("[WARNING] Empty alignment file?! (%s)" % f, file=logfile)
            continue

    # call PriFi for actual primer design
    for f in glob(os.path.join(target_dir, '*.fasta')):
        aln = AlignIO.read(f, 'fasta')
        summary = AlignInfo.SummaryInfo(aln)
        l = aln.get_alignment_length()
        primerpairs = prifipy.findprimers(0, list(aln), summary, l, settings, logfile)
        if not primerpairs:
            print("%s: No valid primer pair found" % f, file=logfile)
        else:
            print('%s: Found %d primer pair suggestions. Writing primer files:' % (f, len(primerpairs)))
            prifipy.writePrimersToFiles(f, primerpairs, 1, logfile)

def design_primers_cl(source_dir, target_dir, prifi, logfile):
    print("\nDesigning primers using PriFi...\n", file=logfile)
    # get rid of previous files
    utils.purge_dir(target_dir)
    aln_files = glob(os.path.join(source_dir, '*.fasta'))
    print("\tChecking for empty alignments...", file=logfile)
    for f in aln_files:
        try:
            align = AlignIO.read(f, 'fasta')
            filename = os.path.basename(f)
            shutil.copyfile(f, os.path.join(target_dir, filename))
        except Exception:
            print("[WARNING] Empty alignment file?! (%s)" % f, file=logfile)
            continue

    # call PriFi for actual primer design
    for f in glob(os.path.join(target_dir, '*.fasta')):
        print(os.getcwd(), file=logfile)
        prifi_params = [prifi, f]
        print(prifi_params, file=logfile)
        sp = subprocess.Popen(prifi_params, stdout=logfile) #, cwd=primer_dir)
        sp.wait()


# export primer-ortholog-reference alignment
def export_primer_alignments(source_dir, orthologs):
    for db_ortho in orthologs:
        primers = db_ortho.primer_sets
        if len(primers) > 0:
            # get ortholog-reference alignment
            aln = AlignIO.read(os.path.join(source_dir, "%s.fasta" % db_ortho.id), 'fasta')
            # format sequences for output
            for rec in aln:
                #rec.id = rec.id if len(rec.id) < 35 else rec.id[:30] + '[...]'
                #rec.description = rec.id
                rec.seq = rec.seq.upper()

            # generate alignment sequence from primers
            pseqs = []
            i = 1
            for ps in primers:
                pos_fw = [int(x) for x in ps.pos_fw.split('-')]
                pos_rv = [int(x) for x in ps.pos_rv.split('-')]
                # generate gapped sequence
                seq = ('-'*pos_fw[0] + ps.seq_fw + '-'*(pos_rv[0]-pos_fw[1]) +
                       Seq(ps.seq_rv).reverse_complement() + '-'*(len(aln[0])-pos_rv[1]))
                rec = SeqRecord(seq, id="%s_%s" % (db_ortho.id, i), description='')
                pseqs.append(rec)
                i += 1
                # determine number of species covered by primer set
                input_seqs = MultipleSeqAlignment([r for r in aln if r.description.find('id_species=')>-1])
                species_ids_fw = set()
                for desc in [r.description for r in input_seqs[:,pos_fw[0]:pos_fw[1]] if r.seq.count('-')==0]:
                    species_ids_fw.update(re.findall('id_species=(\d+)', desc))
                species_ids_rv = set()
                for desc in [r.description for r in input_seqs[:,pos_rv[0]:pos_rv[1]] if r.seq.count('-')==0]:
                    species_ids_rv.update(re.findall('id_species=(\d+)', desc))
                n_species = min(len(species_ids_fw), len(species_ids_rv))
                ps.num_species = n_species
            with open(os.path.join(source_dir, "%s.primer_aln.fasta" % db_ortho.id), 'wt') as f:
                AlignIO.write(MultipleSeqAlignment(pseqs+[r for r in aln]), f, 'fasta')


#############################################################
# 7. primer BLAST                                           #
#    BLAST primers against NCBI nt database for specificity #
#############################################################

# run NCBI BLAST
def blast_primers_online(primer_dir, out_fn, log_fh=sys.stderr):
    primerfile = os.path.join(primer_dir, 'primers.fa')
    print(datetime.datetime.now(), file=log_fh)
    print("Performing remote BLAST search for primers...", file=log_fh)
    handle = NCBIWWW.qblast('blastn', 'refseq_mrna', open(primerfile).read(), entrez_query='txid2[Orgn] OR txid9606[Orgn]')
    print(datetime.datetime.now(), file=log_fh)
    with open(out_fn, 'w') as outfile:
        outfile.write(handle.read())

# run local BLAST
def blast_primers_offline(primer_dir, out_dir):
    print("\nRunning primer BLAST...")
    query_fn = glob(os.path.join(primer_dir, '*.fasta'))[0]
    out_fn = os.path.join(out_dir, 'out.blastn')
    # blast_options = ['-query', query_fn, '-db', genome, '-out', out_fn]
    # cline = ['blastn'] + blast_options
    cline = NcbiblastnCommandline(query=query_fn, db='nt', out=out_fn, outfmt='5') #'"6 std sstrand"')
    print("\t%s\n" % cline)
    cline()

    # parse BLAST hits
    #print("\tLoading BLAST hits...\n")
    #model.load_primer_blast_hits(out_fn)
    #hits = model.get_best_primer_hits()


########################
# final. create report #
########################
def create_report_dir(primer_dir, report_dir):
    if not os.path.exists(report_dir):
        print("\nSetting up report...", file=sys.stderr)
        print("\tCopy report HTML files to %s" % report_dir)
        shutil.copytree(os.path.join('.', 'resources', 'report'), report_dir)

    #print("\nGenerating data for report...\n", file=sys.stderr)
    #utils.generateRecordsJs(primer_dir, os.path.join(report_dir, 'js'))
    #utils.generateAlignmentJs(primer_dir, os.path.join(report_dir, 'js'))
