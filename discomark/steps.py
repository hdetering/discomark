from models import *
import os, shutil, subprocess
from StringIO import StringIO
from glob import glob
from Bio import SeqIO, AlignIO
#from Bio.Align.Applications import MafftCommandline # can only be used with python >=2.7
from Bio.Blast import NCBIWWW
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


################################
# 1. parse predicted orthologs #
################################
def merge_species(input_dir, ortho_dir, orthologs):
    print("Parsing input files...\n")
    # combine species
    print("\nMerging orthologs for all species in folder %s" % ortho_dir)
    for ortho in orthologs:
        f = open(os.path.join(ortho_dir, "%s.fasta" % ortho.id), 'wt')
        for db_seq in ortho.sequences:
            seq = Seq(db_seq.sequence)
            rec = SeqRecord(seq)
            rec.id = db_seq.fasta_id
            rec.description = db_seq.description
            SeqIO.write(rec, f, 'fasta')
        f.close()


###########################
# 2. align ortholog files #
###########################
def align_orthologs(ortho_dir, aligned_dir, orthologs):    
    print("\nAligning ortholog sequences...")    
    # align each ortholog
    for o in orthologs:
        ortho_fn = os.path.join(ortho_dir, "%s.fasta" % o.id)
        align_fn = os.path.join(aligned_dir, '%s.aligned.fasta' % o.id)
        # alignment makes sense only if file contains >1 sequences
        if len(o.sequences) > 1:
            cline = ['mafft','--localpair','--maxiterate','16','--inputorder','--preservecase',ortho_fn]
            print("\t%s " % ' '.join(cline))
            #stdout = subprocess.check_output(cline) # won't work with python <2.7!
            stdout = subprocess.Popen(cline, stdout=subprocess.PIPE).communicate()[0] # this works with python <2.7
            #mafft_cline = MafftCommandline(input=ortho_fn, localpair=True, maxiterate=16)
            #print("\t%s " % mafft_cline)
            # run MAFFT
            #stdout, stderr = mafft_cline()
            with open(align_fn, "wt") as handle:
                handle.write(stdout)
        # otherwise just copy the ortholog file
        else:
            shutil.copyfile(ortho_fn, align_fn)


######################
# 3. trim alignments #
######################
def trim_alignments(aligned_dir, trimmed_dir, trimal):    
    print("\nTrimming alignments...")
    aligned_files = os.walk(aligned_dir).next()[2]
    aligned_files = [os.path.join(aligned_dir, f) for f in os.listdir(aligned_dir) if os.path.isfile(os.path.join(aligned_dir,f))]
    
    for f in aligned_files:
        o_id = os.path.split(f)[1].split('.')[0]
        out = os.path.join(trimmed_dir, "%s.trim.fasta" % o_id)
        trimal_params = [trimal, "-in", f, "-out", out, "-htmlout", "%s.html" % out, "-strictplus"]
        subprocess.call(trimal_params)


# create BLAST database for reference alignment
def makeblastdb(genome):
    path, filename = os.path.split(genome)
    cline = ['makeblastdb', '-in', genome, '-dbtype', 'nucl']
    print("\t%s\n" % cline)
    subprocess.call(cline)

###################################
# 4. map against reference genome #
###################################
def map_to_reference(query_dir, mapped_dir, genome):    
    import os, glob
    from Bio import SeqIO
    from Bio.Blast.Applications import NcbiblastnCommandline
    from Bio.Blast import NCBIXML
    
    # create the BLAST db to map against
    print("\nCreating BLAST database from reference...")
    makeblastdb(genome)
    
    query_files = glob.glob(os.path.join(query_dir, '*.fasta'))
    
    # combine query sequences into a single FASTA file
    query_fn = os.path.join(mapped_dir, 'query.fasta')
    query_file = open(query_fn, 'wt')
    for f in query_files:
        o_id = os.path.split(f)[1].split('.')[0]
        SeqIO.write(SeqIO.parse(open(f, 'rt'), 'fasta'), query_file, 'fasta')
    query_file.close() 

    # run BLAST
    print("\nRunning BLAST...")
    out_fn = os.path.join(mapped_dir, 'blast.out')
    # blast_options = ['-query', query_fn, '-db', genome, '-out', out_fn]
    # cline = ['blastn'] + blast_options
    cline = NcbiblastnCommandline(query=query_fn, db=genome, out=out_fn, outfmt='"6 std sstrand"')
    print("\t%s\n" % cline)
    cline()
    
    return out_fn
        
def add_reference(trimmed_dir, mapped_dir, genome, hits):
    # combine ortholog and reference sequences
    for rec in SeqIO.parse(open(genome, 'rt'), 'fasta'):
        if rec.id in hits:
            rec_hits = hits[rec.id]
            in_fn  = os.path.join(trimmed_dir, "%s.trim.fasta" % rec_hits['ortholog'])
            out_fn = os.path.join(mapped_dir, "%s.ref.fasta" % rec_hits['ortholog'])
            with open(out_fn, 'wt') as out_f:
                # write out ortholog sequences
                for seq in SeqIO.parse(in_fn, 'fasta'):
                    if seq.id in rec_hits['seqs']:
                        # reverse complement sequence if necessary
                        if rec_hits['seqs'][seq.id][1] == 'minus':
                            seq = seq.reverse_complement()
                            seq.id = seq.id + '_rv'
                            seq.description = seq.id
                    else:
                        print("[WARNING] ortholog sequence '%s' not found in Blast hits." % seq.id)
                    SeqIO.write(seq, out_f, 'fasta')
                #out_f.write(in_f.read())
                # write out relevant slice of reference
                start = max(0, rec_hits['range'][0]-100)
                end = min(len(rec), rec_hits['range'][1]+100)
                SeqIO.write(rec[start:end].upper(), out_f, 'fasta')
    
    # align combined files using MAFFT
    print("Realigning Orthologs (including reference)...")
    for f in glob(os.path.join(mapped_dir, '*.ref.fasta')):
        o_id = os.path.split(f)[1].split('.')[0]
        # run MAFFT (preserve input order, so ref seq is last)
        cline = ['mafft','--localpair','--maxiterate','16','--inputorder','--preservecase',f]
        print("\t%s " % ' '.join(cline))
        #stdout = subprocess.check_output(cline) # won't work with python <2.7!
        stdout = subprocess.Popen(cline, stdout=subprocess.PIPE).communicate()[0] # this works with python <2.7
        #mafft_cline = MafftCommandline(input=f, localpair=True, maxiterate=16, inputorder=True, clustalout=True)
        # print("\t%s " % mafft_cline)
        #stdout, stderr = mafft_cline()
        #align = AlignIO.read(StringIO(stdout), "clustal")
        with open(os.path.join(mapped_dir, '%s.mapped.aln' % o_id), 'wb') as handle:
            aln = AlignIO.read(StringIO(stdout), 'fasta')
            AlignIO.write(aln, handle, 'clustal')


#####################
# 5. design primers #
#####################
def design_primers(mapped_dir, primer_dir, prifi):
    import os, glob
    from Bio import AlignIO
    
    print("\nDesigning primers using PriFi...\n")
    
    mapped_files = glob.glob(os.path.join(mapped_dir, '*.mapped.aln'))
    
    print("\tChecking for empty alignments...")
    for f in mapped_files:
        try:
            align = AlignIO.read(f, 'clustal')
        except Exception:
            print("Whoa! Empty alignment file?! (%s)" % f)
            continue
    
        # exchange seq names with numbers
        # (Clustal format truncates to len 30 but need to be unique for PriFi)
        i = 0
        for rec in align:
            rec.id = str(i)
            rec.name = str(i)
            i += 1

        o_id = os.path.split(f)[1].split('.')[0]
        handle = open(os.path.join(primer_dir, o_id+'.prifi.aln'), 'wt')
        AlignIO.write(align, handle, 'clustal')
        handle.close()

    # call PriFi for actual primer design
    import subprocess
    for f in glob.glob(os.path.join(primer_dir, '*.prifi.aln')):
        print(os.getcwd())
        prifi_params = [prifi, f]
        print(prifi_params)
        sp = subprocess.Popen(prifi_params) #, cwd=primer_dir)
        sp.wait()


# export primer-ortholog-reference alignment
def export_primer_alignments(primer_dir, orthologs):
    import os
    from Bio import AlignIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord
    from Bio.Align import MultipleSeqAlignment

    for db_ortho in orthologs:
        primers = db_ortho.primer_sets
        if len(primers) > 0:
            # get ortholog-reference alignment
            aln = AlignIO.read(os.path.join(primer_dir, "%s.prifi.aln" % db_ortho.id), 'clustal')
            # generate alignment sequence from primers
            pseqs = []
            i = 1
            for ps in primers:
                pos_fw = [int(x) for x in ps.pos_fw.split('-')]
                pos_rv = [int(x) for x in ps.pos_rv.split('-')]
                # generate gapped sequence
                seq = ('-'*pos_fw[0] + ps.seq_fw + '-'*(pos_rv[0]-pos_fw[1]) + 
                       Seq(ps.seq_rv).reverse_complement() + '-'*(len(aln[0])-pos_rv[1]))
                rec = SeqRecord(seq, id="%s_%s-%s" % (db_ortho.id, i, i))
                pseqs.append(rec)
                i += 1
            with open(os.path.join(primer_dir, "%s.primer_aln.fasta" % db_ortho.id), 'wt') as f:
                AlignIO.write(MultipleSeqAlignment(pseqs+[r for r in aln]), f, 'fasta')

#############################################################
# 7. primer BLAST                                           #
#    BLAST primers against NCBI nt database for specificity #
#############################################################

# run NCBI BLAST
def blast_primers_online(primer_dir, out_dir):
    print(datetime.datetime.now())
    print("Performing remote primer BLAST...")
    handle = NCBIWWW.qblast('blastn', 'refseq_mrna', open('primers.fa').read(), entrez_query='txid2[Orgn] OR txid9606')
    print(datetime.datetime.now())
    out_fn = os.path.join(out_dir, 'blast_out.xml')
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
def create_report(primer_dir, report_dir):
    import utils

    print("\nGenerating data for report...\n")
    utils.generateRecordsJs()
    utils.generateAlignmentJs()
