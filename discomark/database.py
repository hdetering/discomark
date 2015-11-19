from __future__ import division, print_function
from discomark.models import *
from sqlalchemy import create_engine, desc, distinct, func
from sqlalchemy.orm import sessionmaker
import os, re, sys
from glob import glob
from Bio import SeqIO
from Bio.Blast import NCBIXML

# useful to compare sequences for identity
class DnaSeq:
    id = ''
    dna = ''

    def __init__(self, id, dna):
        self.id = id
        self.dna = dna

    def __eq__(self, other):
        return (isinstance(other, self.__class__) and self.id==other.id and self.dna==other.dna)

    def __hash__(self):
        return hash((self.id, self.dna))

###############################
# DB access abstraction layer #
###############################

class DataBroker():
    """ This class maintains the db session and handles data access. """
    def __init__(self, project_name):
        # choose whether to use an in-memory db or create a db file
        if project_name:
            self.conn_str = 'sqlite:///%s/markers.db' % (project_name)
        else:
            print("\nUsing in-memory database.\n")
            self.conn_str = 'sqlite:///:memory:'

        # connection to database
        self.engine = create_engine(self.conn_str, echo=False)
        # create a database session
        Session = sessionmaker(bind=self.engine)
        self.session = Session()

    def get_session(self):
        return self.session

    # define tables etc...
    def create_schema(self):
        # create database tables
        Base.metadata.create_all(self.engine)

    def get_orthologs(self):
        orthologs = self.session.query(Ortholog).all()
        return orthologs

    # get best (max len) Blast hits from database
    # (result is indexed by reference to facilitate reference fasta processing)
    def get_best_hits(self):
        # retrieve maximal length of Blast hit for each ortholog
        ortho_hits = self.session.query(Ortholog.id, Mapping.refseq, func.max(Mapping.length)).join(Ortholog.sequences).join(Sequence.mappings).group_by(Ortholog).all()
        ref2ortho = {}
        # for each ortholog's best hit, get
        #  - orientation of each sequence hit
        #  - reference sequence range
        for o_id, ref_id, max_len in ortho_hits:
            hits = (
                self.session.query(
                    Ortholog.id,
                    Sequence.fasta_id,
                    Mapping.ref_start,
                    Mapping.ref_end,
                    Mapping.strand
                )
                .join(Ortholog.sequences)
                .join(Sequence.mappings)
                .filter(Ortholog.id==o_id, Mapping.refseq==ref_id)
                .all()
            )
            start = min([x[2] for x in hits])
            end   = max([x[3] for x in hits])
            ref2ortho[ref_id] = {'ortholog':o_id, 'range': (start, end), 'seqs': {x[1]: x[4] for x in hits}}

        return ref2ortho


    # visualizing the table schema
    @staticmethod
    def plot_schema(conn_str):
        from sqlalchemy import MetaData
        from sqlalchemy_schemadisplay import create_schema_graph

        # create the pydot graph object by autoloading all tables via a bound metadata object
        graph = create_schema_graph(metadata=MetaData(conn_str),
                                    show_datatypes=False, # The image would get nasty big if we'd show the datatypes
                                    show_indexes=False, # ditto for indexes
                                    rankdir='LR', # From left to right (instead of top to bottom)
                                    concentrate=False # Don't try to join the relation lines together
                                    )
        graph.write_png('dbschema.png') # write out the file

    #plot_schema(conn_str)

    # helper function to implement unique objects
    def get_or_create(self, model, **kwargs):
        instance = self.session.query(model).filter_by(**kwargs).first()
        if instance:
            return instance
        else:
            instance = model(**kwargs)
            self.session.add(instance)
            return instance

    # load general ortholog info (available orthologs + functional annotations)
    def initialize_db(self, data_dir):
        # create database tables
        self.create_schema()
        session = self.session

        print("\nLoading functional categories ...")
        cat_fn = os.path.join(data_dir, 'fun.txt')
        if not os.path.exists(cat_fn):
            print("\tWarning: file '%s' is missing! Orthologs will not be functionally annotated." % cat_fn)
            return False

        for l in open(cat_fn):
            line = l.rstrip()
            if len(line) == 0:
                continue
            if not line.startswith(' '):
                db_cat = Category(name=line.strip())
                session.add(db_cat)
            else:
                m = re.search(r' \[(\w)\] (.+)', line)
                code = m.group(1)
                desc = m.group(2)
                db_fun = Function(shortcode=code, description=desc)
                db_fun.category = db_cat
                session.add(db_fun)
        session.commit()

        print("\nLoading ortholog annotations ...")
        anno_fn = os.path.join(data_dir, 'ortho2fun.csv')
        if not os.path.exists(anno_fn):
            print("\tWarning: file '%s' is missing! Orthologs will not be functionally annotated." % anno_fn)
            return False

        for line in open(anno_fn):
            o_id, cog, fun, prot = line.strip('\n').split('\t')
            db_ortho = Ortholog(id=int(o_id), prot=prot)
            session.add(db_ortho)
            for f in list(fun):
                db_fun = session.query(Function).filter(Function.shortcode==f).one()
                db_ortho.functions.append(db_fun)
                session.add(db_ortho)
        session.commit()

    # loading input data
    # =====================
    def create_db_from_input(self, input_dir, log_fh=sys.stderr):
        session = self.session

        print("\nLoading data from directory '%s' ..." % input_dir, file=log_fh)
        species = sorted(next(os.walk(input_dir))[1])
        print("\nFound %d species:\n\t%s\n" % (len(species), '\n\t'.join(species)), file=log_fh)

        # traverse species folders
        for sp_name in species:
            db_species = Species(name=sp_name)
            session.add(db_species)

            sp_dir = os.path.join(input_dir, sp_name)
            sp_files = glob(os.path.join(sp_dir, '*.fa'))

            # loop through FASTA files
            for fn in sp_files:
                oid = re.findall("^\d+", os.path.split(fn)[1])[0]
                #oid = os.path.split(fn)[1].split('.')[0]

                db_ortho = self.get_or_create(Ortholog, id=str(oid))
                db_file = File(path=fn)
                db_file.ortholog = db_ortho
                session.add(db_file)

                # make sure sequences are unique
                sequences = set([DnaSeq(r.id, str(r.seq)) for r in SeqIO.parse(open(fn, 'rt'), 'fasta')])
                for seq in sequences:
                    db_seq = Sequence(fasta_id=seq.id, description='', sequence=str(seq.dna.upper()))
                    db_seq.species = db_species
                    db_seq.ortholog = db_ortho
                    session.add(db_seq)

        # save data to database
        session.commit()


    # loading BLAST hits
    # ======================
    def load_blast_hits(self, blast_filename, add=False):
        session = self.session

        # truncate existing table if not in 'add' mode
        if not add:
            tab = Base.metadata.tables['mappings']
            session.execute(tab.delete())
            session.commit()

        # find best reference hits in local alignments
        f = open(blast_filename, 'rt')
        row = f.readline().strip().split()
        seq_id, ref_id, max_len, start, end, strand = (row[0], row[1], int(row[3]), int(row[8]), int(row[9]), row[12])
        # format: qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sstrand
        for line in f:
            row = line.strip().split()
            if seq_id == row[0]:
                # extend reference position range if necessary
                if ref_id == row[1]:
                    reverse = (row[12] == 'minus')
                    hit_start = int(row[8]) if not reverse else int(row[9])
                    hit_end   = int(row[9]) if not reverse else int(row[8])
                    start = min(start, hit_start)
                    end   = max(end, hit_end)
                    max_len = max(max_len, int(row[3]))
                    strand = row[12]
                continue
            else:
                mapping = Mapping(refseq=ref_id, ref_start=start, ref_end=end, length=max_len, strand=strand)
                seq = session.query(Sequence).filter_by(fasta_id=seq_id).first()
                mapping.sequence = seq
                session.add(mapping)
            seq_id, ref_id, max_len, start, end, strand = (row[0], row[1], int(row[3]), int(row[8]), int(row[9]), row[12])

        session.commit()

    def load_blast_hits_xml(blast_filename):
        session = self.session

        # find best (i.e. longest) reference hits in local alignments
        for rec in NCBIXML.parse(open(blast_filename, 'rt')):
            seq_id = rec.descriptions[0].title
            max_len = 0
            ref_id = ''
            strand = ''
            for al in rec.alignments:
                if al.length > max_len:
                    max_len = al.length
                    ref_id = al.hit_def.split()[0]
                    strand = '+' if al.hsps[0].sbjct_start < al.hsps[0].sbjct_end else '-'
            mapping = Mapping(refseq=ref_id, length=max_len, strand=strand)
            seq = session.query(Sequence).filter_by(fasta_id=seq_id).first()
            mapping.sequence = seq
            session.add(mapping)

        session.commit()


    def load_primers(self, primer_dir, add=False):
        session = self.session

        # truncate existing table if not in 'add' mode
        if not add:
            tab = Base.metadata.tables['primer_sets']
            session.execute(tab.delete())
            session.commit()

        regex = '''Primer set (?P<ps_idx>\d+)\s+\((?P<pos_fw>\S+) / (?P<pos_rv>\S+)\)

Fw 5'-(?P<seq_fw>[^\n]+)
Rv 5'-(?P<seq_rv>[^\n]+)

Tm = (?P<Tm_fw>\S+) / (?P<Tm_rv>\S+)
Primer lengths: (?P<len_fw>\d+) / (?P<len_rv>\d+)
Avg\. #sequences in primer alignments: \S+ / \S+
(Estimated p|P)roduct length: (?P<prod_len>\d+)
'''

        primer_files = glob(os.path.join(primer_dir, '*.rep'))
        for pf in primer_files:
            o_id = os.path.split(pf)[1].split('.')[0]
            text = open(pf, 'rtU').read()
            primer_texts = re.findall("(Primer set.+?)---", text, re.DOTALL)
            for t in primer_texts:
                m = re.match(regex, t, re.DOTALL)
                ps = PrimerSet(id_ortholog=o_id)
                ps.ps_idx = int(m.group('ps_idx'))
                ps.prod_len = int(m.group('prod_len'))
                ps.pos_fw = m.group('pos_fw')
                ps.pos_rv = m.group('pos_rv')
                ps.seq_fw = m.group('seq_fw')
                ps.seq_rv = m.group('seq_rv')
                ps.tm_fw = float(m.group('Tm_fw'))
                ps.tm_rv = float(m.group('Tm_rv'))
                session.add(ps)
        session.commit()
        session.close()

    # export primers to file
    def export_primers_to_file(self, filename):
        primer_sets = self.session.query(PrimerSet).all()

        with open(filename, 'wt') as outfile:
            for ps in primer_sets:
                outfile.write(">%d_%d_fw\n%s\n" % (ps.id, ps.id_ortholog, ps.seq_fw))
                outfile.write(">%d_%d_rv\n%s\n" % (ps.id, ps.id_ortholog, ps.seq_rv))

    # load primer BLAST hits from XML file
    def load_primer_blast_hits_xml(self, blast_xml):
        session = self.session
        # find best reference hits in local alignments
        for rec in NCBIXML.parse(open(blast_xml, 'rt')):
            if len(rec.alignments) > 0:
                # query id format: "<id_primerset>_<id_ortholog>_(fw|rv)"
                primer_id = int(rec.query.split('_')[0])
                primer_fwrv = rec.query.split('_')[-1]
                subject_id = rec.alignments[0].accession

                session.query(PrimerSet).filter_by(id=primer_id).update({"blast_%s" % primer_fwrv: subject_id})

        session.commit()

    def primersets_to_records_js(self, target_fn):
        #primer_sets = self.session.query(PrimerSet).all()
        primer_sets = self.session.query(PrimerSet, func.count(distinct(Species.id)).label('n_species')) \
                            .join(Ortholog).join(Sequence).join(Species) \
                            .group_by(PrimerSet.id) \
                            .order_by(desc("n_species"), Ortholog.id) \
                            .all()

        rec_cnt = 0
        js = "var myRecords = [\n"
        for res in primer_sets:
            ps = res[0]
            js += ps.to_json(rec_cnt, res[1]) + ',\n'
            rec_cnt += 1
        js += "\n];"

        with open(target_fn, 'wt') as outfile:
            print(outfile.name)
            outfile.write(js)

    def generateSummaryJs(self, target_fn):
        n_primers = self.session.query(func.count(PrimerSet.id)).one()[0]
        n_markers = (self.session.query(
            func.count(distinct(Ortholog.id)))
            .join(PrimerSet)
            .one()
        )[0]
        categories = (self.session.query(
            Category.name,
            func.count(distinct(Ortholog.id)))
            .join(Function)
            .join(Ortholog, Function.orthologs)
            .join(PrimerSet)
            .group_by(Category.name)
            .all()
        )
        functions = (self.session.query(
            Category.name,
            Function.shortcode,
            func.count(distinct(Ortholog.id)))
                     .join(Function)
            .join(Ortholog, Function.orthologs)
            .join(PrimerSet)
            .group_by(Function.shortcode)
            .order_by(Category.name)
            .all()
        )

        out_str = '''var summary = [{
    'n_primers': %i,
    'n_markers': %i
}];

    var categories = [
    %s
];

var subcats = [
    %s
];''' % (n_primers, n_markers,
             ',\n\t'.join(["['%s', %i]" % x for x in categories]),
             ',\n\t'.join(["['%s', %i]" % (x[1], x[2]) for x in functions]),)

        with open(target_fn, 'wt') as outfile:
            print(outfile.name)
            outfile.write(out_str)
