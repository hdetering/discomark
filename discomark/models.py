################
# data classes #
################

from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import relationship
from sqlalchemy.schema import ForeignKey
Base = declarative_base()

from sqlalchemy import Table, Column, Float, Integer, String, Text

class Category(Base):
    """ Functional categories used in ortholog annotations. """
    __tablename__ = 'categories'

    id        = Column(Integer, primary_key=True)
    name      = Column(String)
    functions = relationship('Function', backref='category')

    def __repr__(self):
        return "<Category('%s')>" % self.name

# many-to-many relationship between Function and Ortholog
fun_orto = Table('function_ortholog', Base.metadata,
                 Column('id_function', Integer, ForeignKey('functions.id')),
                 Column('id_ortholog', Integer, ForeignKey('orthologs.id')))

class Function(Base):
    """ Functional used in ortholog annotations. """
    __tablename__ = 'functions'

    id          = Column(Integer, primary_key=True)
    id_category = Column(Integer, ForeignKey('categories.id'))
    shortcode   = Column(String)
    description = Column(String)
    orthologs   = relationship('Ortholog', secondary=fun_orto, backref='functions')

    def __repr__(self):
        return "<Function('%s','%s')>" % (self.shortcode, self.description)

class Species(Base):
    """ Species' refer to biological species'. """
    __tablename__ = 'species'

    id   = Column(Integer, primary_key=True)
    name = Column(String)

    def __repr__(self):
        return "<Species(name='%s')>" % self.name

class Ortholog(Base):
    """ Orthologs are groups of sequences encoding for the same gene. """
    __tablename__ = 'orthologs'

    id          = Column(Integer, primary_key=True)
    prot        = Column(String)
    sequences   = relationship("Sequence", backref="ortholog")
    files       = relationship("File", backref="ortholog")

    def __repr__(self):
        return "<Ortholog(%d)>" % self.id

class OrthologAnnotation(Base):
    """ OrthologAnnotations specify functional aspects for Orthologs. """
    __tablename__ = 'ortholog_annotations'

    id          = Column(Integer, primary_key=True)
    id_ortholog = Column(Integer, ForeignKey('orthologs.id'))
    funcode     = Column(String)

class Sequence(Base):
    """ Sequences represent generic FASTA sequences. """
    __tablename__ = 'sequences'

    id          = Column(Integer, primary_key=True)
    id_species  = Column(Integer, ForeignKey('species.id'))
    id_ortholog = Column(Integer, ForeignKey('orthologs.id'))
    fasta_id    = Column(String)
    description = Column(String)
    residues    = Column(Text)
    length      = Column(Integer)
    species     = relationship("Species", backref="sequences")

    def __repr__(self):
        return "<Sequence(id='%s')>" % self.fasta_id

class File(Base):
    """ Files refer to input files on the file system. """
    __tablename__ = 'files'

    id          = Column(Integer, primary_key=True)
    id_ortholog = Column(Integer, ForeignKey('orthologs.id'))
    path        = Column(String)

    def __repr__(self):
        t = 'ortholog' if self.id_ortholog != None else '?'
        return "<File(type=%s, path=%s)>" % (t,self.path)

class Mapping(Base):
    """ Mappings are sequence-to-reference alignment hits (from BLAST). """
    __tablename__ = 'mappings'

    id          = Column(Integer, primary_key=True)
    id_sequence = Column(Integer, ForeignKey('sequences.id'))
    refseq      = Column(String)
    ref_start   = Column(Integer)
    ref_end     = Column(Integer)
    length      = Column(Integer)
    strand      = Column(String)
    sequence    = relationship("Sequence", backref="mappings")

    def __repr__(self):
        return "<Mapping(seq='%s', ref='%s', len=%d, strand='%s')>" % (self.sequence.fasta_id, self.refseq, self.length, self.strand)

#class Alignment(Base):
#    """ Alignments are collections of AlignedSequences. """
#    __tablename__ = 'alignments'
#
#    id          = Column(Integer, primary_key=True)
#    id_ortholog = Column(Integer, ForeignKey('orthologs.id'))
#    length      = Column(Integer)
#    ortholog    = relationship("Ortholog", backref="alignments")

#class AlignedSequence(Base):
#    """ AlignedSequences represent the individual sequences in Alignments. """
#    __tablename__ = 'aligned_sequences'
#
#    id          = Column(Integer, primary_key=True)
#    id_sequence = Column(Integer, ForeignKey('sequences.id'))
#    residues    = Column(Text)
#    sequence    = relationship("Sequence")

class PrimerSet(Base):
    """ Primers are templates for PCR amplification of Orthologs. """
    __tablename__ = 'primer_sets'

    id          = Column(Integer, primary_key=True)
    id_ortholog = Column(Integer, ForeignKey('orthologs.id'))
    ortholog    = relationship("Ortholog", backref="primer_sets")
    ps_idx      = Column(Integer)
    prod_len    = Column(Integer)
    pos_fw      = Column(String)
    pos_rv      = Column(String)
    seq_fw      = Column(String)
    seq_rv      = Column(String)
    tm_fw       = Column(Float)
    tm_rv       = Column(Float)
    blast_fw    = Column(String) # NCBI accession
    blast_rv    = Column(String) # NCBI accession
    num_species = Column(Integer)

    def __repr__(self):
        return "<PrimerSet(ortholog='%d', product=%dbp)>" % (self.ortholog.id, self.prod_len)

    def to_json(self, idx, n_spec):
        format_str = '''  {
    "index": "%s",
    "export": "0",
    "markerId": "%s",
    "ps_idx": "%s",
    "no.OfSpecies": "%s",
    "est.ProductLength-(bp)": "%s",
    "fwSequence-(5'-3')": "%s",
    "rvSequence-(5'-3')": "%s",
    "tm": "%s/%s",
    "primerLength": "%s/%s",
    "fwBlastHit": "%s",
    "rvBlastHit": "%s",
    "class": "%s"
  }'''
        return format_str % (idx,
                             self.ortholog.id,
                             self.ps_idx,
                             self.num_species,
                             self.prod_len,
                             self.seq_fw,
                             self.seq_rv,
                             self.tm_fw, self.tm_rv,
                             len(self.seq_fw), len(self.seq_rv),
                             self.blast_fw,
                             self.blast_rv,
                             ','.join([f.shortcode for f in self.ortholog.functions]))

    def to_csv(self, n_spec, sep=','):
        field_values = ["%s_%s" % (self.ortholog.id, self.ps_idx),
                        str(self.ortholog.id),
                        ''.join([f.shortcode for f in self.ortholog.functions]),
                        str(self.num_species),
                        str(self.prod_len),
                        self.seq_fw,
                        self.seq_rv,
                        "%s/%s" % (self.tm_fw, self.tm_rv),
                        "%s/%s" % (len(self.seq_fw), len(self.seq_rv)),
                        self.blast_fw if self.blast_fw else '-',
                        self.blast_rv if self.blast_rv else '-']
        return sep.join(field_values) + '\n'
