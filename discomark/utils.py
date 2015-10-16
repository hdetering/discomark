from __future__ import division, print_function
from glob import glob
import os
import re
import shutil
import sys
from Bio import AlignIO
from Bio.Seq import Seq


class Primer:
    seq = None
    length = 0
    pos = ""
    tm = 0

    def __init__(self, sequence, start_end, temp):
        self.seq = Seq(sequence)
        self.length = len(self.seq)
        self.pos = start_end
        self.tm = temp

class PrimerPair:
    id = ''
    fw = None
    rv = None
    prod_len = 0

    def __init__(self, name, fw_primer, rv_primer, product_len):
        self.id = name
        self.fw = fw_primer
        self.rv = rv_primer
        self.prod_len = product_len

    def toJson(self, idx):
        format_str = '''  {
    "index": "%s",
    "export": "0",
    "orthologId": "%s",
    "product-(bp)": "%s",
    "fwSequence": "%s",
    "rvSequence": "%s",
    "tm": "%s/%s",
    "primerLength": "%s/%s"
  }'''
        return format_str % (idx, self.id, self.prod_len, self.fw.seq, self.rv.seq, self.fw.tm, self.rv.tm, len(self.fw.seq), len(self.rv.seq))

def print_error_and_exit(msg):
    """Prints an error message to STDERR and exits with exit code 1.
    """
    print("Error: %s" % msg, file=sys.stderr)
    sys.exit(1)

def setup_output_folders(out_dir, input_seqs, genome, config):
    """Sets up the output directory structure and copies input data.
    """
    dirs = {
      'ref_dir'     : os.path.join(out_dir, config.get('Data', 'reference_dir')),
      'input_dir'   : os.path.join(out_dir, config.get('Data', 'input_dir')),
      'ortho_dir'   : os.path.join(out_dir, config.get('Data', 'ortho_dir')),
      'aligned_dir' : os.path.join(out_dir, config.get('Data', 'aligned_dir')),
      'trimmed_dir' : os.path.join(out_dir, config.get('Data', 'trimmed_dir')),
      'mapped_dir'  : os.path.join(out_dir, config.get('Data', 'mapped_dir')),
      'primer_dir'  : os.path.join(out_dir, config.get('Data', 'primer_dir')),
      'blast_dir'   : os.path.join(out_dir, config.get('Data', 'blast_dir'))
    }

    # create directories
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    for d in dirs:
        if not os.path.exists(dirs[d]):
            os.makedirs(dirs[d])

    # copy input data
    shutil.copy(genome, os.path.join(dirs['ref_dir'], 'genome.fasta'))
    for d in input_seqs:
        dirname = os.path.basename(os.path.normpath(d))
        shutil.copytree(d, os.path.join(dirs['input_dir'], dirname))


def parsePrifiReport(filename):
    regex = '''Primer set \d+\s+\((?P<pos_fw>\S+) / (?P<pos_rv>\S+)\)

Fw 5'-(?P<seq_fw>[^\n]+)
Rv 5'-(?P<seq_rv>[^\n]+)

Tm = (?P<Tm_fw>\S+) / (?P<Tm_rv>\S+)
Primer lengths: (?P<len_fw>\d+) / (?P<len_rv>\d+)
Avg\. #sequences in primer alignments: \S+ / \S+
(Estimated p|P)roduct length: (?P<prod_len>\d+)
'''

    primer_pairs = []
    o_id = os.path.split(filename)[1].split('.')[0]
    text = open(filename, 'rtU').read()
    primer_texts = re.findall("(Primer set.+?)---", text, re.DOTALL)
    for t in primer_texts:
        m = re.match(regex, t, re.DOTALL)
        prod_len = int(m.group('prod_len'))
        pos_fw = m.group('pos_fw')
        pos_rv = m.group('pos_rv')
        seq_fw = m.group('seq_fw')
        seq_rv = m.group('seq_rv')
        tm_fw = m.group('Tm_fw')
        tm_rv = m.group('Tm_rv')
        p_fw = Primer(seq_fw, pos_fw, tm_fw)
        p_rv = Primer(seq_rv, pos_rv, tm_rv)
        pp = PrimerPair(o_id, p_fw, p_rv, prod_len)

        primer_pairs.append(pp)

    return primer_pairs


def parsePrimers(primer_dir):
    primer_pairs = []
    primer_files = glob(os.path.join(primer_dir, '*.rep'))
    for pf in primer_files:
        primer_pairs += parsePrifiReport(pf)

    return primer_pairs

def generateRecordsJs(primer_dir, target_dir):
    pps = parsePrimers(primer_dir)
    rec_cnt = 0

    js = "var myRecords = [\n"
    for pp in pps:
        rec_cnt += 1
        js += pp.toJson(rec_cnt-1) + ',\n'
    js += "\n];"
    outfile = open(os.path.join(target_dir, 'records.js'), 'w')
    print(outfile.name)
    outfile.write(js)
    outfile.close()

def generateAlignmentJs(primer_dir, target_dir):
    filenames = glob(os.path.join(primer_dir, '*.primer_aln.fasta'))
    out_str = "var alignments = {\n"
    for fn in filenames:
        o_id = os.path.split(fn)[1].split('.')[0]
        aln = AlignIO.read(fn, 'fasta')
        out_str += "'%s': {\n%s\n},\n" % (o_id, ',\n'.join(["'%s': '%s'" % (r.id, r.seq) for r in aln]))

    out_str = out_str[:-2]
    out_str += "\n};"

    f = open(os.path.join(target_dir, 'alignments.js'), 'w')
    print(f.name)
    f.write(out_str)
    f.close()

def generateSummaryJs(project_name, target_dir):
    from database import DataBroker, Category, Function, Ortholog, PrimerSet
    from sqlalchemy import func, distinct

    m = db.DataBroker(project_name)

    n_primers = m.session.query(func.count(PrimerSet.id)).one()[0]
    n_orthologs = (m.session.query(
        func.count(distinct(Ortholog.id)))
        .join(PrimerSet)
        .one()
    )[0]
    categories = (m.session.query(
        Category.name,
        func.count(distinct(Ortholog.id)))
        .join(Function)
        .join(Ortholog, Function.orthologs)
        .join(PrimerSet)
        .group_by(Category.name)
        .all()
    )
    functions = (m.session.query(
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
    'n_orthologs': %i
}];

var categories = [
    %s
];

var subcats = [
    %s
];''' % (n_primers, n_orthologs,
         ',\n\t'.join(["['%s', %i]" % x for x in categories]),
         ',\n\t'.join(["['%s', %i]" % (x[1], x[2]) for x in functions]),)

    f = open(os.path.join(target_dir, 'summary.js'), 'w')
    f.write(out_str)
    f.close()
