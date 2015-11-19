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

    def toJson(self, idx, n_spec):
        format_str = '''  {
    "index": "%s",
    "export": "0",
    "markerId": "%s",
    "species": "%s"
    "product-(bp)": "%s",
    "fwSequence": "%s",
    "rvSequence": "%s",
    "tm": "%s/%s",
    "primerLength": "%s/%s"
  }'''
        return format_str % (idx, self.id, n_spec, self.prod_len, self.fw.seq, self.rv.seq, self.fw.tm, self.rv.tm, len(self.fw.seq), len(self.rv.seq))

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

    # copy resource files
    resource_dir = os.path.join('.', 'resources')
    shutil.copy(os.path.join(resource_dir, 'fun.txt'), out_dir)
    shutil.copy(os.path.join(resource_dir, 'ortho2fun.csv'), out_dir)

    # copy input data
    if genome:
        shutil.copy(genome, os.path.join(dirs['ref_dir'], 'genome.fasta'))
    for d in input_seqs:
        dirname = os.path.basename(os.path.normpath(d))
        shutil.copytree(d, os.path.join(dirs['input_dir'], dirname))


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
