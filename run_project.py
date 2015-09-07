#!/usr/bin/env python
from __future__ import division, print_function

# Run parameters 
# (todo: create command line params)
num_threads = 1

import argparse, ConfigParser, os, sys
from discomark import database, steps, utils

Config = ConfigParser.ConfigParser()
Config.read("discomark.conf")

prifi   = Config.get('Tools', 'prifi')
trimal  = Config.get('Tools', 'trimal')

def parse_args():
    parser = argparse.ArgumentParser(description="Discover phylogenetic markers from orthologous sequences.")

    parser.add_argument('-d', '--dir', help="working directory (where results will be stored)", default='./output')
    parser.add_argument('-s', '--step', help="start from step N", type=int, default=0)
    parser.add_argument('-v', '--verbose', help="increase output verbosity", action='store_true')
    parser.add_argument('--no-primer-blast', help="skip online primer BLAST (use, when running without internet connection", action='store_true')
    args = parser.parse_args()

    # check if project directory exists
    data_dir = os.path.join("./data", args.project)
    if not os.path.isdir(data_dir):
        print("Error: folder '%s' does not exist!" % data_dir)

    return args

if __name__ == '__main__':
    args = parse_args()
    print(args)
    #sys.exit(0)
    data_dir = os.path.join("./data", args.project)
    model = database.DataBroker(args.project)
    model.initialize_db("./data")

    reference   = os.path.join(data_dir, 'genome', 'genome.fasta')
    input_dir   = os.path.join(data_dir, '1_input')
    ortho_dir   = os.path.join(data_dir, '2_orthologs')
    aligned_dir = os.path.join(data_dir, '3_aligned')
    trimmed_dir = os.path.join(data_dir, '4_trimmed')
    mapped_dir  = os.path.join(data_dir, '5_mapped')
    primer_dir  = os.path.join(data_dir, '6_primers')
    blast_dir   = os.path.join(data_dir, '7_primerblast')
    report_dir  = os.path.join(data_dir, 'report')
    
    # 1. parse predicted orthologs
    if args.step <= 0:
        model.create_db_from_input(input_dir)
    orthologs = model.get_orthologs()
    if args.step <= 1:
        steps.merge_species(input_dir, ortho_dir, orthologs)
    # 2. align ortholog files
    if args.step <= 2:
        steps.align_orthologs(ortho_dir, aligned_dir, orthologs)
    # 3. trim alignments
    if args.step <= 3:
        steps.trim_alignments(aligned_dir, trimmed_dir, trimal)
    # 4. map trimmed alignments against reference genome
    if args.step <= 4:
        out_fn = steps.map_to_reference(trimmed_dir, mapped_dir, reference)
        model.load_blast_hits(out_fn)
        hits = model.get_best_hits()
        steps.add_reference(trimmed_dir, mapped_dir, reference, hits)
    # 5. design primers
    if args.step <= 5:
        steps.design_primers(mapped_dir, primer_dir, prifi)
        model.load_primers(primer_dir)
        orthologs = model.get_orthologs()
        steps.export_primer_alignments(primer_dir, orthologs)
    # 6. primer BLAST
    if args.step <= 6 and not args.no_primer_blast:
        steps.blast_primers_online(primer_dir, blast_dir)
    # create report
    steps.create_report(primer_dir, report_dir)
