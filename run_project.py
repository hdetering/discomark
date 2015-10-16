#!/usr/bin/env python
from __future__ import division, print_function

# Run parameters
# (todo: create command line params)
num_threads = 1

import argparse
import os
import sys
# name of configparser module has been changed in Python3
try:
    import configparser # python3
except ImportError:
    import ConfigParser as configparser # python2
from discomark import database, steps, utils

config = configparser.ConfigParser()
config.read("discomark.conf")

prifi   = config.get('Tools', 'prifi')
trimal  = config.get('Tools', 'trimal')

def parse_args():
    parser = argparse.ArgumentParser(description="Discover phylogenetic markers from orthologous sequences.")

    parser.add_argument('-d', '--dir', help="working directory (where results will be stored)", default='./output')
    parser.add_argument('-i', '--input', action='append', help="input folder for sample (at least two '-i' options must be specified)")
    parser.add_argument('-r', '--reference', help="reference genome file (FASTA)")
    parser.add_argument('-s', '--step', help="start from step N", type=int, default=0)
    parser.add_argument('-v', '--verbose', help="increase output verbosity", action='store_true')
    parser.add_argument('--no-primer-blast', help="skip online primer BLAST (use, when running without internet connection", action='store_true')
    args = parser.parse_args()

    # check if output directory exists
    if os.path.isdir(args.dir) and args.step == 0:
        utils.print_error_and_exit("output folder '%s' already exists. Delete it or resume using '-s/--step'." % args.dir)
    elif args.step == 1:
        # check if at least two samples were specified
        if not args.input or not len(args.input) >= 2:
            utils.print_error_and_exit("specify at least two input samples (-i/--input option)")
        else:
            # check if input directories exist and are not empty
            for p in args.input:
                if not os.path.isdir(p):
                    utils.print_error_and_exit("input folder '%s' does not exist." % p)
                elif len([name for name in os.listdir(p) if os.path.isfile(os.path.join(p, name))]) == 0:
                    utils.print_error_and_exit("input folder '%s' does not contain files." % p)

        # check if reference is regular file
        if not os.path.isfile(args.reference):
            utils.print_error_and_exit("reference genome '%s' does not exist or is not a regular file" % args.reference)

    return args

if __name__ == '__main__':
    args = parse_args()
    print(args)
    data_dir = os.path.join(args.dir)

    # 0. initialize folder structure and DB
    reference   = os.path.join(data_dir, 'genome', 'genome.fasta')
    input_dir   = os.path.join(args.dir, config.get('Data', 'input_dir'))
    ortho_dir   = os.path.join(args.dir, config.get('Data', 'ortho_dir'))
    aligned_dir = os.path.join(args.dir, config.get('Data', 'aligned_dir'))
    trimmed_dir = os.path.join(args.dir, config.get('Data', 'trimmed_dir'))
    mapped_dir  = os.path.join(args.dir, config.get('Data', 'mapped_dir'))
    primer_dir  = os.path.join(args.dir, config.get('Data', 'primer_dir'))
    blast_dir   = os.path.join(args.dir, config.get('Data', 'blast_dir'))
    report_dir  = os.path.join(args.dir, config.get('Data', 'report_dir'))

    model = database.DataBroker(args.dir)
    model.initialize_db(args.dir)

    # 1. parse predicted orthologs
    if args.step <= 0:
        utils.setup_output_folders(data_dir, args.input, args.reference, config)
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
        #steps.blast_primers_online(primer_dir, blast_dir)
        pass

    # create report
    steps.create_report(primer_dir, report_dir)
