#!/usr/bin/env python
from __future__ import division, print_function

program = "DiscoMark"
version = "0.9"
# Run parameters
# (todo: create command line params)
num_threads = 1

import argparse
import datetime
import os
import shutil
import sys
try: # name of configparser module has been changed in Python3
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
    parser.add_argument('--no-trim', help="skip alignment trimming step", action='store_true')
    parser.add_argument('--no-primer-blast', help="skip online primer BLAST (use, when running without internet connection", action='store_true')
    args = parser.parse_args()

    # if no arguments were provided, display help and exit
    if len(sys.argv) < 2:
        parser.print_help()
        sys.exit(0)

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
    print("\n%s v%s\n" % (program, version))
    args = parse_args()
    do_ref_map = args.reference and os.path.exists(args.reference)

    # 0. initialize folder structure and DB
    reference   = os.path.join(args.dir, config.get('Data', 'reference_dir'), 'genome.fasta')
    input_dir   = os.path.join(args.dir, config.get('Data', 'input_dir'))
    ortho_dir   = os.path.join(args.dir, config.get('Data', 'ortho_dir'))
    aligned_dir = os.path.join(args.dir, config.get('Data', 'aligned_dir'))
    trimmed_dir = os.path.join(args.dir, config.get('Data', 'trimmed_dir'))
    mapped_dir  = os.path.join(args.dir, config.get('Data', 'mapped_dir'))
    primer_dir  = os.path.join(args.dir, config.get('Data', 'primer_dir'))
    blast_dir   = os.path.join(args.dir, config.get('Data', 'blast_dir'))
    report_dir  = os.path.join(args.dir, config.get('Data', 'report_dir'))

    if args.step <= 0:
        if not os.path.exists(args.dir):
            utils.setup_output_folders(args.dir, args.input, args.reference, config)
        logfile = open(os.path.join(args.dir, 'discomark.log'), 'wt')
        print(datetime.datetime.now(), file=logfile)
        print("Running DiscoMark with the following parameters:\n%s" % args, file=logfile)
        model = database.DataBroker(args.dir)
        model.initialize_db(args.dir)
    else:
        logfile = open(os.path.join(args.dir, 'discomark.log'), 'at')
        print("\n\n\n---%s" % datetime.datetime.now(), file=logfile)
        print("Resuming DiscoMark with the following parameters:\n%s" % args, file=logfile)
        model = database.DataBroker(args.dir)
        # was a reference supplied in this call?
        if do_ref_map and not os.path.exists(reference):
            shutil.copyfile(args.reference, reference)
        # was a reference supplied in an earlier call?
        elif not do_ref_map and os.path.exists(reference):
            do_ref_map = True


    # 1. parse predicted orthologs
    if args.step <= 0:
        model.create_db_from_input(input_dir)
    orthologs = model.get_orthologs()
    if args.step <= 1:
        print("\n[1] Merging orthologs from input folders...")
        steps.merge_species(input_dir, ortho_dir, orthologs, logfile)
    # 2. align ortholog files
    if args.step <= 2:
        print("\n[2] Aligning orthologous sequences...")
        settings = config.items('02_MAFFT_settings')
        steps.align_orthologs(ortho_dir, aligned_dir, orthologs, settings, logfile)
    # 3. trim alignments
    if args.step <= 3 and not args.no_trim:
        print("\n[3] Trimming alignments...")
        settings = config.items('03_TrimAl_settings')
        steps.trim_alignments(aligned_dir, trimmed_dir, trimal, settings, logfile)
    # 4. map trimmed alignments against reference genome
    if args.step <= 4:
        print("\n[4] Mapping alignments to reference...")
        if do_ref_map:
            source_dir = aligned_dir
            settings = config.items('04_BLAST_settings')
            out_fn = steps.map_to_reference(source_dir, mapped_dir, reference, settings, logfile)
            model.load_blast_hits(out_fn)
            hits = model.get_best_hits()
            settings = config.items('04_MAFFT_settings')
            steps.add_reference(source_dir, mapped_dir, reference, hits, settings, logfile)
        else:
            print("\t-> no reference genome provided -> skipping this step...")

    # 5. design primers
    if args.step <= 5:
        print("\n[5] Designing primers based on multiple alignments...")
        source_dir = mapped_dir if do_ref_map else (trimmed_dir if not args.no_trim else aligned_dir)
        steps.design_primers(source_dir, primer_dir, prifi, logfile)
        model.load_primers(primer_dir)
        model.export_primers_to_file(os.path.join(primer_dir, 'primers.fa'))
        orthologs = model.get_orthologs()
        steps.export_primer_alignments(primer_dir, orthologs)
    # 6. primer BLAST
    if args.step <= 6 and not args.no_primer_blast:
        print("\n[6] Searching primer sequences in BLAST database...")
        blast_outfile = os.path.join(blast_dir, 'blast_out.xml')
        steps.blast_primers_online(primer_dir, blast_outfile, logfile)
        model.load_primer_blast_hits_xml(blast_outfile)

    # create report
    steps.create_report_dir(primer_dir, report_dir)
    print("\nGenerating data for report...\n", file=sys.stderr)
    model.primersets_to_records_js(os.path.join(report_dir, 'js', 'records.js'))
    utils.generateAlignmentJs(primer_dir, os.path.join(report_dir, 'js'))
    model.generateSummaryJs(os.path.join(report_dir, 'js', 'summary.js'))

    logfile.close()
