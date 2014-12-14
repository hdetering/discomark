#!/usr/bin/env python
from __future__ import division, print_function

# Run parameters 
# (todo: create command line params)
num_threads = 1

import os, ConfigParser
from discomark import database, steps, utils

Config = ConfigParser.ConfigParser()
Config.read("discomark.conf")

prifi   = Config.get('Tools', 'prifi')
trimal  = Config.get('Tools', 'trimal')
stellar = Config.get('Tools', 'stellar')    

def parse_args():
    import argparse, os

    parser = argparse.ArgumentParser(description="Discover phylogenetic markers from orthologous sequences.")
    parser.add_argument('project', help="Refers to the subfolder in the 'data' directory")
    parser.add_argument('-v', '--verbose', help="increase output verbosity", action='store_true')
    args = parser.parse_args()

    # check if project directory exists
    data_dir = os.path.join("./data", args.project)
    if not os.path.isdir(data_dir):
        print("Error: folder '%s' does not exist!" % data_dir)

    return args

if __name__ == '__main__':
    args = parse_args()
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
    model.create_db_from_input(input_dir)
    orthologs = model.get_orthologs()
    steps.merge_species(input_dir, ortho_dir, orthologs)
    # 2. align ortholog files
    steps.align_orthologs(ortho_dir, aligned_dir, orthologs)
    # 3. trim alignments
    steps.trim_alignments(aligned_dir, trimmed_dir)
    # 4. map trimmed alignments against reference genome
    steps.map_to_reference(trimmed_dir, mapped_dir, reference)
    # 5. design primers
    steps.design_primers(mapped_dir, primer_dir)
    model.load_primers(primer_dir)
    steps.export_primer_alignments(primer_dir)
    # 6. primer BLAST
    steps.blast_primers_online(primer_dir, blast_dir)
    # create report
    steps.create_report(primer_dir, report_dir)
