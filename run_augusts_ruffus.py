#!/usr/bin/env python

from ruffus import *
import ruffus.cmdline as cmdline

import subprocess
import shlex
import sys
import os.path
import logging
import os
from Bio import SeqIO
import re
import argparse

MODULES_KEY = 'MODULESHOME'
if MODULES_KEY in os.environ:
    modules_init = os.path.join(os.environ[MODULES_KEY], 'init/python.py')
    execfile(modules_init)
    # need this for faToTwoBit
    module('load', 'blat/default')

parser = cmdline.get_argparse(description='Use Ruffus to process .out files from genblastA')
parser.add_argument('--working_directory', '-W', default='.')
parser.add_argument('genome_filename', help='FASTA format genome file, must end in .fa or .fasta')
parser.add_argument('hints_filename', help='Augustus hints filename, generate from exonerate')
args = parser.parse_args()

logger, logger_mutex = cmdline.setup_logging(__name__, args.log_file, args.verbose)


FASTA_RE=r'\.(fa|fasta)$'
os.chdir(args.working_directory)

def safe_open(filename, mode='r'):
    try:
        file_obj = open(filename, mode)
    except IOError as e:
        logger.error('Failed to open {}: {}'.format(filename, str(e)))
        sys.exit(1)
    return file_obj

@transform(args.genome_filename,
           regex(FASTA_RE),
           '.idx')
def make_index(input_file, output_file):
    SeqIO.index_db(output_file, input_file, 'fasta')

@transform(args.hints_filename,
           suffix('.hints'),
           '.sorted.hints')
def sort_hints(input_filename, output_filename):
    cmd_string = 'sort -k1,1 -k4,4n -T . -o {} {}'.format(output_filename, input_filename)
    cmd = shlex.split(cmd_string)
    process = subprocess.check_call(cmd) # throws CalledProcessError on non-zero return code

HINTS_SUFFIX='.contig.hints'
CONTIG_SUFFIX='.contig.fasta'
unitig_re = re.compile('@unitig_(\d+)|quiver')
@follows(make_index)
@transform(sort_hints, suffix('.sorted.hints'), ['contig_list.txt','*'+HINTS_SUFFIX,'*'+CONTIG_SUFFIX], args.genome_filename)
def make_contigs_and_split_hints(input_filename, output_filenames, genome_filename):
    prefix = input_filename.replace('.sorted.hints', '')
    genome_filename = genome_filename.replace(prefix, '')
    index_filename = re.sub(FASTA_RE,'.idx', genome_filename)
    hints_file = open(input_filename)
    contig_list_filename = output_filenames[0].replace(prefix, '')
    contig_list_output_file = open(contig_list_filename,'w')
    genome_dict = SeqIO.index_db(index_filename)
    current_contig = None
    contigs_seen = set()
    for line in hints_file:
        if line.startswith('#'):
            continue
        fields = line.split('\t')
        assert len(fields) == 9, 'invalid hints format, expected 9 fields, got this line: {}'.format(line)
        contig_name = fields[0]
        if contig_name != current_contig:
            contigs_seen.add(contig_name)
            if contig_name.startswith('@unitig'):
                match = unitig_re.match(contig_name)
                output_prefix = 'unitig_' + match.group(1)
            else:
                output_prefix = contig_name
            hints_output_file = open(output_prefix + HINTS_SUFFIX,'w')
            contig_output_file = open(output_prefix + CONTIG_SUFFIX, 'w')
            contig_seq = genome_dict[contig_name]
            SeqIO.write(contig_seq, contig_output_file, 'fasta')
            contig_output_file.close()
            contig_list_output_file.write('\t'.join([contig_name, output_prefix + HINTS_SUFFIX, output_prefix + CONTIG_SUFFIX]) + '\n')
            current_contig = contig_name
        hints_output_file.write(line)
    hints_output_file.close()
    # write out all the contigs for which we have no hints, 
    # along with blank hints files
    for contig_name in genome_dict.keys():
        if not contig_name in contigs_seen:
            #TODO: make this into a function, we're re-using code here
            if contig_name.startswith('@unitig'):
                match = unitig_re.match(contig_name)
                output_prefix = 'unitig_' + match.group(1)
            else:
                output_prefix = contig_name
            hints_output_file = open(output_prefix + HINTS_SUFFIX,'w')
            hints_output_file.close() # write a blank hints file
            contig_output_file = open(output_prefix + CONTIG_SUFFIX, 'w')
            contig_seq = genome_dict[contig_name]
            SeqIO.write(contig_seq, contig_output_file, 'fasta')
            contig_output_file.close()
            contig_list_output_file.write('\t'.join([contig_name, output_prefix + HINTS_SUFFIX, output_prefix + CONTIG_SUFFIX]) + '\n')            
    contig_list_output_file.close()

cmdline.run(args)