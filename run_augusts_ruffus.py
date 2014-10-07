#!/usr/bin/env python

from ruffus import *
import ruffus.cmdline as cmdline

from subprocess import Popen, PIPE
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
    cmd_string = 'sort -k1,1 -k4,4n -o {} {}'.format(output_filename, input_filename)
    cmd = shlex.split(cmd_string)
    process = subprocess.check_call(cmd) # throws CalledProcessError on non-zero return code

HINTS_SUFFIX='.contig.hints'
CONTIG_SUFFIX='.contig.fasta'
unitig_re = re.compile('@unitig_(\d+)|quiver')
@split([sort_hints, make_index], ['contig_list.txt','*'+HINTS_SUFFIX,'*'+CONTIG_SUFFIX])
def make_contigs_and_split_hints(input_filenames, output_filenames):
    hints_filename = input_filenames[0]
    hints_file = open(index_filename)
    index_filename = intput_filenames[1]
    contig_list_filename = output_filenames[0]
    contig_list_output_file = open(contig_list_filename,'w')
    genome_dict = SeqIO.index_db(index_filename)
    current_contig = None
    for line in hints_file:
        if line.startswith('#'):
            continue
        fields = line.split('\t')
        assert len(fields) == 9, 'invalid hints format, expected 9 fields, got this line: {}'.format(line)
        contig_name = fields[0]
        if contig_name != current_contig:
            hints_output_file = contig_to_output(contig_name, genome_dict, contig_list_output_file)
            if contig_name.startswith('@unitig'):
                match = unitig_re.match(contig_name)
                output_prefix = 'unitig_' + match.group(1)
            else:
                output_prefix = contig_name
            hints_output_file = open(output_prefix + HINTS_SUFFIX)
            contig_output_file = open(output_prefix + CONTIG_SUFFIX)
            contig_seq = genome_dict[contig_name]
            SeqIO.write(contig_output_file, contig_seq, 'fasta')
            contig_output_file.close()
            contig_list_output_file.write('\t'.join([contig_name, output_prefix + HINTS_SUFFIX, output_prefix + CONTIG_SUFFIX]) + '\n')
            current_contig = contig_name
        hints_output_file.write(line)
    hints_output_file.close()

cmdline.run(args)