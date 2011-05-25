#!/usr/bin/python

import argparse
import os
import sys
import gdc

default_fasta = os.path.join(gdc.get_data_dir(), 'dm3.chr2L.oneline.fa')
default_model = os.path.join(gdc.get_data_dir(), 'example-annotation.txt')

ap = argparse.ArgumentParser()
ap.add_argument('--fasta', default=default_fasta,
                help='FASTA file to pull sequences from')
ap.add_argument('--model', default=default_model,
                help='Model file to create data from')
ap.add_argument('--chrom', default='chr2L',
                help='Chromosome from FASTA to use')
ap.add_argument('--start', default=1, type=int,
                help='Start position of model on chrom. Default is '
                '%(default)s')
ap.add_argument('--scalar', default=1, type=int,
                help='Number of bp each ASCII character represents.  Default '
                'is %(default)s')
ap.add_argument('--readlen', default=10, type=int,
                help='Read length. Default is %(default)s')
ap.add_argument('--debug', action='store_true',
                help='Enable debug mode')
ap.add_argument('--outdir', default='gdc-data',
                help='Output directory to save results in. Default is '
                     '%(default)s')
ap.add_argument('--prefix', default='gdc',
                help='Prefix for each file.  Default is %(default)s')

args = ap.parse_args()

params = dict(chrom_start=args.start,
              chrom=args.chrom,
              scalar=args.scalar,
              read_length=args.readlen,
              debug=args.debug)

g = gdc.GenomeModel(**params)

models = open(args.model).read().splitlines(True)
g.parse(models)

if not os.path.exists(args.outdir):
    os.system('mkdir -p %s' % args.outdir)

prefix = os.path.join(args.outdir, args.prefix)

open(prefix + '.sam', 'w').write(g.reads.to_sam(args.fasta))
open(prefix + '.fastq', 'w').write(g.reads.to_fastq(args.fasta))
open(prefix + '.gff', 'w').write(g.features.to_gff())
open(prefix + '.model', 'w').write(''.join(models))

params.update(args.__dict__)
params_string = ['='.join((key, str(val))) for key, val in params.items()]
params_string = '#' + '\n#'.join(params_string)
open(prefix + '.params', 'w').write(params_string)

# TODO: once GFFutils dependency removed, this should be unnecessary
fasta_out = open(prefix + '.fasta', 'w')
from Bio import SeqIO
p = SeqIO.parse(args.fasta, 'fasta')
fasta_out.write(p.next().format('fasta'))
fasta_out.close()

os.system('samtools view -S -b -T %s %s > %s' % (fasta_out.name,
                                                 prefix + '.sam',
                                                 prefix + '.bam'))
