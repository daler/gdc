#!/usr/bin/python

import argparse
import os
import sys
import gdc

ap = argparse.ArgumentParser()
ap.add_argument('--fasta', help='FASTA file to pull sequences from')
ap.add_argument('--model', help='Model file to create data from')
ap.add_argument('--chrom', help='Chromosome from FASTA to use')
ap.add_argument('--start', default=1, type=int,
                help='Start position of model on chrom. Default is '
                '%(default)s')
ap.add_argument('--scalar', default=1, type=int,
                help='Number of bp each ASCII character represents.  Default '
                'is %(default)s')
ap.add_argument('--readlen', default=10, type=int,
                help='Read length. Default is %(default)s')
ap.add_argument('--debug', action='store_true', help='Enable debug mode')
ap.add_argument('--outdir', default='gdc-data',
                help='Output directory to save results in. Default is '
                     '%(default)s')
ap.add_argument('--prefix', default='gdc',
                help='Prefix for each file.  Default is %(default)s')
ap.add_argument('--demo', action='store_true',
                help='Use an example FASTA file and an example models file. '
                     'All other arguments are honored, with the exception '
                     'of "chrom", which is forced to be "chr2L", the only '
                     'chromosome in the FASTA file.')
args = ap.parse_args()

if args.demo:
    args.fasta = os.path.join(gdc.get_data_dir(), 'dm3.chr2L.oneline.fa')
    args.model = os.path.join(gdc.get_data_dir(), 'example-annotation.txt')
    args.chrom = 'chr2L'

g = gdc.GenomeModel(chrom_start=args.start,
                    chrom=args.chrom,
                    scalar=args.scalar,
                    read_length=args.readlen,
                    debug=args.debug)

models = open(args.model).read().splitlines(True)
g.parse(models)

if not os.path.exists(args.outdir):
    os.system('mkdir -p %s' % args.outdir)

prefix = os.path.join(args.outdir, args.prefix)

open(prefix + '.sam', 'w').write(g.reads.to_sam(args.fasta))
open(prefix + '.fastq', 'w').write(g.reads.to_fastq(args.fasta))
open(prefix + '.gff', 'w').write(g.features.to_gff())
open(prefix + '.model', 'w').write(''.join(models))
open(prefix + '.fasta', 'w').write(open(args.fasta).read())

if args.demo:
    sys.stdout.write('\nSee results, including input model file and fasta '
                     'fasta file, in %s\n' % os.path.abspath(args.outdir))
