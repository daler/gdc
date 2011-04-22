import gdc
import os

here = os.path.dirname(__file__)

# The genome, as a FASTA file with one long line per chrom
fasta_fn = os.path.join(here, 'data/dm3-chr2L.oneline.fa')

# The gene model filename
model_fn = os.path.join(here, 'data/example-annotation.txt')

# First we need to set up a GenomeModel object with the parameters we want
g = gdc.GenomeModel(chrom_start=1, chrom='chr2L', scalar=1, read_length=10, debug=True)

# The .parse() method takes a list of strings (with newlines on the end):
models = open(model_fn).read().splitlines(True)

g.parse(models)

# Now we have a list of features:
g.features

# And a list of reads:
g.reads

