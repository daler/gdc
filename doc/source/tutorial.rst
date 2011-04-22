Tutorial
========

The Genome Dataset Constructor (``gdc``) is a package for making example
genomic data files -- BED, GFF, SAM, FASTQ -- in an easy and intuitive way.
While writing tests for other software, I realized I needed a fairly complex
example data set in order to test all the corner cases.  I quickly realized
that writing such files by hand is itself error-prone.

Enter ``gdc``.  Gene models are constructed using ASCII art.  Reads (e.g., from
Illumina sequencing) can also be specified in the same text file, allowing you
to visually build rich data sets.  ``gdc`` then parses these text files and
converts them into GFF, BED, and -- if you supply a FASTA file -- SAM and FASTQ
as well.

Here is an example of an ASCII gene model that has strand-specific reads in
exons::

    >===|||~~~~~~~||||||==@ #mRNA_t1_g2_+
    $   +++       +++--+    #
    $    ++        +---+--  #
    $   -          --+-++   # 

Running this model through the ``gdc`` will give GFF, SAM, BED, and FASTQ
files, ready for your inspection and incorporation into your tests.  Other
features include support for spliced reads and a rather verbose debug mode so
you can follow exactly how the gene models and reads are being parsed.

Here is a brief tutorial showing how to use ``gdc``, including what the heck
those symbols all mean.

Example data file
-----------------
First let's get the path to an example file that comes with ``gdc``:

.. doctest::
    :options: +NORMALIZE_WHITESPACE

    >>> import gdc
    >>> import os
    >>> datadir = gdc.get_data_dir()
    >>> model_fn = os.path.join(datadir, 'example-annotation.txt')

The file looks like this::


    # file: example-annotation.txt
    #
       >===||||||~~~~|||==@ #mRNA_xs2_g2_+
    $+     -      --     +  #



Symbols used to create models
-----------------------------

The first line contains a gene model, ``>===||||||~~~~|||==@``.  The full set
of symbols (not all are used in this example) are as follows:

======== ============================
symbol   meaning
======== ============================
``>``    left-most end transcript
``=``    3' or 5' UTR
``|``    coding sequence
``~``    intron
``@``    right-most end of transcript
``&``    non-coding exon
======== ============================

Right after the gene model is an annotation for that gene, ``#mRNA_xs2_g2_+``.
It is an underscore-separated definition that is used to construct GFF
features.  The syntax is ``FEATURETYPE_FEATURENAME_PARENT_STRAND``.  So in this
case, it's an ``mRNA`` feature named ``xs2`` (this will eventually go in the ID
field of the GFF line), it's a parent of ``g2``, and it's on the ``+`` strand.

The second line contains the reads.  The full set of symbols for reads (not all
are used in this example) are:

======== =================================
symbol   meaning
======== =================================
``$``    Signals the start of a reads line
``-``    A minus-strand read
``+``    A plus-strand read
``^``    A plus-strand spliced read
``;``    A minus-strand spliced read
======== =================================


Setting up a :class:`GenomeModel` object
----------------------------------------

First, let's set up some parameters for a new :class:`GenomeModel`.  These
parameters configure how the parser will interpret your model; the settings
here will determine the final output.  Here's a description of the arguments
used for :class:`GenomeModel`:

:``chrom_start``:
    The chromosome coordinate that the beginning of the gene model line represents

:``chrom``:
    The chromosome to use.  This will determine the "chrom" field for BED and
    GFF files as well as which chromosome to get the sequence for for FASTQ and
    SAM files.

:``scalar``:
    Each ASCII character in the gene models file represents this many bases.
    This is useful for making very long genes out of a few ASCII characters.
    You can use the ``offset`` parameter to shift gene models and reads by
    values smaller than ``scalar``.

:``read_length``:  The number of bp to make each read.  Spliced reads will have fragments that
    are  ``int(floor(read_length/2.0))`` and ``int(ceil(read_length/2.0))``.

Let's use some straightforward parameters to set up our :class:`GenomeModel`.
These are probably not the values you would want to use for your own data sets
(read sizes are very small, the gene will be less than 20 bp long) but these
values will be good for illustration:

.. doctest::

    >>> g = gdc.GenomeModel(chrom_start=1,
    ...                     chrom='chr2L',
    ...                     scalar=1, 
    ...                     read_length=3)


Now that the parameters are set, we need to pass ``g`` the model.
:class:`GenomeModel` objects take a list of lines as input to their
:meth:`.parse` method.  These are typically created from the input file like this:

.. doctest::

    >>> models = open(model_fn).read().splitlines(True)

Next, we pass them to the :class:`GenomeModel`:

.. doctest::

    >>> g.parse(models)

Now ``g`` contains a list of features and a list of reads. The gene model
``>===||||||~~~~|||==@ #mRNA_xs2_g2_+`` became this list of features:

.. doctest:: 

    >>> for feature in g.features:
    ...     print feature
    Feature UTR 'mRNA:xs2:UTR:5-7': chr2L:5-7 (+)
    Feature CDS 'mRNA:xs2:CDS:8-13': chr2L:8-13 (+)
    Feature intron 'mRNA:xs2:intron:14-17': chr2L:14-17 (+)
    Feature CDS 'mRNA:xs2:CDS:18-20': chr2L:18-20 (+)
    Feature UTR 'mRNA:xs2:UTR:21-22': chr2L:21-22 (+)
    Feature exon 'mRNA:xs2:exon:5-13': chr2L:5-13 (+)
    Feature exon 'mRNA:xs2:exon:18-22': chr2L:18-22 (+)
    Feature mRNA 'mRNA:xs2': chr2L:5-22 (+)
    Feature gene 'g2': chr2L:5-22 (+)


Similarly, the reads defined in the line ``$+     -      --     +  #`` became
this list of reads:

.. doctest::

    >>> for read in g.reads:
    ...     print read
    Feature read 'None': chr2L:2-4 (+)
    Feature read 'None': chr2L:8-10 (-)
    Feature read 'None': chr2L:15-17 (-)
    Feature read 'None': chr2L:16-18 (-)
    Feature read 'None': chr2L:22-24 (+)


Since we specified ``scalar=1``, each character corresponds to exactly one
base.  That means we add some numbers to show the coordinates.  Here are the 1-based coordinates for SAM and GFF::

    # ________ 1-based coords_______________
    #       10        20        30        40
    #
    12345678901234567890123456789012345678901
       >===||||||~~~~|||==@ #mRNA_xs2_g2_+
    $+     -      --     +  #


For BED files, which are 0-based::

    # ________ 0-based coords_______________
    #        10        20        30        40
    #
    01234567890123456789012345678901234567890
       >===||||||~~~~|||==@ #mRNA_xs2_g2_+
    $+     -      --     +  #

.. note::

    The :class:`Feature` objects (i.e. ``feature.start``, ``feature.stop``)
    coordinates are 1-based.

We can compare the output below with these coordinates to make sure the files
are being created correctly.

Making a GFF file
-----------------

It's easy to make a GFF file based on our model (note that the GFF lines
include a newline at the end):

.. doctest::
    :options: +NORMALIZE_WHITESPACE

    >>> print g.features.to_gff()
    chr2L	.	UTR	5	7	0	+	.	ID=mRNA:xs2:UTR:5-7;Parent=mRNA:xs2;
    chr2L	.	CDS	8	13	0	+	.	ID=mRNA:xs2:CDS:8-13;Parent=mRNA:xs2;
    chr2L	.	intron	14	17	0	+	.	ID=mRNA:xs2:intron:14-17;Parent=mRNA:xs2;
    chr2L	.	CDS	18	20	0	+	.	ID=mRNA:xs2:CDS:18-20;Parent=mRNA:xs2;
    chr2L	.	UTR	21	22	0	+	.	ID=mRNA:xs2:UTR:21-22;Parent=mRNA:xs2;
    chr2L	.	exon	5	13	0	+	.	ID=mRNA:xs2:exon:5-13;Parent=mRNA:xs2;
    chr2L	.	exon	18	22	0	+	.	ID=mRNA:xs2:exon:18-22;Parent=mRNA:xs2;
    chr2L	.	mRNA	5	22	0	+	.	ID=mRNA:xs2;Parent=g2;
    chr2L	.	gene	5	22	0	+	.	ID=g2;


Even though we only specified a single mRNA, its parent was indicated in the
annotation, so the gene is implied.  Each feature type was also extracted from
the model and put into its own GFF line.

Making a BED file
-----------------
We can do the same thing to get a BED file.  Note the 0-based coordinates;
features defined by BED files do not include the last base.

.. doctest::
    :options: +NORMALIZE_WHITESPACE

    >>> print g.features.to_bed()
    chr2L	4	7	mRNA:xs2:UTR:5-7	0	+
    chr2L	7	13	mRNA:xs2:CDS:8-13	0	+
    chr2L	13	17	mRNA:xs2:intron:14-17	0	+
    chr2L	17	20	mRNA:xs2:CDS:18-20	0	+
    chr2L	20	22	mRNA:xs2:UTR:21-22	0	+
    chr2L	4	13	mRNA:xs2:exon:5-13	0	+
    chr2L	17	22	mRNA:xs2:exon:18-22	0	+
    chr2L	4	22	mRNA:xs2	0	+
    chr2L	4	22	g2	0	+
    <BLANKLINE>

Individual :class:`Feature` objects also have their own :meth:`.to_bed` and
:meth:`.to_gff` methods for more fine-scale control.


Making a FASTQ file
-------------------

In order to get sequences for SAM and FASTQ files created from :attr:`g.reads`, we need to specify a
genome.  This is currently done by providing a fasta file with the sequence of
each chromosome on its own line.  Here we'll use an example that comes with ``gdc``:

.. doctest::

    >>> fasta_fn = os.path.join(datadir, 'dm3.chr2L.oneline.fa')
    >>> print g.reads.to_fastq(fasta_fn)
    @None
    GAC
    +None
    III
    @None
    TGC
    +None
    III
    @None
    TCT
    +None
    III
    @None
    CTC
    +None
    III
    @None
    GCA
    +None
    III
    <BLANKLINE>

Currently, there's no mechanism for creating quality scores, so here they're
all uniform.


Making a SAM file
-----------------

Similarly, we can get a SAM file:

.. doctest::
    :options: +NORMALIZE_WHITESPACE

    >>> print g.reads.to_sam(fasta_fn)
    None	0	chr2L	2	255	3M	*	0	0	GAC	III	NM:i:0	NH:i:1
    None	16	chr2L	8	255	3M	*	0	0	TGC	III	NM:i:0	NH:i:1
    None	16	chr2L	15	255	3M	*	0	0	TCT	III	NM:i:0	NH:i:1
    None	16	chr2L	16	255	3M	*	0	0	CTC	III	NM:i:0	NH:i:1
    None	0	chr2L	22	255	3M	*	0	0	GCA	III	NM:i:0	NH:i:1
    <BLANKLINE>


Debug mode
----------
If you specify `debug=True` when creating a :class:`GenomeModel`, you'll get
lots of output when parsing.  This is useful to confirm that the parser is
seeing everything you think it should be seeing: 

.. doctest::

    >>> g = gdc.GenomeModel(chrom_start=101,
    ...                     chrom='chr2L',
    ...                     scalar=10, 
    ...                     read_length=36,
    ...                     debug=True)
    >>> g.parse(models)
    counter: 4 chars:     start: 101 stop: 130
    counter: 5 chars: > start: 131 stop: 140
    counter: 8 chars: === start: 141 stop: 170
    counter: 14 chars: |||||| start: 171 stop: 230
    counter: 18 chars: ~~~~ start: 231 stop: 270
    counter: 21 chars: ||| start: 271 stop: 300
    counter: 23 chars: == start: 301 stop: 320
    counter: 24 chars: @ start: 321 stop: 330
    counter: 25 chars:   start: 331 stop: 340
    counter: 38 chars: #mRNA_xs2_g2_+ start: 331 stop: 470
    counter: 5 chars:    > start: 101 stop: 140
    counter: 14 chars: ===|||||| start: 141 stop: 230
    counter: 18 chars: ~~~~ start: 231 stop: 270
    counter: 23 chars: |||== start: 271 stop: 320
    counter: 38 chars: @ #mRNA_xs2_g2_+ start: 311 stop: 470
    counter: 2 chars: $ start: 101 stop: 110
    counter: 3 chars: + start: 111 stop: 120
    counter: 8 chars:       start: 121 stop: 170
    counter: 9 chars: - start: 171 stop: 180
    counter: 15 chars:        start: 181 stop: 240
    counter: 16 chars: - start: 241 stop: 250
    counter: 17 chars: - start: 251 stop: 260
    counter: 22 chars:       start: 261 stop: 310
    counter: 23 chars: + start: 311 stop: 320
    counter: 25 chars:    start: 321 stop: 340
    counter: 25 chars: # start: 331 stop: 340

Advanced usage
--------------
To be written.  Will describe things like ``offset_plus_10``, making spliced
reads, combining multiple :class:`GenomeModel` objects together, noncoding RNA,
and generic features.
