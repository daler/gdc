Genome Dataset Constructor (``gdc``)
------------------------------------

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

Tutorial coming soon.
