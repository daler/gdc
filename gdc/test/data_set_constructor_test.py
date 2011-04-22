
from gdc import *
import os
import GFFutils
import nose.tools as nt
import operator

# numbers in parentheses are GFF coords; 
# other numbers are python string indices.
gene_models = """
>=|~|=@ #mRNA_xs1_g1_+
$+      #
$ ^^^   #
#||
#||
#||
#||
#|CDS 2(3) 2(3)
#UTR 1(2) 1(2)
#
#
"""

# Simple model to test individual reads
gene_models2 = """
   >===||||||~~~~|||==@ #mRNA_xs2_g2_+
$+     -      --     +  #  
#^     ^      ^^     ^
#|     |      ||     |
#|     |      ||     \ 21(22) - 23(24)
#|     |      |\ 15(16) - 17(18)
#|     |      \ 14(15) - 16(17)
#|     \ 7(8) - 9(10) 
#\ 1(2) - 3(4)
#
"""

gene_models3 = """
   >===||||||~~~~|||==@ #mRNA_xs2_g2_+
$           ^^^^^^      #
$           ;;;;;;      #
#           ^
#           |
#           |
#           \ 12(13) - 17(18)
#             (scalar=5 -> 
"""

gene_models4 = """
   >===||||||~~~~|||==@ #mRNA_xs2_g2_+
$     -+    ^^^^^^      #
$           ;;;;;;      #
#           ^
#           |
#           |
#           \ 12(13) - 17(18)
#             (scalar=5 -> 
"""

gene_models5 = """
   >===||||||~~~~|||==@ #mRNA_xs2_g2_+
$     -+    ^^^^^^      #offset_plus_17
$           ;;;;;;      #offset_minus_2
#           ^
#           |
#           |
#           \ 12(13) - 17(18)
#             (scalar=5 -> 
"""

gene_models6 = """
   >===||||||~~~~|||==@ #mRNA_xs4_g000_+
$     -+    ^^^^^^      #offset_plus_17
$           ;;;;;;      #offset_minus_2
#           ^
#           |
#           |
#           \ 12(13) - 17(18)
#             (scalar=5 -> 
"""

gene_models = gene_models.splitlines(True)
g = GenomeModel(chrom_start=1,scalar=1,read_length=3, debug=False)
g.parse(gene_models)
g2 = GenomeModel(chrom_start=1,scalar=5,read_length=3, debug=False)
g2.parse(gene_models)
g3 = GenomeModel(chrom_start=11,scalar=5,read_length=3, debug=False)
g3.parse(gene_models)

gene_models2 = gene_models2.splitlines(True)
g4 = GenomeModel(chrom_start=1, scalar=1, read_length=3, debug=False)
g4.parse(gene_models2)
g5 = GenomeModel(chrom_start=1, scalar=5, read_length=3, debug=False)
g5.parse(gene_models2)
g6 = GenomeModel(chrom_start=11, scalar=5, read_length=3, debug=False)
g6.parse(gene_models2)

gene_models3 = gene_models3.splitlines(True)
g7 = GenomeModel(chrom_start=1, scalar=1, read_length=3, debug=False)
g7.parse(gene_models3)
g8 = GenomeModel(chrom_start=1, scalar=5, read_length=3, debug=False)
g8.parse(gene_models3)
g9 = GenomeModel(chrom_start=11, scalar=5, read_length=3, debug=False)
g9.parse(gene_models3)

gene_models4 = gene_models4.splitlines(True)
g10 = GenomeModel(chrom_start=11, scalar=5, read_length=10, debug=False)
g10.parse(gene_models4)

gene_models5 = gene_models5.splitlines(True)
g11 = GenomeModel(chrom_start=1, scalar=5, read_length=3, debug=False)
g11.parse(gene_models5)

gene_models6 = gene_models6.splitlines(True)
g12 = GenomeModel(chrom='chr3R', chrom_start=101, scalar=5, read_length=3, debug=False)
g12.parse(gene_models6)

here = os.path.dirname(__file__)
genome = GFFutils.Genome(os.path.join(here, 'data/dm3.chr2L.oneline.fa'))

def feature_exists(genome_model_obj, start, stop,
                      featuretype, chrom='chr2L',strand='+'):
    # Checks to see if a feature exists.  Does not check names, only genomic
    # coords and featuretype
    for feature in genome_model_obj.features:
        if (feature.start == start) and (feature.stop == stop) \
           and (feature.chrom == chrom) and (feature.strand == strand) \
           and (feature.featuretype == featuretype):
           return True
    return False    

def read_exists(genome_model_obj, start, stop,
                      featuretype, chrom='chr2L',strand='+'):
    # Checks to see if a feature exists.  Does not check names, only genomic
    # coords and featuretype
    for feature in genome_model_obj.reads:
        if (feature.start == start) and (feature.stop == stop) \
           and (feature.chrom == chrom) and (feature.strand == strand) \
           and (feature.featuretype == featuretype):
           return True
    return False    

def spliced_read_exists(genome_model_obj,start1,stop1,start2,stop2,featuretype, chrom='chr2L', strand='+'):
    for read in genome_model_obj.reads:
        if (read.start == start1) and \
         (read.stop == stop2) and \
         (read.start1 == start1) and \
         (read.start2 == start2) and \
         (read.stop2 == stop2) and \
         (read.featuretype == featuretype) and \
         (read.chrom == chrom) and \
         (read.strand == strand):
            return True
    return False

def test_splice_reads():
    assert spliced_read_exists(g2,start1=15,stop1=15,start2=21,stop2=22, featuretype='spliced_read', strand='+')

    assert spliced_read_exists(g7,start1=13,stop1=13,start2=18,stop2=19, featuretype='spliced_read', strand='+')
    assert spliced_read_exists(g7,start1=13,stop1=13,start2=18,stop2=19, featuretype='spliced_read', strand='-')

    assert spliced_read_exists(g8,start1=65,stop1=65,start2=86,stop2=87, featuretype='spliced_read', strand='-')
    assert spliced_read_exists(g8,start1=65,stop1=65,start2=86,stop2=87, featuretype='spliced_read', strand='+')
     
    assert spliced_read_exists(g9,start1=75,stop1=75,start2=96,stop2=97, featuretype='spliced_read', strand='-')
    assert spliced_read_exists(g9,start1=75,stop1=75,start2=96,stop2=97, featuretype='spliced_read', strand='+')

def test_splices_to_bed():
    expected = "chr2L\t74\t97\tNone\t0\t+\t74\t97\t0,0,0\t2\t1,2\t0,21\n"
    observed = g9.reads[0].to_bed()
    print 'expected:',expected
    print 'observed:',observed
    assert expected == observed

def test_fastq():
    
    # single + read
    read = g10.reads[0]
    print read
    expected =  '@None\n'
    expected += 'ATGAGAGGCA\n'
    expected += '+None\n'
    expected += 'IIIIIIIIII\n'
    observed = read.to_fastq(genome)
    print 'expected:\n',expected
    print 'observed:\n',observed
    assert expected == observed
    
    # single - read
    read = g10.reads[1]
    print read
    expected =  '@None\n'
    expected += 'CTCATTTTCT\n'
    expected += '+None\n'
    expected += 'IIIIIIIIII\n'
    observed = read.to_fastq(genome)
    print 'expected:\n',expected
    print 'observed:\n',observed
    assert expected == observed

    # spliced +
    read = g10.reads[2]
    print read
    expected =  '@None\n'
    expected += 'GAGAAGTAGT\n'
    expected += '+None\n'
    expected += 'IIIIIIIIII\n'
    observed = read.to_fastq(genome)
    print 'expected:\n',expected
    print 'observed:\n',observed
    assert expected == observed

    # spliced -
    read = g10.reads[3]
    print read
    expected =  '@None\n'
    expected += 'ACTACTTCTC\n'
    expected += '+None\n'
    expected += 'IIIIIIIIII\n'
    observed = read.to_fastq(genome)
    print 'expected:\n',expected
    print 'observed:\n',observed
    assert expected == observed

def test_reads():
    assert read_exists(g2, start=6, stop=8, strand='+', featuretype='read')
    
    assert read_exists(g4,start=2, stop=4, strand='+', featuretype='read')
    assert read_exists(g4,start=8, stop=10, strand='-', featuretype='read')
    assert read_exists(g4,start=15, stop=17, strand='-', featuretype='read')
    assert read_exists(g4,start=16, stop=18, strand='-', featuretype='read')
    assert read_exists(g4,start=22, stop=24, strand='+', featuretype='read')
    
    assert read_exists(g5,start=6, stop=8, strand='+', featuretype='read')
    assert read_exists(g5,start=36, stop=38, strand='-', featuretype='read')
    assert read_exists(g5,start=71, stop=73, strand='-', featuretype='read')
    assert read_exists(g5,start=76, stop=78, strand='-', featuretype='read')
    assert read_exists(g5,start=106, stop=108, strand='+', featuretype='read')

    assert read_exists(g6,start=16, stop=18, strand='+', featuretype='read')
    assert read_exists(g6,start=46, stop=48, strand='-', featuretype='read')
    assert read_exists(g6,start=81, stop=83, strand='-', featuretype='read')
    assert read_exists(g6,start=86, stop=88, strand='-', featuretype='read')
    assert read_exists(g6,start=116, stop=118, strand='+', featuretype='read')

def test_num_features():
    assert len(list(g.features)) == 9
    assert len(filter(lambda x: x.featuretype=='gene', g.features)) == 1
    assert len(filter(lambda x: x.featuretype=='exon', g.features)) == 2
    assert len(filter(lambda x: x.featuretype=='UTR', g.features)) == 2
    assert len(filter(lambda x: x.featuretype=='intron', g.features)) == 1
    assert len(filter(lambda x: x.featuretype=='mRNA', g.features)) == 1
    assert len(filter(lambda x: x.featuretype=='CDS', g.features)) == 2

def test_features_exist():
    assert feature_exists(g,start=2, stop=2, featuretype='UTR')
    assert feature_exists(g,start=3, stop=3, featuretype='CDS')
    assert feature_exists(g,start=4, stop=4, featuretype='intron')
    assert feature_exists(g,start=5, stop=5, featuretype='CDS')
    assert feature_exists(g,start=6, stop=6, featuretype='UTR')
    assert feature_exists(g,start=2, stop=6, featuretype='mRNA')
    assert feature_exists(g,start=2, stop=6, featuretype='gene')
    assert feature_exists(g,start=2, stop=3, featuretype='exon')
    assert feature_exists(g,start=5, stop=6, featuretype='exon')
    
    assert feature_exists(g2,start=6, stop=10, featuretype='UTR')
    assert feature_exists(g2,start=11, stop=15, featuretype='CDS')
    assert feature_exists(g2,start=16, stop=20, featuretype='intron')
    assert feature_exists(g2,start=21, stop=25, featuretype='CDS')
    assert feature_exists(g2,start=26, stop=30, featuretype='UTR')
    assert feature_exists(g2,start=6, stop=30, featuretype='mRNA')
    assert feature_exists(g2,start=6, stop=30, featuretype='gene')
    assert feature_exists(g2,start=6, stop=15, featuretype='exon')
    assert feature_exists(g2,start=21, stop=30, featuretype='exon')

    assert feature_exists(g3,start=16, stop=20, featuretype='UTR')
    assert feature_exists(g3,start=21, stop=25, featuretype='CDS')
    assert feature_exists(g3,start=26, stop=30, featuretype='intron')
    assert feature_exists(g3,start=31, stop=35, featuretype='CDS')
    assert feature_exists(g3,start=36, stop=40, featuretype='UTR')
    assert feature_exists(g3,start=16, stop=40, featuretype='mRNA')
    assert feature_exists(g3,start=16, stop=40, featuretype='gene')
    assert feature_exists(g3,start=16, stop=25, featuretype='exon')
    assert feature_exists(g3,start=31, stop=40, featuretype='exon')

def test_add():
    nt.assert_raises(ValueError, operator.add, *(g11, g10))
    g00 = g11 + g12
    print 'still need a good way to test this...'

def test_feature_to_bed():
    feature = Feature('chr2L', start=1, stop=100)
    observed = feature.to_bed()
    expected = 'chr2L\t0\t100\tNone\t0\tNone\n'
    print 'observed:\n',observed
    print 'expected:\n',expected
    assert observed == expected

def test_feature_to_sam():
    
    # single + read
    read = g10.reads[0]
    print read
    expected = [None, 16, 'chr2L', 41, 255, '10M', "*", 0, 0, 'ATGAGAGGCA', 'IIIIIIIIII', 'NM:i:0\tNH:i:1']
    expected = map(str,expected)
    expected = '\t'.join(expected)+'\n'
    observed = read.to_sam(genome)
    print 'expected:\n',expected
    print 'observed:\n',observed
    assert expected == observed

class TestDataList(object):
    def setup(self):
        feature1 = Feature(chrom='chr4', start=1, stop=10, strand='+')
        feature2 = Feature(chrom='chrX', start=50, stop=200, strand='-')
        self.datalist = DataList([feature1, feature2])
   
    def test_to_gff(self):
       expected = ''
       expected += 'chr4\t.\t.\t1\t10\t0\t+\t.\n'
       expected += 'chrX\t.\t.\t50\t200\t0\t-\t.\n'
       observed = self.datalist.to_gff()
       print 'observed:\n',observed
       print 'expected:\n',expected
       assert observed == expected
    
    def test_reprs(self):
        assert str(self.datalist) == '<DataList instance with 2 items>'

    def test_to_bed(self):
        expected = ''
        expected += 'chr4\t0\t10\tNone\t0\t+\n'
        expected += 'chrX\t49\t200\tNone\t0\t-\n'
        observed = self.datalist.to_bed()
        print 'observed:\n',observed
        print 'expected:\n',expected
        assert observed == expected


def test_offset():
    assert spliced_read_exists(g11,start1=63,stop1=63,start2=84,stop2=85, featuretype='spliced_read', strand='-')
    assert spliced_read_exists(g11,start1=82,stop1=82,start2=103,stop2=104, featuretype='spliced_read', strand='+')
     
