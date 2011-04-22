import itertools
import GFFutils
import os
import sys
from math import floor, ceil
"""
Module for creating a data set for testing.  Starts with an ASCII repsentation
of gene models and reads, and creates:

    * GFF files
    * FASTQ files
    * SAM files
    * intermediate data files

Reads are placed at the BEGINNING of a bin on the plus strand and at the END of
a bin on the minus strand

Spliced reads will be placed at the END of the bin on the 5' end and BEGINNING
of the bin on the 3' end.  This makes sure they line up with exon boundaries.

Annotations for features are of the form "#featuretype_name_parentgene".  This
allows you to create tRNA and rRNA features.  If there are CDSs present in the
feature, then exons will be created out of UTRs and CDSs. mRNAs should be
explicitly described as such in the annotation string.

   |||  start=3, stop=6;  len=3
0123456789 # bed coords

   |||  start=4. stop=6;  len=3
1234567890 # gff coords
"""


# Legend for gene models.
legend = {'>':'start',
          '<':'start',
          '=':'UTR',
          '~':'intron',
          '|':'CDS',
          ' ':'spacer',
          '@':'stop',
          '*':'repeat',
          '&':'noncoding_exon',
          '\n':'line end'}

"""
API should be something like

    >>> g = GenomeModel(gene_models)
    >>> g.to_bed()
    >>> g.to_gff()
    >>> g.to_fastq(fasta='test.fa')
    >>> g.to_sam()

"""

def get_data_dir():
    """
    Returns directory of installed data
    """
    here = os.path.abspath(os.path.dirname(__file__))
    return os.path.join(here, 'test/data')

class Feature(GFFutils.GFFFeature):
    def __init__(self,*args,**kwargs):
        """
        start, stop are in GFF coordinate space
        """
        GFFutils.GFFFeature.__init__(self,*args,**kwargs)
        self.value = 0

    def to_bed(self, nfields=6):
        valid_fields = [3, 6, 9, 12]
        if nfields not in valid_fields:
            raise ValueError('Invalid number of fields, must be one of %s' % valid_fields)
        fields = [self.chrom,
                  self.start-1,
                  self.stop,
                  self.id,
                  self.value,
                  self.strand,
                  self.start-1,
                  self.stop,
                  '0,0,0',
                  1,
                  len(self),
                  0]
        return '\t'.join(map(str,fields[:nfields]))+'\n'

    def to_fastq(self, genome):
        """
        *genome* is a GFFutils.Genome instance.  Returns the sequence for this
        feature as a FASTQ, with fake quality scores.
        """
        s = ''
        s += '@%s\n' % self.id
        seq = genome.sequence(chrom=self.chrom, start=self.start, stop=self.stop, strand=self.strand)
        s += seq+'\n'
        s += '+%s\n' % self.id
        s += 'I' * len(seq) + '\n' 
        return s

    def add_offset(self,offset):
        self.start += offset
        self.stop += offset

    def to_sam(self, genome):
        """
        Writes feature as a SAM line.  Needs a genome to get the sequence from;
        *genome* is a GFFutils.Genome instance.  Fake quality scores are
        created.
        """
        s = ''
        seq = genome.sequence(chrom=self.chrom, start=self.start, stop=self.stop, strand=self.strand)
        qual = 'I' * len(seq)
        if self.strand == '-':
            flag = 16
        else:
            flag = 0
        line = [None, flag, self.chrom, self.start, 255, str(len(self))+'M', '*', 0, 0, seq, qual, 'NM:i:0\tNH:i:1']
        line = map(str,line)
        return '\t'.join(line)+'\n'

    def to_gff(self):
        return self.tostring()

class SplicedRead(Feature):
    def __init__(self,*args,**kwargs):
        """
        All the +1/-1 stuff was figured out by testing...
        """
        self.read_length = kwargs.pop('read_length')
        self.scalar = kwargs.pop('scalar')
        Feature.__init__(self,*args,**kwargs)

        # If read length is odd, then first fragment will be one bp smaller.
        len1 = self.read_length/2
        len2 = self.read_length - len1

        self.stop1 = self.start + self.scalar - 1
        self.start1 = self.stop1 - len1 + 1
        self.start2 = self.stop - self.scalar + 1
        self.stop2 = self.start2 + len2 - 1

        self.start = self.start1
        self.stop = self.stop2

        self.featuretype = 'spliced_read'

    def to_fastq(self,genome):
        if self.strand == '+':
            s = ''
            s += '@%s\n'%self.id
            seq1 = genome.sequence(chrom=self.chrom, start=self.start1, stop=self.stop1, strand=self.strand)
            seq2 = genome.sequence(chrom=self.chrom, start=self.start2, stop=self.stop2, strand=self.strand)
            seq = seq1+seq2
            s += seq+'\n'
            s += '+%s\n' % self.id
            s += 'I' * len(seq) + '\n' 
            return s

        if self.strand == '-':
            s = ''
            s += '@%s\n'%self.id
            seq1 = genome.sequence(chrom=self.chrom, start=self.start1, stop=self.stop1, strand=self.strand)
            seq2 = genome.sequence(chrom=self.chrom, start=self.start2, stop=self.stop2, strand=self.strand)
            seq = seq2+seq1
            s += seq+'\n'
            s += '+%s\n' % self.id
            s += 'I' * len(seq) + '\n' 
            return s

    def to_bed(self):
        # only supporting 2 frags for now.
        block_count = 2

        block_sizes = [ self.stop1-self.start1+1, self.stop2-self.start2+1 ] 
        block_sizes = ','.join(map(str,block_sizes))

        block_starts = ','.join(map(str, [self.start1-self.start, self.start2-self.start]))

        fields = [ self.chrom,
                   self.start-1,
                   self.stop,
                   None,
                   0,
                   self.strand,
                   self.start-1,
                   self.stop,
                   '0,0,0',
                   block_count,
                   block_sizes,
                   block_starts]
        return '\t'.join(map(str,fields))+'\n' 

    def add_offset(self,offset):
        self.start1 += offset
        self.start2 += offset
        self.stop1 += offset
        self.stop2 += offset
        self.start += offset
        self.stop += offset

    def __repr__(self):
        return '<SplicedRead %s:[%s-%s, %s-%s](%s)>' % (self.chrom,self.start1, self.stop1, self.start2, self.stop2, self.strand)

class DataList(object):
    def __init__(self,items):
        self.items = items

    def __str__(self):
        return "<DataList instance with %s items>" % len(self.items)

    def __repr__(self):
        s = str(self)+'\n'
        s += '\n'.join('\t'+repr(i) for i in self.items)
        return s

    def to_gff(self):
        """
        Prints items as GFF lines
        """
        s = []
        for item in self.items:
            s.append(item.tostring())
        return ''.join(s)

    def to_bed(self):
        """
        Prints items as BED lines
        """
        s = []
        for item in self.items:
            s.append(item.to_bed())
        return ''.join(s)

    def to_sam(self, fasta):
        """
        Prints items as SAM lines
        """
        s = []
        genome = GFFutils.Genome(fasta)
        for item in self.items:
            s.append(item.to_sam(genome))
        return ''.join(s)

    def __getitem__(self,ind):
        return self.items[ind]

    def to_fastq(self,fasta):
        """
        Creates sequences and fake quality scores.  Sequence names are the same
        as the GFFFeature.id.
        """
        genome = GFFutils.Genome(fasta)
        s = []
        for item in self.items:
            s.append(item.to_fastq(genome))
        return ''.join(s)

class GenomeModel(object):
    def __init__(self,chrom_start=1,chrom='chr2L',scalar=1,read_length=10,debug=False):
        """
        Class to represent the whole genome.  Will have at least one GeneModel.

        chrom_start is in GFF coords.

        """
        self.debug = debug
        self._features = []
        self._reads = []
        self.genes = {}
        self.chrom_start = chrom_start
        self.scalar = scalar
        self.chrom = chrom
        self.read_length = read_length
        self.composite = False

    def __add__(self,other):
        """
        Does not do any checking in overlaps.   So it's possible to overwrite
        gene models
        """
        if not isinstance(other,GenomeModel):
            raise ValueError, "Adding GenomeModel to a non-GenomeModel instance not supported"
        self_ids = set([i.id for i in self.features])
        other_ids = set([i.id for i in other.features])
        if len(self_ids.intersection(other_ids)) > 0:
            raise ValueError, "Feature ID collision -- make sure that features in GenomeModels being added are uniquely named"

        # would be nice to add something here to indicate which features overlap.

        new_GenomeModel = GenomeModel()
        new_GenomeModel._features.extend(self._features)
        new_GenomeModel._features.extend(other._features)
        new_GenomeModel._reads.extend(self._reads)
        new_GenomeModel._reads.extend(other._reads)
        self.composite = True
        return new_GenomeModel

    def parse(self, gene_models):
        if self.composite:
            raise ValueError, "parsing a composite, already-merged GenomeModel is not supported"
        for model in gene_models:
            if len(model.strip()) < 1:
                continue
            if model.startswith('#'):
                continue
            if model.startswith('$'):
                self.parse_reads(model)
            else:
                self.parse_model(model)

    def start_stop(self,chars):
        gff_stop = (self.counter - 2) * self.scalar + self.scalar
        gff_start = gff_stop - ((len(chars) * self.scalar) -1)

        gff_start += self.chrom_start - 1
        gff_stop += self.chrom_start - 1

        if self.debug:
            print 'counter:',self.counter,
            print 'chars:', ''.join(chars),
            print 'start:',gff_start, 'stop:', gff_stop
        sys.stdout.flush()
        return gff_start, gff_stop

    def parse_model(self,model):
        model = model.rstrip()
        featuretypes_seen = set()

        def legend_func(x):
            """
            Function for groupby() that returns the featuretype.  Anything
            that's not in the legend is considered an annotation.
            """
            self.counter += 1
            try:
                return legend[x]
            except KeyError:
                return 'annotation'

        def exon_func(x):
            """
            Returns an 'exon' featuretype.
            """
            self.counter += 1
            try:
                ft = legend[x]
                if ft in ['CDS','UTR','noncoding_exon']:
                    return 'exon'
                else:
                    return None
            except KeyError:
                return None

        featuretypes_to_ignore = [None, 'spacer','annotation','start','stop','noncoding_exon']
        features = []

        # Parse the model.
        self.counter = 0
        for featuretype, chars in itertools.groupby(model,legend_func):

            # Keep track of what we've seen
            featuretypes_seen = featuretypes_seen.union([featuretype])

            # for later manipulation, convert the generator into a list
            chars = list(chars)

            # some output
            start,stop = self.start_stop(chars)

            if featuretype == 'annotation':
                annotation = ''.join(chars)

            # Don't create Features for some types.
            if featuretype in featuretypes_to_ignore:
                continue

            # If you got here, it should be a real Feature to create.
            feature = Feature(chrom=self.chrom,
                              start=start,
                              stop=stop,
                              featuretype=featuretype)
            features.append(feature)

        # Re-parse the model, this time looking for exon features.
        self.counter = 0
        for featuretype, chars in itertools.groupby(model, exon_func):
             chars = list(chars)
             start,stop = self.start_stop(chars)
             if featuretype in featuretypes_to_ignore:
                continue
             feature = Feature(chrom=self.chrom,
                               start=start,
                               stop=stop,
                               featuretype=featuretype)
             features.append(feature)

        # Now parse the annotation.  This will get you:
        #   - the strand (which will be used to update all features)
        #   - the transcript (rRNA, mRNA, whatever)
        #   - the parent gene

        annotation = annotation.replace('#','')
        annotation = annotation.strip().split('_')
        transcript_featuretype = annotation[0]
        transcript_name = annotation[1]
        parent_gene = annotation[2]
        strand = annotation[3]

        # get the transcript start and stop positions.
        transcript_start = min(i.start for i in features)
        transcript_stop = max(i.stop for i in features)
        feature = Feature(chrom=self.chrom,
                          start=transcript_start,
                          stop=transcript_stop,
                          featuretype=transcript_featuretype)
        transcript_id = '%s:%s'%(transcript_featuretype, transcript_name)
        features.append(feature)


        # now edit all the features...
        edited_features = []
        for feature in features:
            feature.strand = strand
            if feature.featuretype == transcript_featuretype:
                feature.add_attribute('ID',transcript_id)
                feature.add_attribute('Parent',parent_gene)
                edited_features.append(feature)
            else:
                feature_id = '%s:%s:%s:%s-%s' % (transcript_featuretype, transcript_name, feature.featuretype, feature.start, feature.stop)
                feature.add_attribute('ID',feature_id)
                feature.add_attribute('Parent',transcript_id)
                edited_features.append(feature)

        self._features.extend(features)


    def parse_reads(self, model):

        if not model.startswith('$'):
            raise ValueError, 'Expected reads model does not start with a "$"'

        # this turned out to be important!!!
        model = model.rstrip() 

        offset = 0
        chrom = self.chrom

        reads_legend = {'+':'plus',
                        '-':'minus',
                        ';':'minus_spliced',
                        '^':'plus_spliced',
                        '$':'start',
                        ' ':'spacer',
                        '\n':'line end'}


        def keyfunc(x):
            """
            Individual reads should not be piled up together but spliced reads
            should.  To detach individual reads, return the counter which will
            be different for each iteration.
            """
            self.counter += 1
            if x in ['-','+']:
                return ' %s read at %s'% (x,self.counter)
            else:
                try:
                    return reads_legend[x]
                except KeyError:
                    return 'annotation'

        self.counter = 0
        chars_to_ignore = ['annotation', 'spacer', 'start']
        reads = []
        for read_type, chars in itertools.groupby(model,keyfunc):
            chars = list(chars)

            start,stop = self.start_stop(chars)

            if read_type == 'annotation':
                annotation = ''.join(chars)

            if read_type in chars_to_ignore:
                continue

            if 'spliced' in read_type:
                if 'minus' in read_type:
                    strand = '-'
                if 'plus' in read_type:
                    strand = '+'
                read = SplicedRead(chrom=chrom,
                                   start=start,
                                   stop=stop,
                                   strand=strand,
                                   read_length=self.read_length,
                                   scalar=self.scalar)
            else:
                strand = ''.join(chars)
                read = Feature(chrom=self.chrom,
                                  start=start,
                                  stop=start+self.read_length-1,
                                  strand=strand,
                                  featuretype='read')
            reads.append(read)

        #try:
        #    print annotation
        #except:
        #    print 'no annoation in line:'
        #    print model
        #    raise
        if 'offset' in annotation:
            _, sign, bp = annotation.split('_')
            if sign == 'plus':
                sign = 1
            if sign == 'minus':
                sign = -1
            bp = int(bp)
            offset = bp*sign

        for read in reads:
            read.add_offset(offset) 
            self._reads.append(read)

    @property
    def reads(self):
        return DataList(self._reads)

    @property
    def features(self):
        """
        Returns a list of features with genes attached with their start/stop
        updated according to the presence of transcripts.
        """
        self.update_gene_features()
        features_to_return = []
        for feature in self._features:
            features_to_return.append(feature)
        for gene in self.genes.values():
            features_to_return.append(gene)
        return DataList(features_to_return)

    def update_gene_features(self):
        for feature in self._features:
            if feature.featuretype in ['mRNA', 'rRNA']:
                parent_gene = feature.attributes.Parent[0]
                self.genes.setdefault(parent_gene, None)
                if self.genes[parent_gene] is None:
                    gene = Feature(chrom=self.chrom,
                                               start=feature.start,
                                               stop=feature.stop,
                                               featuretype='gene',
                                               strand=feature.strand,
                                               value=0,
                                               attributes='ID=%(parent_gene)s;' % locals())
                    self.genes[parent_gene] = gene

                gene = self.genes[parent_gene]
                if feature.start < gene.start:
                    gene.start = feature.start
                if feature.stop > gene.stop:
                    gene.stop = feature.stop
                self.genes[parent_gene] = gene
