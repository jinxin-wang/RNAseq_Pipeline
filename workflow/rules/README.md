# 1. [star.smk](https://github.com/jinxin-wang/RNAseq_Pipeline/blob/main/workflow/rules/README.md#starsmk-rule-star-align-rna-seq-samples-using-star-review-some-options)

## [star.smk] rule star: Align RNA-seq samples using [STAR](https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf), review some options 

#### --sjdbOverhang 99

**--sjdbOverhang** specifies the length of the genomic sequence around the annotated junction
to be used in constructing the splice junctions database. Ideally, this length should be equal
to the ReadLength-1, where ReadLength is the length of the reads. For instance, for Illumina
2x100b paired-end reads, the ideal value is 100-1=99. In case of reads of varying length, the
ideal value is max(ReadLength)-1. In most cases, the default value of 100 will work as
well as the ideal value.

#### --sjdbGTFfile 

[gtf file format](http://genome.ucsc.edu/FAQ/FAQformat#format4)

**--sjdbGTFfile** specifies the path to the file with annotated transcripts in the standard GTF
format. STAR will extract splice junctions from this file and use them to greatly improve
accuracy of the mapping. While this is optional, and STAR can be run without annotations,
using annotations is highly recommended whenever they are available. Starting from 2.4.1a,
the annotations can also be included on the fly at the mapping step.

#### --twopassMode Basic

 Use **--twopassMode** Basic option to run STAR 2-pass mapping for each sample separately. Annotated 
 junctions will be included in both the 1st and 2nd passes. STAR will perform the 1st pass mapping,
then it will automatically extract junctions, insert them into the genome index, and, finally, re-map
all reads in the 2nd mapping pass. This option can be used with annotations, which can be included
either at the run-time (see #1), or at the genome generation step.

--twopass1readsN defines the number of reads to be mapped in the 1st pass. The default and
most sensitive approach is to set it to -1 (or make it bigger than the number of reads in the sample) -
in which case all reads in the input read file(s) are used in the 1st pass. While it can reduce mapping
time by ∼ 40%, it is not recommended to use a small portion of the reads in the 1st step, since it
will significantly reduce sensitivity for the low expressed novel junctions. The idea to use a portion
of the reads in the 1st pass was inspired by Kim, Langmead and Salzberg in Nature Methods 12,
357–360 (2015)

#### --outSAMtype BAM SortedByCoordinate

output sorted by coordinate Aligned.sortedByCoord.out.bam file, similar to samtools sort command. If this option causes problems, it is recommended to reduce --outBAMsortingThreadN from the default 6 to lower values (as low as 1).

#### --outSAMunmapped Within

Unmapped reads can be output into the SAM/BAM Aligned.* file(s) with --outSAMunmapped
Within option. --outSAMunmapped Within KeepPairs will (redundantly) record unmapped mate
for each alignment, and, in case of unsorted output, keep it adjacent to its mapped mate (this only
affects multi-mapping reads). uT SAM tag indicates reason for not mapping:

0 : no acceptable seed/windows, ”Unmapped other” in the Log.final.out

1 : best alignment shorter than min allowed mapped length, ”Unmapped: too short” in the Log.final.out

2 : best alignment has more mismatches than max allowed number of mismatches, ”Unmapped: too many mismatches” in the Log.final.out

3 : read maps to more loci than the max number of multimappng loci, ”Multimapping: mapped to too many loci” in the Log.final.out

4 : unmapped mate of a mapped paired-end read

#### --outSJfilterReads Unique

default: All

string: which reads to consider for collapsed splice junctions output

- All

all reads, unique- and multi-mappers

- Unique

uniquely mapping reads only

#### --quantMode GeneCounts

With **--quantMode** TranscriptomeSAM option STAR will output alignments translated into transcript
coordinates in the Aligned.toTranscriptome.out.bam file (in addition to alignments in genomic 
coordinates in Aligned.*.sam/bam files). These transcriptomic alignments can be used with
various transcript quantification software that require reads to be mapped to transcriptome, such as
RSEM or eXpress. For example, RSEM command line would look as follows:

rsem-calculate-expression ... --bam Aligned.toTranscriptome.out.bam

/path/to/RSEM/reference RSEM

Note, that STAR first aligns reads to entire genome, and only then searches for concordance
between alignments and transcripts.This approach offers certain advantages compared to the alignment to 
transcriptome only, by not forcing the alignments to annotated transcripts. Note that
--outFilterMultimapNmax filter only applies to genomic alignments. If an alignment passes this
filter, it is converted to all possible transcriptomic alignments and all of them are output.

By default, the output satisfies RSEM requirements: soft-clipping or indels are not allowed. Use
--quantTranscriptomeBan Singleend to allow insertions, deletions ans soft-clips in the transcriptomic 
alignments, which can be used by some expression quantification software (e.g. eXpress).


#### --outSAMstrandField intronMotif

For unstranded RNA-seq data, Cufflinks/Cuffdiff require spliced alignments with XS strand attribute,
which STAR will generate with --outSAMstrandField intronMotif option. As required, the XS
strand attribute will be generated for all alignments that contain splice junctions. The spliced
alignments that have undefined strand (i.e. containing only non-canonical unannotated junctions)
will be suppressed.

If you have stranded RNA-seq data, you do not need to use any specific STAR options. Instead,
you need to run Cufflinks with the library option --library-type options. For example, cufflinks
... --library-type fr-firststrand should be used for the “standard” dUTP protocol, including Illumina’s 
stranded Tru-Seq. This option has to be used only for Cufflinks runs and not for STAR
runs.

In addition, it is recommended to remove the non-canonical junctions for Cufflinks runs using
--outFilterIntronMotifs RemoveNoncanonical.


#### --outSAMmultNmax 5

The **--outSAMmultNmax** parameter limits the number of output alignments (SAM lines) for
multimappers. For instance, --outSAMmultNmax 1 will output exactly one SAM line for each 11
mapped read. Note that NH:i: tag in STAR will still report the actual number of loci that
the reads map to, while the the number of reported alignments for a read in the SAM file is
min(NH,--outSAMmultNMax). If --outSAMmultNmax is equal to -1, all the alignments are output 
according to the order specified in --outMultimapperOrder option. If --outSAMmultNmax
is not equal to -1, than top-scoring alignments will always be output first, even for the default
--outMultimapperOrder Old 2.4 option.

## [chk_mapping_metrics.smk] rule samtools_flagstat: mapping metrics with [samtools flagstat](http://www.htslib.org/doc/samtools-flagstat.html)

## [mapping_coverage.smk] rule rseqc_geneBody_coverage: mapping coverage with [rseqc](https://rseqc.sourceforge.net/#genebody-coverage-py)

**-l MIN_MRNA_LENGTH, --minimum_length=MIN_MRNA_LENGTH**

Minimum mRNA length (bp). mRNA smaller than “min_mRNA_length” will be skipped. **default=100, set to 500 in the rule.**

## [transcript_integrity_number.smk] rseqc [tin](https://rseqc.sourceforge.net/#tin-py) and [read_duplication](https://rseqc.sourceforge.net/#read-duplication-py)

rule rseqc_tin:

**-c MINIMUM_COVERAGE, --minCov=MINIMUM_COVERAGE**

Minimum number of read mapped to a transcript. **default=10， set to 30 in the rule.**

## [quantify_reads.smk] quantify reads per annotation

**-r \<order\>, --order=\<order\>**

For paired-end data, the alignment have to be sorted either by read name or by alignment position. If your data is not sorted, use the samtools sort function of samtools to sort it. Use this option, with name or pos for <order> to indicate how the input data has been sorted. The default is name.

If **name** is indicated, htseq-count expects all the alignments for the reads of a given read pair to appear in adjacent records in the input data. For **pos**, this is not expected; rather, read alignments whose mate alignment have not yet been seen are kept in a buffer in memory until the mate is found. While, strictly speaking, the latter will also work with unsorted data, sorting ensures that most alignment mates appear close to each other in the data and hence the buffer is much less likely to overflow.

**-s \<yes/no/reverse\>, --stranded=\<yes/no/reverse\>**

Whether the data is from a strand-specific assay (default: yes)

For **stranded=no**, a read is considered overlapping with a feature regardless of whether it is mapped to the same or the opposite strand as the feature. For **stranded=yes** and single-end reads, the read has to be mapped to the same strand as the feature. For paired-end reads, the first read has to be on the same strand and the second read on the opposite strand. For **stranded=reverse**, these rules are reversed.
 
**-i \<id attribute\>, --idattr=\<id attribute\>**
 
GTF attribute to be used as feature ID. Several [GTF](http://genome.ucsc.edu/FAQ/FAQformat#format4) lines with the same feature ID will be considered as 
parts of the same feature. The feature ID is used to identity the counts in the output table. 
 
The default, suitable for RNA-Seq analysis using an Ensembl GTF file, is **gene_id**

**transcript_id**

In principle, you could instruct htseq-count to count for each of a gene’s transcript individually, by 
specifying --idattr transcript_id. However, all reads mapping to exons shared by several transcripts will 
then be considered ambiguous. (See second question.) Counting them for each transcript that contains the 
exons would be possible but makes little sense for typical use cases. (See first question.) If you want to 
perform differential expression analysis on the level of individual transcripts, maybe have a look at the 
paper on DEXSeq for a discussion on why we prefer performing such analyses on the level of exons instead.

**gene_name**
 
**transcript_name**

## [kallisto.smk] quantifying abundances of transcripts with [kallisto](http://pachterlab.github.io/kallisto/manual.html)
