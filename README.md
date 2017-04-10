# CLAMMS: a scalable algorithm for calling common and rare copy number variants from exome sequencing data

## What CLAMMS is for

As per the title, CLAMMS (Copy number estimation using Lattice-Aligned Mixture Models) is an algorithm for calling copy number variants (CNVs) from exome sequencing read depths. It has two main advantages over previous CNV callers for exome data:

1. CLAMMS is suitable for calling CNVs across the whole allele frequency spectrum, not just rare CNVs. Previous tools require that each sample be compared to a reference panel of samples that are assumed to be diploid in any given region. This assumption does not hold in copy number polymorphic regions (where non-diploid alleles are not rare), leading to improper genotypes.
1. CLAMMS can scale to datasets of tens or hundreds of thousands of samples. Apart from one short processing step (which takes ~30 seconds for 30,000 samples), each sample can be processed in parallel. Unlike previous tools, which have RAM requirements that scale linearly in the number of samples, each CLAMMS process uses a constant amount of RAM regardless of the number of samples.

Please note that CLAMMS is not intended to be used with whole-genome sequencing data or data from cancer samples.

The rest of this README will give instructions on how to use CLAMMS.

## Related Publications

* Methods paper: Packer JS, Maxwell EK, O?~@~YDushlaine C, et al. (2015) CLAMMS: a scalable algorithm for calling common and rare copy number variants from exome sequencing data. Bioinformatics 32 (1): 133-135.) [link](http://bioinformatics.oxfordjournals.org/content/32/1/133) describes the methods of CLAMMS, as well as the results of validation experiments we used to evaluate its performance in comparison to previous tools.
* 50K DiscovEHR Study CNV analysis (pre-print): Maxwell EK, Packer JS, O'Dushlaine C, McCarthy SE, Hare-Harris A, Gonzaga-Jauregui C, et al. (2017) Profiling copy number variation and disease associations from 50,726 DiscovEHR Study exomes. bioRxiv. Survey of CLAMMS CNVs from ~50k DiscovEHR study exomes with paired EHR phenotype associations. See supplemental materials for additional details on CLAMMS validation and quality-control procedures. [http://biorxiv.org/content/early/2017/03/22/119461] (http://biorxiv.org/content/early/2017/03/22/119461)

## Getting Started

First, clone the CLAMMS Github repository and compile the code:

    git clone https://github.com/rgcgithub/clamms.git
    cd clamms
    make

Set the environment variable `CLAMMS_DIR` to the appropriate path using the `export` command

    export $CLAMMS_DIR=/path/to/clamms

Now you will need to generate a file `windows.bed`. This file will list windows along the exome for which CLAMMS will estimate copy numbers, along with metadata for those windows. Most windows will simply be exons from your exome capture design, but large exons (>= 1000 bp) will be split up into equally-sized calling windows of ~500 bp.

To generate `windows.bed`, you will need four input files:

1. targets.bed — a BED file listing your exome capture regions.
1. genome.fa — an indexed FASTA file for the reference genome you are using.
1. mappability.bed — a BED file listing mappability scores across the genome. More detail on this below.
1. clamms_special_regions.bed — provided in the data/ directory with the code distribution (hg19 coordinates).

The chromosome names in the BED files and in the genome index should not have "chr" preceding the number/letter (i.e. "1" instead of "chr1"). The BED files must be sorted using either `bedtools sort` or `sort -k1,1 -k2,2n`.

The FASTA index should be generated from the raw FASTA file using [BWA](http://bio-bwa.sourceforge.net/bwa.shtml): `bwa index genome.fa`.

The mappability score for a given base is one divided by the number of locations in the genome that the k-mer starting at that base aligns to (k = the length of your reads), with up to two mismatches allowed (see [here](http://genome.ucsc.edu/cgi-bin/hgFileUi?db=hg19&g=wgEncodeMapability) for more details). You can download mappability tracks for 75-mers or 100-mers on the GRCh37 human reference genome from the link above and convert them to CLAMMS-ready BED files (requires `bigWigToWig` tool from [UCSC](http://genome.ucsc.edu/goldenpath/help/bigWig.html)):

    wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeMapability/wgEncodeCrgMapabilityAlign75mer.bigWig
    bigWigToWig wgEncodeCrgMapabilityAlign75mer.bigWig wgEncodeCrgMapabilityAlign75mer.wig
    grep -v '^#' wgEncodeCrgMapabilityAlign75mer.wig | sed 's/^chr//g' > mappability.bed

Once you have the input files ready, you can generate `windows.bed` with the following commands. This will take ~5 minutes. Note that the preprocessing script `annotate_windows.sh` requires [Bedtools](http://bedtools.readthedocs.org) to be installed and in your system PATH.

    export INSERT_SIZE=200
    chmod +x $CLAMMS_DIR/annotate_windows.sh
    $CLAMMS_DIR/annotate_windows.sh targets.bed genome.fa mappability.bed $INSERT_SIZE $CLAMMS_DIR/data/clamms_special_regions.bed >windows.bed

The `INSERT_SIZE` variable should be set to a value that is a little bit larger than the average insert size for your sequencing process (so that most reads will come from inserts of size <= this value). For example, we use `INSERT_SIZE = 200` when our mean insert size is ~150 bp. If a window is smaller than `INSERT_SIZE`, it is extended to the length of `INSERT_SIZE` for purposes of calculating it's GC content. This is because according to [Benjamini and Speed (2012)](http://www.ncbi.nlm.nih.gov/pubmed/22323520), GC coverage bias is best estimated based on the GC content of the insert, not the individual reads.

#### Troubleshooting windows.bed file generation
If you have trouble generating the windows.bed file, your input files are likely improperly formatted or have inconsistencies. A few things you should check:

* All BED files are sorted properly, using `sort -k1,1 -k2,2n`. This sorts by chromosome name (string sort) then by start position (numeric sort). You should re-sort all files in the event that your system locale settings differ from those that were used to sort the externally sourced input files.
* Chromosome naming consistency: Make sure that all chromosomes are named consistently (i.e. chromosome 1 is "1", not "chr1"). This must be the case in all input files, including the genome.fa input file.

Here is a simple test you can run on your input files to make sure they are consistent and sorted properly:

    cut -f 1 targets.bed | uniq
    cut -f 1 mappability.bed | uniq
    cut -f 1 clamms_special_regions.bed | uniq

All of these should return the same chromosome names and sort order:

    1
    10
    ...
    19
    2
    20
    21
    22
    3
    4
    ...
    9
    X
    Y

Also check the chromosome names in the genome FASTA file (sort order is not important):

    grep '^>' -m 24 genome.fa
    >1
    >2
    ...
    >22
    >X
    >Y

## Computing depths of coverage

You will need a BED file for each of your samples listng the mean depth of coverage at each of the exact windows listed in `windows.bed`. The coverage files must be named in the following format: `sample_name.coverage.bed`

The sample depth-of-coverage files can be generated from BAM files using [Samtools](http://www.htslib.org/):

    # 30 = minimum mapping quality for a read to be counted
    samtools bedcov -Q 30 windows.bed sample.bam \
    | awk '{ printf "%s\t%d\t%d\t%.6g\n", $1, $2, $3, $NF/($3-$2); }' \
    >sample.coverage.bed

They can also be generated from GATK DepthOfCoverage output files:

    $CLAMMS_DIR/gatk_coverage_to_bed sample.gatk_readDepth_1x_q30.out windows.bed >sample.coverage.bed

This step will almost certainly take longer than the CNV calling process itself. To speed it up, consider processing several samples in parallel using `xargs`:

    cat list.of.samples.txt | xargs -P $NUM_PROCESSES --max-args 1 ./compute_coverage.sh

Where `list.of.samples.txt` lists each sample name (one per line) and `compute_coverage.sh` is a shell script that takes a sample name as its argument and generates its coverage file using one of the two methods shown above.

## Normalizing the coverage data

The first step of CLAMMS is to normalize each individual sample's coverage data to control for GC-bias and the sample's overall average depth of coverage.

    ls *.coverage.bed | cut -d '.' -f 1 | while read SAMPLE
    do
        $CLAMMS_DIR/normalize_coverage $SAMPLE.coverage.bed windows.bed >$SAMPLE.norm.cov.bed
    done

This step can be parallelized using `xargs`:

    echo '$CLAMMS_DIR/normalize_coverage $1.coverage.bed windows.bed >$1.norm.cov.bed' \
    >normalize_coverage.sh && chmod +x normalize_coverage.sh
    cat list.of.samples.txt | xargs -P $NUM_PROCESSES --max-args 1 ./normalize_coverage.sh

## Training the statistical models

To call CNVs for a given sample, CLAMMS compares its coverage data to probability distributions that descibe the expected depth of coverage, conditional on copy number state, at each calling window. These distributions are fit using coverage data from a reference panel of samples that ideally have been sequenced using the same procedures. In this section, we show how to train CLAMMS models under the unrealistic assumption that there are no "batch effects" in your data. Batch effects are systematic variations in coverage due to variability in sample preparation procedures, sequencing procedures, or even input DNA quality. At the end of the tutorial, we will show how to implement our recommended procedure for batch effect correction, which involves selecting a "custom" reference panel for each sample.

A reference panel is specified using a file with two columns: 1) a path to a sample's normalized coverage file, and 2) the sample's sex (optional).

    ls *.norm.cov.bed | while read FILE;
    do
        echo -e -n "$FILE\t"
        grep "^Y" $FILE | awk '{ x += $4; n++; } END { if (x/n >= 0.1) print "M"; else print "F"; }'
    done >ref.panel.files.txt

If you have a reference panel of more than 1000 samples, you may have to increase the limit on the maximum number of open files for your system (i.e. `ulimit -n 20000`). [This link](http://stackoverflow.com/questions/11342167/how-to-increase-ulimit-on-amazon-ec2-instance) explains how to increase max open files on an Amazon EC2 instance.

The CLAMMS models are trained using the `fit_models` program.

    $CLAMMS_DIR/fit_models ref.panel.files.txt windows.bed >models.bed

`models.bed` will have the following columns:

1. chromosome
1. window start coordinate
1. window end coordinate
1. max copy number considered (-1 if window filtered, 6 for known duplication regions, and 3 otherwise; see Supplementary Materials)
1. GC fraction of the window
1. average mappability score of bases in the window
1. homozygous deletion distribution flag (model parameter; see CLAMMS paper for details)
1. lambda (model parameter)
1. mu\_dip (model parameter)
1. sigma\_dip (model parameter)
1. estimated \# samples in the reference panel that have copy number 0 at this window
1. estimated \# samples in the reference panel that have copy number 1 at this window
1. estimated \# samples in the reference panel that have copy number 2 at this window
1. estimated \# samples in the reference panel that have copy number 3 at this window
1. estimated \# samples in the reference panel that have copy number 4 at this window
1. estimated \# samples in the reference panel that have copy number 5 at this window
1. estimated \# samples in the reference panel that have copy number 6 at this window

## Making CNV Calls

Once you have a models file, you can call CNVs using the `call_cnv` program.

    $CLAMMS_DIR/call_cnv sample.norm.cov.bed models.bed --sex $SEX >sample.cnv.bed

The `--sex` argument is optional and takes the values `M` or `F`.

`sample.cnv.bed` will have the following columns.

1. chromosome
1. window start coordinate
1. window end coordinate
1. interval (chr:start-end)
1. sample name/id
1. DEL or DUP
1. most likely integer copy number
1. number of windows in the call
1. Q\_SOME: Phred-scaled quality of any CNV being in this interval.
1. Q\_EXACT: a non-Phred-scaled quality score that measures how closely the coverage profile matches the exact called CNV state and breakpoints. Will document in greater detail later. Any call with Q\_EXACT < 0 is of questionable quality.
1. Q\_LEFT\_EXTEND: Phred-scaled quality of the left breakpoint (based on the likelihood ratio of the stated breakpoint compared to extending the call by 1 window on the left)
1. LEFT\_EXTEND\_COORD: add this to the CNV start coordinate to get the start coordinate of the first window to the left of the called CNV
1. Q\_RIGHT\_EXTEND: phred-scaled quality of the right breakpoint (based on the likelihood ratio of the stated breakpoint compared to extending the call by 1 window on the right)
1. RIGHT\_EXTEND\_COORD: add this to the CNV end coordinate to get the end coordinate of the first window to the right of the called CNV
1. Q\_LEFT\_CONTRACT: phred-scaled quality of the left breakpoint (based on the likelihood ratio of the stated breakpoint compared to shrinking the call by 1 window on the left)
1. LEFT\_CONTRACT\_COORD: add this to the CNV start coordinate to get the start coordinate of the second window of the called CNV
1. Q\_RIGHT\_CONTRACT: phred-scaled quality of the right breakpoint (based on the likelihood ratio of the stated breakpoint compared to shrinking the call by 1 window on the right)
1. RIGHT\_CONTRACT\_COORD: add this to the CNV end coordinate to get the end coordinate of the second-to-last window of the called CNV

## Handling Batch Effects

Calling CNVs on large sample sets is difficult because variability in DNA quality, sample preparation procedures, and sequencing procedures results in systematic biases in coverage data ("batch effects"). Some methods attempt to correct for these biases using dimensionality reduction techniques (i.e. PCA), while others select for each sample a "custom" reference panel of samples that have coverage profiles that are highly correlated to the sample in question. Both approaches become computationally limiting, as they require each sample's coverage profile to be compared to every other sample's coverage profile, resulting in O(*n*<sup>2</sup>) computational complexity.

CLAMMS uses the "custom reference panel" approach to correct batch effects, but instead of examining samples' coverage profile directly, it examines a small number of sequencing quality control (QC) metrics. In practice, we find that samples with similar QC metrics have similar coverage profiles, so the QC metrics can be thought of as a pre-defined dimensionality reduction of the coverage data that to some extent reflects the underlying causes of coverage variance. Since the number of QC metrics we examine is small, we can use a *k*-d tree data structure to efficiently select a reference panel suitable for any given sample. The computational complexity of the CLAMMS CNV calling pipeline is formally O(*n* log *n*), but the reference panel selection procedure is very fast (~30 seconds for 30,000 samples) and all other steps are O(*n*), so in practice CLAMMS achieves linear scalability in the number of samples processed.

Below, we will show how to generate the *k*-d tree and select a reference panel for each sample in your analysis. However, please note that if you are not building CLAMMS into an automated pipeline and are performing a one-time analysis on a dataset of tens or hundreds of samples, it may be sufficient to take a more streamlined approach where a small set of precomputed models are reused within sample batches:

1. Generate a PCA plot from the coverage data from all samples in your small dataset
1. Manually assign samples to batches based on the PCA plot
1. Train a set of CLAMMS models for each batch
1. For each sample, call CNVs using the models for the batch you assigned it to

The `svd` program ([link](http://tedlab.mit.edu/~dr/SVDLIBC/)) makes it easy to compute a PCA, and supports computing only the first *n* principal components. The following code shows an example of computing the first 4 principle components from your coverage data.

    sudo apt-get install gawk
    NUM_SAMPLES=`ls *.norm.cov.bed | wc -l | awk '{print $1}'`
    NUM_WINDOWS=`ls *.norm.cov.bed | head -n 1 | xargs awk '$1 != "X" && $1 != "Y" && $NF == 0 {x++} END {print x}'`
    echo -e "$NUM_SAMPLES\t$NUM_WINDOWS" >matrix.txt
    
    ls *.norm.cov.bed | while read FILE
    do
        awk '$1 != "X" && $1 != "Y" && $NF == 0 { print $4 }' $FILE \
        | gawk -f $CLAMMS_DIR/transpose.gawk >>matrix.txt
    done
    
    svd -d 4 -o svd-output -r dt matrix.txt
    ls *.norm.cov.bed | cut -d '.' -f 1 >sample.names.txt
    tail -n +2 svd-output-Ut | tr ' ' '\t' | gawk -f $CLAMMS_DIR/transpose.gawk \
    | paste sample.names.txt - >pca.coordinates.txt

The output file `pca.coordinates.txt` can be loaded in R:

    coords <- read.table("pca.coordinates.txt", col.names=c("sample", "pc1", "pc2", "pc3", "pc4"), colClasses=c("character", rep("numeric", 4)))

We recommend the library `ggplot2` for plotting.

    library(ggplot2)
    ggplot(coords, aes(x = pc1, y = pc2)) + geom_point()

## Selecting reference panels using the *k*-d tree

To identify a custom reference panel for every sample efficiently, CLAMMS collects seven QC metrics for each sample and performs a fast *k*-nearest neighbors search algorithm (*k=100*) implemented using a *k*-d tree data structure. This is performed with the R package `FNN` ([link](http://cran.r-project.org/web/packages/FNN/index.html)). While different QC metrics and values of *k* can be used, we found the Picard metrics GCDROPOUT, ATDROPOUT, MEANINSERTSIZE, ONBAITVSSELECTED, PCTPFUQREADS, PCTTARGETBASES10X, and PCTTARGETBASES50X to work well in practice. Also note that QC metrics should be normalized to similar scales such that distances between QC metrics are equally weighted.

We have provided an example data set and steps to identify the 20-nearest neighbors in R.

    $CLAMMS_DIR/data/example_qcs.Rdata

The following code assumes that a data frame has been constructed with sample IDs in the first column and raw QC metrics in the subsequent columns, with one sample per row. It also requires the `FNN` package described above.

    # This code requires the FNN (Fast Nearest Neighbors) R package (http://cran.r-project.org/package=FNN)
    require(FNN)
    
    # Load the example data set into data frame 'example.qcs'
    load("example_qcs.Rdata")
    
    # Create a scaled copy of the data frame
    example.qcs.scaled <- example.qcs
    for (i in 2:ncol(example.qcs.scaled)) {
        mini <- min(example.qcs.scaled[,i])
        maxi <- max(example.qcs.scaled[,i])
        example.qcs.scaled[,i] <- apply(example.qcs.scaled, 1, function(row) { 
			row[[i]] <- (as.numeric(row[[i]]) - mini) / (maxi - mini)
		} )
    }
    
    # Get k-nearest neighbors for each sample
    k.param <- 20
    knns <- get.knn(example.qcs.scaled[,c(seq(2,ncol(example.qcs.scaled)))],k=k.param,algorithm="kd_tree")
    
    # Generate a single file for each sample listing its k-nearest neighbor sample IDs
    for (i in 1:nrow(example.qcs.scaled)) {
        fname <- paste(example.qcs.scaled$SAMPLE[i], ".", k.param, "nns.txt", sep="")
        nn.sampleids <- example.qcs.scaled$SAMPLE[ knns$nn.index[i,] ]
        write.table(nn.sampleids, fname, quote=F, row.names=F, col.names=F)
    }

A single file will be generated for each sample with a list of its *k*-nearest neighbor sample IDs (with filenames `<sampleID>.<k>nn.txt`). Mapping this list to a list of normalized coverage BED files will produce the input required for the fit_models command (see sections "Training the Statistical Models" and "Calling CNVs using the selected reference panels").

    sed 's/$/.norm.cov.bed/' <sampleID>.<k>nn.txt > <sampleID>.ref.panel.files.txt
    
To see how well the *k*-nearest neighbors fit for each sample, we compute the distance of each sample to the mean of the cluster corresponding to its *k*-nearest neighbors, then plot the cumulative distribution of this metric over all samples. This should help to identify if there are outlier samples that do not have a good reference panel.

    # To check how well each sample's kNNs fit, compute the distance to its kNN cluster mean
    example.qcs.scaled$DistanceToClusterMean <- sapply(1:nrow(example.qcs.scaled),  function(x) {
        this.knns <- knns$nn.index[x,];
        center <- colMeans(example.qcs.scaled[this.knns, 2:ncol(example.qcs.scaled)]);
        return(as.numeric(dist(rbind(as.numeric(example.qcs.scaled[x, 2:ncol(example.qcs.scaled)]), as.numeric(center)))))
    })

    # Plot distance distribution
    plot(ecdf(example.qcs.scaled$DistanceToClusterMean))

## Calling CNVs using the selected reference panels

    # make a master list of the name, norm.cov.bed filepath, and sex of every sample

    ls *.norm.cov.bed | while read FILE;
    do
        SAMPLE=`echo "$FILE" | cut -d '.' -f 1`
        echo -e -n "$SAMPLE\t$FILE\t"
        grep "^Y" $FILE | awk '{ x += $4; n++; } END { if (x/n >= 0.1) print "M"; else print "F"; }'
    done >sample.file.sex.txt

    # call CNVs (should take ~1 minute/sample)
    # you can run this parallel if you put the inner part of the loop
    # in a shell script and run it using xargs as shown in previous examples

    ls *.norm.cov.bed | cut -d '.' -f 1 | while read SAMPLE;
    do
        SEX=`echo "$SAMPLE" | join - sample.file.sex.txt | tr ' ' '\t' | cut -f 3`
        join $SAMPLE.100nns.txt sample.file.sex.txt | tr ' ' '\t' | cut -f 2- >$SAMPLE.ref.panel.txt
        $CLAMMS_DIR/fit_models $SAMPLE.ref.panel.txt windows.bed >$SAMPLE.models.bed
        $CLAMMS_DIR/call_cnv $SAMPLE.norm.cov.bed $SAMPLE.models.bed --sex $SEX >$SAMPLE.cnv.bed
    done

To alter the sensitivity/specificity profile of CLAMMS CNV calls, modify the `call_cnv --cnv_rate` parameter. The default (3.0e-8) is tuned for specificity, but can be increased to improve sensitivity (particularly for small CNVs) at the cost of an increased FDR.

## Visualizing CNVs for a sample

First, install the R packages `dplyr` and `ggplot2`.
Then, for a sample `$SAMPLE`, run the following command:

    $CLAMMS_DIR/plot_cnv.sh $SAMPLE.cnv.txt $SAMPLE.normalized.coverage.bed $SAMPLE.models.bed

This script will create a directory `clamms_cnv_plots/$SAMPLE/` with PNG images visualizing each CNV called for the sample. An example output image is shown in Section 8 of the Supplementary Materials.
