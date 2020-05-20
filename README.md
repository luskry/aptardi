### aptardi
***a**lternative **p**olyadenylation **t**rascriptome **a**nalysis from **R**NA sequencing and **D**NA sequencing **i**nformation*

### Description
High throughput RNA sequencing (RNA-Seq) is a powerful tool for characterizing and quantitating the expressed transcriptome. Yet constructing the transcriptome from RNA-Seq data alone is a challenging task, particularly at transcript boundaries, i.e. the polyA site. 

As a result, some have utilized the information provided by DNA sequence to more precisely identify polyA sites. However, DNA sequence information alone does not go beyond identifying the genomic sites of polyadenylation, i.e. expressed transcript structures and quantity cannot be determined for specific sampless, information that is crucial for downstream systems genomics studies on health and disease. 

To overcome these limitations, here we introduce aptardi, which combines both RNA-Seq data and DNA sequence information. Namely, aptardi takes as input a transcriptome (gtf/gff format), possibly constructed from RNA-Seq data, and combines RNA-Seq data for the sample with the genome of the sample to identify 3' ends of transcripts using machine learning. The output of aptardi is a new gtf/gff file. Note that aptardi does not evaluate intron-exon junctions but rather only re-annotates 3' ends accordingly.

### Requirements
1. Linux machine
2. Download aptardi (in aptardi folder)

		cd aptardi-1.x
		./configure --prefix=/where/to/install
		make
		make install
		
3. Download the machine learning model and scale in ml_scale folder (unless building your own model)

USAGE

	aptardi {OPTIONS}	

OPTIONS
	
	Required arguments
	
	--o <output directory>			Absolute directory path with read/write permissions 
	--f <fasta file>			Fasta file where headers are chromosomes
	--r <input gtf file or stdin>		Transcript file in gtf/gff format (or standard output from pipe) - this tool was designed to take the output of StringTie, but other formats may work
	--b <bam file>				Sorted bam file of aligned RNA-Seq reads
	
	1. Mode 1: Using pre-built model
	
		Additional required arguments
		
		--n/-n <model file>		Location of model downloaded from ml_scale folder
		--t/-t <scale file>		Location of scale downloaded from ml_scale folder 
		
	2. Mode 2 (advanced usage): Building your own model
	
		--m/-m <machine learning mode>	Enables Mode 2, building your own model, requires reliable genomic locations of polyA sites as the gold standard labels to train model
		
		Additional required arguments
		
		--e/-e <model name>		Name to save custom model in output directory
		--k/-k <scale name>		Name to save custom model's scale in output directory
		--s/-s <polyA sites file>	Tab separated file containing gold standard polyA sites for training model
		
		Additional optional arguments
		
		--c/-c <int>			Set seed for reproducibly building model
		--l/-l <int,int,int>		0-based coordinates of chromosome, strand, and site columns in polyA sites file (comma separated list with no spaces)
		
	Universal optional arguments
		--g/-g <output gtf file>	Name to save output gtf file in output directory
		-h <help>			Prints help/usage
		--version/-v <version>		Prints version
		--d/-d <debugging>		Saves intermediate files to facilitate issues
		-verbose/-vb <verbose>		Prints progress to standard output
		--i/-i <int>			Maximum length analyzed per transcript (default: 300, which is number of 100 base windows analyzed, i.e 300 = 30,000 bases long transcript) 
		--p/-p <float>			Probability threshold, predictions >= threshold are labeled transription stop site (default: 0.5, value must be constrained by (0, 1))
		--a/-a <fr or rf>		Upstream/downstream mate orientations for paired-end alignment against the forward reference strand, fr = firststrand (appropriate for Illumina paired-end library pre, rf = secondstrand (default: fr)

### Generating required input files

1. DNA sequence

	Two options:
	1. Sample specific genome
	2. [Reference genome with headers as chromosomes](https://hgdownload.soe.ucsc.edu/downloads.html)
		
	Ex. [Human:](ftp://hgdownload.soe.ucsc.edu/goldenPath/currentGenomes/Homo_sapiens/bigZips/)
	Choose: hg19.fa.gz



2. Sorted bam file
	
		Ex.1: Standalone
		aptardi --b sorted.bam --f hg38.fa --r stringtie.gtf --g aptardi.gtf --n model.hdf5 --t scale.pk --o output_dir
	


EXAMPLE

	Ex.1: Standalone
	aptardi --b sorted.bam --f hg38.fa --r stringtie.gtf --g aptardi.gtf --n model.hdf5 --t scale.pk --o output_dir
	
	Ex. 2: Pipe standard input to aptardi
	stringtie sorted.bam {OPTIONS} | aptardi --b sorted.bam --f hg38.fa --r - --g aptardi.gtf --n model.hdf5 --t scale.pk --o output_dir
	
	Ex. 3: Write aptardi's gtf to standard output
	aptardi --b sorted.bam --f hg38.fa --r stringtie.gtf --n model.hdf5 --t scale.pk --o output_dir | rsem-prepare-reference --gtf - {OPTIONS}
	
	Ex. 4: Pipe standard input to aptardi and write aptardi's gtf to standard output
	stringtie sorted.bam {OPTIONS} | aptardi --b sorted.bam --f hg38.fa --r - --n model.hdf5 --t scale.pk --o output_dir | rsem-prepare-reference --gtf - {OPTIONS}

### Output
Aptardi analyzes the input gtf file and outputs a new gtf file where transcript ends are re-annotated accordingly. The new gtf file can be used for downstream analyses (i.e. quantitation and systems studies) in the same manner as the input gtf. Note by default aptardi writes to standard output.


## References
1. Pertea, M., Pertea, G., Antonescu, C. et al. StringTie enables improved reconstruction of a transcriptome from RNA-seq reads. Nat Biotechnol 33, 290–295 (2015). https://doi.org/10.1038/nbt.3122
 
Note: If you have any question or suggestions please feel free to email: ryan.lusk@cuanschutz.edu
