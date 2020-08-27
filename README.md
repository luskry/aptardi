### aptardi
***a**lternative **p**olyadenylation **t**rascriptome **a**nalysis from **R**NA sequencing and **D**NA sequencing **i**nformation*

### Description
High throughput RNA sequencing (RNA-Seq) is a powerful tool for characterizing and quantitating the expressed transcriptome. Yet constructing the transcriptome from RNA-Seq data alone is a challenging task, particularly at transcript boundaries, i.e. the polyA site. 

As a result, some have utilized the information provided by DNA sequence to more precisely identify polyA sites. However, DNA sequence information alone does not consider expression of specific samples, information that is crucial for downstream systems genomics studies on health and disease. 

To overcome these limitations, here we introduce aptardi, which combines both RNA-Seq data and DNA sequence information. Namely, aptardi takes as input a transcriptome (gtf/gff format), possibly constructed from RNA-Seq data, and combines RNA-Seq data for the sample with the genome of the sample to identify 3' ends of transcripts using machine learning. The output of aptardi is a new gtf/gff file. Note that aptardi does not evaluate intron-exon junctions but rather only re-annotates 3' ends accordingly.

### Requirements
1. Linux machine

2. [SAMtools (v.1.9 or newer)](http://www.htslib.org/download/)

3. [BEDtools (v.2.29.2 or newer)](https://bedtools.readthedocs.io/en/latest/content/installation.html)

**NOTE: The bioconda samtools and bedtools Python3 libraries will may prevent aptardi from working, please uninstall these Python3 packages if present**

4. Python3 (v.3.7.7 or newer)

5. Besides the standard Python3 libraries, the following libraries (listed version or newer):

tensorflow==2.0.1\
numpy==1.18.1\
pandas==1.0.3\
scikit_learn==0.23.2

These libraries are also available in the python_dependencies folder as requirements.txt

**NOTE: The Python3 version called with #!/usr/bin/env python3 must have all the Python3 dependencies installed**
	
6. aptardi 

	Download aptardi.zip\
	Put in desired folder (we recommend /usr/local/bin)\
	Add executable permission\
	Access PATH file\
	Add aptardi to PATH

		sudo chmod +x aptardi
		vi .bash_profile
		export PATH="/usr/local/bin/aptardi:$PATH"
	
7. Download the machine learning model (model.hdf5) and scale (scale.pk) in ml_scale folder (unless building your own model)

USAGE

	aptardi {OPTIONS}	

OPTIONS
	
	Required arguments
	
	--o <output directory>			Absolute directory path with read/write permissions 
	--f <fasta file>			Fasta file where headers are chromosomes
	--r <input gtf file or stdin>		Transcript file in gtf/gff format (or standard output from pipe) - this tool was designed to take the output of StringTie, but other formats may work
	--b <bam file>				Sorted bam file of aligned RNA-Seq reads
	
	1. Mode 1: Using pre-built model (canonical usage)
	
		Additional required arguments
		
		--n/-n <model file>		Location of model downloaded from ml_scale folder
		--t/-t <scale file>		Location of scale downloaded from ml_scale folder 
		
	2. Mode 2: Building your own model (non-canonical usage)
	
		--m/-m <machine learning mode>	Enables Mode 2, building your own model, requires reliable genomic locations of polyA sites as the gold standard labels to train model
		
		Additional required arguments
		
		--e/-e <model name>		Name to save custom model in output directory
		--k/-k <scale name>		Name to save custom model's scale in output directory
		--s/-s <polyA sites file>	Tab separated file containing gold standard polyA sites for training model
		--l/-l <int,int,int>		0-based coordinates of chromosome, strand, and site columns in polyA sites file (comma separated list with no spaces)
		
		Additional optional arguments
		
		--c/-c <int>			Set seed for reproducibly building model
		
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
	
2. Sorted bam file
	
		Ex. Using HISAT2 with paired end, stranded reads generated from Illumina's protocol (e.g. firststrand)
		(hisat2 -q -p 5 --reorder -t --rna-strandness RF --dta -x <hisat2_index> -1 <myfq1_1.fq,myfq2_1.fq,etc> -2 <myfq1_2.fq,myfq2_2.fq,etc> | samtools view -F 0x4 -bS - | samtools sort - -o sorted.bam 2> sum_sorted_bam.txt
	
3. Reconstruction file [(for reference annotation files (GTF format) click on link, for guide see example below)](https://uswest.ensembl.org/info/data/ftp/index.html)
		
		Ex. Using StringTie with guide and sorted bam file generated above
		stringtie sorted.bam --rf -o stringtie.gtf -G <guide_file>

EXAMPLES 

Demo files (in demo folder, these example files contain data only for chromosome 1):
1. sorted.bam
2. hg38.fa 
3. stringtie.gtf 
4. polya_sites.bed

Any name starting with foo is a file generated by the program and saved. Make sure to specify absolute paths to save these files (i.e. /Users/foo_name/foo_model.hdf5)

	Ex. 1: Standalone using pre-built model
	aptardi --b sorted.bam --f hg38.fa --r stringtie.gtf --g foo_aptardi.gtf --n model.hdf5 --t scale.pk --o foo_output_dir
	
	Ex. 2: Standalone building own model
	aptardi --b sorted.bam --f hg38.fa --r stringtie.gtf --g foo_aptardi.gtf --m --e foo_model.hdf5 --k foo_scale.pk --s polya_sites.bed --o foo_output_dir
	
	Ex. 3: Pipe standard input to aptardi
	stringtie sorted.bam {OPTIONS} | aptardi --b sorted.bam --f hg38.fa --r - --g foo_aptardi.gtf --n model.hdf5 --t scale.pk --o foo_output_dir
	
	Ex. 4: Write aptardi's gtf to standard output
	aptardi --b sorted.bam --f hg38.fa --r stringtie.gtf --n model.hdf5 --t scale.pk --o foo_output_dir | rsem-prepare-reference --gtf - {OPTIONS}
	
	Ex. 5: Pipe standard input to aptardi and write aptardi's gtf to standard output
	stringtie sorted.bam {OPTIONS} | aptardi --b sorted.bam --f hg38.fa --r - --n model.hdf5 --t scale.pk --o foo_output_dir | rsem-prepare-reference --gtf - {OPTIONS}

### Output
Aptardi analyzes the input gtf file and outputs a new gtf file where transcript ends are re-annotated accordingly. The new gtf file can be used for downstream analyses (i.e. quantitation and systems studies) in the same manner as the input gtf. Note by default aptardi writes to standard output.


## References
1. Pertea, M., Pertea, G., Antonescu, C. et al. StringTie enables improved reconstruction of a transcriptome from RNA-seq reads. Nat Biotechnol 33, 290–295 (2015). https://doi.org/10.1038/nbt.3122
2. Kent WJ, Sugnet CW, Furey TS, et al. The human genome browser at UCSC. Genome Res. 2002;12(6):996‐1006. doi:10.1101/gr.229102
3. Kim D, Langmead B, Salzberg SL. HISAT: a fast spliced aligner with low memory requirements. Nat Methods. 2015;12(4):357‐360. doi:10.1038/nmeth.3317
4. Andrew D Yates, Premanand Achuthan, Wasiu Akanni, James Allen, Jamie Allen, Jorge Alvarez-Jarreta, M Ridwan Amode, Irina M Armean, Andrey G Azov, Ruth Bennett, Jyothish Bhai, Konstantinos Billis, Sanjay Boddu, José Carlos Marugán, Carla Cummins, Claire Davidson, Kamalkumar Dodiya, Reham Fatima, Astrid Gall, Carlos Garcia Giron, Laurent Gil, Tiago Grego, Leanne Haggerty, Erin Haskell, Thibaut Hourlier, Osagie G Izuogu, Sophie H Janacek, Thomas Juettemann, Mike Kay, Ilias Lavidas, Tuan Le, Diana Lemos, Jose Gonzalez Martinez, Thomas Maurel, Mark McDowall, Aoife McMahon, Shamika Mohanan, Benjamin Moore, Michael Nuhn, Denye N Oheh, Anne Parker, Andrew Parton, Mateus Patricio, Manoj Pandian Sakthivel, Ahamed Imran Abdul Salam, Bianca M Schmitt, Helen Schuilenburg, Dan Sheppard, Mira Sycheva, Marek Szuba, Kieron Taylor, Anja Thormann, Glen Threadgold, Alessandro Vullo, Brandon Walts, Andrea Winterbottom, Amonida Zadissa, Marc Chakiachvili, Bethany Flint, Adam Frankish, Sarah E Hunt, Garth IIsley, Myrto Kostadima, Nick Langridge, Jane E Loveland, Fergal J Martin, Joannella Morales, Jonathan M Mudge, Matthieu Muffato, Emily Perry, Magali Ruffier, Stephen J Trevanion, Fiona Cunningham, Kevin L Howe, Daniel R Zerbino, Paul Flicek, Ensembl 2020, Nucleic Acids Research, Volume 48, Issue D1, 08 January 2020, Pages D682–D688, https://doi.org/10.1093/nar/gkz966
 
Note: If you have any question or suggestions please feel free to email: ryan.lusk@cuanschutz.edu
