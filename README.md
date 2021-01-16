## aptardi
***a**lternative **p**olyadenylation **t**rascriptome **a**nalysis from **R**NA sequencing and **D**NA sequencing **i**nformation*

## Description
High throughput RNA sequencing (RNA-Seq) is a powerful tool for characterizing and quantitating the expressed transcriptome. Yet constructing the transcriptome from RNA-Seq data alone is a challenging task, particularly at transcript boundaries, i.e. the polyA site. 

As a result, some have utilized the information provided by DNA sequence to more precisely identify polyA sites. However, DNA sequence information alone does not consider expression of specific samples, information that is crucial for downstream systems genomics studies on health and disease. 

To overcome these limitations, here we introduce aptardi, which combines both RNA-Seq data and DNA sequence information. Namely, aptardi takes as input a transcriptome (gtf/gff format), possibly constructed from RNA-Seq data, and combines RNA-Seq data for the sample with the genome (DNA sequence) of the sample to identify 3' ends of transcripts using machine learning. The output of aptardi is a new gtf/gff file that incorporates novel transcripts identified by aptardi. Note that aptardi does not evaluate intron junctions but rather only annotates 3' ends accordingly.

## Requirements
1. Linux machine

2. [SAMtools (v.1.9)](http://www.htslib.org/download/)

3. [BEDtools (v.2.29.2)](https://bedtools.readthedocs.io/en/latest/content/installation.html)

**Note: The bioconda samtools and bedtools Python3 libraries may prevent aptardi from working, please uninstall these Python3 packages if present**

4. Python3 (v.3.7.7) (**PYTHON VERSIONS 3.8x AND ABOVE DO NOT SUPPORT THE TENSORFLOW VERSION REQUIRED**)

5. The following Python3 libraries (listed version), in addition to the standard Python3 libraries:

	tensorflow==2.0.0 (**NEWER VERSIONS ARE CURRENTLY NOT COMPATABLE WITH APTARDI - DO NOT USE**)<br/> 
	numpy==1.18.1\
	pandas==1.0.3\
	scikit_learn==0.22.2.post1 (**v.0.23.2 IS CURRENTLY NOT COMPATABLE WITH APTARDI - DO NOT USE**)

These libraries are also available in the python_dependencies folder as requirements.txt

**Note: The Python3 version called with #!/usr/bin/env python3 must have all the Python3 dependencies installed**
	
6. aptardi (installation time <5 min)

	Download the latest release of aptardi (aptardi.zip under Releases)\
	Put in desired folder (we recommend /usr/local/bin)\
	Add executable permission\
	Access PATH file\
	Add aptardi to PATH (if you put aptardi in /usr/local/bin/ the PATH will almost assuredly be in your profile)

		sudo chmod +x aptardi
		vi .bash_profile
		export PATH="/usr/local/bin/:$PATH"
		
To check if aptardi was properly installed in your PATH, enter the following in your terminal

	which aptardi
	
This command should print to your screen aptardi's PATH
	
7. The machine learning model (model.hdf5) and scale (scale.pk) in ml_scale folder (unless building your own model)

### Usage

	aptardi {OPTIONS}	

### Options
	
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
		--i/-i <int>			Maximum length analyzed per transcript (default: 300, which is number of bins analyzed, i.e 300 and bin size of 100 = 300 * 100 = 30,000 bases long transcript) 
		--p/-p <float>			Probability threshold, predictions >= threshold are labeled transription stop site (default: 0.5, value must be constrained by (0, 1))
		--a/-a <fr or rf>		Upstream/downstream mate orientations for paired-end alignment against the forward reference strand, fr = firststrand (appropriate for Illumina paired-end library pre, rf = secondstrand (default: fr)
		--w/-w <int>			Size of bin to make predictions on (default: 100, choices: 25-200 in 25 base increments, we DO NOT recommend altering this argument using the pre-built model since this model was built using 100 base bins)

### Generating required input files

1. DNA sequence

	Two options:
	1. [Reference genome with headers as chromosomes](https://hgdownload.soe.ucsc.edu/downloads.html)
	2. Sample specific genome (user generated)
	
2. Sorted bam file
	Note: Prior to aligning reads to the genome, we recommend assessing read quality with a tool such as [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) and performing any processing steps based on the results (e.g. trim adapters with [Cutadapt](https://cutadapt.readthedocs.io/en/stable/installation.html)) 
	
		Ex. Using HISAT2 with paired end, stranded reads generated from Illumina's protocol (i.e. firststrand)
		(hisat2 -q -p 5 --reorder -t --rna-strandness RF --dta -x <hisat2_index> -1 <myfq1_1.fq,myfq2_1.fq,etc> -2 <myfq1_2.fq,myfq2_2.fq,etc> | samtools view -F 0x4 -bS - | samtools sort - -o sorted.bam 2> sum_sorted_bam.txt

3. Input gtf file

	Two options:
	1. Reconstruction gtf file that considers expression (e.g. StringTie) - see example below
	2. [Reference annotation that does not consider expression (i.e. Ensembl)](https://uswest.ensembl.org/info/data/ftp/index.html)
		
		```Ex. Using StringTie with guide (i.e. Ensembl reference annotation) and sorted bam file generated above
		stringtie sorted.bam --rf -o stringtie.gtf -G <guide_file>
### Examples

Demo files (in demo folder, these example files contain data only for chromosome 1):
1. sorted.bam
2. hg38.fa.zip (please unzip)
3. stringtie.gtf 
4. polya_sites.bed

Note: 
1. Any command provided to an argument that starts with name (i.e. name_aptardi.gtf) should be given simply a name and not an absolute path; the file will be saved to the output directory path with the name.
2. Any command provided to an argument that does not start with name requires the path; the path can be relative to the current working directory, but when in doubt, specify the absolute path (i.e. /Users/foo_name/model.hdf5).

[StringTie](https://ccb.jhu.edu/software/stringtie/#install) and [RSEM](https://deweylab.github.io/RSEM/) are required to utilize the piping examples.

	Ex. 1: Standalone using pre-built model (~25 min runtime)
	aptardi --b sorted.bam --f hg38.fa --r stringtie.gtf --g name_aptardi.gtf --n model.hdf5 --t scale.pk --o output_dir
	
	Ex. 2: Standalone building own model (~1 hr 15 min runtime)
	aptardi --b sorted.bam --f hg38.fa --r stringtie.gtf --g name_aptardi.gtf --m --e name_model.hdf5 --k name_scale.pk --s polya_sites.bed --l 0,4,1 --o output_dir
	
	Ex. 3: Pipe upstream standard input from StringTie to aptardi
	stringtie sorted.bam {OPTIONS} | aptardi --b sorted.bam --f hg38.fa --r - --g name_aptardi.gtf --n model.hdf5 --t scale.pk --o output_dir
	
	Ex. 4: Pipe aptardi's standard output gtf to downstream RSEM
	aptardi --b sorted.bam --f hg38.fa --r stringtie.gtf --n model.hdf5 --t scale.pk --o output_dir | rsem-prepare-reference --gtf - {OPTIONS}
	
	Ex. 5: Pipe StringTie standard input to aptardi and pipe aptardi's standard output gtf to downstream RSEM
	stringtie sorted.bam {OPTIONS} | aptardi --b sorted.bam --f hg38.fa --r - --n model.hdf5 --t scale.pk --o output_dir | rsem-prepare-reference --gtf - {OPTIONS}

## Output
Aptardi analyzes the input gtf file and outputs a new gtf file with novel aptardi transcripts added. The new gtf file can be used for downstream analyses (i.e. quantitation and systems studies) in the same manner as the input gtf file. Note by default aptardi writes to standard output.


## References
1. Pertea, M., Pertea, G., Antonescu, C. et al. StringTie enables improved reconstruction of a transcriptome from RNA-seq reads. Nat Biotechnol 33, 290–295 (2015). https://doi.org/10.1038/nbt.3122
2. Kent WJ, Sugnet CW, Furey TS, et al. The human genome browser at UCSC. Genome Res. 2002;12(6):996‐1006. doi:10.1101/gr.229102
3. Kim D, Langmead B, Salzberg SL. HISAT: a fast spliced aligner with low memory requirements. Nat Methods. 2015;12(4):357‐360. doi:10.1038/nmeth.3317
4. Andrew D Yates, Premanand Achuthan, Wasiu Akanni, James Allen, Jamie Allen, Jorge Alvarez-Jarreta, M Ridwan Amode, Irina M Armean, Andrey G Azov, Ruth Bennett, Jyothish Bhai, Konstantinos Billis, Sanjay Boddu, José Carlos Marugán, Carla Cummins, Claire Davidson, Kamalkumar Dodiya, Reham Fatima, Astrid Gall, Carlos Garcia Giron, Laurent Gil, Tiago Grego, Leanne Haggerty, Erin Haskell, Thibaut Hourlier, Osagie G Izuogu, Sophie H Janacek, Thomas Juettemann, Mike Kay, Ilias Lavidas, Tuan Le, Diana Lemos, Jose Gonzalez Martinez, Thomas Maurel, Mark McDowall, Aoife McMahon, Shamika Mohanan, Benjamin Moore, Michael Nuhn, Denye N Oheh, Anne Parker, Andrew Parton, Mateus Patricio, Manoj Pandian Sakthivel, Ahamed Imran Abdul Salam, Bianca M Schmitt, Helen Schuilenburg, Dan Sheppard, Mira Sycheva, Marek Szuba, Kieron Taylor, Anja Thormann, Glen Threadgold, Alessandro Vullo, Brandon Walts, Andrea Winterbottom, Amonida Zadissa, Marc Chakiachvili, Bethany Flint, Adam Frankish, Sarah E Hunt, Garth IIsley, Myrto Kostadima, Nick Langridge, Jane E Loveland, Fergal J Martin, Joannella Morales, Jonathan M Mudge, Matthieu Muffato, Emily Perry, Magali Ruffier, Stephen J Trevanion, Fiona Cunningham, Kevin L Howe, Daniel R Zerbino, Paul Flicek, Ensembl 2020, Nucleic Acids Research, Volume 48, Issue D1, 08 January 2020, Pages D682–D688, https://doi.org/10.1093/nar/gkz966
5. Andrews, S. FASTQC. A quality control tool for high throughput sequence data (2010).
6. Martin, M. Cutadapt removes adapter sequences from high-throughput sequencing reads. EMBnet.journal 17, 10, doi:10.14806/ej.17.1.200 (2011).
 
If you have any questions or suggestions please feel free to email: ryan.lusk@cuanschutz.edu
