### aptardi
***a**lternative **p**olyadenylation **t**rascriptome **a**nalysis from **R**NA sequencing and **D**NA sequencing **i**nformation*

### Description
High throughput RNA sequencing (RNA-Seq) is a powerful tool for characterizing and quantitating the expressed transcriptome. Yet constructing the transcriptome from RNA-Seq data alone is a challenging task, particularly at transcript boundaries, i.e. the polyA site. 

As a result, some have utilized the information provided by DNA sequence to more precisely identify polyA sites. However, DNA sequence information alone does not go beyond identifying the genomic sites of polyadenylation, i.e. expressed transcript structures and quantity cannot be determined for specific sampless, information that is crucial for downstream systems genomics studies on health and disease. 

To overcome these limitations, here we introduce aptardi, which combines both RNA-Seq data and DNA sequence information. Namely, aptardi takes as input a transcriptome (gtf/gff format), possibly constructed from RNA-Seq data, and combines RNA-Seq data for the sample with the genome of the sample to identify 3' ends of transcripts using machine learning. The output of aptardi is a new gtf/gff file. Note that aptardi does not evaluate intron-exon junctions but rather only re-annotates 3' ends accordingly.

### Requirements
1. Linux machine
2. Download aptardi
Please run the following command to install aptardi:
		
		./configure
		make
		make install
3. Download the machine learning model and scale in ml_scale folder (unless building your own model)

USAGE

	aptardi {OPTIONS}	

OPTIONS
	
	Required arguments
	--o <output directory>	Absolute directory path to save new gtf file 

	--f <fasta file>	Fasta file where headers are chromosomes
	--r <transcriptome file>	Transcript file in gtf/gff format - this tool was designed to take the output of StringTie, but other formats may work
	--b <bam file>	Sorted bam file of aligned RNA-Seq reads
	1. Mode 1: Using pre-built model
		Additional required arguments
		--n/-n <model file>	Location of model downloaded from ml_scale folder
		--t/-t <scale file>	Location of sclae downloaded from ml_scale folder 
	2. Mode 2: Building your own model
		--m/-m <machine learning mode>	Enables Mode 2, building your own model, requires reliable genomic locations of polyA sites as the gold standard labels to train model
		Additional required arguments
		--e/-e <model name>	Name to save custom model in output directory
		--k/-k <scale name>	Name to save custom model's scale in output directory
		--s/-s <polyA sites fille>	Tab separated file containing gold standard polyA sites for training model
		--l/-l <chrom,strand,site>	0-based coordinates of chromosome, strand, and site columns in polyA sites file
	Optional arguments
		-h <help>	Prints help
		--version/-v <version>	Prints version
		--d/-d <debugging>	Saves intermediate files to facilitate issuess
		-verbose/-vb <verbose>	Prints progress to standard output
	
		
	-testSeq <input_sequence_file>	A FASTA file that contains human genomic sequences of length 200 nts. 

	-testSS <input_RNA_secondary_structure_file>	An input file that contains the RNA secondary structures of the input sequences.
					The tool expects three most energy efficient RNA secondary structures for each input sequence.
					These RNA secondary structures are generated using [RNAshapes](https://academic.oup.com/bioinformatics/article/22/4/500/184565).

	-o <output_file_name>		Output file name is given using this option. If this option is not used then the tool outputs
					AUC and AUPRC values of the prediction. In order to calculate the AUC and AUPRC values the tool 
					needs ground truth data. The ground truth data is added at the end of the title of each
					sequence. E.g. for a positive sequence example, the title is >chr15_100354095_positive_1; on
					the other hand, the title of a negative sequence example is >chr15_100565120_positive_0. 
	2. Building your own model
EXAMPLE

	python DeepPASTA_polyA_site_prediction_testing.py -testSeq sample_sequence_input.hg19.fa -testSS sample_secondary_structure_input.txt  


### Input and output file of the polyA site prediction model
The model takes two files as input: sequence and RNA secondary structure files. The sequence file is a FASTA file that contains two lines per example.
The first line is the title of the sequence and the second line is the 200 nts sequence. The RNA secondary structure has four lines per example.
The first line is the title and the next three lines for three RNA secondary structures. The model outputs AUC and AUPRC values when -o option
is not used. In order to get the AUC and AUPRC values, the user must give the ground truth values using the title. E.g. title_ground_truth_value; 
for a positive sequence example, the title is >chr15_100354095_positive_1; on the other hand, the title of a negative sequence example is >chr15_100565120_positive_0. 
If the user uses -o option, the model will output the predicted likelihood values in an output file. 

	
## Tissue-specific polyA site prediction
In order to predict tissue-specific polyA sites, please use DeepPASTA_tissue_specific_polyA_site_prediction_testing.py of tissue_specific_polyA_site_prediction folder. 
Sample input files are given in the sample_input directory. Commands to run the tissue-specific polyA site prediction tool:

USAGE

	cd tissue_specific_polyA_site_prediction
	python DeepPASTA_tissue_specific_polyA_site_prediction_testing.py {OPTIONS}

OPTIONS

	-test <input_sequence_file>	A FASTA file that contains human genomic sequences of length 200 nts.
	
	-tests <input_RNA_secondary_structure_file>	An input file that contains the RNA secondary structures of the input sequences. 
                                                The tool expects three most energy efficient RNA secondary structures for each input sequence.
						These RNA secondary structures are generated using [RNAshapes](https://academic.oup.com/bioinformatics/article/22/4/500/184565).

	-testl <ground_truth_label_file>	An input file for the ground truth labels. The ground truth labels help the tool to calculate AUC and AUPRC values.

	-o <output_file_name>	Output file name is given using this option. This option prints the result in an output file.

EXAMPLE

	python DeepPASTA_tissue_specific_polyA_site_prediction_testing.py -test sample_sequence_input.hg19.fa -tests sample_secondary_structure_input.txt -testl sample_tissue_specific_label.txt

### Input and output of tissue-specific polyA site prediction model
The model takes two files: sequence and RNA secondary structure files. The sequence file is a FASTA file that contains two lines per example.
The first line is the title of the sequence and the second line is the 200 nts sequence. The RNA secondary structure has four lines per example.
The first line is the title and the next three lines for three RNA secondary structures. The model outputs AUC and AUPRC values when -testl option
is used. Using -testl option the user have to give the ground truth data. For each example, the ground truth data has two lines: these two lines are 
title and read counts (of nine tissues separated by comma), respectively. If -o option is used the model outputs tissue-specific polyA site prediction 
in a file. For each input example, the output file has a line containing title and nine likelihood values (separated by comma) for nine tissues.


## Tissue-specific relatively dominant polyA site prediction
The tool can also predict tissue-specific relatively dominant polyA sites. The files to run the tissue-specific relatively dominant polyA
site prediction model are in tissue_specific_relatively_dominant folder. For an example, if an user wants to run the liver tissue relatively
dominant polyA site prediction model, he/she have to follow the following commands:

USAGE
	
	cd tissue_specific_relatively_dominant/tissue_set_one/liver
	python DeepPASTA_relatively_dominant_liver_testing.py {OPTIONS}

OPTIONS
	
	-test <input_file>	An input file that contains the gene name, two polyA sites, two sequences around the polyA sites, the sequences 
				corresponding RNA secondary structures. If the user wants to caculate the AUC and AUPRC values of the prediction,
				he/she also have to provide the polyA sites corresponding read counts in this file.

	-o <output_file_name>	This option prints the result in an output file.

EXAMPLE

	python DeepPASTA_relatively_dominant_liver_testing.py -test sample_relatively_dominant_input_liver.txt

### Input and output of tissue-specific relatively dominant polyA site prediction model
If the user wants to calculate the AUC and AUPRC values of the prediction, he/she must input a file that contains the gene name, first polyA site location, second
polyA site location, sequence (200 nts) around the first polyA site, sequence (200 nts) around the second polyA site, read count of the first polyA site, 
read count of the second polyA site, RNA secondary structure of the first sequence, and RNA secondary structure of the second sequence. If the user wants 
to output the prediction result to a file using -o option, the input file must contains all the above information except the read counts. The output file contains 
the gene name, the two polyA site locations and the probabilities of relative dominance of these two polyA sites. Sample input files of the tissue-specific relatively 
dominant polyA site prediction model are given in the sample_input folder. 


## Tissue-specific absolutely dominant polyA site prediction
The tool can predict tissue-specific absolutely dominant polyA site when sequence (and RNA secondary structure) around a polyA site of a gene 
is given as input. The files to run the tissue-specific absolutely dominant polyA site prediction model are in tissue_specific_absolutely_dominant
folder. For an example, if an user wants to run the liver tissue absolutely dominant polyA site prediction model, he/she have to follow the following
commands:

USAGE

	cd tissue_specific_absolutely_dominant/liver
	python DeepPASTA_absolutely_dominant_liver_testing.py {OPTIONS}

OPTOINS

	-test <input_file>	An input file that contains the gene name, tissue name, polyA site location, sequence around the polyA site, RNA secondary structure 
				of the sequence. If the user wants to calculte the AUC and AUPRC values of the prediction, he/she have to provide the
				label of the polyA site in this file.

	-o <output_file_name>	This option prints the result in an output file.

EXAMPLE

	python DeepPASTA_absolutely_dominant_liver_testing.py -test sample_absolutely_dominant_input_liver.txt

### Input and output of tissue-specific absolutely dominant polyA site prediction model
If the user wants to calculate the AUC and AUPRC values of the prediction, he/she must input a file that contains the gene name, tissue name, polyA site location, sequence
(200 nts) around the polyA site, RNA secondary structure of the sequence, and absolutely dominant label. If the user wants to output the prediction result to a file using
-o option, the input file must contains all the above information except the label. The output file contains the gene name, polyA site location and likelihood value of the 
prediction. Sample input files of the tissue-specific absolutely dominant polyA site prediction model are given in the sample_input folder.

## Creating environment for DeepPASTA
After installing anaconda2, please run the following commands to create the environment for DeepPASTA:

	conda create -n DeepPASTA_env python=2.7.5
	source activate DeepPASTA_env
	conda install -c anaconda keras
	conda install -c anaconda scikit-learn

## How to generate secondary structures from sequences?
The first step of generating secondary structures uses RNAshapes (Steffen *et al.*, 2005). The output from RNAshapes is then converted to DeepPASTA suitable format
using a sequence of commands. The necessary commands to generate DeepPASTA suitable secondary structures from sequences are given below: 

	./RNAshapes -f <input_file_in_fa_format> -s -c 5 -t 1 -w 100 -W 100 -O 'D{%s\n}' > <output_file_in_txt_format>
	python combining_substructure.py -i <output_from_the_previous_command> -o <output_file_in_txt_format>
	python filtering_number_of_ss.py -n 3 -i <output_from_the_previous_command> -o <output_file_in_txt_format>
	python shape_assign_per_nucleotide.py -c 3 -i <output_from_the_previous_command> -o <output_file_in_txt_format>

RNAshapes (Steffen *et al.*, 2005) is provided in the generating_secondary_structure_from_sequence directory for user convenience. Please follow the manual of RNAshapes 
(Steffen *et al.*, 2005) for the possible options. In filtering_number_of_ss.py and shape_assign_per_nucleotide.py scripts, 3 is used to generate the three most probable 
secondary structreus for a given sequence. 

## Genome-wide polyA site prediction for human
We have used DeepPASTA to perform a genome-wide polyA site prediction for human based on the PolyA-Seq data in (Derti *et al.*, 2012).
The prediction results can be found [**here**](http://www.cs.ucr.edu/~aaref001/genome_wide_prediction/genome_wide_polyA_site_prediction_human.txt).

## Reference
1. Derti, A. *et al.*, (2012) A quantitative atlas of polyadenylation in five mammals. *Genome Research*, **22** (6), 1173-1183 
2. Steffen, P. *et al.*, (2005) RNAshapes: an integrated RNA analysis package based on abstract shapes. *Bioinformatics*, **22** (4), 500-503
 
Note: If you have any question or suggestion please feel free to email: aaref001@ucr.edu or ashraful.arefeen@csebuet.org
