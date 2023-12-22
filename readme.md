**ReporType**

ReporType is an automatic, easy-to-use and flexible pipeline, created with Snakemake, for loci screening and typing. Its application can be particularly useful for rapid genotyping of infectious agents, namely virus and bacteria.

ReporType was designed to accept multiple input formats (from Illumina or ONT reads to Sanger raw files or FASTA files), being suitable for application in a wide variety of pathogens. It relies on multiple software for technology-specific reads QC and de novo assembly, and thus apply ABRicate (https://github.com/tseemann/abricate) for locus screening, culminating in the generation of easy-to-interpret reports towards the identification of pathogen genotypes/subspecies or the screening of loci of interest.

ReporType comes with pre-prepared databases for genotyping of a few virus/bacteria, but can be easily setup to handle custom databases, instructions below. You can also change several analysis parameters, as well as modify parameters of each software used.
The final report consists of a document in table format containing the most relevant results, such as sample name, element found (such as genotype, subspecies, etc), coverage and percentage of identity, the database used and accession number. You will also be able to access detailed ABRIcate output files and intermediate files that are produced by other software (clipped samples, fasta files, etc...).


![alt text](https://github.com/ibigen/ReporType/blob/main/ReporType_workflow.png)



## Instalation
You need to have  [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html) installed.
All the other dependencies will be automatically installed with ReporType.
For installation, you need to:


1. Download this git repository:<br>
`$ git clone https://github.com/ibigen/ReporType/`<br>
`$ cd ReporType`

2. Install running:<br>
`$ chmod +x install.sh`<br>
`$ ./install.sh`<br>

### Databases instalation

Before installing the databases, it is necessary to activate the conda environment created for ReporType to work. You can activate the environment with the activation command: <br>

`$ alias ReporType='conda activate ReporType && snakemake'; conda activate ReporType`<br>

Then install the databases running:<br>
`$ chmod +x databases_install.sh`<br>
`$ ./databases_install.sh`<br>
## Usage

First of all, you need to activate the ReporType environment with the command:<br>
`$ alias ReporType='conda activate ReporType && snakemake'; conda activate ReporType`<br>


Now you must configure your entery params. You have to options, you can open de "config.yaml" file and fill it with your options you configurate them through the command line.<br>

There are some mandatory params for configuration listed below. <br>

**Database input params:** <br>

If you have already install the incorporated databases or created your own: <br>
> **database**: name of the database you wish to use (example: database=my_database).<br>
 
If is the first time using a new database you need to add the path to the formated fasta file (```seq~~~id~~~acession```) contaning the database, a new database will be created with the name of the given fasta file:<br>
> **database**: path to fasta file for new database (example: database=path/to/my_database.fasta).<br>


If you don't have a database file already formatated for abricate, you can provide two files to crate a new database:<br>
Note that, in this case, you should write the name of your new database in the "database" variable.
> **fasta_db**: fasta file with the sequences for your database (example: fasta_db=path/to/sequences.fasta).<br>
> **table_db**:  table (tsv) with three columns: column one (sequence), with the name of each sequence; column two (id), with the identification of each element (genotype, subspecies, etc); and column three (accession), with the acession number for each sequence (example table_db=path/to/table.tsv)<br>
> **database**: name of the database you wish to create (example: database=my_database).<br>

**Samples params:** <br>

> **sample_directory**: path to the folder with the samples you wish to analyse. This folder can contain samples from different technologies, as long as they are all analyzed according to the same database (example: sample_directory=path/to/my_samples_folder/).<br>

> **sample_name**: if you wish to analyse only one sample you must give the sample name, you can provide a list of samples (default=all). Note that in paired end sequences, you must give the sample name without any prefixes.<br>

ReporType optional configuration params includes: <br>

> **output_name**: name of your final csv output file (default: output_name=all_samples)<br>
> **output_directory**: directory for your results (default: output_directory=results/)<br>
> **input_format**: especify the input format you are going to analyse. If you leave it with the default, all samples of the given folder will be analysed. Your opcions are: fasta,nanopore,illumina_single,illumina_paired,sanger, or any. You must separete them with a coma (default: input_format=any)<br>
> **multi_fasta**: if you are going to analyse any multi-fasta files where each sequence corresponds to one sample, give the name of each multi-fasta file. You can chosse "all" if all of your fasta files are multi-fasta(default: multi_fasta=none).<br>
> **threads**: threads you which to use (default: threads=2).<br>
> **prioritize**: in case there is more than one gene detected for the same feature, choose if you want to prioritize greater coverage (cov) or greater identity (id) (default: prioritize=cov).<br>

You can also specify some software params.<br>

**Abricate params:**<br>

> **minid**: minimum DNA %identity (default: minid=1).<br>
> **mincov**: minimum DNA %coverage (default: mincov=1).<br>

**Illumina params:** (for single and paired reads)<br>
> **illuminaclip_single** and **illuminaclip_paired**: Trimmomatic Illuminaclip, directory of your illumina adapters, as well as specific cleaning informations for your file (default: illuminaclip=ILLUMINACLIP:primers/adapters.fasta:3:30:10:6:true)<br>
> **slidingwindow_single** and **slidingwindow_paired**: Trimmomatic Slidingwindow, minimum average quality established for each sequence according to a certain number of bases (default: slidingwindow=SLIDINGWINDOW:4:15).<br>
> **minlen_single** and **minlen_paired**: Trimmomatic Minlen, minimum read size (default: minlen=MINLEN:36).<br>
> **leading_single** and **leading_paired**: Trimmomatic Leading, bases to remove at the beginning of the read (default: leading=LEADING:3)<br>
> **trailing_single** and **trailing_paired**: Trimmomatic Trailing, bases to remove at the end of the read (default: trailing=TRAILING:3)<br>
> **encoding_single** and **encoding_paired**:Trimmomatic encoding: if quality encoding is not specified in fastq file, specify specify the quality encoding (default=in_file, Your options are: 'phred33', 'phred64')

**Nanopore params:**<br>
> **quality**: Nanofilt minimum quality mean for read (default: quality=8).<br>
> **length**: Nanofilt minimun length per read (default: length=50).<br>
> **maxlength**: Nanofilt maximum length per read (default: maxlength=50000)<br>
> **headcrop**: Nanofilt headcrop, bases to remove at the beginning of the read (default: headcrop=30).<br>
> **tailcrop**: Nanofilt tailcrop, bases to remove at the end of the read (default: tailcrop=30).<br>
> **kmer**: Raven k-mer, length of minimizers used to find overlaps (default: kmer=15). <br>
> **polishing**: Raven polishing-rounds, number of times racon is invoked (default: polishing=2).<br>

**Sanger params:**
> **startbase**: Abiview first sequence base to report or display (default: startbase=20).<br>
> **endbase**: Abiview last sequence base to report or display (default: endbase=800). <br>


The optional configuration params also include all the configuration params for Snakemake, that can be consulted [here](https://snakemake.readthedocs.io/en/v5.1.4/executable.html). The most relevant Snakemake executable params are: <br>
> **--cores**: number of CPU to be used, it is mandatory (example: --cores all).<br>
> **-np**: dry-run to verify the jobs you are submiting. <br>
> **--config**: you must use this command before configurate the ReporType params previosly refered (example: --config database=my_database).<br>
> **--snakefile**: you can execute ReporType in any directory using this command to specify the directory for the snakefile of ReporType (example: --snakefile path/to/ReporType/snakefile).<br> 
> **--configfile**: you can execute ReporType in any directory using this command to specify the directory for the config file of ReporType (example: --configfile path/to/ReporType/config.yaml).<br>


## Execution<br>

ReporType is run through the command line, here are some examples, from the simplest to the most complex.

### Configuration with config.yaml file
If you configurate the config.yaml file, you can only run:<br>
`$ ReporType --cores all `<br>

### Configuration with command line<br>

#### Example 1 - Database already used or previosly installed:<br>
`$ ReporType --cores all --config sample_directory=path/to/my_samples_folder/ database=my_database`<br>

#### Example 2 - New database with formatted fasta file: <br>
`$ ReporType --cores all --config sample_directory=path/to/my_samples_folder/ database=path/to/my_database.fasta`<br>

#### Example 3 - New database without formatted fasta file: <br>
`$ ReporType --cores all --config sample_directory=path/to/my_samples_folder/ database=path/to/my_database.fasta fasta_db=path/to/sequences.fasta table_db=path/to/table.tsv`<br>

#### Example 4 - Output params configuration: <br>
`$ ReporType --cores all --config sample_directory=path/to/my_samples_folder/ database=my_database output_name=all_samples output_directory=results`<br>

#### Example 5 - Input format params configuration <br>
##### Example 5.1 - You want to analyze all the samples in your folder and you have two multi fasta files:<br>
`$ ReporType --cores all --config sample_directory=path/to/my_samples_folder/ database=my_database input_format=any multi_fasta=multi_fasta_1,multi_fasta_2`<br>

##### Example 5.2 - You want to analyze all fasta files and samples sequenced with nanopore technology, all your fasta files are multi fasta:<br>
`$ ReporType --cores all --config sample_directory=path/to/my_samples_folder/ database=my_database input_format=fasta,nanopore multi_fasta=all`<br>

#### Example 6 - Configuration of some analysis parameters: <br>
`$ ReporType --cores all --config sample_directory=path/to/my_samples_folder/ database=my_database input_format=fasta,nanopore multi_fasta=all minid=1 mincov=1`<br>

#### Example 7 - To execute a dry run:<br>
`$ ReporType -np --config sample_directory=path/to/my_samples_folder/ database=my_database`<br>

#### Example 8 - To run ReporType out of instalation directory:<br>
`$ ReporType --cores all --snakefile path/to/ReporType/snakefile --configfile path/to/ReporType/config.yaml –-config sample_directory=path/to/my_samples_folder/ database=my_database`<br>


<br>
<br>


When you are donne using ReporType you can deactivate the environment with:<br>
`$ conda deactivate ReporType`<br>


## Uninstall

To uninstall ReporType, you need to delete the conda environment with: <br>
`$ conda env remove --name ReporType`<br>

## Citation
If you run ReporType, please cite this Github page:<br>

Helena Cruz, Miguel Pinheiro, Vítor Borges (2023). ReporType - Flexible bioinformatics tool for targeted loci screening and typing of infectious agents.
https://github.com/insapathogenomics/reportype

 
