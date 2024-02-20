import pandas as pd
import os
import re
from statistics import mean
import subprocess
from Bio import SeqIO
import gzip
import shutil
import tempfile
import argparse
import sys
from datetime import datetime


######################  INPUTS  ###############################

def string_to_list(str):
    str=str.replace(' ','')
    list = str.split(',')
    return list


configfile: "config.yaml"

command_line_args = sys.argv[1:]

snakefile_arg_index = command_line_args.index("--snakefile") + 1



table=command_line_args[snakefile_arg_index].replace("snakefile","table_configuration.py")
create_abricate_db=command_line_args[snakefile_arg_index].replace("snakefile","create_abricate_file.py")

db = str(config["database"])
sample_path = str(config["sample_directory"])


threads = int(config["threads"])
minid=int(config["minid"])
mincov=int(config["mincov"])
output=str(config["output_name"])
output_directory=str(config["output_directory"])
tec_input=str(config["input_format"])
tec_input=string_to_list(tec_input)
multi_fasta=str(config["multi_fasta"])
multi_fasta=string_to_list(multi_fasta)
startbase=int(config["startbase"])
endbase=int(config["endbase"])
fasta_db=str(config["fasta_db"])
table_db=str(config["table_db"])
illuminaclip_single=str(config["illuminaclip_single"])
illuminaclip_paired=str(config["illuminaclip_paired"])
headcrop_single=str(config["headcrop_single"])
headcrop_paired=str(config["headcrop_paired"])
crop_single=str(config["crop_single"])
crop_paired=str(config["crop_paired"])
slidingwindow_single=str(config["slidingwindow_single"])
slidingwindow_paired=str(config["slidingwindow_paired"])
minlen_single=str(config["minlen_single"])
minlen_paired=str(config["minlen_paired"])
leading_single=str(config["leading_single"])
leading_paired=str(config["leading_paired"])
trailing_single=str(config["trailing_single"])
trailing_paired=str(config["trailing_paired"])
quality=int(config["quality"])
length=int(config["length"])
maxlength=int(config["maxlength"])
headcrop=int(config["headcrop"])
tailcrop=int(config["tailcrop"])
kmer=int(config["kmer"])
polishing=int(config["polishing"])
encoding_single=str(config["encoding_single"])
sample_name=str(config["sample_name"])
sample_name=string_to_list(sample_name)


if encoding_single=="in_file":
    encoding_single=str()
else:
    encoding_single="-"+encoding_single

encoding_paired=str(config["encoding_paired"])
if encoding_paired=="in_file":
    encoding_paired=str()
else:
    encoding_paired="-"+encoding_paired
sort=str(config["prioritize"])
if sort=="cov":
    sort="%COVERAGE"
elif sort=="id":
    sort="%IDENTITY"




command_line="ReporType" + ' ' + ' '.join(command_line_args)
data = datetime.now().strftime("%Y-%m-%d")



###################### DB INPUT ###############################


def path_check(path):
    exist = os.path.exists(path)
    return (exist)


def db_check(db):
    for i in shell("abricate --list", iterable=True):
        if db in i:
            check = True
            break
        else:
            check = False
    return (check)


def fasta_import_to_abricate(path):
    path = path.replace('"', "")
    path = str(path)
    return (path)


def abri_default_db():
    x = subprocess.check_output("command abricate --help | grep datadir", shell=True)
    x = str(x)
    x = x.strip()
    x = x.split("]")
    x = x[-2]
    x=x.split("[")
    path=x[-1]
    return (path)


def name_db(path):
    name = path.strip()
    name = name.split("/")
    name = name[-1]
    name = name.split(".")
    name = name[0]
    return (name)


def main_db(db):
    fasta_check = path_check(db)

    if fasta_check == False:
        check = db_check(db)
        if check == False:
            print("Input error: The database name is incorrect or the given path to fasta file does not exist!")
            sys.exit()
        else:
            print('Database available!')
    else:
        fasta_path = fasta_import_to_abricate(db)
        
        join = ['"', '"']
        fasta_path = fasta_path.join(join)
        
        db = name_db(fasta_path)
        db=db.replace('"','')
        check = db_check(db)
        if check == False:
            abricate_db = abri_default_db()
            shell("mkdir {abricate_db}/{db}")
            shell("cp {fasta_path} {abricate_db}/{db}/sequences")
            shell("cd {abricate_db}/{db} ; makeblastdb -in sequences -title {db} -dbtype nucl -hash_index")
            shell("abricate --setupdb")
            check=db_check(db)
            if check == False:
                print("Input error: Failed creating new database! Check the fasta file and respective path.")
            else:            
                print('Database available! ')
        else:
            
            abricate_db = abri_default_db()
            shell("rm {abricate_db}/{db}/*")
            shell("cp {fasta_path} {abricate_db}/{db}/sequences")
            shell("cd {abricate_db}/{db} ; makeblastdb -in sequences -title {db} -dbtype nucl -hash_index")
            shell("abricate --setupdb")
            check=db_check(db)
            if check == False:
                print("Input error: Failed creating new database! Check the fasta file and respective path.")
            else:            
                print('Database available!')
    return(db)


  
if sample_path=="path/to/my_samples_folder/":
    if fasta_db=="path/to/sequences.fasta":
        db=main_db(db)
    else:
        shell('python {create_abricate_db} -f {fasta_db} -s {table_db} -o {db}')
        db=main_db(db)
    sys.exit("Exiting because no samples were submitted")


if fasta_db=="path/to/sequences.fasta":
    db=main_db(db)
else:
    shell('python {create_abricate_db} -f {fasta_db} -s {table_db} -o {db}')
    db=main_db(db)



############## SAMPLES INPUT ##############################


FASTA_EXTENSION = [".fasta",".fsta", ".fasta.gz", ".fa", ".fa.gz", ".fna", ".fna.gz", ".fas"]
FASTQ_EXTENSION = [".fastq", ".fastq.gz", ".fq", ".fq.gz"]
AB1_EXTENSION = [".ab1",".ab",".fsa"]

### size limit of the sequences
MAX_LENGHT_ILLUMINA_FASQC_SEQ = 500     ## because of IONtorrent
MIN_LENGHT_MINION_FASQC_SEQ = 100

################ SEPARATION #######################

def __test_ends_with(file_name, list_ends):
    """ check if the file as this etention """
    for end_str in list_ends:
        if file_name.endswith(end_str): return True
    return False

def __is_nanopore(file_name):
    """
    return True for nanopore, default illumina or IONtorrent
    raise exception if can not detected
    """

    vect_length = []
    with (gzip.open(file_name, mode='rt') if file_name.endswith(".gz") \
        else open(file_name, mode='r')) as handle_read:
        ### read 100 lines
        try:
            count = 0
            for record in SeqIO.parse(handle_read, "fastq"):
                vect_length.append(len(str(record.seq)))
                count += 1
                if (count > 100): break
                #print("mean(vect_length): {} ".format(mean(vect_length)))
        except:
            pass

    ### if read something in last SeqIO.parse
    if (len(vect_length) > 1):
        if (mean(sorted(vect_length, reverse=True)[:5]) <= MAX_LENGHT_ILLUMINA_FASQC_SEQ): return False
        if (mean(vect_length) > MIN_LENGHT_MINION_FASQC_SEQ): return True
        raise Exception("Can not detect file format. Ensure Illumina fastq file.")
    raise Exception("File is not in fastq.gz format.")


def __is_fasta_file(file_name):
    """ test fasta format """
    with (gzip.open(file_name, mode='rt') if file_name.endswith(".gz") \
        else open(file_name, mode='r')) as handle_read:
        try:
            for record in SeqIO.parse(handle_read, "fasta"):
                handle_read.close()
                return True
        except:
            pass
    return False

def remove_extensions_file_name(file_name, vect_to_remove=FASTQ_EXTENSION):
    if os.path.exists (file_name):
        file_name=file_name.strip()
        file_name = file_name.split("/")
        file_name = file_name[-1]
    for to_search in vect_to_remove:
        if (file_name.endswith(to_search)):
            return file_name[:len(file_name) - len(to_search)]
    return file_name


def __get_prefix_file_name(file_name):
    """ returnprefix file name based on patterns """
    m = re.search('[a-zA-Z0-9_\.]+(_[lL]\d+_[rR]\d+_\d+)_[a-zA-Z0-9_\.]+', file_name)
    if (not m is None): return file_name[:m.regs[1][0]]
    m = re.search('[a-zA-Z0-9_\.]+(_[lL]\d+_[rR]\d+)[a-zA-Z0-9_\.]+', file_name)
    if (not m is None): return file_name[:m.regs[1][0]]
    m = re.search('[a-zA-Z0-9_\.]+([_\.][rR]\d+_[lL]\d+)[a-zA-Z0-9_\.]+', file_name)
    if (not m is None): return file_name[:m.regs[1][0]]
    m = re.search('[a-zA-Z0-9_\.]+(_[lL]\d+)[a-zA-Z0-9_\.]+', file_name)
    if (not m is None): return file_name[:m.regs[1][0]]
    m = re.search('[a-zA-Z0-9_\.]+([_\.][rR]\d)[_\.][a-zA-Z0-9_\.]+', file_name)
    if (not m is None): return file_name[:m.regs[1][0]]
    m = re.search('[a-zA-Z0-9_\.]+([_\.][rR]\d+)[a-zA-Z0-9_\.]+', file_name)
    if (not m is None): return file_name[:m.regs[1][0]]
    m = re.search('[a-zA-Z0-9_\.]+([_]reverse)[a-zA-Z0-9_\.]+', file_name.lower())
    if (not m is None): return file_name[:m.regs[1][0]]
    m = re.search('[a-zA-Z0-9_\.]+([_]forward)[a-zA-Z0-9_\.]+', file_name.lower())
    if (not m is None): return file_name[:m.regs[1][0]]
    m = re.search('[a-zA-Z0-9_\.]+(_\d+)[\.][a-zA-Z0-9_\.]+', file_name)
    if (not m is None): return file_name[:m.regs[1][0]]
    m = re.search('[a-zA-Z0-9_\.]+(_\d+)[_\.][a-zA-Z0-9_\.]+', file_name)
    if (not m is None): return file_name[:m.regs[1][0]]
    m = re.search('[a-zA-Z0-9_\.]+(_\d+)[_][a-zA-Z0-9_\.]+', file_name)
    if (not m is None): return file_name[:m.regs[1][0]]
    return remove_extensions_file_name(file_name)


def get_number_file(file_name):
    """ return the number of file
    :out None if not found  """

    m = re.search('[a-zA-Z0-9_\.]+([_\.][rR]\d+)[a-zA-Z0-9_\.]+', file_name)
    if (not m is None): return file_name[m.regs[1][0]:m.regs[1][1]].lower().replace('r', '').replace('_', '').replace('.', '')

    m = re.search('[a-zA-Z0-9_\.]+([_\.]\d+)[_\.][a-zA-Z0-9_\.]+', file_name)
    if (not m is None): return file_name[m.regs[1][0]:m.regs[1][1]].lower().replace('_', '').replace('.', '')

    m = re.search('[a-zA-Z0-9_\.]+([_\.]reverse)[a-zA-Z0-9_\.]+', file_name.lower())
    if (not m is None): return file_name[m.regs[1][0]:m.regs[1][1]].lower().replace('_', '').replace('.', '').replace('reverse', '2')

    m = re.search('[a-zA-Z0-9_\.]+([_\.]forward)[a-zA-Z0-9_\.]+', file_name.lower())
    if (not m is None): return file_name[m.regs[1][0]:m.regs[1][1]].lower().replace('_', '').replace('.', '').replace('forward', '1')
    return None


def get_order_pair(b_test_first_in_pair, vect_names):
    """ from a list retruen the first that belongs to first or second pair """
    for path_and_file in vect_names:
        file_name = os.path.basename(path_and_file)
        if b_test_first_in_pair and get_number_file(file_name) == "1": return path_and_file
        if not b_test_first_in_pair and get_number_file(file_name) == "2": return path_and_file
    return None


def collect_files(*args):
    ### list all possible files
    list_fasta = []
    list_fastq_illumina_single = []
    ## paired end files { 'r1' : [ 'file_1.fastq.gz', 'file_r1.fastq.gz', .. ]
    ##                    'r2' : [ 'file_2.fastq.gz', 'file_r2.fastq.gz', .. ] }
    dt_fastq_illumina_pair = { 'r1' : [], 'r2' : [] }
    list_fastq_nanopore = []
    list_ab1 = []

    ## list all files names
    dt_list_all_files_names = {}

    for item in args:
        for p, _, files in os.walk(os.path.abspath(item)):
            for file in files:

                ## file repeated
                path_file = os.path.join(p, file)
                if file in dt_list_all_files_names:
                    print("File name repeated: " + file)
                    continue

                ## testing file
                if __test_ends_with(path_file, FASTA_EXTENSION):
                    if __is_fasta_file(path_file):
                        dt_list_all_files_names[file] = 1
                        list_fasta.append(path_file)
                    else:
                        print("Error: some problem with this fasta file: " + file)
                    continue
                if __test_ends_with(path_file, FASTQ_EXTENSION):
                    try:
                        if __is_nanopore(path_file):
                            list_fastq_nanopore.append(path_file)
                        else: ### first, all of them are single end
                            list_fastq_illumina_single.append(path_file)
                        dt_list_all_files_names[file] = 1
                    except Exception as e:
                        print("Error: " + str(e) + "\nFile: " + path_file)
                    continue
                if __test_ends_with(file, AB1_EXTENSION):
                    dt_list_all_files_names[file] = 1
                    list_ab1.append(path_file)
                    continue

                ## didn't found any technology
                print("File not recognize: " + path_file)

    ### try to find close
    dt_fastq_illumina_pair_temp = {}
    for path_and_file_name in list_fastq_illumina_single:
        sample_name = __get_prefix_file_name(os.path.basename(path_and_file_name))

        if sample_name in dt_fastq_illumina_pair_temp: dt_fastq_illumina_pair_temp[sample_name].append(path_and_file_name)
        else: dt_fastq_illumina_pair_temp[sample_name] = [path_and_file_name]

    ### get the candidates and join them
    dt_to_remove = {}   ## files to remove because are in pair
    for sample_name in dt_fastq_illumina_pair_temp:
        first_in_pair = get_order_pair(True, dt_fastq_illumina_pair_temp[sample_name])
        second_in_pair = get_order_pair(False, dt_fastq_illumina_pair_temp[sample_name])

        if not first_in_pair is None and not second_in_pair is None and first_in_pair != second_in_pair:
            dt_to_remove[first_in_pair] = 1
            dt_to_remove[second_in_pair] = 1
            dt_fastq_illumina_pair['r1'].append(first_in_pair)
            dt_fastq_illumina_pair['r2'].append(second_in_pair)

    ###  remove the ones in the paired files
    for key in dt_to_remove:
        list_fastq_illumina_single.remove(key)

    ## test files
    if len(list_fasta) + len(list_fastq_illumina_single) +\
           (len(dt_fastq_illumina_pair['r1']) * 2) + len(list_fastq_nanopore) + len(list_ab1)== 0:
        sys.exit("Error, no file founded at '{}' path".format(args[0]))

    return list_fasta, list_fastq_illumina_single, dt_fastq_illumina_pair,\
           list_fastq_nanopore, list_ab1




########################PRE ANALISIS###############################

## test exists
exist = path_check(sample_path)

if exist == False:
    print("Input error: The sample path given is not correct!")
else:
    ## test if windows or not and add the back slash to the end
    if os.name == 'nt':
        sample_path += "\\"
   # elif not sample_path.endswith('/'):
    #    sample_path += '/'

    list_fasta, list_illumina_fastq_single, dt_fastq_illumina_pair,\
        list_fastq_nanopore, list_ab1 = collect_files(sample_path)


### get all samples name, till the end
SAMPLES_NAME_global = []


##############TECNOLOGY CONFIRMATION#################

### get all the detected technologies
tec_app=[]
if len(dt_fastq_illumina_pair['r1']) > 0:
    tec_app.append("illumina_paired")

if len (list_illumina_fastq_single) > 0:
    tec_app.append("illumina_single")

if len (list_fasta)>0:
    tec_app.append("fasta")
    
if len (list_fastq_nanopore)>0:
    tec_app.append("nanopore")

if len (list_ab1)>0:
    tec_app.append("sanger")

### check useres input
for tecnology_input in tec_input:
    if tecnology_input == "any":
        tec_input.remove("any")

###remove non listed tecnoloogies

if len(tec_input)>0:
    if "nanopore" not in tec_input:
        list_fastq_nanopore=[]
    if "fasta" not in tec_input:
        list_fasta=[]
    if "sanger" not in tec_input:
        list_ab1=[]
    if "illumina_single" not in tec_input:
        list_illumina_fastq_single=[]
    if "illumina_paired" not in tec_input:
        dt_fastq_illumina_pair = { 'r1' : [], 'r2' : [] }





########################ANALISIS###############################
rule all:
    input:
        expand("{output_directory}/{output}.tsv",output=output, output_directory=output_directory)

    params:
        tec_input=tec_input,
        tec_app=tec_app,
        data=data,
        command_line=command_line
    run:
        print(f"Inputs for rule all: {input}")
        logfile=output_directory+"/logfile.log"
        with open(logfile, "w") as file:
            file.write("ReporType log file\n\n")
            file.write("Date: " + params.data + "\n")
            file.write("Snakemake execution: " + params.command_line + "\n")
        for tecnology_input in params.tec_input:
            if len(tec_input) > 0:
                if tecnology_input not in params.tec_app:
                    print("Warning: no files detect for", tecnology_input, "tecnology. Verify if you enter the rigth tecnology.")
        if len (params.tec_input) > 0:
            for tecnology_app in params.tec_app:
                if tecnology_app not in params.tec_input:
                    print("Warning: there are samples that are not", ' or '.join(tec_input), "files, those samples will not be analysed. Verify if you enter the rigth tecnologies.")
                    break
        



if len(dt_fastq_illumina_pair['r1']) > 0:

    ## list of the samples
    SAMPLES_NAME_global.extend([__get_prefix_file_name(os.path.basename(file_name)) for file_name in dt_fastq_illumina_pair['r1']])

    dt_extention_r1 = { __get_prefix_file_name(os.path.basename(file_name)) : dt_fastq_illumina_pair['r1'][index] for index, file_name in enumerate(dt_fastq_illumina_pair['r1']) }
    dt_extention_r2 = { __get_prefix_file_name(os.path.basename(file_name)) : dt_fastq_illumina_pair['r2'][index] for index, file_name in enumerate(dt_fastq_illumina_pair['r2']) }

    rule pre_illumina_paired:
        input:
            r1 = lambda wildcards: dt_extention_r1[wildcards.sample],
            r2 = lambda wildcards: dt_extention_r2[wildcards.sample]
        output:
            s1=expand("{output_directory}/intermediate/trimm_paired_sur_1/{{sample}}.fastq.gz", output_directory=output_directory),
            d1=expand("{output_directory}/intermediate/trimm_paired_rem_1/{{sample}}.fastq.gz", output_directory=output_directory),
            s2=expand("{output_directory}/intermediate/trimm_paired_sur_2/{{sample}}.fastq.gz", output_directory=output_directory),
            d2=expand("{output_directory}/intermediate/trimm_paired_rem_2/{{sample}}.fastq.gz", output_directory=output_directory)
        
        params:
            threads = threads,
            illuminaclip_paired=illuminaclip_paired,
            headcrop_paired=headcrop_paired,
            crop_paired=crop_paired,
            slidingwindow_paired=slidingwindow_paired,
            minlen_paired=minlen_paired,
            leading_paired=leading_paired,
            trailing_paired=trailing_paired,    
            encoding_paired=encoding_paired       
        shell:
            """trimmomatic PE -threads {params.threads} {params.encoding_paired} {input.r1} {input.r2} {output.s1} {output.d1} {output.s2} {output.d2} {params.illuminaclip_paired}  {params.headcrop_paired} {params.crop_paired} {params.slidingwindow_paired} {params.leading_paired} {params.trailing_paired} {params.minlen_paired} || (touch {output.s1} {output.s2} && echo Warning: trimmomatic failed, were created empty files)"""

    rule illumina_paired:
        input:
            r1 =expand("{output_directory}/intermediate/trimm_paired_sur_1/{{sample}}.fastq.gz", output_directory=output_directory),
            r2 =expand("{output_directory}/intermediate/trimm_paired_sur_2/{{sample}}.fastq.gz", output_directory=output_directory)
        output:
             expand("{output_directory}/intermediate/fasta_and_fai_files/{{sample}}.fasta", output_directory=output_directory)
        params:
            threads = threads,
            out=output_directory
        shell:
            """spades.py --pe1-1 {input.r1} --pe1-2 {input.r2} --only-assembler -t {params.threads} -o {params.out}/intermediate/spades/{wildcards.sample} ||  echo Warning: spades failed, were created empty files;"""
            """cp {params.out}/intermediate/spades/{wildcards.sample}/contigs.fasta {output} || touch {output}"""#; rm -rf results/intermidiate/spades/{wildcards.sample}"""


if len(list_illumina_fastq_single) > 0:

    ### add sample names for global
    SAMPLES_NAME_global.extend(__get_prefix_file_name(os.path.basename(file_name)) for file_name in list_illumina_fastq_single)

    dt_extention_single = { __get_prefix_file_name(os.path.basename(file_name)) : list_illumina_fastq_single[index] for index, file_name in enumerate(list_illumina_fastq_single) }

    rule pre_illumina_single:
        input:
            r1 = lambda wildcards: dt_extention_single[wildcards.sample]
        output:
            s1=expand("{output_directory}/intermediate/trimm_single_sur_1/{{sample}}.fastq.gz", output_directory=output_directory)
        params:
            threads = threads,
            illuminaclip_single=illuminaclip_single,
            headcrop_single=headcrop_single,
            crop_single=crop_single,
            slidingwindow_single=slidingwindow_single,
            minlen_single=minlen_single,
            leading_single=leading_single,
            trailing_single=trailing_single,
            encoding_single=encoding_single
        shell:
            """trimmomatic SE -threads {params.threads} {params.encoding_single} {input.r1} {output.s1} {params.illuminaclip_single} {params.headcrop_single} {params.crop_single} {params.slidingwindow_single} {params.leading_single} {params.trailing_single} {params.minlen_single} || (touch {output.s1} && echo Warning: trimmomatic failed, were created empty files)"""

    rule illumina_single:
        input:
            r1 = expand("{output_directory}/intermediate/trimm_single_sur_1/{{sample}}.fastq.gz", output_directory=output_directory)
        output:
            expand("{output_directory}/intermediate/fasta_and_fai_files/{{sample}}.fasta", output_directory=output_directory)
        params:
            threads = threads,
            out=output_directory
        shell:
            """spades.py -s {input.r1} --only-assembler -t {params.threads} -o {params.out}/intermediate/spades/{wildcards.sample} || echo Warning: spades failed, were created empty files;"""
            """cp {params.out}/intermediate/spades/{wildcards.sample}/contigs.fasta {output} || touch {output}"""#; rm -rf results/intermidiate/spades/{wildcards.sample}"""


if len(list_fasta) > 0:

    ### sample name for fasta

    SAMPLES_NAME_global.extend(remove_extensions_file_name(file_name,FASTA_EXTENSION) for file_name in list_fasta)

    dt_extention_fasta_file = {remove_extensions_file_name(os.path.basename(file_name),FASTA_EXTENSION) : list_fasta[index] for index, file_name in enumerate(list_fasta)}
 

   
    rule fasta:
        input:
            lambda wildcards: dt_extention_fasta_file[wildcards.sample]
        output:
            expand("{output_directory}/intermediate/fasta_and_fai_files/{{sample}}.fasta", output_directory=output_directory)
        shell:
            "cp {input} {output}"


if len(list_fastq_nanopore) > 0:

    ### sample name for nanopore
    SAMPLES_NAME_global.extend(remove_extensions_file_name(file_name) for file_name in list_fastq_nanopore)

    dt_extention_nano = { remove_extensions_file_name(os.path.basename(file_name)) : list_fastq_nanopore[index] for index, file_name in enumerate( list_fastq_nanopore)}
    
    rule pre_nanopore:
        input:
            lambda wildcards: dt_extention_nano[wildcards.sample]
        params:
            length=length,
            quality=quality,
            maxlength=maxlength,
            headcrop=headcrop,
            tailcrop=tailcrop
        output:
            expand("{output_directory}/intermediate/nanofilt_filtred_files/{{sample}}.fastq.gz", output_directory=output_directory)
        shell:
            """(gunzip -c {input} | NanoFilt -q {params.quality} -l {params.length} --headcrop {params.headcrop} --tailcrop {params.tailcrop} --maxlength {params.maxlength} | gzip > {output}) || (touch {output} && echo Warning: nanofilt failed, were created empty files)"""
 

    rule nanopore:
        input:
            expand("{output_directory}/intermediate/nanofilt_filtred_files/{{sample}}.fastq.gz", output_directory=output_directory)
        output:
            gfa = expand("{output_directory}/intermediate/raven_gfa_files/{{sample}}.fasta", output_directory=output_directory),
            fasta = expand("{output_directory}/intermediate/fasta_and_fai_files/{{sample}}.fasta", output_directory=output_directory)
        params:
            threads=threads,
            kmer=kmer,
            polishing=polishing
        shell:
            """raven -t {params.threads} -k {params.kmer} -p {params.polishing} {input} --graphical-fragment-assembly {output.gfa} --disable-checkpoints || (touch {output.gfa} && echo Warning: raven failed, were created empty files) ; awk \'/^S/{{ printf(">%s\\n%s\\n", $2, $3) }}\' {output.gfa} > {output.fasta} || touch {output.fasta} """


if len(list_ab1) > 0:

    SAMPLES_NAME_global.extend(remove_extensions_file_name(file_name,AB1_EXTENSION) for file_name in list_ab1)

    dt_extention_sanger = { remove_extensions_file_name(os.path.basename(file_name),AB1_EXTENSION) : list_ab1[index] for index, file_name in enumerate(list_ab1)}

    rule sanger:
        input:
            lambda wildcards: dt_extention_sanger[wildcards.sample]
        output:
            expand("{output_directory}/intermediate/fasta_and_fai_files/{{sample}}.fasta", output_directory=output_directory)
        params:
            startbase=startbase,
            endbase=endbase
        shell:
            "abiview {input} -graph none -startbase {params.startbase} -endbase {params.endbase} -osformat2 fasta -outseq {output}  || (touch {output} && echo Warning: abiview failed, were created empty files)"
           



rule abricate:
    input:
        expand("{output_directory}/intermediate/fasta_and_fai_files/{{sample}}.fasta", output_directory=output_directory)
    output:
        expand("{output_directory}/detailed/{{sample}}.tab", output_directory=output_directory)
    params:
        db = db,
        minid=minid,
        mincov=mincov
    shell:
        "abricate --db {params.db} --nopath {input} --minid {params.minid} --mincov {params.mincov} > {output} ;"

if sample_name!=["all"]:
    SAMPLES_NAME_global=sample_name


if len(SAMPLES_NAME_global)==0:
 
    rule no_samples:
        output:
            expand("{output_directory}/{output}.tsv",output=output, output_directory=output_directory)
        shell:
            "touch {output}"


else:
    rule table_match:
        input:
            fasta_file=expand("{output_directory}/intermediate/fasta_and_fai_files/{sample}.fasta",output_directory=output_directory,sample=SAMPLES_NAME_global),
            tab_file=expand("{output_directory}/detailed/{sample}.tab",output_directory=output_directory,sample=SAMPLES_NAME_global)
        output:
            table=expand("{output_directory}/{output}.tsv",output=output, output_directory=output_directory),
            fasta_output=expand("{output_directory}/match_regions/match_region_{sample}.fasta",sample=SAMPLES_NAME_global, output_directory=output_directory)
        params:
            multi=multi_fasta,
            sort=sort,
            table=table  

        run:
            shell("python {table} {input.tab_file} {output.table} {params.multi} {params.sort}")
            for fasta_file, tab, fasta_out in zip(input.fasta_file, input.tab_file, output.fasta_output):
                table_match = pd.read_csv(tab, sep="\t")
                
                if table_match.empty:
                    sample=tab.split('/')[-1].replace('.tab','') 
                    path_file = f"{output_directory}/match_regions/match_region_{sample}.fasta"
                    note = "No match regions found."
                    with open(path_file, 'w') as file:
                        file.write(note)

                else:
                    sequence_id = table_match.loc[:, "SEQUENCE"].tolist()
                    start = table_match.loc[:, "START"].tolist()
                    end = table_match.loc[:, "END"].tolist()
                    gene = table_match.loc[:,"GENE"].tolist()
                    database=table_match.loc[:,"DATABASE"].tolist()
                    sequence_info = {}
                    for seq, st, en, gen, db in zip(sequence_id, start, end, gene, database):
                       if seq not in sequence_info:
                           sequence_info[seq] = {'start': st, 'end': en, 'gene': gen,'database': db }
                       intervals = f"{seq}:{st}-{en}"
                       intervals=intervals.replace('|', '\\|')
                       
                       shell(f"samtools faidx {fasta_file} {intervals} >> {fasta_out}")
                       
                       with open(fasta_out, "r") as infile:
                           lines = infile.readlines()
                       for i in range(len(lines) - 1, -1, -1):
                           if lines[i].startswith('>'):
                               lines[i] = f">{db}_{gen}_{intervals}" + "\n"
                               break
             
                       with open(fasta_out, "w") as outfile:
                           outfile.writelines(lines)


           
