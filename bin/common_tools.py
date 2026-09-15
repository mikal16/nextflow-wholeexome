import os
import re
import subprocess
import csv
import sys
from pyfiglet import Figlet
from sample_sheet import SampleSheet
from xlsxwriter.workbook import Workbook
import pandas as pd
import gzip
import math
import numpy as np
from cyvcf2 import VCF

"""
    This file contains functions and classes with more general purposes that would be reused in multiple files
"""


class Genotype(object):
    """
        Class created to patch a problem with cyvcf2 and SAMPLE_GT field when the variant caller is VARSCAN
        Code was provided by creator of cyvcf2
        see https://github.com/brentp/cyvcf2/issues/58
    """

    __slots__ = ("alleles", "phased")

    def __init__(self, genotype):
        self.alleles = genotype[:-1]
        self.phased = genotype[-1]

    def __str__(self):
        sep = "/|"[int(self.phased)]
        return sep.join("0123."[a] for a in self.alleles)

    __repr__ = __str__


def print_header():
    """
        This function prints the tool's name in ascii characters
    """
    f = Figlet(font="slant")
    header = f.renderText("NGSTK")
    return header


def verify_extension(path, extension):
    """
        Verifies that the file extension of path is the one provided.
        Args:
            path: the path of the file
            extension: the extension to revise
        Returns:
            True if the file path ends with the extension, False otherwise
    """
    status = False
    if path.endswith(extension):
        status = True
    # else:
    #     print("File " + path +  " must be of type " + extension)
    #     status = False
    return status


def verify_file_exists(path, print_message=True):
    """
        Verifies that the file provided in parameters exists.
        Args:
            path: the path of the file
            print_message: boolean
        Returns:
            True if the file exists, False otherwise
    """

    if os.path.exists(get_full_path(path)):
        # print("File path passed as an argument is: " + path)
        status = True

    else:
        if print_message:
            print("This file does not exist: " + path)
            # print("Please enter a valid file path")
        status = False
    return status


def verify_multiple_sample_files_exist(samples, directory, extension):
    """
        Verifies that the file provided in parameters exists for each sample.
        Args:
            samples: a list of sample IDs
            directory: the path to the directory where the files should be
            extension: the extension to complete the file names
        Returns:
            True if all the files exist, False otherwise
    """
    validation = True
    validated_directory = validate_directory(directory)
    for sample in samples:
        file_exists = verify_file_exists(validated_directory + sample + extension)
        if not file_exists:
            validation = False
            break
    return validation


def verify_bed_file(bed_file, output_directory):
    """
        Verifies that the bed file provided is well formated with tabs, if spaces are found instead, reformats the bed file.
        Args:
            bed_file: path to the bed file
            output_directory: path to the directory where the reformated bed file is to be written if necessary

        Return:
            string: the path to the bed file

    """

    bed_file = get_full_path(bed_file)
    if verify_file_exists(bed_file):
        with open(bed_file, "r") as f:
            bed = f.readlines()

        for line in bed:
            line = line.rstrip()
            if " " in line:
                bed_file = reformat_bed_file(bed_file, output_directory)
                break

        return bed_file
    else:
        print("ERROR: This bed file doesn't exist: {bed_file}".format(bed_file=bed_file))
        sys.exit(1)


def reformat_bed_file(bed_file, output_directory):
    """
        Reformat the bed file removing spaces and replacing them with tabs.
        Args:
            bed_file: path to the bed file
            output_directory: path to the directory where the reformated bed file is to be written
        Return:
            string: path to the reformated bed file
    """

    reformated_bed_output_directory = output_directory + "reformated_beds/"
    copy_file_in_directory(bed_file, reformated_bed_output_directory)

    bed_file_name = os.path.basename(bed_file)
    reformated_bed_file = reformated_bed_output_directory + bed_file_name

    # Removes spaces at the end of a line, so they don't get replaced by tabs
    subprocess.run(["sed", "-i", "s/ *$//", reformated_bed_file])
    # Replaces blocks of spaces with a single tab
    subprocess.run(["sed", "-i", "s/ \+/\t/g", reformated_bed_file])

    return reformated_bed_file


def create_directories(path, folders):
    """
        Creates the folders specified in the folders list at the path location, if they don't already exist.
        Args:
            path: Path to output directory
            folders: List of folder names
    """
    valid_path = validate_directory(path)
    for folder in folders:
        try:
            if not os.path.exists(valid_path + folder):
                os.makedirs(valid_path + folder)
        except OSError:
            print("Failed to create the following folder: " + folder)
            raise


def create_directory(path):
    """
        Creates the directory specified by path if it does not already exist.
        Args:
            path: Path of the directory to create
    """
    try:
        if not os.path.exists(path):
            os.makedirs(path)
    except OSError:
        print("Failed to create the following folder: " + path)
        raise


def copy_config_file(config_file, output_directory):
    """
        Copies the config file in the output directory
        Args:
            config_file: path of the config file
            output_directory: path of the output directory where to save the config file
    """
    valid_path = validate_directory(output_directory)
    # create_directory(valid_path + "config/")

    try:
        os.system("cp " + config_file + " " + valid_path)
    except:
        print("Failed to copy " + config_file + " into " + valid_path)
        raise


def copy_file_in_directory(file, output_directory):
    """
        Creates the output directory if it does not exist
        Copies the file in the output directory
        Args:
            file: path of the file
            output_directory: path of the output directory where to copy the file
    """
    valid_path = validate_directory(output_directory)
    file_path = os.path.expandvars(file)
    create_directory(valid_path)

    try:
        subprocess.call(["cp", file_path, valid_path])
    except:
        print("Failed to copy " + file_path + " into " + valid_path)
        raise


def validate_directory(directory):
    """
        Adds a "/" to the end of provided path if needed
        Args:
            path to a directory
        Returns:
            the path to the provided path ending with a "/"
    """
    # current_working_directory = os.path.realpath(os.getcwd())
    # # no slash
    # if len(directory.split("/")) == 1:
    #     output_directory = current_working_directory + "/" + directory
    # # either path/directory or directory/
    # elif len(directory.split("/")) == 2:
    #     if directory.endswith("/"):
    #         output_directory = current_working_directory + "/" + directory
    #     else:
    #         output_directory = current_working_directory + "/" + directory.split("/")[-1]
    # # either path/path/directory or path/path/directory/
    # else:
    #     if directory.endswith("/"):
    #         output_directory = current_working_directory + "/" + directory.split("/")[-2]
    #     else:
    #         output_directory = current_working_directory + "/" + directory.split("/")[-1]
    output_directory = directory
    if not output_directory.endswith("/"):
        output_directory += "/"
    return output_directory


def get_full_path(path):
    """
        Removes the linux environment variable in the path
        Args:
            path containing a linux environment variable at the start (eg: $RESOURCE/path/to/resource)
        Returns:
            full path without variables or same path if the path doesn't start with an environment variable
    """

    if path.startswith("$"):
        path_list = path.split("/")
        environment_variable = path_list[0][1:]
        path_list[0] = os.environ.get(environment_variable)
        if path_list[0] is None:
            print(f"ERROR: the environment variable ${environment_variable} was not found for this path {path}")
            sys.exit(1)
        full_path = "/".join(path_list)

    else:
        full_path = path

    return full_path

def parse_sample_sheet(sample_sheet, verify_empty = False, single_end = False):
    """
        Creates a dictionary containing the sample ID associated with the list of fastq file.
        Args:
            sample_sheet: path to an illumina format sample sheet
            verify_empty: boolean
        Returns:
            a dictionary of fastq files mapped to the sample IDs
    """
    sample_fastq_dict = {}
    parsed_sample_sheet = SampleSheet(sample_sheet)
    # get the directory attribute or returns the directory where the sample sheet is located as a default value.
    directories = getattr(parsed_sample_sheet, "Directories", [os.path.dirname(os.path.abspath(sample_sheet))])
    samples = parsed_sample_sheet.samples
    for sample in samples:
        # check if sample id was duplicated
        if sample["Sample_ID"] in sample_fastq_dict:
            print("ERROR: {sample_id} is duplicated in sample sheet.".format(sample_id=sample["Sample_ID"]))
            exit(1)
        # check if sample id has an underscore that could lead to have another sample id be a substring of this id
        # if "_" in sample["Sample_ID"]:
        #     print("ERROR: {sample_id} has an underscore in its name.".format(sample_id=sample["Sample_ID"]))
        #     exit(1)

        sample_fastq_dict[sample["Sample_ID"]] = {"R1": "", "R2": ""}
        for directory in directories:
            directory_content = os.listdir(directory)
            for file in directory_content:
                if file.startswith(sample["Sample_ID"] + "_") and file.endswith(".fastq.gz"):
                    # check if file is empty
                    if verify_empty:
                        if is_gz_file_empty(directory + "/" + file):
                            print("ERROR: {file} is empty".format(file=directory + "/" + file))
                            # exit(1)
                    if "_R1" in file:
                        # check if R1 was already found
                        if sample_fastq_dict[sample["Sample_ID"]]["R1"] != "":
                            print(
                                "ERROR: Found two R1 fastq file for {sample_id}.\n R1: {R1_file_1} \n R1: {R1_file_2}".format(
                                    sample_id=sample["Sample_ID"],
                                    R1_file_1=sample_fastq_dict[sample["Sample_ID"]]["R1"],
                                    R1_file_2=directory + "/" + file
                                ))
                            exit(1)
                        sample_fastq_dict[sample["Sample_ID"]]["R1"] = directory + "/" + file
                    elif "_R2" in file:
                        # check if R2 was already found
                        if sample_fastq_dict[sample["Sample_ID"]]["R2"] != "":
                            print(
                                "ERROR: Found two R2 fastq file for {sample_id}.\n R2: {R2_file_1} \n R2: {R2_file_2}".format(
                                    sample_id=sample["Sample_ID"],
                                    R2_file_1=sample_fastq_dict[sample["Sample_ID"]]["R2"],
                                    R2_file_2=directory + "/" + file
                                ))
                            exit(1)
                        sample_fastq_dict[sample["Sample_ID"]]["R2"] = directory + "/" + file

            # add break to not loop all directories? should be 2-3 max usually
        if not single_end:
            if sample_fastq_dict[sample["Sample_ID"]]["R1"] == "" or sample_fastq_dict[sample["Sample_ID"]]["R2"] == "":
                message = f"""
    ERROR: Did not find two fastq files for sample: {sample["Sample_ID"]}. 
    R1: {sample_fastq_dict[sample["Sample_ID"]]["R1"]}
    R2: {sample_fastq_dict[sample["Sample_ID"]]["R2"]}
    """
                print(message)
                sample_fastq_dict.pop(sample["Sample_ID"])
                # exit()
            
    return sample_fastq_dict


def is_gz_file_empty(file_name):
    """
        Opens a zipped file and reads the first line to evaluate if the file is empty
        Args:
            file_name: path to the file test
        Returns:
            a boolean stating if the file is empty
    """
    # rb stands for read bytes
    with gzip.open(file_name, "rb") as f:
        data = f.read()
        is_empty = len(data) == 0
    return is_empty


def generate_index_reset():
    """
        Returns bash code to reset the PIDs list index to 0
    """
    code_segment = """
pids=
i=0
"""
    return code_segment


def generate_pid_segment():
    """
        Returns bash code that adds the last process ID into a list
    """
    code_segment = """
pid=$!
pids[$i]=${pid}
i=$((i+1))
"""
    return code_segment


def generate_wait_segment(step):
    """
        Returns a bash code segment to execute a loop to wait on for the current list of PIDs
    """
    code_segment = """
# wait for all pids
for pid in ${{pids[*]}}; do
    wait $pid
    if [ $? -eq 0 ];then
        echo "{step} done"
    else
        echo "ERROR $? {step} failed"
        exit
    fi
done
""".format(step=step)
    return code_segment

def generate_segmentation_marker(job_id=""):
    """
        Returns string for segmentation
    """

    if job_id != "":
        dash_string = "-" * int(60 - ((len(job_id) + 2)/2))
        marker = f"#{dash_string} {job_id} {dash_string}#\n"
    else:
        marker = "#------------------------------------------------------------------------------------------------------------------------#\n\n"

    return marker

def generate_full_command(job_id, command, bgzip="", gzip="", tabix="", verification_segment=True):
    """
        Returns full command with segmentation marker and echo job id. Options to add bgzip, tabix and segment verification steps.
        Args:
            job_id: String with the name of the job
            command: String of the command for the job
            bgzip: String of file to zip with bgzip
            gzip: String of file to zip with gzip
            tabix: String of file to index with tabix
            verification_segment: Boolean to add or not the verification segment step
        Return:
            String containing the full command
             
    """

    start_marker = generate_segmentation_marker(job_id)
    end_marker = generate_segmentation_marker()
    job_id_echo = f"echo \"JOB: {job_id}\"\n"
    
    full_command = start_marker + job_id_echo + command
    
    if verification_segment:
        full_command += generate_step_verification_segment(job_id)
    
    if bgzip != "":
        full_command += apply_bgzip(bgzip)

    if gzip != "":
        full_command += apply_gzip(gzip)
    
    if tabix != "":
        full_command += apply_tabix(tabix)

    full_command += end_marker

    return full_command

def generate_start_time():
    """
        Returns a bash code segment to save a variable with the start time
    """
    code_segment = """
start_time=$(date +%s)
"""
    return code_segment


def generate_end_time():
    """
        Returns a bash code segment to save a variable with the end time and print the elapsed time
    """
    code_segment = """
end_time=$(date +%s)

elapsed=$(( end_time - start_time ))

execution_time=$(date -d@$elapsed -u +%H:%M:%S)

echo "Execution time: $execution_time"
"""
    return code_segment


def create_script_file(script, file_name, output_directory):
    """
        Creates a .sh file with the script in it
        Args:
            script: string containing a bash script.
            file_name: name to give to the script file
            output_directory: path to the directory whgere to save the script
        Returns:
            the path to the created script file
    """
    valid_path = validate_directory(output_directory)
    file_path = valid_path + file_name
    script_file = open(file_path, "w+")
    script_file.write(script)
    script_file.close()
    return file_path


def remove_file_extension(file, extension, remove_path=True):
    """
        Returns the name of a file without the specified extension and without the ending dot
        Args:
            file: the file path to the extension from
            extension: the extension to remove
            remove_path: boolean, if set to false returns the whole path and not just the file name
        Returns:
            the name of a file without the specified extension
    """
    if extension and file.endswith(extension) and "/" not in file:
        if "." not in extension:
            return file[:-len("." + extension)]
        elif "." in extension:
            return file[:-len(extension)]
    if extension and file.endswith(extension) and "/" in file:
        if remove_path:
            file = file.split("/")[-1]
        if "." not in extension:
            return file[:-len("." + extension)]
        elif "." in extension:
            return file[:-len(extension)]
    else:
        print("Failed to remove extension " + extension + "from " + file)
        exit()


def execute_script_file(script_file, log_file):
    """
        Executes the bash script passed in parameters
        Args:
            script_file: the path to the script file
            log_file: the path to the log file
    """
    # give permission to execute the generated script
    subprocess.call(["chmod", "a+x", script_file])
    # write the command to run that includes the log file
    command = f"{script_file} 2>&1 | tee {log_file}"
    # execute the generated script
    subprocess.call(command, shell=True)


def set_permissions_to_output_directory(output_directory):
    """
        Returns a command to change to permissions to the analysis output directory
        Args:
            output_directory: the path to the output directory
    """
    command = f"""
chmod -R go+rw {output_directory}
"""
    return command


def generate_step_verification_segment(job_id):
    """
        Returns a code segment to determine the exit status of the last command
    """
    capital_case_job_id = job_id.upper()
    code_segment = f"""
err=$?
if [ $err -eq 0 ]; then
    echo "{capital_case_job_id} DONE"
else
    echo "ERROR $err {capital_case_job_id} FAILED"
    exit $err
fi
"""
    return code_segment


def apply_tabix(vcf_file):
    """
     Returns a command to apply tabix to a vcf file
    """
    command = f"""
tabix -f -p vcf {vcf_file}
"""
    return command


def apply_bgzip(input_file):
    """
     Returns a command to apply bgzip to a vcf file
    """
    command = f"""
bgzip -f {input_file} 
"""
    return command

def apply_gzip(input_file):
    """
     Returns a command to apply gzip to a file
    """
    command = f"""
gzip {input_file} 
"""
    return command


def extract_file_name_up_to_separator(file_path, separator=""):
    """
        Extracts the file name from the path and remove its extensions up to the separator if a separator is given
        Args:
            file_path: the path to the file from which to extract the name
            separator: Optional. The separator of in the file path
        Returns:
            a file name
    """
    file_name = os.path.split(file_path)[1]
    if separator != "":
        file_name = file_name.split(separator, 1)[0]
    return file_name


def create_report_tools(output_directory, path_to_bin):
    """
        Creates a report with the versions of the used softwares and libraries
        Assumes that all tools are installed in the folder following the same architecture
        Args:
            output_directory: Path to the directory where to save the report
            path_to_bin: Path to the folder were the tools are installed
    """
    tools_list = ["fastp", "bwa", "bcftools", "samtools", "gatk", "gatk3_old", "picard_folder", "vardict", "snpeff", "vt", "vcfanno"]
    tools_dict = {}
    validated_ouput_directory = validate_directory(output_directory)
    validated_bin_folder = validate_directory(path_to_bin)
    for tool in tools_list:
        tools_dict[tool] = []
        symbolic_link = os.readlink(os.path.expandvars(validated_bin_folder) + tool)
        tools_dict[tool].append(os.path.expandvars(validated_bin_folder) + symbolic_link)
        try:
            tools_dict[tool].append(symbolic_link.split("/")[1].split("_")[1])
        except IndexError:
            tools_dict[tool].append(symbolic_link.split("/")[1].split("-")[1])
    python_modules = os.popen("pip freeze")
    # print(python_modules)
    with open(validated_ouput_directory + "report_tools.tsv", "w") as f:
        f.write("Tool\tVersion\tPath\n")
        for tool in tools_dict:
            f.write("%s\t%s\t%s\n" % (tool, tools_dict[tool][1], tools_dict[tool][0]))
        f.write("\nPython libraries\nModule\tVersion\n")
        for module in python_modules:
            if "==" in module:
                module_version = module.split("==")
                f.write("%s\t%s" % (module_version[0], module_version[1]))


###############################
### FILTERING AND REPORTING ###
###############################

def generate_snpeff_coding_annotations():
    """
        Generates and returns the list of annotations considered as coding in snpeff
    """
    annotations = [
        "splice_acceptor_variant",
        "inframe_deletion",
        "conservative_inframe_deletion",
        "protein_protein_contact",
        "initiator_codon_variant",
        "rare_amino_acid_variant",
        "structural_interaction_variant",
        "disruptive_inframe_deletion",
        "disruptive_inframe_insertion",
        "chromosome",
        "stop_retained_variant",
        "exon_loss_variant",
        "splice_donor_variant",
        "feature_ablation",
        "stop_lost",
        "inframe_insertion",
        "conservative_inframe_insertion",
        "splice_region_variant",
        "5_prime_UTR_premature_start_codon_gain_variant",
        "gene_fusion",
        "stop_gained",
        "coding_sequence_variant",
        "inversion",
        "start_lost",
        "frameshift_variant",
        "bidirectional_gene_fusion",
        "missense_variant",
        "duplication",
        "rearranged_at_DNA_level",
    ]
    return annotations


def generate_snpeff_lof_annotations():
    """
        Generates and returns the list of annotations considered as loss of function in snpeff
    """
    annotations = [
        "stop_lost",
        "inversion",
        "splice_donor_variant",
        "splice_acceptor_variant",
        "start_lost",
        "frameshift_variant",
        "bidirectional_gene_fusion",
        "stop_gained",
        "exon_loss_variant",
        "duplication",
        "gene_fusion",
        "chromosome",
        "feature_ablation",
    ]
    return annotations


def get_alternate_allele_frequency(variant, index, variant_caller):
    """
        Calculates the alternate allele frequency of a variant for a sample.
        Args:
            variant: a variant object from cyvcf2.
            index: the index of the sample for which to calculate the alternate allele frequency.
            variant_caller: the vairnat caller used to generate the VCF
        Returns:
            the alternate allele frequency in float format. For example:
            0.2
    """
    depth = get_variant_depth(variant, index, variant_caller)
    alternate_allele_depth = get_variant_alternate_depth(variant, index, variant_caller)
    alternate_allele_freq = calculate_frequency(alternate_allele_depth, depth)
    return alternate_allele_freq


def calculate_frequency(numerator, denominator):
    """
        Calculates a frequency
        Args:
            numerator
            denominator
        Returns: the value of the division or 0
    """
    try:
        frequency = numerator / denominator
    except ZeroDivisionError:
        frequency = 0

    return frequency


def parse_single_annotation(annotation_info_dict: dict, allele_annotation: list) -> dict:
    """
    Parses 1 snpeff annotation field
    Args:
        annotation_info_dict: dictionary of SnpEff fields (length of 15)
        allele_annotation: list of a single SnpEff annotation (16 fields)
    Returns: dictionary with appended annotation information
    """
    annotation_info_dict["annotation"].append(allele_annotation[1])
    annotation_info_dict["impact"].append(allele_annotation[2])
    annotation_info_dict["gene_name"].append(allele_annotation[3])
    annotation_info_dict["gene_id"].append(allele_annotation[4])
    annotation_info_dict["feature_type"].append(allele_annotation[5])
    annotation_info_dict["transcript_biotype"].append(allele_annotation[7])
    if allele_annotation[8] == "":
        annotation_info_dict["rank"].append(f"{allele_annotation[6]}:no rank")
    else:
        annotation_info_dict["rank"].append(f"{allele_annotation[6]}:{allele_annotation[8]}")
    annotation_info_dict["hgvs_c"].append(f"{allele_annotation[6]}:{allele_annotation[9]}")
    if allele_annotation[10] != "":
        annotation_info_dict["hgvs_p"].append(f"{allele_annotation[6]}:{allele_annotation[10]}")
    else:
        annotation_info_dict["hgvs_p"].append(f"{allele_annotation[6]}:")
    annotation_info_dict["cdna_pos"].append(allele_annotation[11])
    annotation_info_dict["cds_pos"].append(allele_annotation[12])
    annotation_info_dict["protein_pos"].append(allele_annotation[13])
    annotation_info_dict["distance"].append(allele_annotation[14])
    annotation_info_dict["errors"].append(allele_annotation[15])

    return annotation_info_dict


def parse_snpeff_annotation(snpeff_annotation_field: str, transcripts_file: str = None) -> dict:
    """
        Parses the snpeff annotation field and returns a dictionary containing the individual annotations.
        Args:
            snpeff_annotation_field:
                string containing the information of snpeff annotation:
                "Allele | Annotation | Annotation_Impact | Gene_Name | Gene_ID | Feature_Type | Feature_ID | Transcript_BioType | Rank | HGVS.c | HGVS.p | cDNA.pos / cDNA.length | CDS.pos / CDS.length | AA.pos / AA.length | Distance | ERRORS / WARNINGS / INFO"
            transcripts_file: txt file with transcript ID and gene name
        Returns:
            record_annotation: a dictionary key represents annotation fields, values are the annotation according to the allele, precisely:
            a dictionary of the following individual annotations:
            {"annotation": <str> Type of change,
            "impact: <str> Annotation impact,
            "gene_name": <str> Gene name,
            "gene_id: <str> Gene ID,
            "feature_type": <str> Feature type,
            "feature_id": <str> Feature ID,
            "transcript_biotype": <str> Transcript biotype,
            "Rank": <str> Rank divided by total (on which exon or intron the change is),
            "hgvs_c": <str> cDNA change (HGVS nomenclature),
            "hgvs_p": <str> Amino acid change (HGVS nomenclature) if variant is coding,
            "cdna_pos": <str> Position in cDNA and trancript's cDNA length (one based),
            "cds_pos": <str> Position and number of coding bases (one based includes START and STOP codons),
            "protein_pos": <str> Position and number of AA (one based, including START, but not STOP),
            "distance": <str> All items in this field are options, so the field could be empty,
            "errors": <str> Add errors, warnings or informative message that can affect annotation accuracy
            "lof": <bool> Loss of function,
            "coding": <bool>
            }
    """
    # read transcripts file (*refseq_name_to_name2.txt)
    transcripts_df = None
    if transcripts_file is not None:
        transcripts_df = pd.read_csv(transcripts_file, sep="\t", names=["TRANSCRIPT", "GENE"])

    record_annotation = {}
    coding_variants_snpeff = generate_snpeff_coding_annotations()
    loss_of_function_snpeff = generate_snpeff_lof_annotations()

    # the case where snpeff annotation is None
    if snpeff_annotation_field is None:
        record_annotation["annotation"] = ""
        record_annotation["impact"] = ""
        record_annotation["gene_name"] = ""
        record_annotation["gene_id"] = ""
        record_annotation["feature_type"] = ""
        record_annotation["feature_id"] = ""
        record_annotation["transcript_biotype"] = ""
        record_annotation["rank"] = ""
        record_annotation["hgvs_c"] = ""
        record_annotation["hgvs_p"] = ""
        record_annotation["cdna_pos"] = ""
        record_annotation["cds_pos"] = ""
        record_annotation["protein_pos"] = ""
        record_annotation["distance"] = ""
        record_annotation["errors"] = ""
        record_annotation["lof"] = False
        record_annotation["coding"] = False

    else:
        snpeff_alleles = snpeff_annotation_field.split(",")
        # create dict of empty lists
        keys = ["annotation", "impact", "gene_name", "gene_id", "feature_type", "transcript_biotype", "rank", "hgvs_c", "hgvs_p",
                "cdna_pos", "cds_pos", "protein_pos", "distance", "errors"]
        annotation_info_dict = dict(zip(keys, ([] for _ in keys)))

        for allele in snpeff_alleles:
            allele = allele.split("|")

            # check if the feature ID of this allele is in the transcripts file of transcripts that are reported in the clinical lab
            if transcripts_df is not None:
                if transcripts_df["TRANSCRIPT"].str.contains(allele[6]).any():  # get all the records matching the query
                    annotation_info_dict = parse_single_annotation(annotation_info_dict=annotation_info_dict, allele_annotation=allele)
            else:
                # no transcripts file provided so we parse every annotation
                annotation_info_dict = parse_single_annotation(annotation_info_dict=annotation_info_dict, allele_annotation=allele)

        # case when feature id in snpeff fields of 1 variant doesn't match any reported transcripts so the annotation is null
        if not any(annotation_info_dict.values()):  # check that all values are empty lists
            # append annotation results to output dictionary
            record_annotation["annotation"] = ""
            record_annotation["impact"] = ""
            record_annotation["gene_name"] = ""
            record_annotation["gene_id"] = ""
            record_annotation["feature_type"] = ""
            record_annotation["feature_id"] = ""
            record_annotation["transcript_biotype"] = ""
            record_annotation["rank"] = ""
            record_annotation["hgvs_c"] = ""
            record_annotation["hgvs_p"] = ""
            record_annotation["cdna_pos"] = ""
            record_annotation["cds_pos"] = ""
            record_annotation["protein_pos"] = ""
            record_annotation["distance"] = ""
            record_annotation["errors"] = ""
            record_annotation["lof"] = False
            record_annotation["coding"] = False

        else:
            # sort ranks (list of transcript id and exon/intron rank. i.e. NM_000123.1:1/10)
            rank_split = [x.split(":") for x in set(annotation_info_dict["rank"])]
            rank_split_str = sorted(rank_split, key=lambda x: (x[0], x[1]))  # sort by transcript id first then by rank

            # transform list of lists to dict to better subset ranks for same transcript
            transcript_rank_dict = {}
            for sublist in rank_split_str:
                key = sublist[0]
                if key in transcript_rank_dict:
                    transcript_rank_dict[key].append(sublist[1])
                else:
                    transcript_rank_dict[key] = [sublist[1]]

            # select the ranks by transcript
            ranks_output = []
            for rank in transcript_rank_dict.values():
                ranks_output.append(";".join(rank))

            # append annotation results to output dictionary
            record_annotation["annotation"] = ",".join(set(annotation_info_dict["annotation"]))
            record_annotation["impact"] = ",".join(set(annotation_info_dict["impact"]))
            record_annotation["gene_name"] = ",".join(set(annotation_info_dict["gene_name"]))
            record_annotation["gene_id"] = ",".join(set(annotation_info_dict["gene_id"]))
            record_annotation["feature_type"] = ",".join(set(annotation_info_dict["feature_type"]))
            record_annotation["feature_id"] = ",".join(transcript_rank_dict.keys())
            record_annotation["transcript_biotype"] = ",".join(set(annotation_info_dict["transcript_biotype"]))
            record_annotation["rank"] = ",".join(ranks_output)
            record_annotation["hgvs_c"] = ",".join(set(annotation_info_dict["hgvs_c"]))
            record_annotation["hgvs_p"] = ",".join(set(annotation_info_dict["hgvs_p"]))
            record_annotation["cdna_pos"] = ",".join(set(annotation_info_dict["cdna_pos"]))
            record_annotation["cds_pos"] = ",".join(set(annotation_info_dict["cds_pos"]))
            record_annotation["protein_pos"] = ",".join(set(annotation_info_dict["protein_pos"]))
            record_annotation["distance"] = ",".join(set(annotation_info_dict["distance"]))
            record_annotation["errors"] = ",".join(set(annotation_info_dict["errors"]))
            record_annotation["lof"] = False
            record_annotation["coding"] = False

            annotations = record_annotation["annotation"].split(",")
            for annotation in annotations:
                if "&" in annotation:  # Fix problem with double annotations (e.g. "missense_variant&splice_region_variant")
                    annotation = annotation.split("&")
                    for sub_annotation in annotation:
                        if not record_annotation["coding"]:
                            if sub_annotation in coding_variants_snpeff:
                                record_annotation["coding"] = True
                        if not record_annotation["lof"]:
                            if sub_annotation in loss_of_function_snpeff:
                                record_annotation["lof"] = True
                    if record_annotation["coding"] and record_annotation["lof"]:
                        break
                else:
                    if not record_annotation["coding"]:
                        if annotation in coding_variants_snpeff:
                            record_annotation["coding"] = True
                    if not record_annotation["lof"]:
                        if annotation in loss_of_function_snpeff:
                            record_annotation["lof"] = True
                    if record_annotation["coding"] and record_annotation["lof"]:
                        break

    return record_annotation


def calculate_gnomad_allele_frequency(variant):
    """
        Calculates the allele frequency in gnomad for the specified variant.
        Args:
            variant: <object> cyvcf2.cyvcf2.Variant object
        Returns:
            A float of the gnomad allele frequency of this variant. For example:
            0.34
    """
    try:
        gnomad_exome_alternate_allele_count = variant.INFO.get("GNDEXAC")
        gnomad_genome_alternate_allele_count = variant.INFO.get("GNDWGAC")
        gnomad_exome_total_allele_count = variant.INFO.get("GNDEXAN")
        gnomad_genome_total_allele_count = variant.INFO.get("GNDWGAN")

        if gnomad_exome_alternate_allele_count is None:
            gnomad_exome_alternate_allele_count = 0

        if gnomad_genome_alternate_allele_count is None:
            gnomad_genome_alternate_allele_count = 0

        if gnomad_exome_total_allele_count is None:
            gnomad_exome_total_allele_count = 0

        if gnomad_genome_total_allele_count is None:
            gnomad_genome_total_allele_count = 0

        return (gnomad_exome_alternate_allele_count + gnomad_genome_alternate_allele_count) / (
                gnomad_exome_total_allele_count + gnomad_genome_total_allele_count)

    except ZeroDivisionError:
        return 0


def get_variant_depth(variant, index, variant_caller):
    """
        Returns the variant depth for a specific sample
        Args:
            variant: a variant object from cyvcf2.
            index: the index of the sample
            variant_caller: the vairnat caller used to generate the VCF
        Returns: the variant depth
    """
    if variant_caller.upper() == "GATK" or variant_caller.upper() == "VARDICT":
        variant_depth = variant.gt_depths[index]
        # cyvcf2 replaces missing values by -1 when using gt_depths
        if variant_depth is None or variant_depth == -1:
            variant_depth = 0
    elif variant_caller.upper() == "VARSCAN":
        variant_depth = variant.format("DP")[index][0]
        # cyvcf2 replaces missing values by -2147483648 when using format. See Github issue: https://github.com/brentp/cyvcf2/issues/172
        if variant_depth is None or variant_depth == -2147483648:
            variant_depth = 0
    else:
        print("ERROR: the variant caller is not recognized when validating variants")
        exit()

    return variant_depth


def get_variant_alternate_depth(variant, index, variant_caller):
    """
        Returns the variant alternate allele depth for a specific sample
        Args:
            variant: a variant object from cyvcf2.
            index: the index of the sample
            variant_caller: the variant caller used to generate the VCF
        Returns: the variant alternate allele depth
    """
    if variant_caller.upper() == "GATK" or variant_caller.upper() == "VARDICT":
        variant_alternate_depth = variant.gt_alt_depths[index]
        # cyvcf2 replaces missing values by -1 when using variant_alternate_depth
        if variant_alternate_depth is None or variant_alternate_depth == -1:
            variant_alternate_depth = 0
    elif variant_caller.upper() == "VARSCAN":
        variant_alternate_depth = variant.format("AD")[index][0]
        # cyvcf2 replaces missing values by -2147483648 when using format. See Github issue: https://github.com/brentp/cyvcf2/issues/172
        if variant_alternate_depth is None or variant_alternate_depth == -2147483648:
            variant_alternate_depth = 0
    else:
        print("ERROR: the variant caller is not recognized when validating variants")
        exit()

    return variant_alternate_depth


def validate_variant(variant, index, depth_min, alternate_depth_min, alternate_frequency_min, variant_caller):
    """
        Validates if a variant passes the filtering options
        Args:
            variant: the cyvcf2 Variant object to test
            index: the index of the sample to get the variant data from
            depth_min: the minimum depth to keep the variant
            alternate_depth_min: the minimum depth of the alternate allele to keep the variant
            alternate_frequency_min: the minimum alternate allele frequency to keep the variant
            variant_caller: variant caller tool
        Returns:
            a boolean stating if the variant passes the filter
    """
    validation = False
    undetermined_genotypes = ["0/.", "./.", "./0", "0|.", ".|.", ".|0"]

    variant_depth = get_variant_depth(variant, index, variant_caller)
    variant_alternate_depth = get_variant_alternate_depth(variant, index, variant_caller)
    variant_alternate_frequency = get_alternate_allele_frequency(variant, index, variant_caller)

    if (variant_depth >= int(depth_min) and
            variant_alternate_depth >= int(alternate_depth_min) and
            variant_alternate_frequency >= float(alternate_frequency_min)):
        validation = True

    if (variant_caller.upper() == "VARSCAN" and
            not (str([Genotype(genotype) for genotype in variant.genotypes][index]) not in undetermined_genotypes)):
        validation = False
        # the true statement for reference
        # if str([Genotype(genotype) for genotype in variant.genotypes][index]) not in undetermined_genotypes:

    return validation


def get_gnomad_format(variant):
    """
        Builds the gnomad notation of a variant.
        Assumes that the VCF was decomposed -> only one allele per variant.
        Args:
            cyvcf2 Variant object
        Returns:
            the gnomad notation of the variant
    """
    try:
        if len(variant.ALT) > 1:
            # # this could be applied, but it should not be necessary because the vcf was decomposed
            # for alternate in variant.ALT:
            #     alleles_alternate += alternate + ","
            #     alleles_alternate = alleles_alternate[:-1]
            print(
                "Please make sure that the VCF has been decomposed. More then one alternative alleles were found in the following variant:")
            print(variant)
            sys.exit()
        elif len(variant.ALT) == 0:
            gnomad_format = str(variant.CHROM) + "-" + str(variant.start + 1) + "-" + variant.REF + "-."
        else:
            gnomad_format = str(variant.CHROM) + "-" + str(variant.start + 1) + "-" + variant.REF + "-" + variant.ALT[0]
        if gnomad_format.startswith("chr"):
            gnomad_format = re.sub("^%s" % "chr", "", gnomad_format)
    except:
        print(f"""Failed to get the gnomad notation of the following variant:
CHROM: {str(variant.CHROM)}
START: {str(variant.start + 1)}
REF: {variant.REF}
ALT: {variant.ALT}""")
        sys.exit()

    return gnomad_format


# TODO Convert all to only parse_zeroes_for_report if needed
def parse_percentage_for_report(dataframe, percentage_columns):
    """
        Convert percentage column to type string and replaces all strings 100.0 to the string 100 and 0.0 to 0 to write a report
        Args:
            dataframe: a pandas DataFrame containing the columns to parse
            percentage_columns: the percentage columns to parse
    """
    parsed_dataframe = pd.DataFrame()
    for column in dataframe:
        if column in percentage_columns:
            parsed_dataframe[column] = dataframe[column].astype(str).str.replace("100.0", "100",
                                                                                 regex=False).str.replace("\.0+$", "",
                                                                                                          regex=True)
        else:
            parsed_dataframe[column] = dataframe[column]
    parsed_dataframe.index = dataframe.index
    return parsed_dataframe


def parse_zeroes_for_report(dataframe):
    """
        Convert all float columns to type string and replaces all strings ".0" to an empty string to write a report
        Args:
            dataframe: a pandas DataFrame containing the columns to parse
    """
    parsed_dataframe = pd.DataFrame()
    for column in dataframe:
        if dataframe[column].dtype == float:
            # This regex avoids removing the 0.002. The regex matches a .0 at the end of a string
            parsed_dataframe[column] = dataframe[column].astype(str).str.replace("\.0+$", "", regex=True).replace("nan",
                                                                                                                  "",
                                                                                                                  regex=False)
        else:
            parsed_dataframe[column] = dataframe[column]
    parsed_dataframe.index = dataframe.index
    return parsed_dataframe


def generate_excel_file(input_file_path, separator="\t"):
    """
    Creates an excel file from a tsv file
    Args:
        input_file_path: The path to a tsv or csv file
        separator: string that seperate the columns in the tsv or csv file, default: "\t"
    Note:
        - Excel cell can have a maximum of 32767 characters, when the data is too big for a cell it is truncated.
        see https://support.microsoft.com/en-us/office/excel-specifications-and-limits-1672b34d-7043-467e-8e27-269d656771c3#ID0EDBD=Newer_versions
    """

    if separator == "\t" and input_file_path.endswith(".tsv"):
        output_file = os.path.dirname(input_file_path) + "/" + remove_file_extension(input_file_path, ".tsv") + ".xlsx"
    elif separator == "," and input_file_path.endswith(".csv"):
        output_file = os.path.dirname(input_file_path) + "/" + remove_file_extension(input_file_path, ".csv") + ".xlsx"
    else:
        print(
            "ERROR: Unsupported seperator: \"{separator}\" for generate_excel_file, only tabs (tsv) and comma (csv) are supported")
        output_file = None

    if output_file is not None:
        # Create an XlsxWriter workbook object and add a worksheet.
        workbook = Workbook(output_file)
        worksheet = workbook.add_worksheet()

        # Reads the csv or tsv file and writes in the excel file line by line
        with open(input_file_path) as input_file:
            csv_reader = csv.reader(input_file, delimiter=separator)

            for row, data in enumerate(csv_reader):
                for column_number in range(len(data)):
                    worksheet.write(row, column_number, data[column_number])

        # Close the excel file
        workbook.close()


def apply_generate_excel_file(ngstk_path, input_file_path, separator="\t"):
    """
        Creates the command to create an excel file from a tsv or csv file
        Args:
            ngstk_path: The path to the ngstk folder
            input_file_path: The path to a tsv or csv file
            separator: string that seperate the columns in the tsv or csv file, default: "\t"

    """
    command = """
python -c "import sys;
sys.path.append(\\"{ngstk_path}\\");
import common_tools
common_tools.generate_excel_file(\\"{input_file_path}\\", \\"{separator}\\")"
""".format(
        ngstk_path=ngstk_path,
        input_file_path=input_file_path,
        separator=separator,
    )

    return command


def find_dict_file_from_reference(reference_genome):
    """
        Returns the reference with a .dict extension instead of the expected .fa
        Args:
            reference_genome: Path to a reference fasta file
        Aussumes that the .dict file is in the same location as the fasta file.
    """
    reference_dict = remove_file_extension(reference_genome, ".fa", False) + ".dict"
    if not verify_file_exists(reference_dict):
        print(
            "ERROR: Expected to find {reference_dict}. There must be a .dict file in the same folder and with the same name as the .fa reference")
        exit()
    return reference_dict


def format_number_for_report(number, decimals=2):
    """
        Rounds the number if it's provided or returns None.
        Args:
            number: <float> or empty or nan
            decimals: the number of decimals to keep. Default is 2.
        Returns:
            "None" in case there is no number provided
            the rounded number to the input decimals
            0 if near zero
    """
    if number == "None" or number is None:
        return "None"
    else:
        # Making sure the number is a float
        number = float(number)
        if abs(number) >= 1:
            return round(number, decimals)
        elif np.isclose(number, 0, rtol=1e-08, atol=1e-08):
            return 0
        else:
            if decimals == 2:
                round_value = -int(math.floor(math.log10(abs(number))) - 1)
            else:
                round_value = -int(math.floor(math.log10(abs(number))) - 1) + decimals - 2
            return round(number, round_value)

def format_del_and_dup(hgvs_string: str) -> str:
    """
        Formats the HGVS string
        :param hgvs_string: HGVS nomenclature
        :return: formatted HGVS nomenclature
    """
    if "del" in hgvs_string:
        # this should have split the string in 2
        del_segments = hgvs_string.split("del")
        # if there is also an insertion, we want to keep that information
        if "ins" in del_segments[1]:
            # this should have split the string in 2
            ins_segments = del_segments[1].split("ins")
            formatted_hgvs_string = del_segments[0] + "delins" + ins_segments[1]
        else:
            formatted_hgvs_string = del_segments[0] + "del"
    elif "dup" in hgvs_string:
        # this should have split the string in 2
        dup_segments = hgvs_string.split("dup")
        # there are no dupins
        if "ins" in dup_segments[1]:
            sys.exit(f"A dupins was found: {hgvs_string}")
        else:
            formatted_hgvs_string = dup_segments[0] + "dup"
    else:
        formatted_hgvs_string = hgvs_string
    return formatted_hgvs_string


def correct_hgvs_transcripts(transcripts_file: str, report_data: pd.DataFrame) -> pd.DataFrame:
    """
        Filters the annotation based on the transcripts we're reporting. These transcripts can be found in the
        $RESOURCES files (refseq_name_to_name2.txt) where the transcript ID and gene names can be found.
        We create new columns in the report where only reported transcripts and their HGVS nomenclature are found.
        Value is set to None is no transcript is found in the TSV provided.

        Columns "SNPEFF_GENE_NAME", "SNPEFF_HGVS_C" and "SNPEFF_HGVS_P" are required for a report for the CNV pipeline or
        "VARIANT_GENE_NAME", "VARIANT_HGVS_C", and "VARIANT_HGVS_P" are required for a report for the Ampliseq or Enrichment pipeline.

        :param transcripts_file: txt file with transcript ID and gene name
        :param report_data: dataframe of the report we want to filter.
        :return: dataframe of the corrected report with 2 new columns HGVS_C and HGVS_P
    """
    # read transcripts file (refseq_name_to_name2.txt)
    transcripts_df = pd.read_csv(transcripts_file, sep="\t", names=["TRANSCRIPT", "GENE"])
    transcripts_dict = {}
    for row in transcripts_df.itertuples():
        if row.GENE not in transcripts_dict:
            transcripts_dict[row.GENE] = []
        transcripts_dict[row.GENE].append(row.TRANSCRIPT)

    # create empty dataframes with columns names and the HGVS columns
    validated_report_df = pd.DataFrame(columns=report_data.columns.tolist() + ["HGVS_C", "HGVS_P"])

    for row in report_data.itertuples():
        # append row to dataframe
        
        validated_report_df = pd.concat([validated_report_df, report_data.loc[[row.Index]]], ignore_index=True)
        # change dtypes explicitly of columns because it puts float dtype by default
        validated_report_df = validated_report_df.astype(dtype={"HGVS_C": "object", "HGVS_P": "object"})
        if ("VARIANT_HGVS_C" in report_data.columns and row.VARIANT_HGVS_C == "" or
                "SNPEFF_HGVS_C" in report_data.columns and row.SNPEFF_HGVS_C == ""):
            # put None in the columns and skip to next row because no hgvs was found
            validated_report_df.at[row.Index, "HGVS_C"] = "None"
            validated_report_df.at[row.Index, "HGVS_P"] = "None"
            continue

        # check if gene in current row is found in the file of reported genes and transcripts
        genes = row.SNPEFF_GENE_NAME.split(",") if "SNPEFF_GENE_NAME" in report_data.columns else row.VARIANT_GENE_NAME.split(",")
        gene_not_found = True
        input_transcripts = []
        for gene in genes:
            try:
                input_transcripts = transcripts_dict[gene]
                gene_not_found = False
                break
            except:
                continue

        if gene_not_found:
            column_to_use = row.SNPEFF_GENE_NAME if "SNPEFF_GENE_NAME" in report_data.columns else row.VARIANT_GENE_NAME
            print(
                f"WARNING: None of the genes in the value {column_to_use} are in the transcripts file provided.")
            validated_report_df.at[row.Index, "HGVS_C"] = "None"
            validated_report_df.at[row.Index, "HGVS_P"] = "None"
            continue

        transcripts = row.SNPEFF_HGVS_C.split(",") if "SNPEFF_HGVS_C" in report_data.columns else row.VARIANT_HGVS_C.split(",")

        # it is slower to loop twice through the transcripts but a little more comprehensible
        valid_transcripts = []
        for transcript in transcripts:
            transcript_value = transcript.split(":")[0]
            transcript_number = transcript_value.split(".")[0]
            cdna_change = transcript.split(":")[1]
            for input_transcript in input_transcripts:
                if transcript_value == input_transcript:
                    valid_transcripts.append(format_del_and_dup(hgvs_string=transcript))
                    break
                elif input_transcript.split(".")[0] == transcript_number:
                    new_transcript = ":".join([input_transcript, cdna_change])
                    valid_transcripts.append(format_del_and_dup(hgvs_string=new_transcript))
                    break
        # append None if empty list
        if not valid_transcripts:
            valid_transcripts.append("None")

        validated_report_df.at[row.Index, "HGVS_C"] = ",".join(valid_transcripts)

        # no del/dup format correction needed for the amino acid changes
        valid_transcripts_p = []
        if ("VARIANT_HGVS_P" in report_data.columns and row.VARIANT_HGVS_P == "" or
                "SNPEFF_HGVS_P" in report_data.columns and row.SNPEFF_HGVS_P == ""):
            validated_report_df.at[row.Index, "HGVS_P"] = "None"
        else:
            transcripts_p = row.SNPEFF_HGVS_P.split(",") if "SNPEFF_HGVS_P" in report_data.columns else row.VARIANT_HGVS_P.split(",")
            for transcript in transcripts_p:
                transcript_value = transcript.split(":")[0]
                transcript_number = transcript_value.split(".")[0]
                cdna_change = transcript.split(":")[1]
                for input_transcript in input_transcripts:
                    if transcript_value == input_transcript:
                        valid_transcripts_p.append(transcript)
                        break
                    elif input_transcript.split(".")[0] == transcript_number:
                        new_transcript = ":".join([input_transcript, cdna_change])
                        valid_transcripts_p.append(new_transcript)
                        break
        # append None if empty list
        if not valid_transcripts_p:
            valid_transcripts_p.append("None")

        validated_report_df.at[row.Index, "HGVS_P"] = ",".join(valid_transcripts_p)

    # reorder HGVS columns if the report is for ampliseq or enrichment pipeline
    if "VARIANT_HGVS_P" in report_data.columns:
        hgvs_column = validated_report_df.pop("HGVS_C")
        index = validated_report_df.columns.tolist().index("VARIANT_HGVS_C") + 2
        validated_report_df.insert(index, "HGVS_C", hgvs_column)

        hgvs_column = validated_report_df.pop("HGVS_P")
        index = validated_report_df.columns.tolist().index("HGVS_C") + 1
        validated_report_df.insert(index, "HGVS_P", hgvs_column)

    return validated_report_df

def get_number_of_samples(vcf_file):
    parsed_vcf = VCF(vcf_file)
    number_of_samples = len(parsed_vcf.samples)
    return number_of_samples

def resolve_path_with_dollar_sign(path):
    parts = path.split("/")
    resolved_path_parts = []
    for part in parts:
        if "$" in part:
            variable_name = part.split("$")[1]
            variable_value = os.environ.get(variable_name)
            resolved_path_parts.append(variable_value)
        else:
            resolved_path_parts.append(part)
    resolved_path = "/".join(resolved_path_parts)
    return resolved_path
