import pandas as pd
import numpy as np
import os
import cyvcf2
import common_tools
from common_tools import Genotype
import sys
import configparser

class ReportVariants:
    """
    This class contains functions to generate report files in .tsv format from a VCF file.
    
    Attributes:
        vcf_file: a string indicating the path of a VCF file.
        variant_caller: a string indicating the variant caller.
        reporting_file: a string indicating the path to the CSV file describing the reports to generate.
        output_directory: a string indicating the path of the output directory.
        denovo_parameters: a dict containing the denovo analysis parameters: 
            minimum depth of coverage, proband alternate depth, parent alternate depth, proband alternate frequency, parent alternate frequency
        report_parameters: a dictionary where the keys are the names of the report columns and the values are booleans confirming if the column must be included.
        hgvs_transcripts: a string to the path of a file with gene names and transcripts that are reported
        coding_variants_snpeff: a list of annotations considered as coding in snpeff
        loss_of_function_snpeff: a list of annotations considered as loss of function in snpeff
        reports_list: a list of Sample ID in reports files
        ngstk_path: a string of the path to the ngstk code


    Methods:
        parse_reference_resources(report_columns)
            Parse the columns which are resources such as canonical RefSeq, GnomAD z scores, hgnc and HPO
        validate_gene_list(gene_list_files)
            Verifies if there is a gene lists were provided by the client and if the genes are in HUGO Gene Nomenclature Committee (HGNC) symbols. 
            Generates a file containing the genes that are not validated from HGNC and returns a dictionary file names pertaining to the gene list.
        add_gene_list(gene_names)
            Verifies if the gene names are in the gene lists provided by the client and returns name of lists file(s) pertaining to the gene(s).
        get_canonical_transcript(hgvs_value, gene_names)
            Retrieves the canonical transcript corresponding to the gene(s).
        get_refseq_select_transcript(hgvs_value, refseq_select_transcript_ids)
            Retrieves the hgvs value corresponding to the refseq select transcript
        return_string_format_boolean(field)
            Returns "True" or "False" in string format.
        hpo(gene_name)
            Parses the Human Phenotype Ontology data and returns diseases associated to a gene and their inheritance mode.
        get_gnomad_z_score(gene_name)
            Parses the gnomad metric z scores file, appends missense, synonymous and loss of function z scores then it returns three z scores.
        get_reports_list(reporting_file)
            Parses the reporting file and returns the lists of samples IDs in each report.
        gene_count(report_data, gene_count_dict)
            Generates the gene count column for the report_data dictionary.
        get_sample_recessive_variants(proband_id, mother_id="", father_id="", depth_min=10, alternate_depth_min=3, alternate_frequency_min=0.05)
            Generates a list of recessive variants for a specific sample(proband_id).
        get_phasing_recessive_variants(mother_index, father_index, variant_1, variant_2)
            Determines the phasing of 2 variants in genes where the heterozygous count equals 2.
        denovo(proband_id, father_id, mother_id, depth_coverage_proband, depth_coverage_parents, alternate_depth_proband, alternate_depth_parents, alternate_frequency_proband, alternate_frequency_parents)
            Searches for denovo variants in the proband sample in the vcf file.
        get_zygosity(alternate_allele_freq)
            Determines the zygosity of a variant according to the alternate allele frequency.
        get_variant_cluster(sample_id, variants_list)
            Generates a list of variants that are clustered.
        generate_sample_variants_dictionary(depth_min, alternate_depth_min, alternate_frequency_min)
            Generates a dictionary containing sample IDs as keys with relevant variants and variant clusters.
        get_variant_genes_snpeff(variant)
            Gets a list of genes from the snpeff annotation of a variant
        get_report_output_file_name(samples_to_report)
            Name the report file depending on the number of samples to report
        create_report_summary(output_directory)
            Creates a file to save the reporting information
        format_meta_scores_columns(meta_score)
            Format the META scores from the VCF file for the report.
        generate_report(sample_variants_dictionary, report_parameters)
            Generates the reports in TSV format.
        generate_report_variants_command(config_file)
            Generates a bash command to execute run_report_variants.py
        apply_labkey_filter_report_variants(labkey_filter_list, report_variants_file)
            Generates a bash command to apply the filter_report_variants.py script located in the labkey folder
        generate_report_mei(vcf_file, gnomad_data_string)
            Generates a TSV report from the merged MEI VCF
        generate_report_mei_command(vcf_file, config_file)
            Generates a bash command to execute the generate_report_mei function
    """

    # constructor for class
    def __init__(self, config_file, vcf_file, variant_caller, output_directory, annotation_tool, reporting_file="", report_parameters=""):
        self.config = configparser.ConfigParser()
        self.config.read(config_file)
        self.vcf_file = vcf_file
        self.variant_caller = variant_caller
        self.reporting_file = reporting_file
        self.output_directory = output_directory
        self.annotation_tool = annotation_tool
        self.report_parameters = report_parameters
        self.hgvs_transcripts = False #hgvs_transcripts
        self.coding_variants_snpeff = common_tools.generate_snpeff_coding_annotations()
        self.loss_of_function_snpeff = common_tools.generate_snpeff_lof_annotations()
        # Generate a lis of the reports to make
        self.reports_list = self.get_reports_list(self.reporting_file)
        # sys.argv[0] refers to the executed pipeline script
        script_path = os.path.dirname(os.path.abspath(sys.argv[0]))
        # the -1 removes the last element corresponding to the analysis folder
        self.ngstk_path = os.sep.join(script_path.split(os.sep)[:-1])

    def parse_reference_resources(self, report_columns):
        """
        Parse the columns which are resources such as canonical RefSeq, GnomAD z scores, hgnc and HPO
        Args:
            report_columns: list of column names in report file
        Returns: -
        """

        if "CANONICAL_HGVS_C" in report_columns or "CANONICAL_HGVS_P" in report_columns:
            # Parse canonical refseq file
            canonical_refseq = common_tools.get_full_path(self.config.get("REPORT_VARIANTS", "CANONICAL_REFSEQ"))
            canonical_refseq_file = pd.read_csv(canonical_refseq, sep="\t", names=["chr", "start", "end", "clusterId", "canonical.transcript", "gene", "refseq", "Ensembl"])
            
            # gene_canonical_transcript_dict: a dictionary key represents gene name, value represents canonical transcript according to the gene. For example:
            # {"MASP2": "NM_006610.3:c.395C>T"}
            self.gene_canonical_transcript_dict = pd.Series(canonical_refseq_file.refseq.values, index=canonical_refseq_file.gene).to_dict()

        if "GNOMAD_SYN_Z" in report_columns or "GNOMAD_MIS_Z" in report_columns or "GNOMAD_LOF_Z" in report_columns:
            # Parse gnomad metrics z scores
            gnomad_metrics = common_tools.get_full_path(self.config.get("REPORT_VARIANTS", "GNOMAD_METRICS"))
            self.gnomad_file = pd.read_csv(gnomad_metrics, sep="\t")
        
        if "GENE_LIST" in report_columns:

            gene_list_files = self.config.get("REPORT_VARIANTS", "GENE_PANEL")
            # Parse provided gene list directories and hgnc file
            provided_gene_lists_dict_temporary = {}
            if (gene_list_files):
                self.report_parameters["GENE_LIST"] = True
                hgnc_genes = common_tools.get_full_path(self.config.get("REPORT_VARIANTS", "HGNC_GENES"))
                self.hgnc_file = pd.read_csv(hgnc_genes, sep="\t", dtype="object")
                provided_gene_lists_dict_temporary = self.validate_gene_list(gene_list_files)
            self.provided_gene_lists_dict = provided_gene_lists_dict_temporary

        # self.denovo_depth_coverage_proband, self.denovo_depth_coverage_parents, self.denovo_alternate_depth_proband,
        # self.denovo_alternate_depth_parents, self.denovo_alternate_frequency_proband, self.denovo_alternate_frequency_parents = denovo_parameters

        if "HPO_DISEASE" in report_columns:
            # Parse HPO files
            hpo_genes = common_tools.get_full_path(self.config.get("REPORT_VARIANTS", "HPO_GENES"))
            gene_diseases_data = pd.read_csv(hpo_genes, sep="\t")

            hpo_diseases = common_tools.get_full_path(self.config.get("REPORT_VARIANTS", "HPO_DISEASES"))
            disease_terms_data = pd.read_csv(hpo_diseases, sep="\t")

            hpo_disease_inheritances = common_tools.get_full_path(self.config.get("REPORT_VARIANTS", "HPO_DISEASE_INHERITANCES"))
            disease_inheritance_data = pd.read_csv(hpo_disease_inheritances)
            inheritance_mode = list(disease_inheritance_data)

            hpo_gene_diseases = common_tools.get_full_path(self.config.get("REPORT_VARIANTS", "HPO_GENES_DISEASES"))
            gene_diseases_parsed_data = pd.read_csv(hpo_gene_diseases, sep="\t", names=["GENES_SYMBOL", "DISEASE"])
            gene_diseases_parsed_data.set_index("GENES_SYMBOL", inplace=True)
            self.gene_diseases_parsed = gene_diseases_parsed_data.to_dict()["DISEASE"]

            # gene_diseases key represents gene id, value is a list of disease ids. For example:
            # {"ABCA5":  ["OMIM:135400","ORPHA:2026"]}
            gene_diseases_data.set_index("GENE_SYMBOL", inplace=True)
            self.gene_diseases = gene_diseases_data.to_dict()["DISEASE_IDS"]

            # disease_to_inheritance key represents disease id, value is mode of inheritance. For example:
            # {"ORPHA:199318": "autosomal_dominant_inheritance"}
            disease_to_inheritance_temporary = {}
            for inheritance_term in inheritance_mode:
                for disease_id in disease_inheritance_data[inheritance_term]:
                    if disease_id != "nan":
                        if disease_id not in disease_to_inheritance_temporary.keys():
                            disease_to_inheritance_temporary[disease_id] = [inheritance_term]
                        elif disease_id in disease_to_inheritance_temporary.keys():
                            disease_to_inheritance_temporary[disease_id].append(inheritance_term)
            self.disease_to_inheritance = disease_to_inheritance_temporary

            # disease_id_to_name key represents disease id, value is the disease name correspond to the id. For example:
            # {"ORPHA:1228": "Banki Syndrome"}
            disease_terms_data.set_index("DISEASE_ID", inplace=True)
            self.disease_id_to_name = disease_terms_data.to_dict()["DISEASE_NAME"]

    def validate_gene_list(self, gene_list_files):
        """
            Verifies if there is a gene lists were provided by the client and if the genes are in HUGO Gene Nomenclature Committee (HGNC) symbols. 
            Generates a file containing the genes that are not validated from HGNC and returns a dictionary file names pertaining to the gene list.
            Args:
                gene_list_files: <str> one or more paths to a provided gene list seperated by comma. For example:
                "/project/6033480/resources/gene_list/pediatric_cancer.txt,/project/6033480/resources/gene_list/intellectual_disability_genomics_england_panelapp.tsv"
            Returns:
                A dictionary if the gene lists are provided by the client or False.
                provided_gene_lists_dict: a dictionary with file names as keys and a list of gene names as value. For example:
                {"signature_cancer_susceptibility_genes":  ["AIP","ALK","APC","ARMC5","ATM","ATR,"AXIN2"]} 
        """
        provided_gene_lists_dict = {}
        genes_not_validated_dict = {}

        if len(gene_list_files) != 0:
            # Get list of genes validated by hgnc
            hgnc_genes_list = self.hgnc_file["symbol"].tolist()

            for genes_list_file in gene_list_files.split(","):
                current_file = genes_list_file.rstrip().lstrip()
                # Remove extension file
                extension_file = genes_list_file.split(".")[1]
                file_name = common_tools.remove_file_extension(genes_list_file, extension_file)
                # Get genes list 
                provided_genes_list = np.array([line.rstrip().lstrip() for line in open(os.path.expandvars(current_file))])
                # Check validates genes and not validates genes from HGNC
                common_genes = np.isin(provided_genes_list, hgnc_genes_list)
                validate_genes = provided_genes_list[common_genes]
                not_validate_genes = provided_genes_list[np.logical_not(common_genes)]

                # Exception in one gene which is in previous gene symbol
                not_validate_genes = not_validate_genes[not_validate_genes != "FAM175A"]
                # delete "" in the lists
                not_validate_genes = not_validate_genes[not_validate_genes != ""]
                validate_genes = validate_genes[validate_genes != ""]
                provided_genes_list = provided_genes_list[provided_genes_list != ""]

                # Create dictionaries one for file name corresponding to the genes list and the second file name corresponding to the genes that are not validated by hgnc
                provided_gene_lists_dict[file_name] = provided_genes_list.tolist()
                genes_not_validated_dict[file_name] = not_validate_genes.tolist()

            reporting_genes_not_validated_file = open(self.output_directory + "list_genes_not_validated_hgnc.txt", "a")
            for file_name, list_genes_not_validated in genes_not_validated_dict.items():
                if len(list_genes_not_validated) > 1:
                    reporting_genes_not_validated_file.write("Those genes: " + str(list_genes_not_validated).strip("[]").replace("'",
                                                                                                                                 "") + " are not validated from " + file_name + " file.")
                elif len(list_genes_not_validated) == 1:
                    reporting_genes_not_validated_file.write(
                        "This gene: " + str(list_genes_not_validated).strip("[]").replace("'", "") + " is not validated from " + file_name + " file.")
            reporting_genes_not_validated_file.close()

            return provided_gene_lists_dict

        else:
            return False

    def add_gene_list(self, gene_names):
        """
            Verifies if the gene names are in the gene lists provided by the client and returns name of lists file(s) pertaining to the gene(s).
            Args:
                gene_names: <str> one or more gene names separated by comma. For example:
                "BRCA1,SAMD11"
            Returns:
                A string of file names(s). For example:
                "signature_cancer_susceptibility_genes,pediatric_cancer,cancer_childhood_genomics_england_pannelapp"
        """
        file_name_lists = []
        gene_names = gene_names.split(",")

        for file_name, genes_list in self.provided_gene_lists_dict.items():
            for gene in gene_names:
                if gene in genes_list:
                    file_name_lists.append(file_name)

        if len(file_name_lists) == 0:
            return None
        elif len(file_name_lists) > 0:
            return ",".join(set(file_name_lists))

    def get_canonical_transcript(self, hgvs_value, gene_names):
        """
            Retrieves the canonical transcript corresponding to the gene(s).
            Args:
                hgvs_value: <str> one or more hgvs coding dna seperated by comma. For example:
                "NM_139208.2:c.395C>T,NM_006610.3:c.395C>T"
                gene_names: <str> one or more gene names seperated by comma. For example:
                "BRCA1,SAMD11"
            Returns:
                canonical_transcript_list: a list of canonical transcripts. For example:
                ["NM_007375.3:c.*4561T>G", "NM_006610.3:c.1731A>C"]
        """
        canonical_transcript_list = []
        gene_names = gene_names.split(",")
        hgvs_value = hgvs_value.split(",")

        for variant in hgvs_value:
            refseq_id = variant.split(":")[0].split(".")[0]
            for gene in gene_names:
                try:
                    if refseq_id in self.gene_canonical_transcript_dict[gene]:
                        canonical_transcript_list.append(variant)
                except KeyError:
                    pass

        return ",".join(set(canonical_transcript_list))

    def get_refseq_select_transcript(self, hgvs_value, refseq_select_transcript_ids):
        """
            Retrieves the hgvs value corresponding to the refseq select transcript.
            Args:
                hgvs_value: <str> one or more hgvs coding dna seperated by comma. For example:
                "NM_139208.2:c.395C>T,NM_006610.3:c.395C>T"
                refseq_select_transcript_ids: <str> one or more refseq select transcript id seperated by a comma. For example:
                "NM_139208.2, NM_006610.3"
            Returns:
                refseq_select_transcript_list: <str> on or more hgvs values seperated by comma. For example:
                "NM_000075.3:c.*174T>C,NM_005981.3:c.*1061A>G"
        """

        hgvs_transcript_list = []
        hgvs_value = hgvs_value.split(",")
        refseq_transcript_ids_list = [refseq_select_transcript_id.split(".")[0] for refseq_select_transcript_id in
                                      refseq_select_transcript_ids.split(",")]

        for variant in hgvs_value:
            refseq_id = variant.split(":")[0].split(".")[0]

            if refseq_id in refseq_transcript_ids_list:
                hgvs_transcript_list.append(variant)

        return ",".join(set(hgvs_transcript_list))

    def return_string_format_boolean(self, field):
        """
            Returns "True" or "False" in string format.
            Args:
                field: Which is a None or True
            Returns:
                "False" in case in None
                "True" in case it is not None.
        """
        if field == None:
            return "False"
        elif field:
            return "True"

    def hpo(self, gene_names):
        """
            Parses the Human Phenotype Ontology data and returns diseases associated to a gene and their inheritance mode.
            Args:
                gene_names: <str> name of genes seperated by comma. For example:
                "BRCA,ABCA5"
            Returns:
                A string composed of the disease name, the disease ID and the inheritance mode.the hpo term in hpo column or None for no gene that has hpo term. For example:
                Behçet Disease(ORPHA:117;multifactorial_inheritance)
        """
        all_genes_hpo_terms = []
        genes_list = gene_names.split(",")

        for gene in genes_list:
            hpo_term = []
            if gene in self.gene_diseases.keys():
                disease_id = self.gene_diseases[gene]
                diseases_list = disease_id.split(",")
                for disease in diseases_list:
                    if disease in self.disease_to_inheritance.keys():
                        hpo_term.append(
                            [str(self.disease_id_to_name[disease]) + "(" + str(disease) + ";" + ",".join(self.disease_to_inheritance[disease]) + ")"])
                    else:
                        hpo_term.append([str(self.disease_id_to_name[disease]) + "(" + str(disease) + ")"])
            if hpo_term:
                all_genes_hpo_terms.append(gene + ":" + "|".join(str(term).strip("'[]").replace("\"", "") for term in hpo_term))
            else:
                all_genes_hpo_terms.append(gene + ":" + "None")

        return "/".join(all_genes_hpo_terms)

    def get_hpo_diseases(self, gene_names):
        """
            Retrieve diseases from the hpo gene diseases file based on gene name
            Args:
                gene_names: <str> name of genes seperated by comma. For example:
                "BRCA,ABCA5"
            Returns:
                A string composed of the disease name, the disease ID and the inheritance mode, the hpo term or None for gene that has no hpo term. For example:
                BRCA:Behçet Disease(ORPHA:117;Multifactorial inheritance)
        """
        all_genes_hpo_terms = []
        genes_list = gene_names.split(",")

        for gene in genes_list:
            if gene in self.gene_diseases_parsed.keys():
                all_genes_hpo_terms.append(gene + ":" + self.gene_diseases_parsed[gene])             
            else:
                all_genes_hpo_terms.append(gene + ":" + "None")
        return "/".join(all_genes_hpo_terms)
        
    def get_gnomad_z_score(self, gene_names):
        """
            Parses the gnomad metric z scores file, appends missense, synonymous and loss of function z scores then it returns three z scores. 
            Args:
                gene_names : <str> Gene name. For example: "BRCA"
            Returns: 
                dictionary with the following values missense_z_score, synonymous_z_score, lof_z_score. For example:
                {"missense": -3.46, "synonymous": -2.79, "lof": -0.5}
        """

        genes_list = gene_names.split(",")
        z_score_dict = {}
        missense_z_score = []
        synonymous_z_score = []
        lof_z_score = []

        for gene in genes_list:
            if gene in self.gnomad_file.gene.values:
                missense_z_score.append(str(round(self.gnomad_file["mis_z"][self.gnomad_file["gene"] == gene].values[0], 2)))
                synonymous_z_score.append(str(round(self.gnomad_file["syn_z"][self.gnomad_file["gene"] == gene].values[0], 2)))
                lof_z_score.append(str(round(self.gnomad_file["lof_z"][self.gnomad_file["gene"] == gene].values[0], 2)))
            else:
                missense_z_score.append("None")
                synonymous_z_score.append("None")
                lof_z_score.append("None")
        z_score_dict["missense"] = ",".join(missense_z_score)
        z_score_dict["synonymous"] = ",".join(synonymous_z_score)
        z_score_dict["lof"] = ",".join(lof_z_score)

        return z_score_dict

    def get_reports_list(self, reporting_file):
        """
            Parses the reporting file and returns the lists of samples IDs in each report. 
            Args:
                reporting_file : <str> file location of the samples order file. For example: 
                "/home/user/projects/analysis/data/report_200614_3.csv"
            Returns: 
                reports_list : a list containing one or more lists of sample IDs to include in each report file. For example:
                reports_list = {["1823912", "1829468"], ["1829469"], ["1829470"]}
        """
        reports_list = []

        # if no file provided: print all samples in one report
        if reporting_file == "":
            reports_list.append("all")
        else:
            with open(reporting_file) as samples_id_file:
                for line in samples_id_file:
                    line = line.strip()
                    if line != "SAMPLE_ID":
                        samples = line.split(",")
                        reports_list.append(samples)

        return reports_list

    def gene_count(self, report_data, gene_count_dict):
        """
            Generates the gene count column for the report_data dictionary.
            Args:
                the current report data, but really only the VARIANT_GENE_NAME column is needed
                the gene_count_dict: A dictionary with gene ID as keys and values are the count of each gene
            Returns:
                a list of the gene counts that will make the GENE_COUNT column in the report data. For example:
                ["ATAD3B: 3", "ATAD3A: 2", "CEP104: 4", "CEP104: 4", "PIK3CD: 6", "PIK3CD: 6", "TARDBP: 3,MASP2: 3"]
        """
        gene_count_list = []
        for genes in report_data["VARIANT_GENE_NAME"]:  # .iteritems(): index,
            gene_count_strings_list = []
            for gene in genes.split(","):
                gene_count_strings_list.append(gene + ":" + str(gene_count_dict[gene]))
                # if gene in gene_count_dict.keys():
                # if report_data["GENE_COUNT"].iloc[index] == None:
                #     report_data.at[index, "GENE_COUNT"] = gene + ":" + str(gene_count_dict[gene])
                # else:
                #     report_data.at[index, "GENE_COUNT"] = report_data["GENE_COUNT"].iloc[index] + "," + gene + ":" + str(gene_count_dict[gene])
            gene_count_string = ",".join(gene_count_strings_list)
            gene_count_list.append(gene_count_string)

        return gene_count_list  # report_data

    def get_sample_recessive_variants(self, proband_id, mother_id="", father_id="", depth_min=10, alternate_depth_min=3,
                                      alternate_frequency_min=0.05):
        """
            Generates a list of recessive variants for a specific sample (proband_id).
            If the parent IDs are given conducts the phasing analysis
            A variant is considered recessive if it has at least 2 heterozygous variants in the same gene
            or at least 1 homozygous variant in a gene
            Args:
                proband_id: the ID of the sample on which to do the analysis
                mother_id: the ID of the proband's mother
                father_id: the ID of the proband's father
                depth_min: depth filter used to filter the variants
                alternate_depth_min: alternate allele depth filter used to filter the variants
                alternate_frequency_min: alternate allele frequency filter used to filter the variants
            Returns:
                a dictionary with variants as keys and the phasing as value
        """
        vcf = cyvcf2.VCF(self.vcf_file)
        samples = vcf.samples
        proband_index = samples.index(proband_id)

        # Record genes and genotypes
        genes_genotypes_counts = {}
        for variant in vcf:
            # if (self.variant_caller.upper() == "GATK" or self.variant_caller.upper() == "VARDICT"):
            #     validation = variant.gt_depths[proband_index] >= int(depth_min) and variant.gt_alt_depths[proband_index] >= int(alternate_depth_min) and (variant.gt_alt_depths[proband_index]/variant.gt_depths[proband_index]) >= float(alternate_frequency_min)
            # elif self.variant_caller.upper() == "VARSCAN":
            #     validation = str([Genotype(genotype) for genotype in variant.genotypes][proband_index]) != "0/." and str([Genotype(genotype) for genotype in variant.genotypes][proband_index]) != "./." and variant.format("DP")[proband_index][0] >= int(depth_min) and variant.format("AD")[proband_index][0] >= int(alternate_depth_min) and (variant.format("AD")[proband_index][0]/variant.format("DP")[proband_index][0]) >= float(alternate_frequency_min)
            # else:
            #     validation = False

            validation = common_tools.validate_variant(variant, proband_index, depth_min, alternate_depth_min, alternate_frequency_min,
                                                       self.variant_caller)
            if validation:
                gnomad_format = common_tools.get_gnomad_format(variant)

                if self.annotation_tool == "snpeff":
                    variant_annotation_info_field_name = "ANN"
                if variant.INFO.get(variant_annotation_info_field_name) is not None:
                    # calculate variant frequency
                    alternate_allele_freq = common_tools.get_alternate_allele_frequency(variant, proband_index, self.variant_caller)
                    # determine zygosity
                    zygosity = self.get_zygosity(alternate_allele_freq)
                    # count genotypes according to the alternate allele frequency
                    # if alternate_allele_freq >= 0.1 and alternate_allele_freq < 0.8:
                    #     genotypes["heterozygous"] += 1
                    # elif alternate_allele_freq >= 0.8:
                    #     genotypes["homozygous"] += 1
                    # # elif alternate_allele_freq < 0.1:
                    # #     continue
                    # else: 
                    if zygosity not in ["heterozygous", "homozygous"]:
                        continue

                    # get the genes
                    if self.annotation_tool == "snpeff":
                        genes = self.get_variant_genes_snpeff(variant)
                    for gene in genes:
                        if gene is not None:
                            if gene not in genes_genotypes_counts.keys():
                                genes_genotypes_counts[gene] = {}
                                genes_genotypes_counts[gene]["genotypes"] = {"heterozygous": 0, "homozygous": 0}
                                genes_genotypes_counts[gene]["variants"] = [variant]
                                genes_genotypes_counts[gene]["heterozygous_variants"] = []
                            else:
                                genes_genotypes_counts[gene]["variants"].append(variant)
                            genes_genotypes_counts[gene]["genotypes"][zygosity] += 1
                            if zygosity == "heterozygous":
                                genes_genotypes_counts[gene]["heterozygous_variants"].append(variant)

        # Record recessive genes
        recessive_genes = []
        heterozygous_2_count_genes = []
        for gene in genes_genotypes_counts.keys():
            if genes_genotypes_counts[gene]["genotypes"]["heterozygous"] > 1 or genes_genotypes_counts[gene]["genotypes"]["homozygous"] > 0:
                recessive_genes.append(gene)
            # if the parent information is present, record the genes with a heterozygous count of 2
            if mother_id and father_id:
                if genes_genotypes_counts[gene]["genotypes"]["heterozygous"] == 2:
                    if len(genes_genotypes_counts[gene]["heterozygous_variants"]) == 2:
                        heterozygous_2_count_genes.append(gene)
        print("recessive genes total: " + str(len(recessive_genes)))
        print("genes with two heterozygous variants total: " + str(len(heterozygous_2_count_genes)))

        # select variants
        recessive_variants = {}
        for gene in recessive_genes:
            # recessive_variants[gene] = {}
            for variant in genes_genotypes_counts[gene]["variants"]:
                gnomad_format = common_tools.get_gnomad_format(variant)
                # if variant in recessive_variants.keys():
                if gnomad_format in recessive_variants.keys():
                    recessive_variants[gnomad_format]["recessive"] += ",{gene}:True".format(gene=gene)
                else:
                    recessive_variants[gnomad_format] = {}
                    recessive_variants[gnomad_format]["phasing"] = None
                    recessive_variants[gnomad_format]["recessive"] = "{gene}:True".format(gene=gene)

                if mother_id and father_id:
                    if gene in heterozygous_2_count_genes:
                        if variant in genes_genotypes_counts[gene]["heterozygous_variants"]:
                            if len(genes_genotypes_counts[gene]["heterozygous_variants"]) == 2:
                                variant_index = genes_genotypes_counts[gene]["heterozygous_variants"].index(variant)
                                other_variant_index = 1 - variant_index
                                mother_index = samples.index(mother_id)
                                father_index = samples.index(father_id)
                                phasing = self.get_phasing_recessive_variants(mother_index, father_index,
                                                                              genes_genotypes_counts[gene]["heterozygous_variants"][variant_index],
                                                                              genes_genotypes_counts[gene]["heterozygous_variants"][
                                                                                  other_variant_index])
                                if recessive_variants[gnomad_format]["phasing"] is not None:
                                    recessive_variants[gnomad_format]["phasing"] += ",{gene}:{phasing}".format(gene=gene, phasing=phasing)
                                else:
                                    recessive_variants[gnomad_format]["phasing"] = "{gene}:{phasing}".format(gene=gene, phasing=phasing)
                            else:
                                print("Could not determine phasing for gene {gene} because we could not find 2 variants. Found {count}".format(
                                    gene=gene, count=len(genes_genotypes_counts[gene]["heterozygous_variants"])))
                                if recessive_variants[gnomad_format]["phasing"] is not None:
                                    recessive_variants[gnomad_format]["phasing"] += ",{gene}:NOT_2_COUNT".format(gene=gene)
                                else:
                                    recessive_variants[gnomad_format]["phasing"] = "{gene}:NOT_2_COUNT".format(gene=gene)
                    else:
                        # if recessive_variants[gnomad_format]["phasing"] is not None:
                        #     recessive_variants[gnomad_format]["phasing"] += ",{gene}:NOT_2_COUNT".format(gene=gene)
                        # else:
                        #     recessive_variants[gnomad_format]["phasing"] = "{gene}:NOT_2_COUNT".format(gene=gene)
                        recessive_variants[gnomad_format]["phasing"] = "NA"

                else:
                    # if recessive_variants[gnomad_format]["phasing"] is not None:
                    #     recessive_variants[gnomad_format]["phasing"] += "NA"
                    # else:
                    recessive_variants[gnomad_format]["phasing"] = "NA"

        print("recessive variants: " + str(len(recessive_variants.keys())))
        return recessive_variants, str(len(heterozygous_2_count_genes))

    def get_phasing_recessive_variants(self, mother_index, father_index, variant_1, variant_2):
        """
            Determines the phasing of 2 variants in genes where the heterozygous count equals 2.
            Args:
                mother_index: the index of the mother's sample
                father_index: the index of the father's sample
                variant_1: cyvcf2 Variant object of the currently analyzed variant
                variant_2: cyvcf2 Variant object corresponding to the second variant in the gene
            Returns:
                NO_CALL if the zygosity could not be determined in one of the parents
                REF if both parents are considered reference for variant_1
                CIS-FATHER or CIS-MOTHER if both variants are in a parent and reference in the other
                FATHER or MOTHER if variant_1 is in a parent and reference in the other and variant_2 is reference in both parents
                TRANS if each parent carries the opposite variant and is reference for the other variant
                ERROR if no case covers the calculated values
        """
        variant_statuses = ["heterozygous", "homozygous"]
        zigosity_code = {
            "heterozygous": "h",
            "homozygous": "H",
            "reference": "r",
        }

        # calculate the alternate allele frequencies
        mother_alternate_allele_freq_variant_1 = common_tools.get_alternate_allele_frequency(variant_1, mother_index, self.variant_caller)
        father_alternate_allele_freq_variant_1 = common_tools.get_alternate_allele_frequency(variant_1, father_index, self.variant_caller)
        mother_alternate_allele_freq_variant_2 = common_tools.get_alternate_allele_frequency(variant_2, mother_index, self.variant_caller)
        father_alternate_allele_freq_variant_2 = common_tools.get_alternate_allele_frequency(variant_2, father_index, self.variant_caller)

        # get zygosity for each parent for each variant
        mother_zygosity_variant_1 = self.get_zygosity(mother_alternate_allele_freq_variant_1)
        mother_zygosity_variant_2 = self.get_zygosity(mother_alternate_allele_freq_variant_2)
        father_zygosity_variant_1 = self.get_zygosity(father_alternate_allele_freq_variant_1)
        father_zygosity_variant_2 = self.get_zygosity(father_alternate_allele_freq_variant_2)

        # make decision on phasing
        # if the zygosity call was inconclusive
        if "no_zygosity_call" in [mother_zygosity_variant_1, mother_zygosity_variant_2, father_zygosity_variant_1, father_zygosity_variant_2]:
            phasing = "NO_CALL"

        # if the analyzed variant is reference in both parents
        elif mother_zygosity_variant_1 == "reference" and father_zygosity_variant_1 == "reference":
            phasing = "REF"

        # if the variant is absent in the mother
        elif mother_zygosity_variant_1 == "reference" and mother_zygosity_variant_2 == "reference":
            # if both variants are present in the father
            if father_zygosity_variant_1 in variant_statuses and father_zygosity_variant_2 in variant_statuses:
                phasing = "CIS-FATHER"
            # if the analyzed variant is present in the father
            elif father_zygosity_variant_1 in variant_statuses and father_zygosity_variant_2 == "reference":
                phasing = "FATHER"
            else:
                phasing = "ERROR"

        # if the variant is absent in the father
        elif father_zygosity_variant_1 == "reference" and father_zygosity_variant_2 == "reference":
            # if both variants are present in the mother
            if mother_zygosity_variant_1 in variant_statuses and mother_zygosity_variant_2 in variant_statuses:
                phasing = "CIS-MOTHER"
            # if the analyzed variant is present in the mother
            elif mother_zygosity_variant_1 in variant_statuses and mother_zygosity_variant_2 == "reference":
                phasing = "MOTHER"
            else:
                phasing = "ERROR"

        # if each variant comes from a different parent
        elif (mother_zygosity_variant_1 in variant_statuses and mother_zygosity_variant_2 == "reference" and
              father_zygosity_variant_1 == "reference" and father_zygosity_variant_2 in variant_statuses):
            phasing = "TRANS"

        elif (mother_zygosity_variant_2 in variant_statuses and mother_zygosity_variant_1 == "reference" and
              father_zygosity_variant_2 == "reference" and father_zygosity_variant_1 in variant_statuses):
            phasing = "TRANS"

        # if none of the cases is caught
        else:
            print("PHASING UNKNOWN for {variant}: ".format(variant=common_tools.get_gnomad_format(variant_1)) +
                  ",".join([mother_zygosity_variant_1, mother_zygosity_variant_2, father_zygosity_variant_1, father_zygosity_variant_2]))
            phasing = "{mother_zygosity_variant_1}{father_zygosity_variant_1}/{mother_zygosity_variant_2}{father_zygosity_variant_2}".format(
                mother_zygosity_variant_1=zigosity_code[mother_zygosity_variant_1],
                mother_zygosity_variant_2=zigosity_code[mother_zygosity_variant_2],
                father_zygosity_variant_1=zigosity_code[father_zygosity_variant_1],
                father_zygosity_variant_2=zigosity_code[father_zygosity_variant_2],
            )

        return phasing
        
    def denovo(self, proband_id, father_id, mother_id, depth_min=10, alternate_depth_min=3, alternate_frequency_min=0.05):
        """
            Searches for denovo variants in the proband sample in the vcf file.
            Args:
                proband_id: <str> ID of proband. For example:
                "PRO-001"
                father_id: <str> ID of father. For example:
                "PRO-002"
                mother_id: <str> ID of mother. For example:
                "PRO-003"
                denovo_parameters: dict containing the denovo analysis parameters
                depth_min : <int> depth coverage of proband. For example:10
                alternate_depth_min : <int> alternate depth. For example:3
                alternate_frequency_min : <float> maximum alternate frequency. For example: 0.05
            Returns:
                gnomad_format_denovo: a list containing the gnomad notation of the denovo variants. For example:
                ["1-48539539", "2-34355334", "2-322424242"]
        """
        # Declaring list
        gnomad_format_denovo = []

        denovo_parameters = {
            "proband_depth_coverage": self.config.get("REPORT_VARIANTS", "DENOVO_MIN_DEPTH_PROBAND"),
            "parents_depth_coverage": self.config.get("REPORT_VARIANTS", "DENOVO_MIN_DEPTH_PARENTS"),
            "proband_alternate_depth": self.config.get("REPORT_VARIANTS", "DENOVO_MIN_ALT_DEPTH_PROBAND"),
            "parents_alternate_depth": self.config.get("REPORT_VARIANTS", "DENOVO_MAX_ALT_DEPTH_PARENTS"),
            "proband_alternate_frequency": self.config.get("REPORT_VARIANTS", "DENOVO_MIN_ALT_AF_PROBAND"),
            "parents_alternate_frequency": self.config.get("REPORT_VARIANTS", "DENOVO_MAX_ALT_AF_PARENTS")
        }

        parsed_vcf = cyvcf2.VCF(self.vcf_file)
        # Get family indexes from vcf
        proband_index = parsed_vcf.samples.index(proband_id)
        father_index = parsed_vcf.samples.index(father_id)
        mother_index = parsed_vcf.samples.index(mother_id)

        for variant in parsed_vcf:
            validation = common_tools.validate_variant(variant, proband_index, depth_min, alternate_depth_min, alternate_frequency_min,
                                                       self.variant_caller)
            if validation:
                gnomad_format = common_tools.get_gnomad_format(variant)
                # Depth of coverage >= 10 in trio 
                if (variant.gt_depths[proband_index] >= int(denovo_parameters["proband_depth_coverage"]) and
                        variant.gt_depths[father_index] >= int(denovo_parameters["parents_depth_coverage"]) and
                        variant.gt_depths[mother_index] >= int(denovo_parameters["parents_depth_coverage"])):
                    # Alt allele supported by >= 4 reads for proband and <= 2 for parents
                    if (variant.gt_alt_depths[proband_index] >= int(denovo_parameters["proband_alternate_depth"]) and
                            variant.gt_alt_depths[father_index] <= int(denovo_parameters["parents_alternate_depth"]) and
                            variant.gt_alt_depths[mother_index] <= int(denovo_parameters["parents_alternate_depth"])):
                        # Get  alternate frequency for family 
                        alt_frequency_proband = common_tools.get_alternate_allele_frequency(variant, proband_index, self.variant_caller)
                        alt_frequency_father = common_tools.get_alternate_allele_frequency(variant, father_index, self.variant_caller)
                        alt_frequency_mother = common_tools.get_alternate_allele_frequency(variant, mother_index, self.variant_caller)
                        # Confirm parents are reference homozygous
                        # get zygosity for each parent for each variant
                        mother_zygosity = self.get_zygosity(alt_frequency_mother)
                        father_zygosity = self.get_zygosity(alt_frequency_father)
                        if mother_zygosity == "reference" and father_zygosity == "reference":
                            # Alt allele in >= 20% of total reads in proband and <= 5% in parents
                            if (float(alt_frequency_proband) >= float(denovo_parameters["proband_alternate_frequency"]) and
                                    float(alt_frequency_father) <= float(denovo_parameters["parents_alternate_frequency"]) and
                                    float(alt_frequency_mother) <= float(denovo_parameters["parents_alternate_frequency"])):
                                gnomad_format_denovo.append(gnomad_format)
        print("denovo variants total: " + str(len(gnomad_format_denovo)))
        return gnomad_format_denovo

    def get_zygosity(self, alternate_allele_freq):
        """
            Determines the zygosity of a variant according to the alternate allele frequency.
            Args:
                alternate allele frequency
            Returns:
                heterozygous if the frequency is over 0.1 and less than 0.8
                homozygous if the frequency is over 0.8
                reference if the frequency is less than 0.1
                no_zygosity_call in any other case
        """
        if 0.1 <= alternate_allele_freq < 0.8:
            zygosity = "heterozygous"
        elif alternate_allele_freq >= 0.8:
            zygosity = "homozygous"
        elif alternate_allele_freq < 0.1:
            zygosity = "reference"
        else:
            zygosity = "no_zygosity_call"
        return zygosity

    def get_variant_cluster(self, sample_id, variants_list):
        """
            Generates a list of variants that are clustered.
            Args:
                sample_id : <str> ID of sample. For example:
                "PRO-001"
                variants_list: list of variants present in the sample.
            Returns:
                a list of variant that are in clusters.
        """
        clustered_variants = []
        vcf_specific_sample = cyvcf2.VCF(self.vcf_file, samples=sample_id)
        for variant in variants_list:
            chr = variant.CHROM
            window = 5  # means window of 10 (+- 5 of position)
            # Previously the ID was set this way
            # "-".join([variant.CHROM, str(variant.start), variant.REF, "".join(variant.ALT)])
            variant_id = common_tools.get_gnomad_format(variant)
            for specific_variant in vcf_specific_sample(chr + ":" + str(int(variant.start - window)) + "-" + str(int(variant.start + window))):
                # Previously the ID was set this way
                # "-".join([specific_variant.CHROM, str(specific_variant.start), specific_variant.REF, "".join(specific_variant.ALT)])
                specific_variant_id = common_tools.get_gnomad_format(specific_variant)
                if specific_variant_id == variant_id:
                    continue
                else:
                    if specific_variant.genotypes:
                        current_genotype = Genotype(specific_variant.genotypes[0])
                        if str(current_genotype) in ["0/1", "1/1", "0|1", "1|1"]:
                            clustered_variants.append(variant)
        return clustered_variants

    def generate_sample_variants_dictionary(self, vcf_file, depth_min, alternate_depth_min, alternate_frequency_min, clustering=True):
        """
            Generates a dictionary containing sample IDs as keys with relevant variants and variant clusters.
            Args:
                vcf_file: vcf file to parse
                depth_min: <int> minimum depth of analysis. For example: 10
                alternate_depth_min: <int> minimum alternate allele depth. For example: 5
                alternate_frequency_min: <float> minimum alternate allele frequency. For example: 0.05
                clustering: boolean to add variant clusters
            Returns: 
                sample_variants_dictionary: a dictionary key represents sample ID and value represents a dictionary with two keys, s key of a list of variants and key of an index sample. For example : 
                sample_variants_dictionary = {"PRO-001": {variants : [cyvcf2.Variant], "sample_index": 3, "variant_cluster": [cyvcf2.Variant]}}
        """
        sample_variants_dictionary = {}
        vcf = cyvcf2.VCF(vcf_file)
        samples = vcf.samples
        for variant in vcf:
            for index, sample in enumerate(samples):
                if sample not in sample_variants_dictionary:
                    sample_variants_dictionary[sample] = {"variants": [], "sample_index": index, "var_cluster": []}
                validation = common_tools.validate_variant(variant, index, depth_min, alternate_depth_min, alternate_frequency_min,
                                                           self.variant_caller)
                if validation:
                    sample_variants_dictionary[sample]["variants"].append(variant)
                    # sample_variants_dictionary[sample]["sample_index"] = index

        for sample_id, sample_data_dict in sample_variants_dictionary.items():
            variants_list = sample_data_dict["variants"]
            if clustering:
                sample_variants_dictionary[sample_id]["var_cluster"] = self.get_variant_cluster(sample_id, variants_list)

        return sample_variants_dictionary

    def get_variant_genes_snpeff(self, variant):
        """
            Gets a list of genes from the snpeff annotation of a variant
            Args:
                variant: a cyvcf2 variant object
            Returns:
                a list of genes
        """
        # get the genes
        snpeff_alleles = variant.INFO.get("ANN").split(",")
        genes = []
        for allele in snpeff_alleles:
            # the gene name is the 4th column of the allele annotation. gene id is 5th
            gene = allele.split("|")[3]
            if gene not in genes:
                genes.append(gene)
        return genes

    def get_report_output_file_name(self, samples_to_report):
        """
            Name the report file depending on the number of samples to report
            Args:
                samples_to_report: list from self.reports_list
            Returns:
                string: the file name of the report file to be created
        """
        if str(samples_to_report) == "all":
            output_file = "report_variants.tsv"
        else:
            if len(samples_to_report) == 1:
                output_file = "report_variants_" + samples_to_report[0] + ".tsv"
            elif len(samples_to_report) == 2:
                output_file = "report_variants_" + samples_to_report[0] + "_duo.tsv"
            elif len(samples_to_report) == 3:
                output_file = "report_variants_" + samples_to_report[0] + "_trio.tsv"
            elif len(samples_to_report) == 4 and samples_to_report[3] == "family_trio":
                output_file = "report_variants_" + samples_to_report[0] + "_trio.tsv"
            else:
                output_file = "report_variants_" + samples_to_report[0] + "_multi_samples.tsv"
        return output_file

    def create_report_summary(self, output_directory):
        """
            Creates a file to save the reporting information
            Args:
                output_directory: where the file will be saved
            Returns:
                a file object
        """
        # Change report summary file name in generate_report_variants_command as well
        output_file = output_directory + "report_summary.tsv"
        report_summary = open(output_file, "w")
        report_summary.write(
            "\t".join(["SAMPLE", "FINAL_VARIANT_COUNT", "RECESSIVE_VARIANTS", "DENOVO", "2_HET_GENES"]) + "\n")  # , "INITIAL_VARIANT_COUNT"
        return report_summary

    def format_meta_scores_columns(self, meta_score: str):
        """
        Format the META scores from the VCF file for the report. Remove duplicated scores and dots.
        Args:
            meta_score: string from the VCF INFO field
        Returns: formatted meta score

        """
        if meta_score is None or meta_score == ".":
            meta_score = "None"
        else:
            score_set = set(meta_score.split(","))
            # remove dots, works even if the set doesn't contain it
            score_set.discard(".")
            # if there is still something in the set
            if score_set:
                filter_list = []
                for score in score_set:
                    val = score.split(",")
                    if "." in val:
                        val.remove(".")
                    filter_list.append(val)
                meta_score = ";".join({i for sublist in filter_list for i in sublist})
            else:
                meta_score = "None"
        return meta_score

    def generate_report(self, minimum_depth, minimum_depth_alternate, minimum_alternate_frequency):
        """
        Generates the reports in TSV format.
        Args:
            minimum_depth: minimum ref depth
            minimum_depth_alternate: minimum alt depth
            minimum_alternate_frequency: minimum alt freq
        Returns:
        """

        report_columns = [column_name for (column_name, boolean) in self.report_parameters.items() if boolean]
        self.parse_reference_resources(report_columns)

        sample_variants_dictionary = self.generate_sample_variants_dictionary(self.vcf_file, minimum_depth, minimum_depth_alternate,
                                                                              minimum_alternate_frequency)
        report_summary = self.create_report_summary(self.output_directory)

        for samples_to_report in self.reports_list:
            gene_count_dict = {}

            # Removing DENOVO column to add it only if the samples_to_report is a family trio
            if "DENOVO" in report_columns:
                report_columns.remove("DENOVO")

            # defining these as empty to avoid undefined error when comparing
            proband_id = ""
            mother_id = ""
            father_id = ""

            output_file = self.get_report_output_file_name(samples_to_report)

            if str(samples_to_report) == "all":
                samples_to_report = sample_variants_dictionary.keys()
            elif len(samples_to_report) == 4 and samples_to_report[3] == "family_trio":
                proband_id = samples_to_report[0]
                mother_id = samples_to_report[1]
                father_id = samples_to_report[2]
                # remove family_trio from list
                samples_to_report = samples_to_report[:-1]
                if "DENOVO" not in report_columns:
                    report_columns.insert(-1, "DENOVO")
                if "RECESSIVE" not in report_columns:
                    report_columns.insert(-1, "RECESSIVE")
                if "PHASING" not in report_columns:
                    report_columns.insert(-1, "PHASING")

            report_data = pd.DataFrame(columns=report_columns)
            # for sample_id, variants_data in sample_variants_dictionary.items():
            #     if sample_id not in samples_to_report:
            #         continue

            for sample_id in samples_to_report:
                try:
                    variants_data = sample_variants_dictionary[sample_id]
                except KeyError:
                    print("WARNING: Did not find any variant data for sample {sample_id}".format(sample_id=sample_id))
                if len(variants_data["variants"]) > 0:
                    print("Sample {sample} has {variants} variants".format(sample=sample_id, variants=str(len(variants_data["variants"]))))
                    if sample_id == proband_id and "DENOVO" in report_columns:
                        gnomad_position_denovo = self.denovo(proband_id, father_id, mother_id)
                        # self.denovo_depth_coverage_proband, self.denovo_depth_coverage_parents,
                        # self.denovo_alternate_depth_proband, self.denovo_alternate_depth_parents, self.denovo_alternate_frequency_proband, self.denovo_alternate_frequency_parents)
                    else:
                        gnomad_position_denovo = []
                    if "RECESSIVE" in report_columns:
                        if sample_id == proband_id:
                            recessive_variants, heterozygous_genes = self.get_sample_recessive_variants(sample_id, mother_id, father_id)
                        else:
                            recessive_variants, heterozygous_genes = self.get_sample_recessive_variants(sample_id)
                    else:
                        recessive_variants = []
                        heterozygous_genes = 0
                    # pre_filter_count = str(len(pre_filtered_sample_variants_dictionary[sample_id]["variants"])) , pre_filter_count
                    sample_summary = "\t".join(
                        [sample_id, str(len(variants_data["variants"])), str(len(recessive_variants)), str(len(gnomad_position_denovo)),
                         str(heterozygous_genes)])
                    report_summary.write(sample_summary + "\n")
                    for variant in variants_data["variants"]:
                        report_data_dict = {}

                        if "SAMPLE_ID" in report_columns:
                            report_data_dict["SAMPLE_ID"] = sample_id

                        if "POSITION" in report_columns:
                            if str(variant.CHROM).startswith("chr"):
                                report_data_dict["POSITION"] = str(variant.CHROM) + ":" + str(variant.start + 1)
                            else:
                                report_data_dict["POSITION"] = "chr" + str(variant.CHROM) + ":" + str(variant.start + 1)

                        if "ALLELES" in report_columns:
                            alleles = variant.REF
                            alleles += "/"
                            if len(variant.ALT) > 0:
                                for alt in variant.ALT:
                                    alleles += alt + ","
                            # removing the trailing coma
                            alleles = alleles[:-1]
                            report_data_dict["ALLELES"] = alleles

                        if "FILTER" in report_columns:
                            filter_field = ""
                            if variant in sample_variants_dictionary[sample_id]["var_cluster"]:
                                if variant.FILTER is None:
                                    filter_field = "varCluster"
                                else:
                                    filter_field = variant.FILTER + ", varCluster"
                            else:
                                filter_field = variant.FILTER
                            if filter_field == "":
                                # why as a string and not as None?
                                report_data_dict["FILTER"] = "None"
                            else:
                                report_data_dict["FILTER"] = filter_field

                        sample_index = variants_data["sample_index"]
                        # see https://github.com/brentp/cyvcf2/issues/58
                        genotype_list = [Genotype(one_genotype) for one_genotype in variant.genotypes]

                        if "SAMPLE_GT" in report_columns:
                            report_data_dict["SAMPLE_GT"] = genotype_list[sample_index]

                        if "SAMPLE_GT_QUAL" in report_columns:
                            if self.variant_caller.upper() == "VARDICT":
                                report_data_dict["SAMPLE_GT_QUAL"] = None
                            else:
                                try:
                                    report_data_dict["SAMPLE_GT_QUAL"] = common_tools.format_number_for_report(variant.format("GQ")[sample_index][0])
                                except TypeError:
                                    report_data_dict["SAMPLE_GT_QUAL"] = None

                        if self.variant_caller.upper() == "GATK" or self.variant_caller.upper() == "VARDICT":
                            if "SAMPLE_DEPTH" in report_columns:
                                report_data_dict["SAMPLE_DEPTH"] = variant.gt_depths[sample_index]
                            if "SAMPLE_DEPTH_REF" in report_columns:
                                report_data_dict["SAMPLE_DEPTH_REF"] = variant.gt_ref_depths[sample_index]
                            if "SAMPLE_DEPTH_ALT" in report_columns:
                                report_data_dict["SAMPLE_DEPTH_ALT"] = variant.gt_alt_depths[sample_index]
                            # report_data_dict["SAMPLE_ALT_AF"] = common_tools.format_number_for_report(variant.gt_alt_depths[sample_index]/variant.gt_depths[sample_index])
                        if self.variant_caller.upper() == "VARSCAN":
                            if "SAMPLE_DEPTH" in report_columns:
                                report_data_dict["SAMPLE_DEPTH"] = variant.format("DP")[sample_index][0]
                            if "SAMPLE_DEPTH_REF" in report_columns:
                                report_data_dict["SAMPLE_DEPTH_REF"] = variant.format("RD")[sample_index][0]
                            if "SAMPLE_DEPTH_ALT" in report_columns:
                                report_data_dict["SAMPLE_DEPTH_ALT"] = variant.format("AD")[sample_index][0]
                            # report_data_dict["SAMPLE_ALT_AF"] = common_tools.format_number_for_report(variant.format("AD")[sample_index][0]/variant.format("DP")[sample_index][0])

                        if "SAMPLE_ALT_AF" in report_columns:
                            report_data_dict["SAMPLE_ALT_AF"] = common_tools.format_number_for_report(
                                common_tools.get_alternate_allele_frequency(variant, sample_index, self.variant_caller))

                        if self.annotation_tool == "snpeff":
                            annotation_data = common_tools.parse_snpeff_annotation(variant.INFO.get("ANN"))
                        genes = annotation_data["gene_name"].split(",")
                        for gene in genes:
                            if gene in gene_count_dict.keys():
                                gene_count_dict[gene] = int(gene_count_dict[gene]) + 1
                            elif gene not in gene_count_dict.keys():
                                gene_count_dict[gene] = 1

                        if "VARIANT_GENE_NAME" in report_columns:
                            report_data_dict["VARIANT_GENE_NAME"] = annotation_data["gene_name"]
                        if "VARIANT_ANNOTATION" in report_columns:
                            report_data_dict["VARIANT_ANNOTATION"] = annotation_data["annotation"]
                        if "VARIANT_HGVS_C" in report_columns:
                            report_data_dict["VARIANT_HGVS_C"] = annotation_data["hgvs_c"]
                        if "VARIANT_HGVS_P" in report_columns:
                            report_data_dict["VARIANT_HGVS_P"] = annotation_data["hgvs_p"]
                        if "VARIANT_EXON_INTRON" in report_columns:
                            report_data_dict["VARIANT_EXON_INTRON"] = annotation_data["exon"]

                        if "CANONICAL_HGVS_C" in report_columns:
                            canonical_cdna = self.get_canonical_transcript(annotation_data["hgvs_c"], annotation_data["gene_name"])
                            if canonical_cdna != "":
                                report_data_dict["CANONICAL_HGVS_C"] = canonical_cdna
                            # in case that there is no canonical transcript for the gene name - add None
                            else:
                                report_data_dict["CANONICAL_HGVS_C"] = "None"

                        if "CANONICAL_HGVS_P" in report_columns:
                            canonical_protein = self.get_canonical_transcript(annotation_data["hgvs_p"], annotation_data["gene_name"])
                            if canonical_protein != "":
                                report_data_dict["CANONICAL_HGVS_P"] = canonical_protein
                            # in case that there is no canonical transcript for the gene ID - add None
                            else:
                                report_data_dict["CANONICAL_HGVS_P"] = "None"

                        if "SELECT_HGVS_C" in report_columns:
                            ref_select_value = variant.INFO.get("REF_SELECT")
                            if ref_select_value is None or ref_select_value == "":
                                report_data_dict["SELECT_HGVS_C"] = "None"
                            else:
                                refseq_select_hgvs = self.get_refseq_select_transcript(annotation_data["hgvs_c"], ref_select_value)
                                if refseq_select_hgvs == "":
                                    refseq_select_hgvs = "None"
                                report_data_dict["SELECT_HGVS_C"] = refseq_select_hgvs

                        if "SELECT_HGVS_P" in report_columns:
                            ref_select_value = variant.INFO.get("REF_SELECT")
                            if ref_select_value is None or ref_select_value == "":
                                report_data_dict["SELECT_HGVS_C"] = "None"
                            else:
                                refseq_select_hgvs = self.get_refseq_select_transcript(annotation_data["hgvs_p"], ref_select_value)
                                if refseq_select_hgvs == "":
                                    refseq_select_hgvs = "None"
                                report_data_dict["SELECT_HGVS_P"] = refseq_select_hgvs

                        if "VARIANT_CODING" in report_columns:
                            report_data_dict["VARIANT_CODING"] = annotation_data["coding"]
                        if "VARIANT_LOF" in report_columns:
                            report_data_dict["VARIANT_LOF"] = annotation_data["lof"]
                        if "DBSNP_ID" in report_columns:
                            report_data_dict["DBSNP_ID"] = variant.INFO.get("DBSID")
                        if "CLINVAR_ID" in report_columns:
                            report_data_dict["CLINVAR_ID"] = variant.INFO.get("CLNVID")
                        if "CLINVAR_SIG" in report_columns:
                            report_data_dict["CLINVAR_SIG"] = variant.INFO.get("CLNSIG")
                        if "CLINVAR_SIG_CONF" in report_columns:
                            report_data_dict["CLINVAR_SIG_CONF"] = variant.INFO.get("CLNSIGCONF")
                        if "CLINVAR_DN" in report_columns:
                            report_data_dict["CLINVAR_DN"] = variant.INFO.get("CLNDN")
                        if "CLINVAR_REVSTAT" in report_columns:
                            report_data_dict["CLINVAR_REVSTAT"] = variant.INFO.get("CLNREVSTAT")
                        if "COSMIC_ID" in report_columns:
                            report_data_dict["COSMIC_ID"] = variant.INFO.get("CMCID")
                        if "COSMIC_GENE" in report_columns:
                            report_data_dict["COSMIC_GENE"] = variant.INFO.get("CMCGENE")
                        if "COSMIC_CDS" in report_columns:
                            report_data_dict["COSMIC_CDS"] = variant.INFO.get("CMCCDS")
                        if "COSMIC_AA" in report_columns:
                            report_data_dict["COSMIC_AA"] = variant.INFO.get("CMCAA")
                        # This is necessary because one variant can have multiple cosmic entries with different values 
                        # and in this case it is returned as a tuple and not a string like for the previous cosmic values
                        if "COSMIC_CNT" in report_columns:
                            # Replace the parentheses from casting a tuple to a string and spaces after the comma
                            report_data_dict["COSMIC_CNT"] = str(variant.INFO.get("CMCCNT")).replace("(", "").replace(")", "").replace(" ", "")

                        if "GNOMAD_POS" in report_columns:
                            gnomad_format = common_tools.get_gnomad_format(variant)
                            report_data_dict["GNOMAD_POS"] = gnomad_format

                        if "GNOMAD_EX_AC" in report_columns:
                            gnomad_ex_ac = variant.INFO.get("GNDEXAC")
                            if gnomad_ex_ac is None:
                                gnomad_ex_ac = "None"
                            report_data_dict["GNOMAD_EX_AC"] = gnomad_ex_ac

                        if "GNOMAD_EX_AN" in report_columns:
                            gnomad_ex_an = variant.INFO.get("GNDEXAN")
                            if gnomad_ex_an is None:
                                gnomad_ex_an = "None"
                            report_data_dict["GNOMAD_EX_AN"] = gnomad_ex_an

                        if "GNOMAD_EX_AF" in report_columns:
                            gnomad_ex_af = variant.INFO.get("GNDEXAF")
                            if gnomad_ex_af is None:
                                gnomad_ex_af = "None"
                            report_data_dict["GNOMAD_EX_AF"] = common_tools.format_number_for_report(gnomad_ex_af)

                        if "GNOMAD_EX_FILT" in report_columns:
                            report_data_dict["GNOMAD_EX_FILT"] = variant.INFO.get("GNDEXFILT")

                        if "GNOMAD_EX_HOM" in report_columns:
                            report_data_dict["GNOMAD_EX_HOM"] = variant.INFO.get("GNDEXHOM")

                        if "GNOMAD_WG_AC" in report_columns:
                            gnomad_wg_ac = variant.INFO.get("GNDWGAC")
                            if gnomad_wg_ac is None:
                                gnomad_wg_ac = "None"
                            report_data_dict["GNOMAD_WG_AC"] = gnomad_wg_ac

                        if "GNOMAD_WG_AN" in report_columns:
                            gnomad_wg_an = variant.INFO.get("GNDWGAN")
                            if gnomad_wg_an is None:
                                gnomad_wg_an = "None"
                            report_data_dict["GNOMAD_WG_AN"] = gnomad_wg_an

                        if "GNOMAD_WG_AF" in report_columns:
                            gnomad_wg_af = variant.INFO.get("GNDWGAF")
                            if gnomad_wg_af is None:
                                gnomad_wg_af = "None"
                            report_data_dict["GNOMAD_WG_AF"] = common_tools.format_number_for_report(gnomad_wg_af)

                        if "GNOMAD_WG_FILT" in report_columns:
                            report_data_dict["GNOMAD_WG_FILT"] = variant.INFO.get("GNDWGFILT")

                        if "GNOMAD_WG_HOM" in report_columns:
                            report_data_dict["GNOMAD_WG_HOM"] = variant.INFO.get("GNDWGHOM")

                        if "GNOMAD_AC" in report_columns:
                            gnomad_ex_ac = variant.INFO.get("GNDEXAC")
                            gnomad_wg_ac = variant.INFO.get("GNDWGAC")
                            if gnomad_ex_ac is None and gnomad_wg_ac is None:
                                gnomad_ac = "None"
                            else:
                                if gnomad_ex_ac is None:
                                    gnomad_ex_ac = 0
                                if gnomad_wg_ac is None:
                                    gnomad_wg_ac = 0
                                gnomad_ac = gnomad_ex_ac + gnomad_wg_ac
                            report_data_dict["GNOMAD_AC"] = gnomad_ac

                        if "GNOMAD_AN" in report_columns:
                            gnomad_ex_an = variant.INFO.get("GNDEXAN")
                            gnomad_wg_an = variant.INFO.get("GNDWGAN")
                            if gnomad_wg_an is None and gnomad_ex_an is None:
                                gnomad_an = "None"
                            else:
                                if gnomad_ex_an is None:
                                    gnomad_ex_an = 0
                                if gnomad_wg_an is None:
                                    gnomad_wg_an = 0
                                gnomad_an = gnomad_ex_an + gnomad_wg_an

                            report_data_dict["GNOMAD_AN"] = gnomad_an

                        if "GNOMAD_AF" in report_columns:
                            report_data_dict["GNOMAD_AF"] = common_tools.format_number_for_report(
                                common_tools.calculate_gnomad_allele_frequency(variant))

                        if "GNOMAD_HOM" in report_columns:
                            gnomad_ex_hom = variant.INFO.get("GNDEXHOM")
                            gnomad_wg_hom = variant.INFO.get("GNDWGHOM")
                            if gnomad_ex_hom is None and gnomad_wg_hom is None:
                                gnomad_hom = "None"
                            else:
                                if gnomad_ex_hom is None:
                                    gnomad_ex_hom = 0
                                if gnomad_wg_hom is None:
                                    gnomad_wg_hom = 0
                                gnomad_hom = gnomad_ex_hom + gnomad_wg_hom
                            report_data_dict["GNOMAD_HOM"] = gnomad_hom

                        if "GNOMAD_FILT" in report_columns:
                            gnomad_ex_filter = variant.INFO.get("GNDEXFILT")
                            gnomad_wg_filter = variant.INFO.get("GNDWGFILT")
                            if gnomad_ex_filter is None and gnomad_wg_filter is None:
                                gnomad_filter = "None"
                            else:
                                if gnomad_ex_filter is None:
                                    gnomad_ex_filter = ""
                                if gnomad_wg_filter is None:
                                    gnomad_wg_filter = ""
                                gnomad_filter = gnomad_ex_filter + gnomad_wg_filter
                            report_data_dict["GNOMAD_FILT"] = gnomad_filter

                        if "GNOMAD_WG_POPMAX" in report_columns:
                            gnomad_wg_pop_max = variant.INFO.get("GNDWGPOPMAX")
                            if gnomad_wg_pop_max is None:
                                gnomad_wg_pop_max = "None"
                            report_data_dict["GNOMAD_WG_POPMAX"] = gnomad_wg_pop_max

                        if "GNOMAD_WG_AN_POPMAX" in report_columns:
                            gnomad_wg_an_pop_max = variant.INFO.get("GNDWGAN_POPMAX")
                            if gnomad_wg_an_pop_max is None:
                                gnomad_wg_an_pop_max = "None"
                            report_data_dict["GNOMAD_WG_AN_POPMAX"] = gnomad_wg_an_pop_max

                        if "GNOMAD_WG_AC_POPMAX" in report_columns:
                            gnomad_wg_ac_pop_max = variant.INFO.get("GNDWGAC_POPMAX")
                            if gnomad_wg_ac_pop_max is None:
                                gnomad_wg_ac_pop_max = "None"
                            report_data_dict["GNOMAD_WG_AC_POPMAX"] = gnomad_wg_ac_pop_max

                        if "GNOMAD_WG_AF_POPMAX" in report_columns:
                            gnomad_wg_af_pop_max = variant.INFO.get("GNDWGAF_POPMAX")
                            if gnomad_wg_af_pop_max is None:
                                gnomad_wg_af_pop_max = "None"
                            report_data_dict["GNOMAD_WG_AF_POPMAX"] = common_tools.format_number_for_report(gnomad_wg_af_pop_max)

                        if "GNOMAD_WG_NHOM_POPMAX" in report_columns:
                            gnomad_wg_hom_pop_max = variant.INFO.get("GNDWGNHOM_POPMAX")
                            if gnomad_wg_hom_pop_max is None:
                                gnomad_wg_hom_pop_max = "None"
                            report_data_dict["GNOMAD_WG_NHOM_POPMAX"] = gnomad_wg_hom_pop_max

                        if "GNOMAD_EX_POPMAX" in report_columns:
                            gnomad_ex_pop_max = variant.INFO.get("GNDEXPOPMAX")
                            if gnomad_ex_pop_max is None:
                                gnomad_ex_pop_max = "None"
                            report_data_dict["GNOMAD_EX_POPMAX"] = gnomad_ex_pop_max

                        if "GNOMAD_EX_AN_POPMAX" in report_columns:
                            gnomad_ex_an_pop_max = variant.INFO.get("GNDEXAN_POPMAX")
                            if gnomad_ex_an_pop_max is None:
                                gnomad_ex_an_pop_max = "None"
                            report_data_dict["GNOMAD_EX_AN_POPMAX"] = gnomad_ex_an_pop_max

                        if "GNOMAD_EX_AC_POPMAX" in report_columns:
                            gnomad_ex_ac_pop_max = variant.INFO.get("GNDEXAC_POPMAX")
                            if gnomad_ex_ac_pop_max is None:
                                gnomad_ex_ac_pop_max = "None"
                            report_data_dict["GNOMAD_EX_AC_POPMAX"] = gnomad_ex_ac_pop_max

                        if "GNOMAD_EX_AF_POPMAX" in report_columns:
                            gnomad_ex_af_pop_max = variant.INFO.get("GNDEXAF_POPMAX")
                            if gnomad_ex_af_pop_max is None:
                                gnomad_ex_af_pop_max = "None"
                            report_data_dict["GNOMAD_EX_AF_POPMAX"] = common_tools.format_number_for_report(gnomad_ex_af_pop_max)

                        if "GNOMAD_EX_NHOM_POPMAX" in report_columns:
                            gnomad_ex_hom_pop_max = variant.INFO.get("GNDEXNHOM_POPMAX")
                            if gnomad_ex_hom_pop_max is None:
                                gnomad_ex_hom_pop_max = "None"
                            report_data_dict["GNOMAD_EX_NHOM_POPMAX"] = gnomad_ex_hom_pop_max

                        if "GNOMAD_EX_NON_CANCER_AF" in report_columns:
                            gnomad_ex_non_cancer_af = variant.INFO.get("GNDEXNON_CANCER_AF")
                            if gnomad_ex_non_cancer_af is None:
                                gnomad_ex_non_cancer_af = "None"
                            report_data_dict["GNOMAD_EX_NON_CANCER_AF"] = gnomad_ex_non_cancer_af

                        if "SEGDUP" in report_columns:
                            report_data_dict["SEGDUP"] = self.return_string_format_boolean(variant.INFO.get("SEGDUP"))

                        if "REP_MASK" in report_columns:
                            report_data_dict["REP_MASK"] = self.return_string_format_boolean(variant.INFO.get("REP_MASK"))

                        if "REP_SIMPLE" in report_columns:
                            report_data_dict["REP_SIMPLE"] = self.return_string_format_boolean(variant.INFO.get("REP_SIMPLE"))

                        # if "CADD_RAW" in report_columns :
                        #     cadd_raw = variant.INFO.get("CADD_RAW")
                        #     if cadd_raw is None:
                        #         cadd_raw = 0
                        #     report_data_dict["CADD_RAW"] = common_tools.format_number_for_report(cadd_raw)

                        if "CADD_PHRED_SNV" in report_columns:
                            cadd_phred = variant.INFO.get("CADD_PHRED_SNV")
                            if cadd_phred is None:
                                cadd_phred = "None"
                            report_data_dict["CADD_PHRED_SNV"] = common_tools.format_number_for_report(cadd_phred)

                        if "CADD_PHRED_INDEL" in report_columns:
                            cadd_phred = variant.INFO.get("CADD_PHRED_INDEL")
                            if cadd_phred is None:
                                cadd_phred = "None"
                            report_data_dict["CADD_PHRED_INDEL"] = common_tools.format_number_for_report(cadd_phred)

                        if "CADD_PHRED" in report_columns:
                            cadd_phred = variant.INFO.get("CADD_PHRED")
                            if cadd_phred is None:
                                cadd_phred = "None"
                            report_data_dict["CADD_PHRED"] = common_tools.format_number_for_report(cadd_phred)

                        if "REMM" in report_columns:
                            remm = variant.INFO.get("REMM")
                            if remm is None:
                                remm = "None"
                            report_data_dict["REMM"] = common_tools.format_number_for_report(remm)

                        if "REVEL" in report_columns:
                            revel = variant.INFO.get("REVEL")
                            if revel is None:
                                revel_score = "None"
                            elif isinstance(revel, str) and len(revel.split(",")) > 1:
                                revel_multiple_scores_list = [str(common_tools.format_number_for_report(score, 3)) for score in revel.split(",")]
                                revel_score = ",".join(revel_multiple_scores_list)
                            else:
                                revel_score = common_tools.format_number_for_report(revel, 3)

                            report_data_dict["REVEL"] = revel_score

                        spliceai_columns = ["SPLICEAI_GENE", "SPLICEAI_DS_MAX", "SPLICEAI_TYPE", "SPLICEAI_POS"]
                        if any(column in report_columns for column in spliceai_columns):
                            spliceai_info = variant.INFO.get("SPLICEAI")
                            if spliceai_info is None:
                                if "SPLICEAI_GENE" in report_columns:
                                    report_data_dict["SPLICEAI_GENE"] = "None"
                                if "SPLICEAI_DS_MAX" in report_columns:
                                    report_data_dict["SPLICEAI_DS_MAX"] = "None"
                                if "SPLICEAI_TYPE" in report_columns:
                                    report_data_dict["SPLICEAI_TYPE"] = "None"
                                if "SPLICEAI_POS" in report_columns:
                                    report_data_dict["SPLICEAI_POS"] = "None"
                            else:
                                spliceai_gene = []
                                spliceai_ds_max = []
                                spliceai_type = []
                                spliceai_pos = []
                                # Split SPLICEAI field by gene:
                                for spliceai in spliceai_info.split(","):
                                    spliceai_keys = ["ALLELE", "SYMBOL", "DS_AG", "DS_AL", "DS_DG", "DS_DL", "DP_AG", "DP_AL", "DP_DG", "DP_DL"]
                                    spliceai_values = spliceai.split("|")
                                    # Build dictionary with all fields:
                                    spliceai_dict = dict(zip(spliceai_keys, spliceai_values))
                                    # Record gene:
                                    gene = spliceai_dict["SYMBOL"]
                                    spliceai_gene.append(gene)
                                    # Manage maximum delta score, types and positions:
                                    scores_dict = {key.split("_")[1]: float(value) for (key, value) in spliceai_dict.items() if key.startswith("DS_")}
                                    score_max = max(scores_dict.values())
                                    score_max_types = [key for key, value in scores_dict.items() if value == score_max]
                                    pos_dict = {key.split("_")[1]: value for (key, value) in spliceai_dict.items() if key.startswith("DP_")}
                                    # Recod maximum delta score:
                                    spliceai_ds_max.append(str(score_max))
                                    # Record types and positions
                                    if score_max == 0:
                                        spliceai_type.append("no_consequence")
                                        spliceai_pos.append("no_consequence")
                                    else:
                                        spliceai_types_dict = {"AG": "acceptor_gain", "AL": "acceptor_loss", "DG": "donor_gain", "DL": "donor_loss"}
                                        spliceai_max_types = [spliceai_types_dict[item] for item in score_max_types]
                                        spliceai_type.append(",".join(spliceai_max_types))
                                        spliceai_max_pos = [pos_dict[item] for item in score_max_types]
                                        spliceai_pos.append(",".join(spliceai_max_pos))
                                # Record report columns:
                                if "SPLICEAI_GENE" in report_columns:
                                    report_data_dict["SPLICEAI_GENE"] = "|".join(spliceai_gene)
                                if "SPLICEAI_DS_MAX" in report_columns:
                                    report_data_dict["SPLICEAI_DS_MAX"] = "|".join(spliceai_ds_max)
                                if "SPLICEAI_TYPE" in report_columns:
                                    report_data_dict["SPLICEAI_TYPE"] = "|".join(spliceai_type)
                                if "SPLICEAI_POS" in report_columns:
                                    report_data_dict["SPLICEAI_POS"] = "|".join(spliceai_pos)

                        if "PRIMATEAI" in report_columns:
                            primateai = variant.INFO.get("PRIMATEAI")
                            if primateai is None:
                                primateai_score = "None"
                            elif isinstance(primateai, str) and len(primateai.split(",")) > 1:
                                primateai_multiple_scores_list = [str(common_tools.format_number_for_report(score)) for score in primateai.split(",")]
                                primateai_score = ",".join(primateai_multiple_scores_list)
                            else:
                                primateai_score = common_tools.format_number_for_report(primateai)
                            report_data_dict["PRIMATEAI"] = primateai_score

                        if "AM_SCORE" in report_columns:
                            am_score = variant.INFO.get("AM_SCORE")
                            if am_score is None:
                                am_score = "None"
                            report_data_dict["AM_SCORE"] = am_score

                        if "AM_CLASS" in report_columns:
                            am_class = variant.INFO.get("AM_CLASS")
                            if am_class is None:
                                am_class = "None"
                            report_data_dict["AM_CLASS"] = am_class

                        if "DBSCSNV_ADA" in report_columns:
                            dbscsnv_ada = variant.INFO.get("DBSCSNV_ADA")
                            if dbscsnv_ada is None:
                                dbscsnv_ada = "None"
                            report_data_dict["DBSCSNV_ADA"] = common_tools.format_number_for_report(dbscsnv_ada)

                        if "DBSCSNV_RF" in report_columns:
                            dbscsnv_rf = variant.INFO.get("DBSCSNV_RF")
                            if dbscsnv_rf is None:
                                dbscsnv_rf = "None"
                            report_data_dict["DBSCSNV_RF"] = common_tools.format_number_for_report(dbscsnv_rf)

                        if "GNOMAD_SYN_Z" in report_columns or "GNOMAD_MIS_Z" in report_columns or "GNOMAD_LOF_Z" in report_columns:
                            gnomad_z_score_dict = self.get_gnomad_z_score(annotation_data["gene_name"])

                        if "GNOMAD_SYN_Z" in report_columns:
                            synonymous_z_score = gnomad_z_score_dict["synonymous"]
                            report_data_dict["GNOMAD_SYN_Z"] = synonymous_z_score

                        if "GNOMAD_MIS_Z" in report_columns:
                            missense_z_score = gnomad_z_score_dict["missense"]
                            report_data_dict["GNOMAD_MIS_Z"] = missense_z_score

                        if "GNOMAD_LOF_Z" in report_columns:
                            lof_z_score = gnomad_z_score_dict["lof"]
                            report_data_dict["GNOMAD_LOF_Z"] = lof_z_score
                        
                        if "HPO_DISEASE" in report_columns: 
                            list_hpo_old = self.hpo(annotation_data["gene_name"])
                            report_data_dict["HPO_DISEASE_OLD"] = list_hpo_old
                            list_hpo = self.get_hpo_diseases(annotation_data["gene_name"])
                            report_data_dict["HPO_DISEASE"] = list_hpo

                        if "VARIANT_NUM_CALLED" in report_columns: 
                            report_data_dict["VARIANT_NUM_CALLED"] = variant.num_called

                        if "VARIANT_NUM_HOM_REF" in report_columns:
                            report_data_dict["VARIANT_NUM_HOM_REF"] = variant.num_hom_ref

                        if "VARIANT_NUM_HET" in report_columns:
                            report_data_dict["VARIANT_NUM_HET"] = variant.num_het

                        if "VARIANT_NUM_HOM_ALT" in report_columns:
                            report_data_dict["VARIANT_NUM_HOM_ALT"] = variant.num_hom_alt

                        if "VARIANT_NUM_UNKNOWN" in report_columns:
                            report_data_dict["VARIANT_NUM_UNKNOWN"] = variant.num_unknown

                        if "GENE_LIST" in report_columns:
                            file_names_list = self.add_gene_list(annotation_data["gene_name"])
                            report_data_dict["GENE_LIST"] = file_names_list

                        if "RECESSIVE" in report_columns and "PHASING" in report_columns:
                            # assumes the variant is in the list of recessive variants
                            try:
                                report_data_dict["RECESSIVE"] = recessive_variants[gnomad_format]["recessive"]
                                report_data_dict["PHASING"] = recessive_variants[gnomad_format]["phasing"]
                            # catches error if it is not
                            except KeyError:
                                report_data_dict["RECESSIVE"] = "False"
                                report_data_dict["PHASING"] = "NA"

                        if "DENOVO" in report_columns:
                            if sample_id == proband_id:
                                if report_data_dict["GNOMAD_POS"] in gnomad_position_denovo and report_data_dict["SEGDUP"] == "False":
                                    report_data_dict["DENOVO"] = "True"
                                else:
                                    report_data_dict["DENOVO"] = "False"
                            elif sample_id == mother_id or sample_id == father_id:
                                report_data_dict["DENOVO"] = "None"

                        if "BAYESDEL_NO_AF_SCORE" in report_columns:
                            bayesdel_no_af = variant.INFO.get("BAYESDEL_NOAF")
                            if bayesdel_no_af is None:
                                bayesdel_no_af = "None"
                            report_data_dict["BAYESDEL_NO_AF_SCORE"] = bayesdel_no_af

                        if "BAYESDEL_ADD_AF_SCORE" in report_columns:
                            bayesdel_add_af = variant.INFO.get("BAYESDEL_ADDAF")
                            if bayesdel_add_af is None:
                                bayesdel_add_af = "None"
                            report_data_dict["BAYESDEL_ADD_AF_SCORE"] = bayesdel_add_af

                        if "METASVM_SCORE" in report_columns:
                            metasvm_score = variant.INFO.get("METASVM_SCORE")
                            report_data_dict["METASVM_SCORE"] = self.format_meta_scores_columns(meta_score=metasvm_score)

                        if "METALR_SCORE" in report_columns:
                            metalr_score = variant.INFO.get("METALR_SCORE")
                            report_data_dict["METALR_SCORE"] = self.format_meta_scores_columns(meta_score=metalr_score)

                        if "METARNN_SCORE" in report_columns:
                            metarnn_score = variant.INFO.get("METARNN_SCORE")
                            report_data_dict["METARNN_SCORE"] = self.format_meta_scores_columns(meta_score=metarnn_score)

                        report_data = pd.concat([report_data, pd.DataFrame([report_data_dict])],ignore_index=True)

            if "GENE_COUNT" in report_columns:
                report_data["GENE_COUNT"] = self.gene_count(report_data, gene_count_dict)

            # HGVS correction if transcripts file is present in the config file
            if self.hgvs_transcripts:
                report_data_final = common_tools.correct_hgvs_transcripts(transcripts_file=self.hgvs_transcripts, report_data=report_data)
            else:
                report_data_final = report_data

            report_data_final = report_data_final.fillna(value="None")
            report_data_final = common_tools.parse_zeroes_for_report(report_data_final)
            report_data_final.to_csv(self.output_directory + output_file, index=False, sep="\t")
        report_summary.close()

    def generate_report_variants_command(self, config_file):
        """
            Generates a bash command to execute run_report_variants.py
            Args:
                config_file: the path to the config file used for the analysis
            Returns:
                a dictionary containing the command and the job ID and a list of output files created
        """
        job_id = "report_variants"

        # Create list of output_files
        summary_report_output_file = self.output_directory + "report_summary.tsv"
        output_files = [summary_report_output_file]

        for samples_to_report in self.reports_list:
            output_files.append(self.output_directory + self.get_report_output_file_name(samples_to_report))

        command = f"""
python {self.ngstk_path}/run_report_variants.py -vcf {self.vcf_file} -cf {config_file} -o {self.output_directory} -vc {self.variant_caller} -at {self.annotation_tool}"""

        if self.reporting_file != "":
           command += f" -rf {self.reporting_file}"
        else:
            command += "\n"

        full_command = common_tools.generate_full_command(job_id, command)
        
        step_data = {
            "command": full_command,
            "job_id": job_id,
            "output": output_files
        }
        return step_data

    def apply_labkey_filter_report_variants(self, labkey_filter_list, report_variants_file):
        """
            Generates a bash command to apply the filter_report_variants.py script located in the labkey folder
            This script reads a list of variants on the CMDL LabKey instance and filters variants according to their alternate allele frequency
            Args:
                labkey_filter_list: the suffix of the list to use to filter the report. Lists are located here: CMDL/MolGenet/Bioinfo and have the prefix filter_variants_
                report_variants_file: the report to filter. Must have a GNOMAD_POS and SAMPLE_ALT_AF column
            Returns:
                a dictionary containing the command, the job ID and the output file created
        """
        job_id = "labkey_filter_report_variants"
        output_file = self.output_directory + "report_variants_filtered.tsv"
        command = f"""
echo "JOB: {job_id}"
python {self.ngstk_path}/labkey/filter_report_variants.py -p {labkey_filter_list} -r {report_variants_file} -o {output_file}
"""
        command += common_tools.generate_step_verification_segment(job_id)
        step_data = {
            "command": command,
            "job_id": job_id,
            "output": output_file
        }
        return step_data

    def generate_report_mei(self, vcf_file, gnomad_data_string):
        """
            Generates a TSV report from the merged MEI VCF
            Args:
                vcf_file: path to the VCF file
                gnomad_data_string: a string that is either True or False that asserts if the GNOMAD file for the z scores if found in the config file
        """
        report_columns = [
            "SAMPLE_ID",
            "POSITION",
            "MEI_FAMILY",
            "VARIANT_GENE_NAME",
            "VARIANT_ANNOTATION",
            "VARIANT_HGVS_C",
            "VARIANT_HGVS_P",
            # "VARIANT_EXON_INTRON",
            "VARIANT_CODING",
            "VARIANT_LOF",
            "GNOMAD_SYN_Z",
            "GNOMAD_MIS_Z",
            "GNOMAD_LOF_Z",
            "VARIANT_NUM_CALLED",
            "INSERTION_DIRECTION",
            "CLIPPED_READS_IN_CLUSTER",
            "ALIGNMENT_SCORE",
            "ALIGNEMENT_PERCENT_LENGTH",
            "ALIGNEMENT_PERCENT_IDENTITY",
            "CLIPPED_SEQUENCE",
            "CLIPPED_SIDE",
            "START_IN_MEI",
            "STOP_IN_MEI",
            "POLYA_POSITION",
            "POLYA_SEQ",
            "POLYA_SUPPORTING_READS",
            "TSD",
            "TSD_LENGTH",
        ]
        self.parse_reference_resources(report_columns)
        report_data = pd.DataFrame(columns=report_columns)
        vcf = cyvcf2.VCF(vcf_file)
        gnomad_data = eval(gnomad_data_string)

        for variant in vcf:
            me_info = variant.INFO.get("MEINFO").split(",")[0].split("_")
            position = me_info[0]
            family = me_info[1]
            direction = me_info[2]

            samples = variant.INFO.get("SAMPLE").split(",")
            for sample in samples:
                alignment_info = variant.INFO.get("ALIGN").split("_")
                alignment_score = alignment_info[0]
                alignment_length = alignment_info[1]
                alignment_identity = alignment_info[2]
                start_in_mei = alignment_info[3]
                stop_in_mei = alignment_info[4]

                clip_info = variant.INFO.get("CLIP").split("_")
                clipped_reads = clip_info[0]
                clipped_side = clip_info[1]

                polya_info = variant.INFO.get("POLYA").split("_")
                polya_position = polya_info[0]
                polya_seq = polya_info[1]
                polya_reads = polya_info[2]

                tsd_info = variant.INFO.get("TSD").split("_")
                tsd = tsd_info[0]
                tsd_length = tsd_info[1]

                sample_row = {"SAMPLE_ID": sample, "POSITION": "chr" + position, "MEI_FAMILY": family}
                if self.annotation_tool == "snpeff":
                    annotation_data = common_tools.parse_snpeff_annotation(variant.INFO.get("ANN"))
                    sample_row["VARIANT_GENE_NAME"] = annotation_data["gene_name"]
                    sample_row["VARIANT_ANNOTATION"] = annotation_data["annotation"]
                    sample_row["VARIANT_HGVS_C"] = annotation_data["hgvs_c"]
                    sample_row["VARIANT_HGVS_P"] = annotation_data["hgvs_p"]
                    # sample_row["VARIANT_EXON_INTRON"] = annotation_data["exon"]
                    sample_row["VARIANT_CODING"] = annotation_data["coding"]
                    sample_row["VARIANT_LOF"] = annotation_data["lof"]
                    if gnomad_data:
                        gnomad_z_score_dict = self.get_gnomad_z_score(annotation_data["gene_name"])
                        sample_row["GNOMAD_SYN_Z"] = gnomad_z_score_dict["synonymous"]
                        sample_row["GNOMAD_MIS_Z"] = gnomad_z_score_dict["missense"]
                        sample_row["GNOMAD_LOF_Z"] = gnomad_z_score_dict["lof"]
                sample_row["VARIANT_NUM_CALLED"] = len(samples)
                sample_row["INSERTION_DIRECTION"] = direction
                sample_row["CLIPPED_READS_IN_CLUSTER"] = clipped_reads
                sample_row["ALIGNMENT_SCORE"] = common_tools.format_number_for_report(alignment_score)
                sample_row["ALIGNEMENT_PERCENT_LENGTH"] = common_tools.format_number_for_report(alignment_length)
                sample_row["ALIGNEMENT_PERCENT_IDENTITY"] = common_tools.format_number_for_report(alignment_identity)
                sample_row["CLIPPED_SEQUENCE"] = variant.ALT[0]
                sample_row["CLIPPED_SIDE"] = clipped_side
                sample_row["START_IN_MEI"] = start_in_mei
                sample_row["STOP_IN_MEI"] = stop_in_mei
                sample_row["POLYA_POSITION"] = polya_position
                sample_row["POLYA_SEQ"] = polya_seq
                sample_row["POLYA_SUPPORTING_READS"] = polya_reads
                sample_row["TSD"] = tsd
                sample_row["TSD_LENGTH"] = tsd_length

                report_data = pd.concat([report_data, pd.DataFrame([sample_row])], ignore_index=True)

        # HGVS correction if transcripts file is present in the config file
        if self.hgvs_transcripts:
            report_data_final = common_tools.correct_hgvs_transcripts(transcripts_file=self.hgvs_transcripts, report_data=report_data)
        else:
            report_data_final = report_data

        report_data_final = report_data_final.fillna(value="None")
        report_data_final = common_tools.parse_zeroes_for_report(report_data_final)
        report_data_final.to_csv(self.output_directory + "report_mei.tsv", index=False, sep="\t")

    def generate_report_mei_command(self, vcf_file, config_file):
        """
            Generates a bash command to execute the generate_report_mei function
            Args:
                vcf_file: the path to the VCF file to report
                config_file: the path to the config file
            Returns:
                a dictionary containing the command, the job ID and the output file created
        """

        job_id = "report_mei"
        output_file = self.output_directory + "report_mei.tsv"
        config = configparser.ConfigParser()
        config.read(config_file)
        try:
            reference_resources = {
                "gnomad_metrics": common_tools.get_full_path(config.get("REPORT_VARIANTS", "GNOMAD_METRICS")),
            }
            gnomad_data = True
        except:
            print("WARNING: Gnomad Z scores will be omitted from the MEI report")
            reference_resources = {}
            gnomad_data = False

        try:
            hgvs_transcripts = config.get("REPORT_VARIANTS", "HGVS_TRANSCRIPTS")
        except:
            hgvs_transcripts = ""

        command = f"""
echo "JOB: {job_id}"
python -c "import sys;
sys.path.append(\\"{self.ngstk_path}\\")
import reportVariants;
command_generator = reportVariants.ReportVariants(\\"{self.vcf_file}\\", \\"{self.variant_caller}\\", \\"{self.output_directory}\\", \\"{self.annotation_tool}\\", reference_resources=\\"{reference_resources}\\", hgvs_transcripts=\\"{hgvs_transcripts}\\")
command_generator.generate_report_mei(\\"{vcf_file}\\", \\"{gnomad_data}\\")"
"""
        command += common_tools.generate_step_verification_segment(job_id)

        step_data = {
            "command": command,
            "job_id": job_id,
            "output": output_file
        }
        return step_data
