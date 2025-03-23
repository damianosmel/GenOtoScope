from genotoscope.utils import create_dir, is_gsvar_file_extension, update_assignment, extract_all_variant_types, \
	split_df4threads, is_var_type4PVS1, sec2hour_min_sec, extract_rule_result, is_var2exclude, format_dict2str, \
	deformat_str2dict, load_hugo_genes_df, deactivate_lower_strength_rules, parse_maxentscan_pred, parse_dbscsnv_pred, \
	is_variant_pathogenic, is_variant_conserved, aggregate_patho_predictions, aggregate_splicing_predictions
from genotoscope.variants.RefineLossOfFunction import RefineLossOfFunction
from genotoscope.variants.ProcessClinVar import ProcessClinVar
from genotoscope.variants.AssignBP3 import AssignBP3
from genotoscope.variants.AssignPM4 import AssignPM4
from genotoscope.variants.AssignPM1 import AssignPM1
from genotoscope.variants.AssignACMG import AssignACMG
from genotoscope.variants.VariantInfo import VariantInfo

from os.path import join, exists, basename, dirname
from os import scandir, remove
from pandas import read_csv, read_excel, unique, concat
import logging
from datetime import datetime
from timeit import default_timer as timer

from math import isnan
from multiprocessing.dummy import Pool
from threading import get_ident

from pyensembl import EnsemblRelease
from pybedtools import BedTool
import vcf


class ClassifyVariants:
	"""
	Class to process GSVar file and assign classification type per variant on each GSvar file,
	following the ACMG/AMP guidelines doi: 10.1002/humu.23630
	"""

	def __init__(self, data_path, gsvar_folder, genes_inheritance_file, thresholds_inheritance_file,
	             clinvar_root, clinvar_version, clinvar_file, clinvar_stars_file, min_review_stars,
	             beds_root, critical_prot_version, critical_prot_regions_file, critical_prot_regions_no_benign_file,
	             gnomAD_version, clinical_exons_file,
	             uniprot_root, uniprot_version, uniprot_domains_file, uniprot_repeat_file,
	             pheno_annot_root, phenotype_relevant_exons_file, pm1_regions_file,
	             hugo_genes_file, pheno_high_af_variants_file, temp_path,
	             out_path, filter_hl_genes, min_pathogenicity, num_threads,
	             logging_root, logging_level):
		"""
		ClassifyVariants constructor

		Parameters
		----------
		data_path : str
			absolute path containing annotation data to classify variant
		gsvar_folder : str
			folder containing GSvar files
		genes_inheritance_file : str
			genes inheritance mode filename
		thresholds_inheritance_file : str
			rules thresholds by inheritance mode
		clinvar_root : str
			ClinVar folder root
		clinvar_version : str
			ClinVar version
		clinvar_file : str
			clinvar file name
		clinvar_stars_file : str
			clinvar review star tsv filename
		min_review_stars : int
			minimum review stars for a ClinVar variant to be used
		beds_root : str
			beds files root
		critical_prot_version : str
			critical protein version (clinvar + uniprot)
		critical_prot_regions_file : str
			critical protein regions bed filename
		critical_prot_regions_no_benign_file : str
			critical protein regions without benign variants (used in PM1)
		gnomAD_version : str
			gnomAD version
		clinical_exons_file : str
			clinical significant exons bed filename
		uniprot_root : str
			uniprot annotations root
		uniprot_version : str
			uniprot version
		uniprot_domains_file : str
			uniprot domain annotation bed filename
		uniprot_repeat_file : str
			uniprot repeat annotation bed filename
		pheno_annot_root : str
			phenotype relevant annotation root
		phenotype_relevant_exons_file : str
			phenotype relevant exons bed filename
		pm1_regions_file : str
			protein regions firing PM1 as specified by DOI: 10.1002/humu.23630
		hugo_genes_file : str
			hugo genes filename
		pheno_high_af_variants_file : str
			file containing hearing loss variants that occur in high AF
		temp_path : str
			absolute path where temporary folder will be placed
		out_path : str
			absolute path where output will be placed
		filter_hl_genes : bool
			classify only known HL genes (True), otherwise classify all presented genes
		min_pathogenicity : float
			pathogenicity probability threshold to extract variants with ACMG class equal to 3 (alongside with variants with ACMG class of 4/5)
		num_threads : int
			number of threads for parallel execution of rules
		logging_root : logging
			initialized root by main
		logging_level : str
			logging level
		Returns
		-------
		None
		"""
		### input paths ###
		self.data_path = data_path
		# change on 21st Sept.
		# self.gsvar_path = join(data_path, gsvar_folder)
		self.gsvar_path = gsvar_folder
		self.genes_inheritance_file = genes_inheritance_file
		self.thresholds_inheritance_file = thresholds_inheritance_file
		self.current_sample_name = None

		### output paths ###
		self.out_path = out_path
		create_dir(out_path)
		self.temp_path = temp_path
		create_dir(temp_path)

		### PyEnsembl ###
		# release 75 uses human reference genome GRCh37
		self.ensembl_data = EnsemblRelease(75)

		### ClinVar ###
		self.clinvar_root = join(self.data_path, clinvar_root)
		self.clinvar_path = join(self.clinvar_root, clinvar_version)
		self.clinvar_file = clinvar_file
		self.min_review_stars = min_review_stars
		self.clinvar_stars_df = self.load_clinvar_review_stars_df(clinvar_stars_file)

		### Annotation beds ###
		self.beds_root = join(self.data_path, beds_root)
		# Protein critical regions (clinvar + uniprot)
		critical_prot_path = join(self.beds_root, critical_prot_version)
		self.critical_prot_regions_file = join(critical_prot_path, critical_prot_regions_file)
		self.critical_prot_regions_no_benign_file = join(critical_prot_path, critical_prot_regions_no_benign_file)
		# Clinical significant exons (gnomAD)
		clinical_exons_path = join(self.beds_root, gnomAD_version)
		self.clinical_exons_file = join(clinical_exons_path, clinical_exons_file)

		### Phenotype relevant ###
		pheno_annot_path = join(self.beds_root, pheno_annot_root)
		self.phenotype_exons_file = join(pheno_annot_path, phenotype_relevant_exons_file)
		self.pm1_regions_file = join(pheno_annot_path, pm1_regions_file)

		### UniProt annotations ###
		uniprot_root = join(self.data_path, uniprot_root)
		self.uniprot_path = join(uniprot_root, uniprot_version)
		self.uniprot_domains_file = join(self.uniprot_path, uniprot_domains_file)
		self.uniprot_repeat_file = join(self.uniprot_path, uniprot_repeat_file)
		self.repeats_no_domains_file = None

		### Inheritance info ###
		self.genes_inheritance_df = None
		self.digram2inheritance = {"AR": "autosomal recessive", "AD": "autosomal dominant", "XL": "X-linked",
		                           "MT": "mitochondrial"}
		self.inheritance2digram = {"autosomal recessive": "AR", "autosomal dominant": "AD", "X-linked": "XL",
		                           "mitochondrial": "MT"}

		### gnomAD ###
		self.position2gnomAD_sub = {0: "AFR", 1: "AMR", 2: "EAS", 3: "NFE", 4: "SAS"}

		### Hugo genes ###
		self.hugo_genes_df = load_hugo_genes_df(self.data_path, hugo_genes_file)

		### Phenotype variants that occur in high AF ###
		if exists(join(self.data_path, pheno_high_af_variants_file)):
			self.pheno_high_af_variants = vcf.Reader(filename=join(self.data_path, pheno_high_af_variants_file),
			                                         compressed=True, encoding='utf-8')
		else:
			self.pheno_high_af_variants = None

		### Filter genes of variants ###
		self.classify_only_hl_genes = filter_hl_genes

		### Variant summarization ###
		self.min_pathogenicity = min_pathogenicity

		### Number of threads ###
		self.num_threads = num_threads

		### Logging ###
		self.init_logger(logging_root, "GenOtoScope_Classify", logging_level)

	def init_logger(self, logging_root, logger_name, logging_level):
		"""
		Initiliaze logger
		Credits: Credits: https://stackoverflow.com/questions/20240464/python-logging-file-is-not-working-when-using-logging-basicconfig

		Parameters
		----------
		logging_root : logging
			logging root object
		logger_name : str
			logger name
		logging_level : str
			logging level

		Returns
		-------
		None
		"""
		self.logger = logging_root.getLogger(logger_name)
		self.logger.setLevel(logging_level)

		# create file handler
		fh = logging_root.FileHandler(join(self.out_path,
		                                   'genotoscope_classify_' + datetime.today().strftime("%d_%m_%Y") + '.log'))
		fh.setLevel(logging_level)
		# create console handler
		ch = logging.StreamHandler()
		ch.setLevel(logging.ERROR)

		# create format
		formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
		fh.setFormatter(formatter)
		ch.setFormatter(formatter)
		# add handler to the logger
		self.logger.addHandler(fh)
		self.logger.addHandler(ch)

	def prepare_annotation_files(self):
		"""
		Prepare bed annotation files needed for ACMG rules

		Returns
		-------
		str
			repeats with no overlapping domain file name
		"""
		# prepare repeats without known function
		self.logger.info("Preparing repeats without known function")
		[proteome_id, species_taxid, annot_type, ref_version, _] = basename(self.uniprot_domains_file).strip(
			".bed").split("_")
		repeats_no_domains_filename = "_".join([proteome_id, species_taxid, "repeats_no_domains", ref_version]) + ".bed"
		if not exists(join(self.uniprot_path, repeats_no_domains_filename)):
			self.logger.info("Preparing annotation bed file: Uniprot repeats with no overlapping domain")
			uniprot_repeats = BedTool(self.uniprot_repeat_file)
			uniprot_domains = BedTool(self.uniprot_domains_file)
			repeats_no_domains = uniprot_repeats.subtract(uniprot_domains, s=True).sort()
			repeats_no_domains.saveas(join(self.uniprot_path, repeats_no_domains_filename))
		return join(self.uniprot_path, repeats_no_domains_filename)

	def load_clinvar_review_stars_df(self, clinvar_review_stars_file):
		"""
		Load clinvar review stars into dataframe

		Parameters
		----------
		clinvar_stars_file : str
			clinvar stars file name

		Returns
		-------
		pandas.DataFrame
			loaded clinvar stars data frame
		"""
		clinvar_stars_df = read_csv(join(self.clinvar_root, clinvar_review_stars_file), sep="\t", header=0)
		clinvar_stars_df = clinvar_stars_df.rename(lambda x: x.strip().replace(' ', '_'), axis='columns')
		clinvar_stars_df['Review_status'] = clinvar_stars_df['Review_status'].str.strip()
		clinvar_stars_df['Number_of_gold_stars'] = clinvar_stars_df['Number_of_gold_stars'].str.strip()
		return clinvar_stars_df

	def process_column_line(self, gsvar_file, is_single_file_analysis=False):
		"""
		Process only the column line in order to be able to parsed by pandas

		Feature: ignore no Unicode character
		Credits: https://stackoverflow.com/questions/22216076/unicodedecodeerror-utf8-codec-cant-decode-byte-0xa5-in-position-0-invalid-s

		Parameters
		----------
		gsvar_file : str
			GSvar file name
		is_single_file_analysis : bool
			A single GSvar file will be classified (True), otherwise a folder with GSvar files will be

		Returns
		-------
		str
			updated GSvar file name
		"""
		# Read lines and convert the #chr to chr
		self.logger.info("Convert column name from #chr to chr")
		num_changed_lines = 0

		if is_single_file_analysis:
			gsvar_path = dirname(self.gsvar_path)
			gsvar_file = basename(gsvar_file)
		else:
			gsvar_path = self.gsvar_path
		updated_gsvar_file = gsvar_file.split(".GSvar")[0] + "_column_line.GSvar"
		with open(join(gsvar_path, gsvar_file), 'rb') as gsvar_in, open(join(self.temp_path, updated_gsvar_file),
		                                                                'w') as gsvar_out:
			for line in gsvar_in:
				line = line.decode(errors='ignore')
				if line[:4] == "#chr":
					line = line[1:]
					num_changed_lines += 1
				gsvar_out.write(line.strip() + "\n")
		assert num_changed_lines == 1, "AssertionError: only the column line should be changed"  # assert that only one line has been changed
		self.logger.info("Save converted GSvar file at %s", self.temp_path)
		return updated_gsvar_file

	def find_strictest_inheritance_mode(self, rule_row):
		"""
		Function to find the strictest inheritance mode between the possible ones

		Parameters
		----------
		rule_row : pandas.Series
			rule row

		Returns
		-------
		str
			initials of strictest inheritance mode ("AD" or "AR")
		"""
		if str(rule_row["threshold_direction"]) == "greater" or str(
				rule_row["threshold_direction"]) == "greater_than_or_equal":
			# rule threshold direction is MAF > or >=
			# thus strictest inheritance mode is the mode with maximum threshold
			if float(rule_row["AD"]) >= float(rule_row["AR"]):
				return "AD"
			else:
				return "AR"
		elif str(rule_row["threshold_direction"]) == "lower" or str(
				rule_row["threshold_direction"]) == "lower_than_or_equal":
			# rule threshold direction is MAF < or <=
			# thus strictest inheritance mode is the mode with lowest threshold
			if float(rule_row["AD"]) <= float(rule_row["AR"]):
				return "AD"
			else:
				return "AR"
		else:
			# rule threshold direction is MAF = or unknown direction
			return "unknown"

	def load_genes_inheritance_df(self):
		"""
		Load genes inheritance

		Returns
		-------
		None
		"""
		self.logger.info(
			"Load genes inheritance mode from {}".format(join(self.data_path, self.genes_inheritance_file)))
		self.genes_inheritance_df = read_excel(join(self.data_path, self.genes_inheritance_file), index_col=0, header=0)

	def load_thresholds_inheritance_df(self):
		"""
		Load thresholds per inheritance mode dataframe

		Returns
		-------
		None
		"""
		self.logger.info("Load rule thresholds per inheritance mode from {}".format(
			join(self.data_path, self.thresholds_inheritance_file)))
		self.inheritance_thresholds_df = read_csv(join(self.data_path, self.thresholds_inheritance_file), index_col=0)

		self.inheritance_thresholds_df["strictest_inheritance"] = self.inheritance_thresholds_df.apply(
			self.find_strictest_inheritance_mode, axis=1)
		self.logger.debug("Inheritance thresholds:\n {}".format(self.inheritance_thresholds_df.head()))

	def load_gsvar_df(self, gsvar_dir_parent, gsvar_file):
		"""
		Load gsvar dataframe
		Credits: https://stackoverflow.com/questions/14158868/python-skip-comment-lines-marked-with-in-csv-dictreader

		Parameters
		----------
		gsvar_dir_parent : str
			directory parent of GSvar
		gsvar_file : str
			GSvar file name

		Returns
		-------
		pandas.DataFrame
			loaded gsvar data frame
		bool
			loaded gsvar contains 0 variant rows: True, otherwise False
		"""
		self.logger.info("Load converted GSvar file")
		### ### ###
		# get the directory parent of the GSvar file
		### ### ###
		if gsvar_dir_parent == "temp":
			dir_parent = self.temp_path
		elif gsvar_dir_parent == "output":
			dir_parent = self.out_path
		else:
			dir_parent = None
		gsvar_df = read_csv(join(dir_parent, gsvar_file),
		                    sep='\t',  # field separator
		                    comment='#',  # comment
		                    header=0,  # get header
		                    # index_col=[0, 1, 2],  # set as index the chr, start stop
		                    skipinitialspace=True,
		                    skip_blank_lines=True,
		                    error_bad_lines=False,
		                    warn_bad_lines=True
		                    )
		is_gsvar_empty = False
		if gsvar_df.shape[0] != 0:
			# add the 'chr' prefix in front of each coordinate, if the vcf does not contain it
			if str(gsvar_df.loc[0, "chr"])[0:3] != "chr":
				self.logger.debug("Adding chr prefix on chromosome column")
				gsvar_df["chr"] = gsvar_df["chr"].apply(lambda chr_num: "chr" + str(chr_num))
			self.logger.debug("Loaded GSvar head:\n {}".format(gsvar_df.head(5)))
			# keep the ordered list of input GSvar file
			self.current_gsvar_columns = list(gsvar_df.columns)
		else:
			self.logger.info("Loaded GSvar file {} contains 0 variant rows.".format(gsvar_file))
			is_gsvar_empty = True
		return gsvar_df, is_gsvar_empty

	def examine_allele_freq(self, variants):
		"""
		Examine the allele frequency of the variants from gnomAD

		Parameters
		----------
		variants : pandas.Series
			variant row to be examined for allele freq

		Returns
		-------
		pandas.Series
			variant row with assigned class based on allele freq
		"""
		# examine gnomAD MAF, then decide for discarding or not the variant
		if isnan(float(variants["gnomAD"])):
			af_class = "AF: keep variant"
			af_comment = "AF: variant not present in gnomAD, assume AF=0"
		elif float(variants["gnomAD"]) > 0.01:
			af_class = "AF: discard variant"
			af_comment = "AF: high MAF"
		else:
			af_class = "AF: keep variant"
			af_comment = "AF: MAF<= 0.01"

		if not isnan(float(variants["gnomAD"])):
			# examine if variant should be excluded from MAF criterion
			if is_var2exclude(self.pheno_high_af_variants,
			                  VariantInfo(self.current_sample_name, str(variants["gene"]).split(",")[0],
			                              str(variants["variant_type"]),
			                              None,
			                              str(variants["chr"]),
			                              str(variants["start"]), str(variants["end"]), str(variants["ref"]),
			                              str(variants["obs"]))):
				af_class = "AF: keep variant"
				af_comment = "AF: variant in list of phenotype relevant and pathogenic variants with high AF"
		# update assigned class and comment from gnomAD MAF criterion
		variants["ACMG_rules"] = update_assignment(variants["ACMG_rules"], af_class)
		variants["ACMG_rules_comment"] = update_assignment(variants["ACMG_rules_comment"], af_comment)
		return variants

	@staticmethod
	def var_pos2string(variant_row):
		"""
		Convert variant position to string

		Parameters
		----------
		variant_row : pandas.Series
			variant row

		Returns
		-------
		str
			variant position
		"""
		return "variant position: " + " ".join([variant_row["chr"], str(variant_row["start"]), str(variant_row["end"])])

	def examine_classified_class(self, variants):
		"""
		Examine classified class (column: classification)

		Parameters
		----------
		variants : pandas.Series
			variant row to be examined for the classified class

		Returns
		-------
		pandas.Series
			variant row with assigned class based on classified class
		"""
		# examine if the variant is already classified by hospital db (NGSD)

		if "classification" not in variants.index.values:
			# check if NGSD classification is not found in the column
			ngsd_class = "NGSD: not classified category"
			ngsd_comment = "NGSD: not available field"
		elif isinstance(variants["classification"],float):
			if float(variants["classification"]) < 3:
				ngsd_class = "NGSD: discard variant"
				ngsd_comment = "NGSD: classification <3"
			else:
				ngsd_class = "NGSD: already classified category={}".format(variants["classification"])
				ngsd_comment = "NGSD: already classified variant"
		elif isinstance(variants["classification"],str):
			ngsd_class = "NGSD: already classified with string class={}".format(variants["classification"])
			ngsd_comment = "NGSD: string class as classification"
		else:
			ngsd_class = "NGSD: not classified category"
			ngsd_comment = "NGSD: database does not contain classified category for variant"
		# update assigned class and comment
		variants["ACMG_rules"] = ngsd_class
		variants["ACMG_rules_comment"] = ngsd_comment
		return variants

	def examine_pp3_bp4(self, variants):
		"""
		Examine REVEL and MaxEntScan scores to assess PP3 and BP4:
		PP3: "Multiple lines of computational evidence support a deleterious effect on the gene or gene product
		(conservation, evolutionary, splicing impact, etc.)"
		BP4: "Multiple lines of computational evidence suggest no impact on gene or gene product (conservation, evolutionary,
		splicing impact, etc.)"

		MaxEntScan_ratio = (observed_score - reference_score) / reference_score
		is_pathogenic = Revel or CADD score
		if variant_type == any exonic other than 'synonymous':
			if is_pathogenic == False and PhyloP >= 1.6:
				BP4 = True
			else if is_pathogenic == True and PhyloP < 1.6:
				PP3 = True
		else if variant_type == synonymous:
			BP4 = NA
			if is_pathogenic == True and PhyloP < 1.6:
				PP3 = True
		else if variant_type == splice region or splice acceptor or splice donor:
			if MaxEntScan_ratio >= 0.15 and PhyloP >= 1.6:
				PP3 = True
			BP4 = False, BP4 is used for missense variants (BP7 for splicing variants)
		else:
			# not applicable for this variant type
			BP4 = False
			PP3 = False

		Parameters
		----------
		variants :  pandas.Series
			variant row to examine for PP3 and BP4

		Returns
		-------
		pandas.Series
			variant row with assigned PP3 and BP4 rules
		"""
		self.logger.debug("Examine PP3 and BP4 rules for variant: {}".format(ClassifyVariants.var_pos2string(variants)))
		current_variant_types = extract_all_variant_types(str(variants["variant_type"]))
		intronic_types = {"splice_region", "splice_acceptor", "splice_donor"}
		exonic_types = {"missense", "start_lost", "stop_gained", "stop_lost", "frameshift", "inframe_insertion",
		                "inframe_deletion"}
		synonymous_type = {"synonymous"}
		if len(exonic_types.intersection(set(current_variant_types))) > 0:
			### ### ### ### ### ### ### ### ### ###
			# get REVEL and CADD predictions      #
			# calculate the average of all existing scores #
			### ### ### ### ### ### ### ### ### ###
			patho_predictions = {'REVEL': 0, 'CADD': 0}
			for predictor in patho_predictions.keys():
				patho_predictions[predictor], _ = is_variant_pathogenic(variants, predictor)

			is_var_pathogenic = aggregate_patho_predictions(patho_predictions)
			is_var_conserved = is_variant_conserved(variants)
			self.logger.debug("Variant is predicted as: {}".format(is_var_pathogenic))
			self.logger.debug("Variant is conserved?: {}".format(is_var_conserved))
			if is_var_conserved == 1 and is_var_pathogenic == "pathogenic":
				PP3_class = "PP3: True"
				PP3_comment = "PP3: " + " aggregation of REVEL and CADD predicts pathogenic mutation and PhyloP predicts conserved site"
				BP4_class = "BP4: False"
				BP4_comment = "BP4: " + " aggregation of REVEL and CADD predicts pathogenic mutation and PhyloP predicts conserved site"
			elif is_var_conserved == 0 and is_var_pathogenic == "benign":
				PP3_class = "PP3: False"
				PP3_comment = "PP3: " + " aggregation of REVEL and CADD predicts benign mutation and PhyloP predicts no conserved site"
				BP4_class = "BP4: True"
				BP4_comment = "BP4: " + " aggregation of REVEL and CADD predicts benign mutation and PhyloP predicts no conserved site"
			else:
				BP4_class, PP3_class = "BP4: False", "PP3: False"
				if is_var_pathogenic == "pathogenic":
					patho_comment = "aggregation of REVEL and CADD predicts pathogenic mutation"
				elif is_var_pathogenic == "benign":
					patho_comment = "aggregation of REVEL and CADD predicts benign mutation"
				else:
					patho_comment = "no available REVEL nor CADD predictions to support evidence"
				if is_var_conserved == 1:
					conserved_comment = "PhyloP predicts conserved site"
				elif is_var_conserved == 0:
					conserved_comment = "PhyloP predicts no conserved site"
				else:
					conserved_comment = "no available conservation prediction to support evidence"
				PP3_comment = " ".join(["PP3:", patho_comment, "and", conserved_comment])
				BP4_comment = " ".join(["BP4:", patho_comment, "and", conserved_comment])

		elif len(intronic_types.intersection(set(current_variant_types))) > 0:
			### ### ### ###
			# splice region#
			### ### ### ###
			splicing_predictions = {"MaxEntScan": None, "dbscSNV": None}
			if "MaxEntScan" in variants:
				splicing_predictions["MaxEntScan"] = parse_maxentscan_pred(variants)
			if "dbscSNV" in variants:
				splicing_predictions["dbscSNV"] = parse_dbscsnv_pred(variants)

			### ### ### ### ### ### ###
			# assign PP3 and BP4 rules#
			### ### ### ### ### ### ###
			self.logger.debug("splicing predictions: {}".format(splicing_predictions))
			if (splicing_predictions["MaxEntScan"] is None and splicing_predictions["dbscSNV"] is None):
				# no splicing predictor available
				PP3_class, BP4_class = "PP3: False", "BP4: False"
				PP3_comment = "PP3: no available splicing predictor to support evidence"
				BP4_comment = "BP4: no available splicing predictor to support evidence"
			else:
				### ### ### ### ### ### ### ### ### ###
				# get REVEL and CADD predictions      #
				# calculate the average of all existing scores #
				### ### ### ### ### ### ### ### ### ###
				patho_predictions = {'REVEL': 0, 'CADD': 0}
				for predictor in patho_predictions.keys():
					patho_predictions[predictor], _ = is_variant_pathogenic(variants, predictor)
				is_var_pathogenic = aggregate_patho_predictions(patho_predictions)

				### ### ###
				# base conservation #
				### ### ###
				is_var_conserved = is_variant_conserved(variants)

				### ### ###
				# aggregation for splicing impact
				### ### ###
				impacts_splicing = aggregate_splicing_predictions(splicing_predictions)

				self.logger.debug("Variant is classified as: {}".format(is_var_pathogenic))
				self.logger.debug("Is variant conserved?: {}".format(is_var_conserved))
				self.logger.debug("Does variant affect splicing?: {}".format(impacts_splicing))

				if impacts_splicing and is_var_conserved == 1 and is_var_pathogenic == "pathogenic":
					# variant affects splicing and it's conserved
					PP3_class = "PP3: True"
					PP3_comment = "PP3: " + "aggregation of MaxEntScan and dbscSNV predicts splicing impact, " + "aggregation of REVEL and CADD predicts pathogenic impact" + " and PhyloP predicts conserved site"
					BP4_class = "BP4: False"
					BP4_comment = "BP4: " + "aggregation of MaxEntScan and dbscSNV predicts splicing impact, " + "aggregation of REVEL and CADD predicts pathogenic impact" + " and PhyloP predicts conserved site"
				elif not impacts_splicing and is_var_conserved == 0 and is_var_pathogenic == "benign":
					PP3_class = "PP3: False"
					PP3_comment = "PP3: " + "aggregation of MaxEntScan and dbscSNV predicts no splicing impact, " + "aggregation of REVEL and CADD predicts begign impact" + " and PhyloP predicts no conserved site"
					BP4_class = "BP4: True"
					BP4_comment = "BP4: " + "aggregation of MaxEntScan and dbscSNV predicts no splicing impact, " + "aggregation of REVEL and CADD predicts begign impact" + " and PhyloP predicts no conserved site"
				else:
					# variant either does not affect splicing or it's not a conserved base
					PP3_class, BP4_class = "PP3: False", "BP4: False"
					# assess impact
					if impacts_splicing:
						splicing_comment = "aggregation of MaxEntScan and dbscSNV predicts splicing impact"
					else:
						splicing_comment = "aggregation of MaxEntScan and dbscSNV predicts no splicing impact"
					# assess pathogenicity
					if is_var_pathogenic == "pathogenic":
						patho_comment = "aggregation of REVEL and CADD predicts pathogenic impact"
					elif is_var_pathogenic == "benign":
						patho_comment = "aggregation of REVEL and CADD predicts benign impact"
					else:
						patho_comment = "no available REVEL nor CADD predictions to support evidence"
					# assess conservation
					if is_var_conserved == 1:
						conservation_comment = "PhyloP predicts conserved site"
					elif is_var_conserved == 0:
						conservation_comment = "PhyloP predicts no conserved site"
					else:
						conservation_comment = "no evidence for conservation"
					PP3_comment = " ".join(
						["PP3:", splicing_comment, "and", patho_comment, "and", conservation_comment])
					BP4_comment = " ".join(
						["BP4:", splicing_comment, "and", patho_comment, "and", conservation_comment])
		elif len(synonymous_type.intersection(set(current_variant_types))) > 0:
			### ##### ###
			# synonymous #
			### ##### ###

			### ### ###
			# examine only PP3 rule for synonymous variants
			### #### ###
			patho_predictions = {'REVEL': 0, 'CADD': 0}
			for predictor in patho_predictions.keys():
				patho_predictions[predictor], _ = is_variant_pathogenic(variants, predictor)

			is_var_pathogenic = aggregate_patho_predictions(patho_predictions)
			is_var_conserved = is_variant_conserved(variants)
			self.logger.debug("Variant is predicted as: {}".format(is_var_pathogenic))
			self.logger.debug("Variant is conserved?: {}".format(is_var_conserved))
			if is_var_conserved == 1 and is_var_pathogenic == "pathogenic":
				PP3_class = "PP3: True"
				PP3_comment = "PP3: " + " aggregation of REVEL and CADD predicts pathogenic mutation and PhyloP predicts conserved site"
			elif is_var_conserved == 0 and is_var_pathogenic == "benign":
				PP3_class = "PP3: False"
				PP3_comment = "PP3: " + " aggregation of REVEL and CADD predicts benign mutation and PhyloP predicts no conserved site"
			else:
				PP3_class = "PP3: False"
				if is_var_pathogenic == "pathogenic":
					patho_comment = "aggregation of REVEL and CADD predicts pathogenic mutation"
				elif is_var_pathogenic == "benign":
					patho_comment = "aggregation of REVEL and CADD predicts benign mutation"
				else:
					patho_comment = "no available REVEL nor CADD predictions to support evidence"
				if is_var_conserved == 1:
					conserved_comment = "PhyloP predicts conserved site"
				elif is_var_conserved == 0:
					conserved_comment = "PhyloP predicts no conserved site"
				else:
					conserved_comment = "no available conservation prediction to support evidence"
				PP3_comment = " ".join(["PP3:", patho_comment, "and", conserved_comment])
			BP4_class = "BP4: NA"
			BP4_comment = "BP4 rule not applicable"
		else:
			PP3_class = "PP3: NA"
			PP3_comment = "PP3: rule not applicable"
			BP4_class = "BP4: NA"
			BP4_comment = "BP4: rule not applicable"

		### ### ###
		# update assignment class and comment
		### ### ###
		# add the PP3 rule
		variants["ACMG_rules"] = update_assignment(variants["ACMG_rules"], PP3_class)
		variants["ACMG_rules_comment"] = update_assignment(variants["ACMG_rules_comment"], PP3_comment)
		# add the BP4 rule
		variants["ACMG_rules"] = update_assignment(variants["ACMG_rules"], BP4_class)
		variants["ACMG_rules_comment"] = update_assignment(variants["ACMG_rules_comment"], BP4_comment)
		return variants

	def examine_bp7(self, variants):
		"""
		Examine BP7 rule:
		"A synonymous (silent) variant for which splicing prediction algorithms predict no impact to splice consensus
		sequence nor the creation of a new splice site AND the nucleotide is not highly conserved"

		Parameters
		----------
		variants :  pandas.Series
			variant row to examine BP7 rule on

		Returns
		-------
		pandas.Series
			variant row with assigned BP7
		"""
		self.logger.debug("Examine BP7 rule for variant: {}".format(ClassifyVariants.var_pos2string(variants)))
		current_variant_types = extract_all_variant_types(str(variants["variant_type"]))
		intronic_types = {"splice_region", "intron", "splice_acceptor"}

		if "synonymous" in current_variant_types:
			### ### ###
			# Parse available splicing predictor
			### ### ####
			splicing_predictions = {"MaxEntScan": None, "dbscSNV": None}
			if "MaxEntScan" in variants:
				splicing_predictions["MaxEntScan"] = parse_maxentscan_pred(variants)
			if "dbscSNV" in variants:
				splicing_predictions["dbscSNV"] = parse_dbscsnv_pred(variants)

			### ### ###
			# assign BP7 rule
			### ### ###
			if (splicing_predictions["MaxEntScan"] is None and splicing_predictions["dbscSNV"] is None):
				# no available splicing predictor
				BP7_class = "BP7: False"
				BP7_comment = "BP7: no available splicing predictor to support evidence"
			else:
				### ### ### ### ### ### ### ### ### ###
				# get REVEL and CADD predictions      #
				# calculate the average of all existing scores #
				### ### ### ### ### ### ### ### ### ###
				patho_predictions = {'REVEL': 0, 'CADD': 0}
				for predictor in patho_predictions.keys():
					patho_predictions[predictor], _ = is_variant_pathogenic(variants, predictor)
				is_var_pathogenic = aggregate_patho_predictions(patho_predictions)

				### ### ###
				# base conservation
				### ### ###
				is_var_conserved = is_variant_conserved(variants)

				### ### ###
				# aggregation for splicing impact
				### ### ###
				impacts_splicing = aggregate_splicing_predictions(splicing_predictions)

				self.logger.debug("Variant is predicted as: {}".format(is_var_pathogenic))
				self.logger.debug("Is variant conserved?: {}".format(is_var_conserved))
				self.logger.debug("Does variant affect splicing?: {}".format(impacts_splicing))
				if is_var_conserved == 0 and not impacts_splicing and is_var_pathogenic == "benign":
					BP7_class = "BP7: True"
					BP7_comment = "BP7: " + "aggregation of MaxEntScan and dbscSNV predicts no splicing impact, " + " aggregation of REVEL and CADD predicts benign impact" + " and PhyloP predicts no conserved site"
				else:
					BP7_class = "BP7: False"
					if impacts_splicing:
						splicing_comment = "aggregation of MaxEntScan and dbscSNV predicts splicing impact"
					else:
						splicing_comment = "aggregation of MaxEntScan and dbscSNV predicts no splicing impact"
					if is_var_conserved == 1:
						conserved_comment = "PhyloP predicts conserved site"
					elif is_var_conserved == 0:
						conserved_comment = "PhyloP predicts no conserved site"
					else:
						conserved_comment = "no available conservation prediction to support evidence"
					if is_var_pathogenic == "pathogenic":
						patho_comment = "aggregation of REVEL and CADD predicts pathogenic impact"
					elif is_var_pathogenic == "benign":
						patho_comment = "aggregation of REVEL and CADD predicts benign impact"
					else:
						patho_comment = "no available REVEL nor CADD predictions to support evidence"
					BP7_comment = " ".join(["BP7:", splicing_comment, "and", conserved_comment, "and", patho_comment])
		else:
			# BP7 is for intronic variant types only
			BP7_class = "BP7: NA"
			BP7_comment = "BP7: not applicable for this variant type"

		### ### ###
		# update assignment class and comment for BP7 rule
		### ### ###
		variants["ACMG_rules"] = update_assignment(variants["ACMG_rules"], BP7_class)
		variants["ACMG_rules_comment"] = update_assignment(variants["ACMG_rules_comment"], BP7_comment)
		return variants

	def examine_quality_parameters(self, variants):
		"""
		Examine quality parameters (quality column)

		Parameters
		----------
		variants : pandas.Series
			variant row to be examined for quality parameters

		Returns
		-------
		pandas.Series
			variant row with assigned class based on quality parameters
		"""
		# get qual parameters dictionary, e.g: {'QUAL': 1568.0, 'DP': 79.0, 'AF': 0.94, 'MQM': 27}
		quality_parameters = {}
		for parameter in variants["quality"].split(";"):
			if parameter.split("=")[1].isdigit():
				quality_parameters[parameter.split("=")[0]] = float(parameter.split("=")[1])

		# examine quality, depth and MQM
		if "DP" in quality_parameters and "QUAL" in quality_parameters and "MQM" in quality_parameters:
			if quality_parameters["DP"] < 15:
				self.logger.debug("Assign: re-sequencing, because DP=%f <15x", quality_parameters["DP"])
				self.logger.debug("=> " + ClassifyVariants.var_pos2string(variants))
				quality_class = "QUAL: re-sequencing"
				quality_comment = "QUAL: low depth"
			elif quality_parameters["DP"] < 20 or quality_parameters["QUAL"] < 250 or quality_parameters["MQM"] < 60:
				quality_class = "QUAL: check IGV"
				comment, log = [], []
				if quality_parameters["DP"] < 20:
					comment.append("low depth")
					log.append("DP=" + str(quality_parameters["DP"]) + "<20")
				if quality_parameters["QUAL"] < 250:
					comment.append("low quality")
					log.append("QUAL=" + str(quality_parameters["QUAL"]) + "<250")
				if quality_parameters["MQM"] < 60:
					comment.append("low mapping quality")
					log.append("MQM=" + str(quality_parameters["MQM"]) + "<60")
				self.logger.debug("Assign: check IGV, because %s", ",".join(log))
				self.logger.debug("=> " + ClassifyVariants.var_pos2string(variants))
				quality_comment = "QUAL:" + ",".join(comment)
			else:
				quality_class = "QUAL: True"
				quality_comment = "QUAL: passing as DP>=20,QUAL>=250,MQM>=60"
		else:
			self.logger.debug("Not all quality parameters are found")
			quality_class = "QUAL: NA"
			quality_comment = "QUAL: not all quality parameters are found"
		# update assigned class and comment by quality rule result
		variants["ACMG_rules"] = update_assignment(variants["ACMG_rules"], quality_class)
		variants["ACMG_rules_comment"] = update_assignment(variants["ACMG_rules_comment"], quality_comment)
		return variants

	def extract_inheritance_modes(self, gene_inheritance_info):
		"""
		Extract inheritance modes from gene information row

		Parameters
		----------
		gene_inheritance_info: pandas.Series
			gene inheritance information row

		Returns
		-------
		str, str
			preferred inheritance, alternative inheritance
		"""
		preferred_mode = gene_inheritance_info["Mode of inheritance"]
		alternative_mode = None

		if sum(gene_inheritance_info[["autosomal dominant", "autosomal recessive", "X-linked", "mitochondrial"]]) == 1:
			# if all possible inheritance flags sum to 1 there is only preferred inheritance
			alternative_mode = None
		else:  # find the alternative inheritance mode
			preferred_column = self.digram2inheritance[preferred_mode]
			assert gene_inheritance_info[
				       preferred_column] == 1, "AssertionError: value of preferred inheritance, {}, should be 1".format(
				preferred_column)
			possible_alternatives = ["autosomal dominant", "autosomal recessive", "X-linked", "mitochondrial"]
			possible_alternatives.remove(preferred_column)

			for possible_alternative in possible_alternatives:
				# assuming that only one alternative exists,
				# loop through them and find which one is,
				# save it and finish looping
				if gene_inheritance_info[possible_alternative] == 1:
					alternative_mode = self.inheritance2digram[possible_alternative]

		return preferred_mode, alternative_mode

	def extract_inheritance_modes_acmg(self, gene_name, gene_inheritance_infos, hpo_term):
		"""
		Extract inheritance mode for gene from ACMG HL expert panel, based on phenotype (HPO term)

		Paper with HL-relevant gene inheritance modes:
		DiStefano, Marina T., et al. "ClinGen expert clinical validity curation of 164 hearing loss gene–disease pairs."
		Genetics in Medicine 21.10 (2019): 2239-2247.
		DOI: doi: 10.1038/s41436-019-0487-0.

		Parameters
		----------
		gene_name : str
			gene name
		gene_inheritance_infos: pandas.Series
			gene inheritance information row(s)
		hpo_term : string
			human phenotype term describing the patient's disease

		Returns
		-------
		str, str
			preferred inheritance, alternative inheritance
		"""
		self.logger.debug(
			"Extract preferred and alternative inheritance modes of gene={}, known by ACMG HL expert panel, by HPO term={}".format(
				gene_name, hpo_term))

		### ### ###
		# loop through the inheritance rows of the gene
		# and extract map of disease name to inheritance mode
		### ### ###
		disease2inheritance = {}
		if len(gene_inheritance_infos.shape) == 1:
			# single entry of gene in inheritance mode file
			disease2inheritance[str(gene_inheritance_infos["Disease association"]).strip()] = [
				str(gene_inheritance_infos["Inheritance Pattern"]).strip()]
		else:
			# multiple entries of gene in inheritance mode file
			for _, gene_inheritance_row in gene_inheritance_infos.iterrows():
				if str(gene_inheritance_row["Disease association"]).strip() not in disease2inheritance:
					disease2inheritance[str(gene_inheritance_row["Disease association"]).strip()] = [
						str(gene_inheritance_row["Inheritance Pattern"]).strip()]
				else:
					disease2inheritance[str(gene_inheritance_row["Disease association"]).strip()].append(
						str(gene_inheritance_row["Inheritance Pattern"]).strip())
		self.logger.debug("Disease to inheritance modes: {}".format(disease2inheritance))

		### ### ###
		# select as preferred inheritance the one
		# that is known to be associated with the phenotype term
		### ### ###
		if hpo_term:
			if hpo_term in disease2inheritance:
				if len(disease2inheritance[hpo_term]) == 1:
					# gene contains a single known inheritance mode for input HPO
					self.logger.debug("The gene contains a single known inheritance mode for input HPO term")
					preferred_inheritance, alternative_inheritance = disease2inheritance[hpo_term][0], None
				else:
					# gene contains more than one known inheritance modes for input HPO
					self.logger.debug("The gene contains more than one known inheritance mode for input HPO term")
					preferred_inheritance, alternative_inheritance = disease2inheritance[hpo_term][0], \
					                                                 disease2inheritance[hpo_term][1]
			else:
				# HPO term in not matching known disease associations
				self.logger.debug(
					"HPO term={} is not matching known disease associations of the gene={}".format(hpo_term))
				preferred_inheritance, alternative_inheritance = None, None
		else:
			### ### ###
			# HPO term is not supplied
			# if single known disease => assign preferred inheritance
			# else if all inheritance of each disease are equal to each other => use them as preferred and alternative
			# else => assign unknown inheritance modes
			### ### ###
			self.logger.debug("HPO term is not supplied")
			if len(disease2inheritance) == 1:
				single_disease = next(iter(disease2inheritance), None)
				known_inheritances = disease2inheritance[single_disease]
				if len(known_inheritances) == 1:
					# single disease association, with single inheritance
					self.logger.debug("Single disease association with single inheritance")
					preferred_inheritance, alternative_inheritance = known_inheritances[0], None
				else:
					# single disease association, with multiple inheritances
					self.logger.debug("Single disease association with more than one known inheritances")
					preferred_inheritance, alternative_inheritance = known_inheritances[0], known_inheritances[1]
			else:
				### ### ###
				# check if all inheritances
				# of all associated diseases are equal
				### ### ###
				all_inheritances_equal = True
				first_disease_inheritances = next(iter(disease2inheritance.values()))
				for disease_inheritances in disease2inheritance.values():
					if disease_inheritances != first_disease_inheritances:
						all_inheritances_equal = False
						break
				if all_inheritances_equal:
					self.logger.debug("All associated diseases contain equal inheritance modes")
					if len(first_disease_inheritances) == 2:
						preferred_inheritance, alternative_inheritance = first_disease_inheritances[0], \
						                                                 first_disease_inheritances[1]
					else:
						preferred_inheritance, alternative_inheritance = first_disease_inheritances[0], None
				else:
					# all diseases do not have the same inheritance modes
					preferred_inheritance, alternative_inheritance = None, None
		self.logger.debug("Final preferred: {}, alternative: {}".format(preferred_inheritance, alternative_inheritance))
		return preferred_inheritance, alternative_inheritance

	def is_threshold_rule_triggered(self, rule_name, inheritance_mode, pop_MAF):
		"""
		Examine if threshold rule is triggered

		Parameters
		----------
		rule_name : str
			name of examined rule
		inheritance_mode : str
			digram of inheritance mode to be used for rule threshold
		pop_MAF : float
			sub-population MAF to be used for the rule examination

		Returns
		-------
		bool
			rule is triggered (True), not triggered (False)
		"""
		self.logger.debug("Examine if rule: {} is triggered, with pop MAF={}".format(rule_name, pop_MAF))
		threshold_direction = str(self.inheritance_thresholds_df.loc[rule_name, "threshold_direction"])
		self.logger.debug("threshold direction: {}".format(threshold_direction))

		if inheritance_mode == "unknown":
			# if inheritance mode is unknown, use the strictest inheritance as testing inheritance
			inheritance_mode = str(self.inheritance_thresholds_df.loc[rule_name, "strictest_inheritance"])
		if threshold_direction == "greater":
			if pop_MAF > float(self.inheritance_thresholds_df.loc[rule_name, inheritance_mode]):
				return True
			else:
				return False
		elif threshold_direction == "greater_than_or_equal":
			self.logger.debug(
				"threshold: {}".format(float(self.inheritance_thresholds_df.loc[rule_name, inheritance_mode])))
			if pop_MAF >= float(self.inheritance_thresholds_df.loc[rule_name, inheritance_mode]):
				return True
			else:
				return False
		elif threshold_direction == "lower":
			if pop_MAF < float(self.inheritance_thresholds_df.loc[rule_name, inheritance_mode]):
				return True
			else:
				return False
		elif threshold_direction == "lower_than_or_equal":
			if pop_MAF <= float(self.inheritance_thresholds_df.loc[rule_name, inheritance_mode]):
				return True
			else:
				return False

	def examine_autosomal_recessive_MAF(self, gnomAD_subpopulations, exclude_var):
		"""
		Examine MAF for autosomal recessive gene for hearing loss
		based on Oza et al. 2018, PMID: 30311386

		Parameters
		----------
		gnomAD_subpopulations: dict of str: float
			gnomAD subpopulation MAF
		exclude_var : bool
			to exclude variant from BA1 and BS1 rules (True), otherwise False

		Returns
		-------
		str
			rules examination result
		str
			comments for each rule examination result
		dict of str : list of int
			examined rules
		"""
		### ### ###
		# a) examine each rule per sub-population
		# b) aggregate rule assignment per sub-population
		# c) deactivate lower strength version of rules if the stronger versions are triggered
		### ### ###
		self.logger.debug("Examine for autosomal recessive")
		examined_rules = {"PM2": [0], "PM2_Supporting": [0], "BA1": [0], "BS1": [0], "BS1_Supporting": [0]}
		examined_rules_comment = {"PM2": [], "PM2_Supporting": [], "BA1": [], "BS1": [], "BS1_Supporting": []}

		# a) examine how many times a rule is triggered by a MAF of a subpopulation
		for rule_name in examined_rules.keys():
			for subpop, subpop_MAF in gnomAD_subpopulations.items():
				if self.is_threshold_rule_triggered(rule_name, "AR", subpop_MAF):
					examined_rules[rule_name].append(1)
					examined_rules_comment[rule_name].append("{}: {} {} {}".format(subpop, subpop_MAF,
					                                                               self.inheritance_thresholds_df.loc[
						                                                               rule_name, "threshold_direction"].replace(
						                                                               "_", " "),
					                                                               self.inheritance_thresholds_df.loc[
						                                                               rule_name, "AR"]))

		# b) aggregate MAF examination per sub-population
		ACMG_rules, ACMG_rules_comment = [], []
		for rule_name, variant_class_presence in examined_rules.items():
			if (rule_name == "BA1" or rule_name == "BS1" or rule_name == "BS1_Supporting") and exclude_var:
				ACMG_rules.append(rule_name + ": False")
				ACMG_rules_comment.append(
					rule_name + ": variant in list of phenotype relevant and pathogenic variants with high AF")
				# set all sub-populations not to trigger rule
				examined_rules[rule_name] = [0]
			elif sum(variant_class_presence) > 0:
				ACMG_rules.append(rule_name + ": True")
				ACMG_rules_comment.append(rule_name + ": " + ",".join(examined_rules_comment[rule_name]))
			else:
				ACMG_rules.append(rule_name + ": False")
				ACMG_rules_comment.append(rule_name + ": MAF not triggering rule")

		# c) deactivate lower strength versions of PM2 and BS1 (if needed)
		ACMG_rules_valid, ACMG_rules_valid_comment, examined_rules_valid = deactivate_lower_strength_rules(
			examined_rules, ACMG_rules,
			ACMG_rules_comment)

		if len(ACMG_rules_valid) > 0:
			return "|".join(ACMG_rules_valid), "|".join(ACMG_rules_valid_comment), examined_rules_valid

	def examine_autosomal_dominant_MAF(self, gnomAD_subpopulations, exclude_var):
		"""
		Examine MAF for autosomal dominant gene for hearing loss
		based on Oza et al. 2018, PMID: 30311386

		Parameters
		----------
		gnomAD_subpopulations: dict of str: float
			gnomAD subpopulation MAF
		exclude_var : bool
			to exclude variant from BA1 and BS1 rules (True), otherwise False

		Returns
		-------
		str
			rules examination result
		str
			comments for each rule examination result
		dict of str : list of int
			examined rules
		"""
		### ### ###
		# a) examine each rule per sub-population
		# b) aggregate rule assignment over all sub-populations
		# c) deactivate PM2_Supporting if PM2 is triggered equally deactivate BS1_Supporting upon BS1 trigger
		### ### ###
		self.logger.debug("Examine thresholds for autosomal dominant")
		examined_rules = {"PM2": [0], "BA1": [0], "BS1": [0]}
		examined_rules_comment = {"PM2": [], "BA1": [], "BS1": []}

		# a) examine how many times a rule is triggered by a MAF of a sub-population
		for rule_name in examined_rules.keys():
			for subpop, subpop_MAF in gnomAD_subpopulations.items():
				if self.is_threshold_rule_triggered(rule_name, "AD", subpop_MAF):
					examined_rules[rule_name].append(1)
					examined_rules_comment[rule_name].append("{}: {} {} {}".format(subpop, subpop_MAF,
					                                                               self.inheritance_thresholds_df.loc[
						                                                               rule_name, "threshold_direction"].replace(
						                                                               "_", " "),
					                                                               self.inheritance_thresholds_df.loc[
						                                                               rule_name, "AD"]))
		# b) aggregate MAF examination per sub-population
		ACMG_rules, ACMG_rules_comment = [], []
		for rule_name, variant_class_presence in examined_rules.items():
			if (rule_name == "BA1" or rule_name == "BS1") and exclude_var:
				ACMG_rules.append(rule_name + ": False")
				ACMG_rules_comment.append(
					rule_name + ": variant in list of phenotype relevant and pathogenic variants with high AF")
				# set all sub-populations not to trigger rule
				examined_rules[rule_name] = [0]
			elif sum(variant_class_presence) > 0:
				ACMG_rules.append(rule_name + ": True")
				ACMG_rules_comment.append(
					rule_name + ": " + ",".join(examined_rules_comment[rule_name]))
			else:
				ACMG_rules.append(rule_name + ": False")
				ACMG_rules_comment.append(rule_name + ": MAF not triggering rule")
		self.logger.debug("Rules after AD: {}".format(ACMG_rules))
		# c) deactivate PM2 supporting and/or BS1_Supporting if PM2 and/or BS1 is triggered
		ACMG_rules_valid, ACMG_rules_valid_comment, examined_rules_valid = deactivate_lower_strength_rules(
			examined_rules, ACMG_rules, ACMG_rules_comment)
		self.logger.debug("Rules after correction:{}".format(ACMG_rules_valid))
		if len(ACMG_rules_valid) > 0:
			return "|".join(ACMG_rules_valid), "|".join(ACMG_rules_valid_comment), examined_rules_valid

	def examine_unknown_inheritance_MAF(self, gnomAD_subpopulations, exclude_var):
		"""
		Examine MAF for unknown inheritance gene for hearing loss
		based on Oza et al. 2018, PMID: 30311386

		Parameters
		----------
		gnomAD_subpopulations : dict of str: float
			gnomAD subpopulation MAF
		exclude_var : bool
			exclude variant for BA1 and BS1 rule (True), otherwise False

		Returns
		-------
		str
			rules examination result
		str
			comments for each rule examination result
		dict of str : list of int
			examined rules
		"""
		### ### ###
		# a) find all candidate rules (e.g. a rule that does contain a threshold for at least one inheritance mode)
		# b) examine and find assignment for each candidate rule using the strictest inheritance mode
		# c) aggregate per sub-populations
		# d) deactivate PM2 and BS1 supporting versions if the stronger versions are triggered
		self.logger.debug("Examine rules for a gene with unknown inheritance")

		# a) find all candidate rules
		candidate_rules = []
		for rule_name, rule_row in self.inheritance_thresholds_df.iterrows():
			if not (isnan(float(rule_row["AR"])) and isnan(float(rule_row["AD"]))):
				candidate_rules.append(rule_name)
			else:
				self.logger.debug(
					"Rule: {} contains unknown MAF thresholds for AR and AD inheritance modes".format(rule_name))
		# b) use strictest inheritance mode for each candidate rule
		# initialize examined rules and comment
		examined_rules, examined_rules_comment = {}, {}
		for rule_name in candidate_rules:
			examined_rules[rule_name] = [0]
			examined_rules_comment[rule_name] = []
		for rule_name in candidate_rules:
			for subpop, subpop_MAF in gnomAD_subpopulations.items():
				if self.is_threshold_rule_triggered(rule_name, "unknown", subpop_MAF):
					examined_rules[rule_name].append(1)
					examined_rules_comment[rule_name].append("{}: {} {} {}".format(subpop, subpop_MAF,
					                                                               self.inheritance_thresholds_df.loc[
						                                                               rule_name, "threshold_direction"].replace(
						                                                               "_", " "),
					                                                               self.inheritance_thresholds_df.loc[
						                                                               rule_name, "strictest_inheritance"]))

		# c) aggregate sub-population assignments
		ACMG_rules, ACMG_rules_comment = [], []
		for rule_name, variant_class_presence in examined_rules.items():
			if (rule_name == "BA1" or rule_name == "BS1" or rule_name == "BS1_Supporting") and exclude_var:
				ACMG_rules.append(rule_name + ": False")
				ACMG_rules_comment.append(
					rule_name + ": variant in list of phenotype relevant and pathogenic variants with high AF")
				# set all sub-populations not to trigger rule
				examined_rules[rule_name] = [0]
			elif sum(variant_class_presence) > 0:
				ACMG_rules.append(rule_name + ": True")
				ACMG_rules_comment.append(
					rule_name + ": " + ",".join(examined_rules_comment[rule_name]))
			else:
				ACMG_rules.append(rule_name + ": False")
				ACMG_rules_comment.append(rule_name + ": MAF not triggering rule")

		# d) deactivate PM2 supporting and/or BS1_Supporting if PM2 and/or BS1 is triggered
		ACMG_rules_valid, ACMG_rules_valid_comment, examined_rules_valid = deactivate_lower_strength_rules(
			examined_rules,
			ACMG_rules, ACMG_rules_comment)
		if len(ACMG_rules_valid) > 0:
			return "|".join(ACMG_rules_valid), "|".join(ACMG_rules_valid_comment), examined_rules_valid

	def extract_gene_inheritance_mode(self, gene_name, use_phenotype_inheritance):
		"""
		Extract gene inheritance mode from known gene-inheritance mode relations

		Parameters
		----------
		gene_name : str
			gene name
		use_phenotype_inheritance : bool
			use ACMG inheritance modes of genes based on HPO term (True), otherwise use MHH HG department known inheritance modes for genes

		Returns
		-------
		str
			preferred inheritance mode
		str
			alternative inheritance mode
		"""
		# todo: add HPO as input argument
		hpo_term = None
		# examine if gene has an known inheritance mode
		if gene_name in self.genes_inheritance_df.index:
			if use_phenotype_inheritance:
				return self.extract_inheritance_modes_acmg(gene_name, self.genes_inheritance_df.loc[gene_name],
				                                           hpo_term)
			else:
				return self.extract_inheritance_modes(
					self.genes_inheritance_df.loc[gene_name])
		else:  # if gene has unknown inheritance mode, return both values as None
			return None, None

	def aggregate_preferred_alternative_rules_assignments(self, examined_rules_preferred, inheritance_preferred,
	                                                      examined_rules_alternative, inheritance_alternative):
		"""
		Aggregate rules assignments for preferred and alternative inheritance modes

		Parameters
		----------
		examined_rules_preferred : dict of str : int
			assignments of examined rules for preferred inheritance mode
		inheritance_preferred : str
			preferred inheritance mode
		examined_rules_alternative : dict of str : int
			assignments of examined rules for alternative inheritance mode
		inheritance_alternative : str
			alternative inheritance mode

		Returns
		-------
		str
			aggregated assignment
		str
			aggregated flags
		"""
		### ### ###
		# a) aggregate rules that exist for both inheritance modes
		# b) aggregate rules that exist for preferred mode only
		# c) aggregate rules that exist for alternative mode only
		### ### ###
		self.logger.debug("Aggregate assignments for preferred and alternative inheritance mode")
		self.logger.debug("Preferred: {}, mode: {}".format(examined_rules_preferred, inheritance_preferred))
		self.logger.debug("Alternative: {}, mode: {}".format(examined_rules_alternative, inheritance_alternative))
		merged_assignments = {}
		merging_flags = []

		# a) aggregate assignments that exist for both preferred and alternative inheritance mode
		for rule_name in examined_rules_preferred:
			if rule_name in examined_rules_alternative:

				if sum(examined_rules_preferred[rule_name]) > 0 and sum(examined_rules_alternative[rule_name]) > 0:
					# variant passes threshold for both inheritance modes
					merged_assignments[rule_name] = "True"
				elif sum(examined_rules_preferred[rule_name]) > 0 or sum(examined_rules_alternative[rule_name]) > 0:
					# variant passes threshold for only one inheritance mode
					self.logger.debug(
						"For rule: {}, strictest_inheritance= {}".format(rule_name, self.inheritance_thresholds_df.loc[
							rule_name, "strictest_inheritance"]))
					if self.inheritance_thresholds_df.loc[
						rule_name, "strictest_inheritance"] == inheritance_preferred:
						if sum(examined_rules_preferred[rule_name]) > 0:
							merged_assignments[rule_name] = "True"
						else:
							merged_assignments[rule_name] = "False"
							merging_flags.append(rule_name + ": triggered not for the strictest inheritance mode")
					elif self.inheritance_thresholds_df.loc[
						rule_name, "strictest_inheritance"] == inheritance_alternative:
						if sum(examined_rules_alternative[rule_name]) > 0:
							merged_assignments[rule_name] = "True"
						else:
							merged_assignments[rule_name] = "False"
							merging_flags.append(rule_name + ": triggered not for the strictest inheritance mode")
				else:
					# variant does not pass thresholds for both inheritance modes
					merged_assignments[rule_name] = "False"

		# b) add the preferred-specific rules to the aggregated results
		for rule_name in examined_rules_preferred:
			if rule_name not in merged_assignments:
				if sum(examined_rules_preferred[rule_name]) > 0:
					merged_assignments[rule_name] = "True"
				else:
					merged_assignments[rule_name] = "False"

		# c) add the alternative-specific rules to the aggregated results
		for rule_name in examined_rules_alternative:
			if rule_name not in merged_assignments:
				if sum(examined_rules_alternative[rule_name]) > 0:
					merged_assignments[rule_name] = "True"
				else:
					merged_assignments[rule_name] = "False"
		self.logger.debug("merged_assignments: {}".format(merged_assignments))
		self.logger.debug("merging_flags: {}".format(merging_flags))
		aggregated_assignment = "||".join(
			[rule_name + ": " + rule_assignment for rule_name, rule_assignment in merged_assignments.items()])
		assignment_flags = "||".join(merging_flags)
		self.logger.debug("Aggregated assignment: {}".format(aggregated_assignment))
		self.logger.debug("Aggregated merging flags: {}".format(assignment_flags))
		return aggregated_assignment, assignment_flags

	def examine_inheritance_specific_MAF(self, variants):
		"""
		Examine MAF based on gene inheritance mode

		Parameters
		----------
		variants : pandas.Series
			variants row to be examined for inheritance specific MAF

		Returns
		-------
		pandas.Series
			variant row with assigned class based on gene-specific inheritance mode
		"""
		self.logger.debug(
			"Examine inheritance-specific rule for: ".format(ClassifyVariants.var_pos2string(variants)))
		gene = str(variants["gene"]).split(",")[0]
		current_variant = VariantInfo(self.current_sample_name, gene, str(variants["variant_type"]),
		                              None,
		                              str(variants["chr"]),
		                              str(variants["start"]), str(variants["end"]), str(variants["ref"]),
		                              str(variants["obs"]))
		### ### ###
		# parse gnomAD subpopulations MAF
		### ### ###
		comment, log = ["Inheritance_specific_rules:"], []
		ACMG_flags = None

		### ### ###
		# examine if gene has an known inheritance mode
		### ### ###
		use_phenotype_inheritances = False
		preferred_inheritance, alternative_inheritance = self.extract_gene_inheritance_mode(gene,
		                                                                                    use_phenotype_inheritances)
		if not (preferred_inheritance or alternative_inheritance):
			comment.append("unknown inheritance mode")
			preferred_inheritance = "unknown"
		self.logger.debug("Resulted inheritances: preferred: {}, alternative: {}".format(preferred_inheritance,
		                                                                                 alternative_inheritance))
		### ### ###
		# extract gnomAD sub-populations
		### ### ###
		gnomAD_subpopulations = {}
		if "gnomAD_sub" in variants:
			# GSvar file contains gnomAD subpopulation information
			if "," not in str(variants["gnomAD_sub"]):
				# if gnomAD subpopulation is not known, fall back to gnomAD
				if isnan(float(variants["gnomAD"])):
					# not present in gnomAD, not present in gnomAD subpopulations
					self.logger.debug("Not present in BOTH gnomAD subpopulations and general gnomAD, assume MAF=0")
					comment.append("Not present in BOTH gnomAD subpopulations and general gnomAD, assume MAF=0")
					gnomAD_subpopulations["general"] = 0
				else:
					# not present in gnomAD subpopulations, but present in gnomAD
					gnomAD_subpopulations["general"] = float(variants["gnomAD"])
					self.logger.debug("Not present in gnomAD, fall back to general gnomAD: {}".format(
						gnomAD_subpopulations["general"]))
					comment.append("Not present in gnomAD, fall back to general gnomAD: {}".format(
						gnomAD_subpopulations["general"]))
			else:
				# parse gnomAD subpopulations
				max_MAF_subpop, max_MAF_subpop_index, subpop_MAFs = 0.0, 0, []
				for subpop_index, subpop_MAF in enumerate(str(variants["gnomAD_sub"]).split(",")):
					if subpop_MAF == '':
						# MAF is not found for this sub-population make it equal to 0
						subpop_MAF = 0.0
					# update max MAF
					if float(subpop_MAF) > max_MAF_subpop:
						max_MAF_subpop = float(subpop_MAF)
						max_MAF_subpop_index = subpop_index
					subpop_MAFs.append(float(subpop_MAF))

				gnomAD_subpopulations[self.position2gnomAD_sub[max_MAF_subpop_index]] = max_MAF_subpop
		else:
			# GSvar does not contain gnomAD subpopulation, so use general gnomAD AF
			if isnan(float(variants["gnomAD"])):
				# not present in gnomAD, not present in gnomAD subpopulations
				gnomAD_subpopulations["general"] = 0
				self.logger.debug(
					"gnomAD subpopulation column is not given and variant not present in general gnomAD, assume MAF=0")
				comment.append(
					"gnomAD subpopulation column is not given and variant not present in general gnomAD, assume MAF=0")
			else:
				# not present in gnomAD subpopulations, but present in gnomAD
				gnomAD_subpopulations["general"] = float(variants["gnomAD"])
				self.logger.debug("gnomAD subpopulation column is not given, fall back to general gnomAD: {}".format(
					gnomAD_subpopulations["general"]))
				comment.append("gnomAD subpopulation column is not given, fall back to general gnomAD: {}".format(
					gnomAD_subpopulations["general"]))
		self.logger.debug("Parsed gnomAD subpopulations: {}".format(gnomAD_subpopulations))

		### ### ###
		# a) check if variant is contained in the known list of pathogenic variants with high AF
		# b) for each variant check the preferred inheritance at least
		# c) if variant-affecting gene contains two inheritance modes, concatenate assignments
		### ### ###

		# a) check if variants is contained in the known list of pathogenic variants with high AF
		exclude_var_benign_criteria = is_var2exclude(self.pheno_high_af_variants, current_variant)
		self.logger.debug("Exclude variant: {}".format(exclude_var_benign_criteria))

		# b) examine preferred inheritance mode
		self.logger.debug("Examine preferred inheritance mode: {}".format(preferred_inheritance))
		if preferred_inheritance == "AD":  # autosomal dominant
			ACMG_assignments_preferred, ACMG_assignment_comment_preferred, examined_rules_preferred = self.examine_autosomal_dominant_MAF(
				gnomAD_subpopulations, exclude_var_benign_criteria)
		elif preferred_inheritance == "AR" or preferred_inheritance == "XL":  # autosomal recessive or X-linked
			ACMG_assignments_preferred, ACMG_assignment_comment_preferred, examined_rules_preferred = self.examine_autosomal_recessive_MAF(
				gnomAD_subpopulations, exclude_var_benign_criteria)
		elif preferred_inheritance == "unknown" or preferred_inheritance == "MT":
			ACMG_assignments_preferred, ACMG_assignment_comment_preferred, examined_rules_preferred = self.examine_unknown_inheritance_MAF(
				gnomAD_subpopulations, exclude_var_benign_criteria)

		# c) concatenate MAF filters for preferred and alternative
		if alternative_inheritance:
			self.logger.debug("Examine alternative inheritance mode:")
			if alternative_inheritance == "AD":  # autosomal dominant
				ACMG_assignments_alternative, ACMG_assignment_comment_alternative, examined_rules_alternative = self.examine_autosomal_dominant_MAF(
					gnomAD_subpopulations, exclude_var_benign_criteria)
			elif alternative_inheritance == "AR" or alternative_inheritance == "XL":  # autosomal recessive or X-linked
				ACMG_assignments_alternative, ACMG_assignment_comment_alternative, examined_rules_alternative = self.examine_autosomal_recessive_MAF(
					gnomAD_subpopulations, exclude_var_benign_criteria)
			elif alternative_inheritance == "MT":
				ACMG_assignments_alternative, ACMG_assignment_comment_alternative, examined_rules_alternative = self.examine_unknown_inheritance_MAF(
					gnomAD_subpopulations, exclude_var_benign_criteria)
			# concatenate assignments for two inheritance modes
			ACMG_assignments = preferred_inheritance + "=" + ACMG_assignments_preferred + "||" + alternative_inheritance + "=" + ACMG_assignments_alternative
			ACMG_assignment_comment = " ".join(
				comment) + " preferred: " + preferred_inheritance + "," + ACMG_assignment_comment_preferred + "||" + "alternative: " + alternative_inheritance + "," + ACMG_assignment_comment_alternative
		else:
			# only one inheritance mode is known for current gene
			ACMG_assignments = preferred_inheritance + "=" + ACMG_assignments_preferred
			ACMG_assignment_comment = " ".join(
				comment) + " single_inheritance: " + preferred_inheritance + "," + ACMG_assignment_comment_preferred
		# update assigned class and assignment class columns
		variants["ACMG_rules"] = update_assignment(variants["ACMG_rules"], ACMG_assignments)
		variants["ACMG_rules_comment"] = update_assignment(variants["ACMG_rules_comment"],
		                                                   ACMG_assignment_comment)
		if ACMG_flags:
			if variants["ACMG_rule_flags"] == "-":
				variants["ACMG_rule_flags"] = "Warning: " + ACMG_flags
			else:
				self.logger.info(
					"Inheritance-specific rule adds to warning flag for a gene not related to HL: ".format(
						ClassifyVariants.var_pos2string(variants)))
				variants["ACMG_rule_flags"] = update_assignment(variants["ACMG_rule_flags"], ACMG_flags)
		self.logger.debug("Inheritance specific rules = {}".format(ACMG_assignments))
		return variants

	def refine_pvs1_rule(self, variants):
		"""
		Refine the possible PVS1 rule for each variant

		Parameters
		----------
		variants : pandas.dataframe
			chunk of variants dataframe to run PVS1

		Returns
		-------
		pandas.dataframe
			updated variants chunk with assigned class based on PVS1 rule
		"""
		self.logger.info("Refine PVS1 rule with thread id: {}\n".format(get_ident()))
		# initialize a new RefineLoF object for each thread
		RefineLoF = RefineLossOfFunction(self.data_path, self.clinvar_path,
		                                 self.clinvar_file, self.clinvar_stars_df, self.min_review_stars,
		                                 self.beds_root, self.critical_prot_regions_file,
		                                 self.clinical_exons_file, self.pheno_high_af_variants,
		                                 self.phenotype_exons_file, self.uniprot_domains_file, self.hugo_genes_df)

		# refine PVS1 for each variant row
		# then save assigned class and comment to corresponding columns of the row
		variants_num, processed_var_counter = variants.shape[0], 0
		logged_processed_percentage = 0
		for var_idx, var_row in variants.iterrows():
			if is_var_type4PVS1(str(var_row["variant_type"])):
				# self.logger.debug("Refine PVS1 for variant: {}".format(AnnotateVariants.var_pos2string(var_row)))
				ACMG_rules, ACMG_rules_comment, var_coding_seq = RefineLoF.run(self.current_sample_name,
				                                                               str(var_row["chr"]),
				                                                               int(var_row["start"]),
				                                                               int(var_row["end"]),
				                                                               str(var_row["gene"]),
				                                                               str(var_row["ref"]), str(var_row["obs"]),
				                                                               str(var_row["variant_type"]),
				                                                               str(var_row["coding_and_splicing"]))
			else:
				ACMG_rules, ACMG_rules_comment = "PVS1: NA", "PVS1: not applicable for this variant type"
				var_coding_seq = None
			variants.at[var_idx, "ACMG_rules"] = update_assignment(var_row["ACMG_rules"], ACMG_rules)
			variants.at[var_idx, "ACMG_rules_comment"] = update_assignment(var_row["ACMG_rules_comment"],
			                                                               ACMG_rules_comment)
			variants.at[var_idx, "var_coding_seq"] = format_dict2str(var_coding_seq)
			processed_var_counter = processed_var_counter + 1
			current_processed_percentage = (processed_var_counter / variants_num) * 100
			if current_processed_percentage - logged_processed_percentage >= 25.0:
				logged_processed_percentage = current_processed_percentage
				self.logger.info(
					"Thread id: {}, {}, Refine PVS1 progress: {:.2f}\n ".format(
						get_ident(), ClassifyVariants.var_pos2string(var_row), logged_processed_percentage))

		return variants

	def examine_pm4(self, variants):
		"""
		Run PM4 assessment

		Parameters
		----------
		variants : pandas.DataFrame
			chunk of varaints dataframe to run PM4 for

		Returns
		-------
		pandas.DataFrame
			updated variants chunk with assigned PM4 rule
		"""
		self.logger.info("Examine PM4 rule with thread id: {}\n".format(get_ident()))
		# initialize a new PM4 object for current thread
		assign_PM4 = AssignPM4(self.data_path, self.beds_root, self.uniprot_repeat_file)

		# examine PM4 rule for each variant row
		# then save assigned class and commetn to correspoding columns of row
		variants_num, processed_var_counter = variants.shape[0], 0
		logged_processed_percentage = 0
		for var_idx, var_row in variants.iterrows():
			if extract_rule_result(var_row["ACMG_rules"], "PVS1") == "False" and (
					"inframe" in str(var_row["variant_type"]) or "stop_lost" in str(var_row["variant_type"])):
				ACMG_rules, ACMG_rules_comment = assign_PM4.run(self.current_sample_name, str(var_row["chr"]),
				                                                int(var_row["start"]),
				                                                int(var_row["end"]),
				                                                str(var_row["gene"]), str(var_row["ref"]),
				                                                str(var_row["obs"]),
				                                                str(var_row["variant_type"]),
				                                                str(var_row["coding_and_splicing"]),
				                                                deformat_str2dict(str(var_row["var_coding_seq"])))
			else:
				ACMG_rules, ACMG_rules_comment = "PM4: NA", "PM4: not applicable for this variant type"
			variants.at[var_idx, "ACMG_rules"] = update_assignment(var_row["ACMG_rules"], ACMG_rules)
			variants.at[var_idx, "ACMG_rules_comment"] = update_assignment(var_row["ACMG_rules_comment"],
			                                                               ACMG_rules_comment)
			# update process progress
			processed_var_counter = processed_var_counter + 1
			current_processed_percentage = (processed_var_counter / variants_num) * 100
			if current_processed_percentage - logged_processed_percentage >= 25.0:
				logged_processed_percentage = current_processed_percentage
				self.logger.info(
					"Thread id: {}, {}, Examine PM4 progress: {:.2f}\n ".format(
						get_ident(), ClassifyVariants.var_pos2string(var_row), logged_processed_percentage))
		return variants

	def examine_bp3(self, variants):
		"""
		Run BP3 assessment

		Parameters
		----------
		variants : pandas.DataFrame
			chunk of variants dataframe to run BP3 for

		Returns
		-------
		pandas.DataFrame
			updated variants chunk with assigned BP3 rule
		"""
		self.logger.info("Examine BP3 rule with thread id: {}\n".format(get_ident()))
		# initialize a new BP3 object for current thread
		assign_BP3 = AssignBP3(self.data_path, self.beds_root, self.repeats_no_domains_file)

		# examine BP3 rule for each variant row
		# then save assigned class and comment to corresponding columns of the row
		variants_num, processed_var_counter = variants.shape[0], 0
		logged_processed_percentage = 0
		for var_idx, var_row in variants.iterrows():
			if "inframe" in str(var_row["variant_type"]) or (
					extract_rule_result(var_row["ACMG_rules"], "PVS1") == "FAIL" and "stop_gained" in str(
				var_row["variant_type"])):
				ACMG_rules, ACMG_rules_comment = assign_BP3.run(self.current_sample_name, str(var_row["chr"]),
				                                                int(var_row["start"]),
				                                                int(var_row["end"]),
				                                                str(var_row["gene"]), str(var_row["ref"]),
				                                                str(var_row["obs"]),
				                                                str(var_row["variant_type"]),
				                                                str(var_row["coding_and_splicing"]))
			else:
				ACMG_rules, ACMG_rules_comment = "BP3: NA", "BP3: not applicable for this variant type"
			variants.at[var_idx, "ACMG_rules"] = update_assignment(var_row["ACMG_rules"], ACMG_rules)
			variants.at[var_idx, "ACMG_rules_comment"] = update_assignment(var_row["ACMG_rules_comment"],
			                                                               ACMG_rules_comment)
			# update process progress
			processed_var_counter = processed_var_counter + 1
			current_processed_percentage = (processed_var_counter / variants_num) * 100
			if current_processed_percentage - logged_processed_percentage >= 25.0:
				logged_processed_percentage = current_processed_percentage
				self.logger.info(
					"Thread id: {}, {}, Examine BP3 progress: {:.2f}\n ".format(
						get_ident(), ClassifyVariants.var_pos2string(var_row), logged_processed_percentage))
		return variants

	def examine_pm1(self, variants):
		self.logger.info("Examine PM1 rule with process id: {}\n".format(get_ident()))
		# initialize PM1 rule for current thread
		assign_PM1 = AssignPM1(self.beds_root, self.critical_prot_regions_no_benign_file, self.pm1_regions_file)

		# examine PM1 rule per variant row
		# then save assigned class and comment to corresponding columns of the row
		variants_num, processed_var_counter = variants.shape[0], 0
		logged_processed_percentage = 0
		for var_idx, var_row in variants.iterrows():
			if str(var_row["variant_type"]) == "missense":
				ACMG_rules, ACMG_rules_comment = assign_PM1.run(self.current_sample_name, str(var_row["chr"]),
				                                                int(var_row["start"]), int(var_row["end"]),
				                                                str(var_row["gene"]), str(var_row["ref"]),
				                                                str(var_row["obs"]), str(var_row["variant_type"]),
				                                                str(var_row["coding_and_splicing"]))

			else:
				ACMG_rules, ACMG_rules_comment = "PM1: NA", "PM1: not applicable for this variant type"
			variants.at[var_idx, "ACMG_rules"] = update_assignment(var_row["ACMG_rules"], ACMG_rules)
			variants.at[var_idx, "ACMG_rules_comment"] = update_assignment(var_row["ACMG_rules_comment"],
			                                                               ACMG_rules_comment)
			# update process progress
			processed_var_counter = processed_var_counter + 1
			current_processed_percentage = (processed_var_counter / variants_num) * 100
			if current_processed_percentage - logged_processed_percentage >= 25.0:
				logged_processed_percentage = current_processed_percentage
				self.logger.info(
					"Thread id: {}, {}, Examine PM1 progress: {:.2f}\n ".format(
						get_ident(), ClassifyVariants.var_pos2string(var_row), logged_processed_percentage))
		return variants

	def process_clinvar_rules(self, variants):
		"""
		Examine ClinVar related rules

		Parameters
		----------
		variants : pandas.DataFrame
			chunk of variants dataframe to examine ClinVar rules for

		Returns
		-------
		pandas.dataframe
			updated variants chunk with assigned class based on ClinVar evidence
		"""
		self.logger.info("Process ClinVar rules with process id: {}\n".format(get_ident()))
		processClinVar = ProcessClinVar(self.data_path, self.clinvar_path, self.clinvar_file, self.clinvar_stars_df,
		                                self.min_review_stars, self.hugo_genes_df)

		# examine ClinVar rules per variant row
		# then save assigned class and comment to corresponding columns of the row
		variants_num, processed_var_counter = variants.shape[0], 0
		logged_processed_percentage = 0
		for var_idx, var_row in variants.iterrows():
			if "missense" in str(var_row["variant_type"]):
				# self.logger.debug("Run PS1 and PM5 for variant: {}".format(AnnotateVariants.var_pos2string(var_row)))
				ACMG_rules, ACMG_rules_comment = processClinVar.run(self.current_sample_name, str(var_row["chr"]),
				                                                    int(var_row["start"]),
				                                                    int(var_row["end"]), str(var_row["gene"]),
				                                                    str(var_row["ref"]), str(var_row["obs"]),
				                                                    str(var_row["variant_type"]),
				                                                    str(var_row["coding_and_splicing"]))
				variants.at[var_idx, "ACMG_rules"] = update_assignment(var_row["ACMG_rules"], ACMG_rules)
				variants.at[var_idx, "ACMG_rules_comment"] = update_assignment(var_row["ACMG_rules_comment"],
				                                                               ACMG_rules_comment)
			else:
				aggregated_PS1_PM5_class = update_assignment("PS1: NA", "PM5: NA")
				aggregated_PS1_PM5_comment = update_assignment("PS1: not applicable for this variant type",
				                                               "PM5: not applicable for this variant type")
				variants.at[var_idx, "ACMG_rules"] = update_assignment(var_row["ACMG_rules"],
				                                                       aggregated_PS1_PM5_class)
				variants.at[var_idx, "ACMG_rules_comment"] = update_assignment(var_row["ACMG_rules_comment"],
				                                                               aggregated_PS1_PM5_comment)
			processed_var_counter = processed_var_counter + 1
			current_processed_percentage = (processed_var_counter / variants_num) * 100
			if current_processed_percentage - logged_processed_percentage >= 25.0:
				logged_processed_percentage = current_processed_percentage
				self.logger.info(
					"Thread id: {}, {}, Examined ClinVar rules progress: {:.2f}\n ".format(
						get_ident(), ClassifyVariants.var_pos2string(var_row), logged_processed_percentage))
		return variants

	def extract_variant_types(self, variants_df, gsvar_filename):
		"""
		Extract variant types and their associated coding and splicing

		Parameters
		----------
		variants_df : pandas.DataFrame
			variants dataframe to print unique variant types for
		gsvar_filename : str
			GSvar file name

		Returns
		-------
		None
		"""
		unique_type_filename = gsvar_filename.split(".GSvar")[0] + "_uniq_var_types.txt"
		print("Extract unique variant types into {}".format(unique_type_filename))
		with open(join(self.out_path, unique_type_filename), 'w') as unique_type_file:
			for unique_type in unique(variants_df['variant_type']):
				unique_type_file.write("uniq type: " + unique_type + "\n")
				unique_type_file.write("associated coding and splicing info: \n")
				for _, coding_and_splicing in variants_df.loc[
					variants_df["variant_type"] == unique_type, "coding_and_splicing"].iteritems():
					unique_type_file.write(str(coding_and_splicing) + "\n")
				unique_type_file.write("~~~\n")

	def filter_hl_genes(self, variants_df):
		"""
		Subset variant dataframe to only variants affecting known HL genes,
		if respective option of command line is on

		Parameters
		----------
		variants_df : pandas.DataFrame
			whole set of input variants

		Returns
		-------
		pandas.DataFrame
			filtered dataframe
		"""
		### ### ###
		# 1. First mark the HL genes
		# 2. Filter variants to the ones affecting HL gene (if applicable)
		### ### ###

		variants_df = variants_df.apply(self.mark_hl_genes, axis=1)
		if self.classify_only_hl_genes:
			self.logger.info("Filter variants affecting known HL genes")
			filtered_variants_df = variants_df[variants_df["known_HL_gene"]]
			self.logger.info(
				"Number of filtered variants affecting known HL genes: {}".format(filtered_variants_df.shape[0]))
			return filtered_variants_df
		else:
			self.logger.info(
				"Mark variants affecting genes which are not related to HL phenotype")
			return variants_df

	def mark_hl_genes(self, variant_row):
		"""
		Mark variant row for the presense of a HL-related gene

		Parameters
		----------
		variant_row : pandas.Series
			variant row to be examined for HL-related gene

		Returns
		-------
		pandas.Series
			variant row with marked HL-related gene
		"""
		gene = str(variant_row["gene"]).split(",")[0]
		if gene not in self.genes_inheritance_df.index:
			variant_row["ACMG_rule_flags"] = "Warning: Gene not in HL genes list, provided by ClinGen"
			variant_row["known_HL_gene"] = False
		else:
			variant_row["known_HL_gene"] = True
		return variant_row

	def mark_mt_genes(self, variant_row):
		"""
		Mark variants affecting mitochondrial genes

		Parameters
		----------
		variant_row : pandas.Series
			input variant row

		Returns
		-------
		pandas.Series
			variant row with warning mark
		"""
		chromosome_name = str(variant_row["chr"])
		if chromosome_name == "chrMT":
			if variant_row["ACMG_rule_flags"] == "-":
				variant_row[
					"ACMG_rule_flags"] = "Warning: ACMG guidelines for HL may not be applicable for mitochondrial variants"
			else:
				variant_row["ACMG_rule_flags"] = update_assignment(variant_row["ACMG_rule_flags"],
				                                                   "ACMG guidelines for HL may not be applicable for mitochondrial variants")
		return variant_row

	def examine_acmg_rules(self, variants_df):
		"""
		Examine each implemented ACMG rules,
		in a sequence similar to the one shown in flowchart 12.05.2020

		Parameters
		----------
		variants_df : pandas.DataFrame
			variants dataframe to examine rules for
		Returns
		-------
		pandas.DataFrame
			variants with assigned class after passing flowchart rules
		"""

		self.logger.info("Passing each variants from each flowchart rule")

		self.logger.info("Rule 1. Examine classification column")
		variants_df = variants_df.apply(self.examine_classified_class, axis=1)

		self.logger.info("Rule 2. Examine allele frequency (gnomAD column)")
		variants_df = variants_df.apply(self.examine_allele_freq, axis=1)

		self.logger.info("Rule 3. Examine quality parameters")
		variants_df = variants_df.apply(self.examine_quality_parameters, axis=1)

		self.logger.info("Rule 4. Examine MAF based on gene inheritance mode")
		time_start = timer()
		variants_df = variants_df.apply(self.examine_inheritance_specific_MAF, axis=1)
		time_end = timer()
		self.logger.info("PM2, BA1, BS1, PM2_Supporting and BS1_Supporting (if applicable): {}".format(
			sec2hour_min_sec(time_end - time_start)))

		### ### ### ### ###
		##  PS1 and PM5  ##
		### ### ### ### ###
		self.logger.info("Rule(s) 5. Pathogenic rules related to ClinVar")
		time_start = timer()
		variants_df = self.process_clinvar_rules(variants_df)
		time_end = timer()
		self.logger.info("PS1 and PM5: Elapsed time: {}".format(sec2hour_min_sec(time_end - time_start)))

		'''
		self.logger.info("Rule(s) 5. Pathogenic rules related to ClinVar")
		time_start = timer()
		# split variants dataframe into chunks for parallel processing
		ps1_pm5_bool_rows = variants_df.apply(self.select_ps1_pm5_variant_rows, axis=1)
		no_ps1_pm5_bool_rows = [not bool_row for bool_row in ps1_pm5_bool_rows]
		ps1_pm5_df = variants_df.loc[ps1_pm5_bool_rows]

		# use threads for PS1 and PM5 applicable variants
		variants_chunks = split_df4threads(variants_df, ps1_pm5_df, "PS1 and PM5", self.num_threads)
		# create our pool with `num_threads` processes
		pool = Pool(processes=self.num_threads)

		# apply assignment of ClinVar related rules to each dataframe chunk
		assigned_ps1_pm5_chunked = pool.map(self.process_clinvar_rules, variants_chunks)
		# save assigned rows of chunks to the whole dataframe
		for i in range(len(assigned_ps1_pm5_chunked)):
			variants_df.iloc[assigned_ps1_pm5_chunked[i].index] = assigned_ps1_pm5_chunked[i]

		# pool.close()
		# pool.join()
		# use main thread for variants not applicable for pvs1
		self.logger.info("Run PS1 and PM5 for variants with not applicable type")
		variants_df.loc[no_ps1_pm5_bool_rows] = self.process_clinvar_rules(variants_df.loc[no_ps1_pm5_bool_rows])
		time_end = timer()
		self.logger.info("PS1 and PM5: Elapsed time: {}".format(sec2hour_min_sec(time_end - time_start)))
		'''

		### ### ### ### ###
		## Refined PVS1  ##
		### ### ### ### ###
		self.logger.info("Rule 6. Very strong pathogenic, refine PVS1 rule based on Tayoun et al. Hum. Mutat. 2018")
		time_start = timer()
		variants_df = self.refine_pvs1_rule(variants_df)
		time_end = timer()
		self.logger.info("PVS1: Elapsed time: {}".format(sec2hour_min_sec(time_end - time_start)))

		'''
		## refine PVS1 based on bed files: protein critical regions and clinical significant exons - Oct 2020 ##
		self.logger.info("Rule 6. Very strong pathogenic, refine PVS1 rule based on Tayoun et al. Hum. Mutat. 2018")
		time_start = timer()

		# split dataframe into variants with type applicable for rule and the rest of variants
		pvs1_bool_rows = variants_df.apply(self.select_pvs1_variant_rows, axis=1)
		no_pvs1_bool_rows = [not bool_row for bool_row in pvs1_bool_rows]
		pvs1_df = variants_df.loc[pvs1_bool_rows]

		# use threads for PVS1 applicable variants
		# split variants dataframe into chunks for parallel processing
		variants_chunks = split_df4threads(variants_df, pvs1_df, "PVS1", self.num_threads)
		# create our pool with `num_threads` processes
		pool = Pool(processes=self.num_threads)

		# apply assignment of PVS1 rule to each dataframe chunk
		assigned_PVS1_chunked = pool.map(self.refine_pvs1_rule, variants_chunks)
		# save assigned rows of chunks to the whole dataframe
		for i in range(len(assigned_PVS1_chunked)):
			variants_df.iloc[assigned_PVS1_chunked[i].index] = assigned_PVS1_chunked[i]

		pool.close()
		pool.join()
		# use main thread for variants not applicable for PVS1
		self.logger.info("Run PVS1 for variants with not applicable type")
		variants_df.loc[no_pvs1_bool_rows] = self.refine_pvs1_rule(variants_df.loc[no_pvs1_bool_rows])
		time_end = timer()
		self.logger.info("PVS1: Elapsed time: {}".format(sec2hour_min_sec(time_end - time_start)))
		'''

		### ### ### ###
		#   PM1 rule  #
		### ### ### ###
		self.logger.info("Rule 7. Pathogenic moderate rule PM1")
		time_start = timer()
		variants_df = self.examine_pm1(variants_df)
		time_end = timer()
		self.logger.info("PM1: Elapsed time: {}".format(sec2hour_min_sec(time_end - time_start)))

		### ### ### ###
		#   PM4 rule  #
		### ### ### ###
		self.logger.info("Rule 8. Pathogenic moderate rule PM4")
		time_start = timer()
		variants_df = self.examine_pm4(variants_df)
		time_end = timer()
		self.logger.info("PM4: Elapsed time: {}".format(sec2hour_min_sec(time_end - time_start)))

		'''
		self.logger.info("Rule 8. Pathogenic moderate rule PM4")
		time_start = timer()
		# split dataframe into variants with type applicable for the rule and the rest of variants
		pm4_bool_rows = variants_df.apply(self.select_pm4_variant_rows, axis=1)
		no_pm4_bool_rows = [not pm4_bool_row for pm4_bool_row in pm4_bool_rows]
		pm4_df = variants_df.loc[pm4_bool_rows]

		# use threads for PM4 applicable variants
		# split variants dataframe into chunks for parallel processing
		variants_chunks = split_df4threads(variants_df, pm4_df, "PM4", self.num_threads)
		assigned_pm4_chunked = pool.map(self.examine_pm4, variants_chunks)
		# save assigned rows of chunks to the whole dataframe
		for i in range(len(assigned_pm4_chunked)):
			variants_df.iloc[assigned_pm4_chunked[i].index] = assigned_pm4_chunked[i]

		pool.close()
		pool.join()
		# use main thread for variants not applicable for PM4
		self.logger.info("Run PM4 for variant with not applicable type")
		variants_df.loc[no_pm4_bool_rows] = self.examine_pm4(variants_df.loc[no_pm4_bool_rows])
		time_end = timer()
		self.logger.info("PM4: Elapsed time: {}".format(sec2hour_min_sec(time_end - time_start)))
		'''

		### ### ### ### ###
		#     BP3 rule    #
		### ### ### ### ###
		self.logger.info("Rule 9. Benign supporting rule BP3")
		time_start = timer()
		variants_df = self.examine_bp3(variants_df)
		time_end = timer()
		self.logger.info("BP3: Elapsed time: {}".format(sec2hour_min_sec(time_end - time_start)))

		'''
		self.logger.info("Rule 9. Benign supporting rule BP3")
		time_start = timer()

		# split dataframe into variants with type applicable for rule and the rest of variants
		bp3_bool_rows = variants_df.apply(self.select_bp3_variant_rows, axis=1)
		no_bp3_bool_rows = [not bp3_bool_row for bp3_bool_row in bp3_bool_rows]
		bp3_df = variants_df.loc[bp3_bool_rows]

		# use threads for BP3 applicable variants
		# split variants dataframe into chunks for parallel processing
		variants_chunks = split_df4threads(variants_df, bp3_df, "BP3", self.num_threads)
		assigned_bp3_chunked = pool.map(self.examine_bp3, variants_chunks)
		# save assigned rows of chunks to the whole dataframe
		for i in range(len(assigned_bp3_chunked)):
			variants_df.iloc[assigned_bp3_chunked[i].index] = assigned_bp3_chunked[i]

		# pool.close()
		# pool.join()
		# use main thread for variants not applicable for BP3
		self.logger.info("Run BP3 for variants with not applicable type")
		variants_df.loc[no_bp3_bool_rows] = self.examine_bp3(variants_df.loc[no_bp3_bool_rows])
		time_end = timer()
		self.logger.info("BP3: Elapsed time: {}".format(sec2hour_min_sec(time_end - time_start)))
		'''

		### ### ### ###
		# PP3 and BP4 #
		### ### ### ###
		self.logger.info("Rules 10-11. Supporting rules, PP3 and BP4")
		variants_df = variants_df.apply(self.examine_pp3_bp4, axis=1)

		### ### ### ###
		#     BP7     #
		### ### ### ###
		self.logger.info("Rule 12. Benign supporting rule, BP7")
		variants_df = variants_df.apply(self.examine_bp7, axis=1)

		return variants_df

	def acmg_classify(self, variants):
		"""
		Link variant to a pathogenic class following ACMG recommendations

		Parameters
		----------
		variants : pandas.Series
			variant row

		Returns
		-------
		pandas.Series
			variant row with assigned pathogenicity class based examined ACMG rules
		"""
		variant_info = VariantInfo(self.current_sample_name, str(variants["gene"]).split(",")[0],
		                           str(variants["variant_type"]), None,
		                           str(variants["chr"]),
		                           str(variants["start"]), str(variants["end"]), str(variants["ref"]),
		                           str(variants["obs"]))
		odds_pathogenicity_very_strong, x, pathogenicity_prior = 350, 2, 0.1
		assignACMG = AssignACMG(variant_info, odds_pathogenicity_very_strong, x, pathogenicity_prior)
		preferred_inheritance, alternative_inheritance = assignACMG.prepare4classify(variants["ACMG_rules"],
		                                                                             variants["ACMG_rules_comment"])

		### ### #### #### ### ###
		# calculate ACMG class
		# and pathogenicity probability for available inheritance modes
		### ### #### #### ### ###
		variants["ACMG_class_preferred"], variants["ACMG_class_numeric_preferred"] = assignACMG.compute_categorical(
			"preferred", preferred_inheritance)
		variants["ACMG_class_prob_preferred"] = "{:1.5f}".format(
			assignACMG.compute_probability("preferred", preferred_inheritance))
		if alternative_inheritance:
			variants["ACMG_class_alternative"], variants[
				"ACMG_class_numeric_alternative"] = assignACMG.compute_categorical("alternative",
			                                                                       alternative_inheritance)
			variants["ACMG_class_prob_alternative"] = "{:1.5f}".format(
				assignACMG.compute_probability("alternative", alternative_inheritance))
		else:
			variants["ACMG_class_alternative"], variants["ACMG_class_numeric_alternative"], variants[
				"ACMG_class_prob_alternative"] = "NA", "NA", "NaN"
		return variants

	def select_pvs1_variant_rows(self, variants):
		"""
		Return if variant is applicable for PVS1 rule

		Parameters
		----------
		variants : pandas.Series
			variant row

		Returns
		-------
		bool
			variant is applicable for PVS1 rule (True), otherwise is not (False)
		"""

		return is_var_type4PVS1(str(variants["variant_type"]))

	def select_ps1_pm5_variant_rows(self, variants):
		"""
		Return if variant is applicable for PS1 and PM5 rules

		Parameters
		----------
		variants : pandas.Series
			variant row

		Returns
		-------
		bool
			variant is applicable for PS1 and PM5 (True), otherwise is not (False)
		"""
		if str(variants["variant_type"]) == "missense":
			return True
		else:
			return False

	def select_bp3_variant_rows(self, variants):
		"""
		Return if input variant row is applicable for BP3 rule

		Parameters
		----------
		variants : pandas.Series
			variant row

		Returns
		-------
		bool
			variant is applicable for BP3 (True), otherwise is not (False)
		"""
		if "inframe" in str(variants["variant_type"]) or (
				extract_rule_result(variants["ACMG_rules"], "PVS1") == "False" and "stop_gained" in str(
			variants["variant_type"])):
			return True
		else:
			return False

	def select_pm4_variant_rows(self, variants):
		"""
		Return if input variant row is applicable for PM4 rule

		Parameters
		----------
		variants : pandas.Series
			variant row

		Returns
		-------
		bool
			variant is applicable for PM4 (True), otherwise is not (False)
		"""
		if extract_rule_result(variants["ACMG_rules"], "PVS1") == "False" and (
				"inframe" in str(variants["variant_type"]) or "stop_lost" in str(variants["variant_type"])):
			return True
		else:
			return False

	def filter_pathogenic_variants(self, classified_variant_row, pathogenic_prob_threshold):
		"""
		Filter pathogenic variants from all classified variants

		Parameters
		----------
		classified_variant_row: pandas.Series
			classified variant row

		pathogenic_prob_threshold: float
			threshold on the pathogenicity probability to assume that a VUS variant is probably pathogenic

		Returns
		-------
		bool
			classified variant should be filtered (True), otherwise False
		"""
		### ### ### ###
		# filter pathogenic variants by
		# the classified class to be likely pathogenic or pathogenic
		# or the classified class to be VUS and the pathogenicity probability >= threshold (argument)
		### ### ### ###
		if classified_variant_row["ACMG_class_numeric_preferred"] == 4 or classified_variant_row[
			"ACMG_class_numeric_preferred"] == 5:
			return True
		elif classified_variant_row["ACMG_class_numeric_alternative"] == 4 or classified_variant_row[
			"ACMG_class_numeric_alternative"] == 5:
			return True
		elif classified_variant_row["ACMG_class_numeric_preferred"] == 3 and classified_variant_row[
			"ACMG_class_prob_preferred"] >= pathogenic_prob_threshold:
			return True
		elif classified_variant_row["ACMG_class_numeric_alternative"] == 3 and classified_variant_row[
			"ACMG_class_prob_alternative"] >= pathogenic_prob_threshold:
			return True
		else:
			return False

	def summarize_classification(self, classified_variants_gsvar_file):
		"""
		Summarize classification of variants by:
		Extracting pathogenic variants into a separate GSvar file

		Parameters
		----------
		classified_variants_gsvar_file : str
			classified variants gsvar file name

		Returns
		-------
		None
		"""
		self.logger.info("Summarize classification by first loading all classified variants")
		classified_variants_df, is_gsvar_empty = self.load_gsvar_df("output", classified_variants_gsvar_file)
		if not is_gsvar_empty and self.min_pathogenicity != -1:
			self.logger.info(
				"Extract all variants classified with the 4 or 5 ACMG class or variants of classified as VUS but their pathogenicity probability >={}".format(
					self.min_pathogenicity))
			filtered_pathogenic_variants_df = classified_variants_df[
				classified_variants_df.apply(self.filter_pathogenic_variants,
				                             pathogenic_prob_threshold=self.min_pathogenicity, axis=1)]
			self.logger.info(
				"Number of extracted pathogenic variant(s): {}".format(filtered_pathogenic_variants_df.shape[0]))
			if filtered_pathogenic_variants_df.shape[0] > 0:
				filtered_pathogenic_variants_filename = classified_variants_gsvar_file.split("_genotoscope.GSvar")[
					                                        0] + "_pathogenic_vars_genotoscope.GSvar"
				filtered_pathogenic_variants_df.to_csv(join(self.out_path, filtered_pathogenic_variants_filename),
				                                       sep="\t",
				                                       header=True, index=False)
				self.logger.info("Saving extracted pathogenic variants in %s",
				                 join(self.out_path, filtered_pathogenic_variants_filename))

	def assign_variants_class(self, gsvar_file_name, is_single_file_analysis=False):
		"""
		Assign variants class for input gsvar file

		Parameters
		----------
		gsvar_file_name : str
			GSvar file name
		is_single_file_analysis : bool
			A single GSvar file will be classified (True), otherwise a folder with GSvar files will be

		Returns
		-------
		str
			output gsvar file with classified variants
		"""
		self.current_sample_name = gsvar_file_name.split(".")[0].strip()
		updated_gsvar_file_name = self.process_column_line(gsvar_file_name, is_single_file_analysis)
		gsvar_df, is_gsvar_empty = self.load_gsvar_df("temp", updated_gsvar_file_name)

		if not is_gsvar_empty:
			### ### ### ### ### ### ### ### ### ### ### ###
			# if gsvar file is not empty
			# => continue with classification of variants
			### ### ### ### ### ### ### ### ### ### ### ###
			# init the two new columns
			gsvar_df["ACMG_rules"] = ""
			gsvar_df["ACMG_rules_comment"] = ""
			gsvar_df["ACMG_rule_flags"] = "-"
			gsvar_df["var_coding_seq"] = ""
			gsvar_df["known_HL_gene"] = ""

			### ### ###
			# Filter variants affecting known HL genes
			# Mark variants affecting mitochondrial genome
			### ### ###
			filtered_gsvar_df = self.filter_hl_genes(gsvar_df)
			self.logger.info("Mark variants affecting mitochondrial genome")
			filtered_gsvar_df = filtered_gsvar_df.apply(self.mark_mt_genes, axis=1)

			if filtered_gsvar_df.shape[0] != 0:
				### ### ###
				# if filtered-in variants exist,
				# continue with ACMG criteria examination and classification
				### ### ###

				### ### ###
				# examine ACMG rules
				### ### ###
				examined_gsvar_df = self.examine_acmg_rules(filtered_gsvar_df)

				### ### ###
				# remove not needed columns
				### ### ###
				examined_gsvar_df.drop("var_coding_seq", axis=1, inplace=True)
				examined_gsvar_df.drop("known_HL_gene", axis=1, inplace=True)

				### ### ###
				# classify variant by ACMG category
				### ### ###
				self.logger.info("ACMG class: Apply ACMG recommendations to combine met criteria and classify variants")
				classified_gsvar_df = examined_gsvar_df.apply(self.acmg_classify, axis=1)

				# preserve the initial order of the columns of the input GSvar file
				ordered_columns = self.current_gsvar_columns + ["ACMG_rules", "ACMG_rules_comment", "ACMG_rule_flags",
				                                                "ACMG_class_preferred", "ACMG_class_numeric_preferred",
				                                                "ACMG_class_prob_preferred", "ACMG_class_alternative",
				                                                "ACMG_class_numeric_alternative",
				                                                "ACMG_class_prob_alternative"]
				classified_gsvar_df = classified_gsvar_df[ordered_columns]

				# save csv file
				gsvar_out_filename = basename(gsvar_file_name).split(".GSvar")[0] + "_genotoscope.GSvar"
				classified_gsvar_df.to_csv(join(self.out_path, gsvar_out_filename), sep="\t", header=True, index=False)
				self.logger.info("Saving assigned classes for variants in %s", join(self.out_path, gsvar_out_filename))
				# delete processed duplicate file or original GSvar from temp folder
				remove(join(self.temp_path, updated_gsvar_file_name))
				self.logger.info("Removing processed duplicate file of original GSvar from temp folder")
				return gsvar_out_filename
			else:
				### ### ###
				# no variants are filtered-in,
				# therefore save an empty classification file
				### ### ###
				self.logger.info("Filtering for variants affecting HL-relate genes, resulted to no variants")
				gsvar_out_filename = basename(gsvar_file_name).split(".GSvar")[0] + "_genotoscope.GSvar"
				filtered_gsvar_df.to_csv(join(self.out_path, gsvar_out_filename), sep="\t", header=True, index=False)
				self.logger.info("Saving empty classification file in %s", join(self.out_path, gsvar_out_filename))
				# even if filtered-in variants are 0, delete processed duplicate file of original GSvar from temp folder
				remove(join(self.temp_path, updated_gsvar_file_name))
				self.logger.info("Removing processed duplicate file of original GSvar from temp folder")
				return gsvar_out_filename

		else:
			### ### ### ### ### ### ###
			# is gsvar is empty
			# => create empty genotoscope output file
			### ### ### ### ### ### ###
			gsvar_out_filename = basename(gsvar_file_name).split(".GSvar")[0] + "_genotoscope.GSvar"
			gsvar_df.to_csv(join(self.out_path, gsvar_out_filename), sep="\t", header=True, index=False)
			self.logger.info("Error - Input GSvar file was empty, saving empty GenOtoScope classification file in %s",
			                 join(self.out_path, gsvar_out_filename))
			return None

	def prepare_public_data_files(self):
		"""
		Prepare public data file

		Returns
		-------
		None
		"""
		# compute repeat regions that are not intersecting domain regions
		self.repeats_no_domains_file = self.prepare_annotation_files()
		# load inheritance mode of genes
		self.load_genes_inheritance_df()
		# load thresholds per inheritance mode
		self.load_thresholds_inheritance_df()

	def run_single_sample(self):
		self.logger.info("Classify variants and assign ACMG class for single GSvar file (sample)")
		self.prepare_public_data_files()
		is_single_file_analysis = True
		self.logger.info("==== ==== Start ==== ====")
		if is_gsvar_file_extension(self.gsvar_path):
			self.logger.info("######################### ~ #########################")
			self.logger.info("### Start:  %s ###", basename(self.gsvar_path))
			time_start = timer()
			output_gsvar_file = self.assign_variants_class(self.gsvar_path, is_single_file_analysis)
			time_end = timer()
			self.logger.info(
				"Elapsed time to classify GSvar: {}".format(sec2hour_min_sec(time_end - time_start)))
			self.summarize_classification(output_gsvar_file)
			self.logger.info("### End:  %s  ###", basename(self.gsvar_path))
		self.logger.info("==== ====  End  ==== ====")

	def run_multiple_samples(self):
		"""
		Run variant assignment for each sample/GSvar file

		Returns
		-------
		None
		"""
		self.logger.info("Classify variants and assign ACMG class for multiple samples")
		excluded_file_names = []
		self.prepare_public_data_files()

		self.logger.info("==== ==== Start ==== ====")
		with scandir(self.gsvar_path) as gsvar_dir:
			for content in gsvar_dir:
				if is_gsvar_file_extension(content.name) and content.is_file():
					self.logger.info("######################### ~ #########################")
					self.logger.info("### Start:  %s ###", content.name)
					time_start = timer()
					output_gsvar_file = self.assign_variants_class(content.name)
					time_end = timer()
					self.logger.info(
						"Elapsed time to classify GSvar: {}".format(sec2hour_min_sec(time_end - time_start)))
					self.summarize_classification(output_gsvar_file)
					self.logger.info("### End:  %s  ###", content.name)
				else:
					excluded_file_names.append(content.name)
		self.logger.info("==== ====  End  ==== ====")

		if len(excluded_file_names) > 0:
			self.logger.info("=== Excluded ===")
			self.logger.info("# Excluding following contents as they don't have GSvar extension:")
			self.logger.info("\n".join(excluded_file_names))
