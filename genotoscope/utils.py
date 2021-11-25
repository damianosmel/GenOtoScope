from os import makedirs
from genotoscope.variants.VariantInfo import VariantInfo
import argparse
import multiprocessing
import logging
from pandas import read_csv
from os.path import join, basename, splitext
from math import isnan
from numpy import mean
from Bio.Seq import Seq
from Bio.Data import IUPACData
import vcf
import yaml

logger = logging.getLogger("GenOtoScope.utils")


def create_dir(base_path):
	"""
	Create directory

	Parameters
	----------
	base_path : str
		full path of directory to be created

	Returns
	-------
	None
	"""
	makedirs(base_path, exist_ok=True)


def is_gsvar_file_extension(file_name):
	"""
	Check if file has GSvar extension

	Parameters
	----------
	file_name : str
		file name

	Returns
	-------
	bool
		True if file has gsvar extension, False otherwise
	"""

	if file_name[-6:] == ".GSvar":
		return True
	else:
		return False


def split_name_from_extension(file_name):
	"""
	Split full file name in name and extension

	Parameters
	----------
	file_name : str
		input file name to

	Returns
	-------
	list of str
		extracted file name and file extension
	"""
	return splitext(basename(file_name))


def trim_variant_pos_edit(var_pos_edit, process_left2right):
	"""
	Trim variant pos edit signature to remove not characters in the start or end of the variants

	Parameters
	----------
	var_pos_edit : str
		variant position edit signature
	process_left2right : bool
		trim from left to right or right to left

	Returns
	-------
	str
		trimmed variant pos edit signature
	"""
	if process_left2right:
		check_point = 0
		step = +1
	else:
		check_point = -1
		step = -1
	character = var_pos_edit[check_point]
	while not character.isalpha():
		check_point = check_point + step
		character = var_pos_edit[check_point]

	var_pos_edit_trimmed = None
	if process_left2right:
		var_pos_edit_trimmed = var_pos_edit[check_point:len(var_pos_edit)]
	else:
		var_pos_edit_trimmed = var_pos_edit[0:check_point + 1]
	var_pos_edit_trimmed = var_pos_edit_trimmed.strip()
	return var_pos_edit_trimmed


def parse_coding_splicing_column(transcript_splicing_info, gene_name, hgvs_parser):
	"""
	Parse coding and splicing column from GSvar file

	Parameters
	----------
	transcript_splicing_info : str
		transcript and splicing information (column) from GSVar
	gene_name : str
		name of gene containing current transcripts
	hgvs_parser : hgvs.parser
		initialized HGVS parser object

	Returns
	-------
	list of dict of str : str
		dictionary mapping transcript id -> transcript information
	"""
	logger.debug("Parse coding & splicing column")
	transcripts_info = []
	for transcript_info in transcript_splicing_info.split(","):
		transcript_input = transcript_info.split(":")
		transcript_info = {}
		if len(transcript_input) == 7 or len(transcript_input) == 8:
			transcript_info["gene_name"] = transcript_input[0]
			# remove the version of the transcript id
			transcript_info["transcript_id"] = transcript_input[1].split(".")[0]
			transcript_info["type_variant"] = set(extract_all_variant_types(transcript_input[2]))
			transcript_info["exon"] = transcript_input[4]
			# currently, exclude variant that do not contain any type of position of their transcript
			# or contain non coding position
			if "c." in transcript_input[5] and "*" not in transcript_input[5]:
				# parse coding positions
				# current transcript contains coding position of variant
				transcript_info["var_coding"] = hgvs_parser.parse_c_posedit(transcript_input[5].split("c.")[1])
				if "delins" not in str(transcript_info["var_coding"]):
					# parse single edit sequence
					var_edit = str(transcript_info["var_coding"].edit)
					if "del" in var_edit[0:3] or "dup" in var_edit[0:3]:
						transcript_info["var_seq"] = [
							transcript_input[5].split("c.")[1].split(str(transcript_info["var_coding"]))[1]]
					elif "ins" in var_edit[0:3]:
						transcript_info["var_seq"] = [var_edit.split("ins")[1]]
				else:
					# parse double edit (deletion and then insertion) sequences
					logger.debug("Parse double edit (deletion and then insertion)")
					coding_info = transcript_input[5].split("c.")[1]
					del_seq, ins_seq = coding_info.split('del')[1].split('ins')
					logger.debug('del: {}, ins: {}'.format(del_seq, ins_seq))
					transcript_info["var_seq"] = [del_seq, ins_seq]
				if transcript_input[6] == '':
					transcript_info["var_protein"] = None
				else:
					# transcript_info["var_protein"] = hgvs_parser.parse_p_posedit(transcript_input[6].split("p.")[1])
					transcript_info["var_protein"] = transcript_input[6].split("p.")[1]
			else:
				# exclude transcript if:
				# it does not contain any kind of position
				# it does not contain coding positions
				# exclude affected transcript containing the * character
				continue
		else:
			# print(
			# 	"Current version supports parsing of 8 pieces of information from the coding_and_transcript column")
			# print("A dictionary with mostly None values will be initialized and returned")
			transcript_info["gene_id"] = gene_name
			transcript_info["transcript_id"] = None
			transcript_info["type_variant"] = None
			transcript_info["var_coding"] = None
			transcript_info["exon"] = None
		transcripts_info.append(transcript_info)
	logger.debug("Parsed affected transcripts info: {}".format(transcripts_info))
	return transcripts_info


def extract_all_variant_types(variant_type_column):
	"""
	Extract all mentioned variant types from variant type column

	Parameters
	----------
	variant_type_column : str
		variant type column

	Returns
	-------
	list of str
		all variant types mentioned in the variant type column
	"""
	logger.debug("Extract all variant types")
	var_types_all = []
	if variant_type_column == "nan":
		# not available variant type
		var_types_all = []
	elif ',' in variant_type_column and '&' in variant_type_column:
		var_types = variant_type_column.split(",")
		for var_type in var_types:
			var_types_all += var_type.split("&")
	elif ',' in variant_type_column:
		var_types_all = variant_type_column.split(",")
	elif '&' in variant_type_column:
		var_types_all = variant_type_column.split("&")
	else:
		var_types_all = variant_type_column

	if type(var_types_all) is list:
		return var_types_all
	else:
		return [var_types_all]


def arrange_transcripts_per_variant_type(variant_types, transcripts_info, rule_specific_var_types):
	"""
	Arrange transcripts per variant type

	Parameters
	----------
	variant_types: list of str
		list of variant types
	transcripts_info : list of dict of str : str
		input transcripts information
	rule_specific_var_types : set of str
		variant types specific for currently running rule

	Returns
	-------
	dict of str : dict of str : str
		dictionary: key: variant type, value: all transcripts containing that variant type
	"""
	# logger.debug("Arrange transcripts per variant type")
	# filter variant types to keep only the ones that are important for the currently running rule
	var_types = set(variant_types)
	# start lost should be select before any other variant type so sort list in descending lexicographical order
	rule_types = list(rule_specific_var_types.intersection(var_types))

	if "start_lost" in rule_types:
		rule_types.remove("start_lost")
		rule_types = ["start_lost"] + sorted(rule_types)
	else:
		rule_types = sorted(rule_types)

	# arrange affected transcripts to variant type uniquely
	transcripts_per_var_type = {}
	previous_trans_matched_indices = []
	for var_type in rule_types:
		var_type_full = var_type.replace("'", "_prime_")
		if len(previous_trans_matched_indices) > 0:
			# remove from all transcripts the already matched ones
			matched_infos = [transcripts_info[matched_index] for matched_index in previous_trans_matched_indices]
			for matched_info in matched_infos:
				transcripts_info.remove(matched_info)

		previous_trans_matched_indices = []
		for transcript_idx, transcript_info in enumerate(transcripts_info):
			# check if transcript variant type matches current variant type
			trans_matches_var_type, trans_var_type_idx = False, 0
			# for trans_var_type in list(transcript_info['type_variant']):
			while trans_var_type_idx < len(transcript_info['type_variant']) and not trans_matches_var_type:
				trans_var_type = list(transcript_info['type_variant'])[trans_var_type_idx]
				if var_type_full in trans_var_type:
					trans_matches_var_type = True
				trans_var_type_idx += 1

			if trans_matches_var_type:
				# update matched transcripts
				previous_trans_matched_indices.append(transcript_idx)
				# arrange transcript on variant type and
				# remove transcript for next variant types
				if var_type not in transcripts_per_var_type:
					transcripts_per_var_type[var_type] = [transcript_info]
				else:
					transcripts_per_var_type[var_type].append(transcript_info)
	# logger.debug("~~~ ~~~")
	return transcripts_per_var_type


def extract_codons(sequence):
	"""
	Extract codons from given sequence, which start where the frame start offset is

	Parameters
	----------
	sequence : str
		sequence to extract codons from

	Returns
	-------
	list of str
		input sequence codons
	"""
	codons = []
	for pos in range(0, len(sequence), 3):
		codons.append(sequence[pos:pos + 3])
	return codons


def select_max_strength(assigned_evidences):
	"""
	Select assigned evidence with highest evidence for all assigned

	Parameters
	----------
	assigned_evidences : list of str
		assigned evidences

	Returns
	-------
	str
		evidence with maximum strength
	"""
	strength2weight = {"PVS1": 4, "PVS1_Strong": 3, "PVS1_Moderate": 2, "PVS1_Supporting": 1, "PM5_Strong": 3,
	                   "PM5_Moderate": 2}
	evidence_weights = [strength2weight[evid] for evid in assigned_evidences]
	max_strength_idx = max((v, i) for i, v in enumerate(evidence_weights))[1]
	return assigned_evidences[max_strength_idx]


def aggregate_examined_rules(assigned_classes, assignment_comments, rule_name):
	"""
	Aggregate assigned classes and assignment comments

	Parameters
	----------
	assigned_classes: list of str
		assigned classes per affected transcript
	assignment_comments : list of str
		comment per affected transcript assignment
	rule_name : str
		rule name corresponding to the assigned classes

	Returns
	-------
	str
		aggregated class
	str
		aggregated assignment comment
	"""
	logger.debug("Aggregate assigned classes and comments")
	logger.debug("Assigned classes: {}".format(assigned_classes))
	# logger.debug("Assignment comments: {}".format(assignment_comments))
	if len(assigned_classes) > 1:
		unique_assigned_classes = set(assigned_classes)
		if len(unique_assigned_classes.difference(set(["False"]))) == 0:
			aggregated_class = "False" + "|" + "|".join(assigned_classes)
			aggregated_comment = "|".join(assignment_comments)
		else:
			triggered_evidences = list(unique_assigned_classes.difference(set(["False"])))
			if len(triggered_evidences) > 1:
				max_strength_evid = select_max_strength(triggered_evidences)
			else:
				max_strength_evid = triggered_evidences[0]
			aggregated_class = max_strength_evid + "|" + "|".join(assigned_classes)
			aggregated_comment = "|".join(assignment_comments)
		return rule_name + ": " + aggregated_class, rule_name + ": " + aggregated_comment
	else:
		return rule_name + ": " + assigned_classes[0], rule_name + ": " + assignment_comments[0]


def update_assignment(current_assignment, new_assignment):
	"""
	Update assignment column value

	Parameters
	----------
	current_assignment : str
		current assignment
	new_assignment : str
		new assignment

	Returns
	-------
	str
		updated assignment value
	"""
	return current_assignment + "||" + new_assignment


def var_region_intersects_annotation(var_region_start, var_region_end, annot_start, annot_end):
	"""
	Examine if variant-affected region intersects protein annotation

	Parameters
	----------
	var_region_start : int
		variant-affected region start position
	var_region_end: int
		variant-affected region end position
	annot_start : int
		annotation region start position
	annot_end : int
		annotation region end position

	Returns
	-------
	bool
		variant-affected region intersects protein annotation (True), otherwise (False)
	"""
	# print("Examine if variant-affected region (start:{}, end:{}) intersects annotation region (start:{}, end:{})".format(var_region_start,var_region_end,annot_start,annot_end))

	### ### ### ### ### ### ### ###
	# 2 cases for variant intersects annotation:
	# 1) variant-affected region smaller than annotation region
	# |=== variant region ===|
	# |=== === annotation region === ===|
	# 2) annotation region larger or equal to variant-affected region
	# |=== == variant region == ===|
	# |=== annotation region ===|
	### ### ### ### ### ### ### ###
	is_var_region_in_annotation = False
	var_region_length = var_region_end - var_region_start + 1
	annot_length = annot_end - annot_start + 1

	if var_region_length < annot_length:
		# 1) variant-affected region smaller than annotation region
		if var_region_start <= annot_end:
			# if variant-affected region starts before annotation end,
			# check if intersecting region is at least 25% of the variant-affected region
			#                               |=== variant region ===|
			#   |=== === annotation region === ===|
			intersecting_length = annot_end - var_region_start + 1
			if intersecting_length / var_region_length >= 0.25:
				is_var_region_in_annotation = True
	else:
		# 2) annotation region larger or equal to variant-affected region
		if annot_start <= var_region_end:
			# if annotation starts before variant-affected region,
			# check if the intersecting region is at least 25% of the annotation region
			# |=== === === variant region === === ===|
			#                           |=== annotation region ===|
			intersecting_length = var_region_end - annot_start + 1
			if intersecting_length / annot_length >= 0.25:
				is_var_region_in_annotation = True
	# if is_var_region_in_annotation:
	# 	print("variant region ({},{}) intersects annotation region ({},{})".format(var_region_start,var_region_end,annot_start,anot_end))
	return is_var_region_in_annotation


def split_df4threads(all_variants_df, rule_variants_df, rule_name, num_threads):
	"""
	Split rule-specific variants sub-dataframe into chunks.
	The number of chunks equals the number of threads

	Credits: https://stackoverflow.com/questions/40357434/pandas-df-iterrows-parallelization

	Parameters
	----------
	all_variants_df: pandas.DataFrame
		all variants dataframe
	rule_variants_df : pandas.DataFrame
		subset dataframe of variants applicable for rule
	rule_name : str
		ACMG rule name
	num_threads : int
		number of threads to used for splitting

	Returns
	-------
	List of list of int
		list of chunks, each chunk is a list of the corresponding indices of the dataframe
	"""
	if num_threads == 0:
		# create as many threads as there are machine CPUs minus one
		num_processes = multiprocessing.cpu_count() - 1
	# calculate the chunk size as an integer
	chunk_size = int(rule_variants_df.shape[0] / num_threads)
	logger.info(
		"Rule: {}, for {} threads, split dataframe with {} variant(s), into chunks of {} variant(s) per thread".format(
			rule_name,
			num_threads, rule_variants_df.shape[0],
			chunk_size))
	# this solution was reworked from the above link.
	# will work even if the length of the dataframe is not evenly divisible by num_threads
	chunks = [all_variants_df.iloc[rule_variants_df.index[i:i + chunk_size]] for i in
	          range(0, rule_variants_df.shape[0], chunk_size)]
	return chunks


def is_var_type4PVS1(var_type):
	"""
	Check if variant type is appropriate for PVS1 assessment

	Parameters
	----------
	var_type : str
		variant type column

	Returns
	-------
	bool
		variant type is for PVS1 assessment (True), False otherwise
	"""
	if "start_lost" in var_type:
		return True
	elif "splice_acceptor" in var_type or "splice_donor" in var_type:
		return True
	elif "stop_gained" in var_type or "stop_lost" in var_type or "frameshift" in var_type or "inframe_insertion" in var_type or "inframe_deletion" in var_type:
		return True
	else:
		return False


def sec2hour_min_sec(seconds):
	"""
	Convert elapsed seconds to hours, minutes and seconds
	#Credits: https://codereview.stackexchange.com/questions/174796/convert-seconds-to-hours-minutes-seconds-and-pretty-print

	Parameters
	----------
	seconds : long int
		elapsed seconds

	Returns
	-------
	str
		string with converted hours, minutes, seconds and microseconds
	"""
	microseconds = int(seconds * 1000000)

	if microseconds != 0:
		seconds, microseconds = divmod(microseconds, 1000000)
		minutes, seconds = divmod(seconds, 60)
		hours, minutes = divmod(minutes, 60)
		periods = [('hour(s)', hours), ('minute(s)', minutes), ('second(s)', seconds), ('microsecond(s)', microseconds)]
		return ', '.join('{} {}'.format(value, name) for name, value in periods if value)
	else:
		return str(microseconds) + ' microseconds'


def extract_rule_result(examined_rules, rule_name):
	"""
	Parameters
	----------
	examined_rules : str
		examined ACMG rules
	rule_name : str
		rule to extract its assigned result

	Returns
	-------
	str
		assignment result of specific rule
	"""
	[_, rules_result] = examined_rules.split(rule_name + ":")
	rule_result = rules_result.split("||")[0]
	return rule_result.split("|")[0].strip()


def get_codon_index(seq_codons, target_codon):
	"""
	Search target codon on sequence codons
	and return its index

	Parameters
	----------
	seq_codons : list of str
		sequence codons
	target_codon : str
		codon to search for

	Returns
	-------
	int
		target codon index (1-based)
	"""
	try:
		codon_index = seq_codons.index(target_codon) + 1
	except ValueError:
		codon_index = len(seq_codons)
	return codon_index


def compose_clinvar_change_signature(ref, pos, alt):
	"""
	Create (DNA)change signature of clinvar record
	composed by its reference base, position and observed base
	Parameters
	----------
	ref: str
	reference base of clinvar record
	pos: int
	change position of clinvar record
	alt: str
	observed base of clinvar record

	Returns
	-------
	str
		DNA change signature
	"""
	return "_".join([str(ref), str(pos), str(alt)])


def is_var2exclude(pheno_high_af_variants, current_variant):
	"""
	Check if variant should be excluded by BA1 and BS1 rules
	based on ACMG specifications for hearing loss, Oza et al. 2018

	Parameters
	----------
	current_variant : VariantInfo
		currently investigated variant
	pheno_high_af_variants: vcf.Reader
		phenotype relevant high AF variants

	Returns
	-------
	bool
		variant included in the exclusion list (True), otherwise is not included (False)
	"""
	logger.debug("Check for exclusion of current variant: {}".format(current_variant.to_string()))
	exclude_var = False
	if pheno_high_af_variants:
		# use 0-based half-open coordinate system
		# to extract variants from the vcf file (exclusion list)
		try:
			matching_pos_high_af_vars = list(
				pheno_high_af_variants.fetch(current_variant.chrom.split("chr")[1],
				                             int(current_variant.genomic_start) - 1,
				                             int(current_variant.genomic_end)))
		except ValueError:
			matching_pos_high_af_vars = []

		# check exclusion variants matching by position
		# if they match also the ref and obs sequence of the current variant
		for candidate in matching_pos_high_af_vars:
			logger.debug("candidate: {}".format(candidate.POS))
			candidate_var = VariantInfo("", None, None,
			                            None, "chr" + str(candidate.CHROM), str(candidate.POS), str(candidate.POS),
			                            str(candidate.REF), str(candidate.ALT[0]))
			logger.debug("candidate exclusion variant: {}".format(candidate_var.to_string()))
			if current_variant.is_same_var(candidate_var):
				logger.debug("Current variant matches known pathogenic variant with high AF")
				exclude_var = True
	return exclude_var


def format_dict2str(map):
	"""
	Convert dictionary to string

	Parameters
	----------
	map : dict of str: str
		dictionary to be converted

	Returns
	-------
	str
		string formatted dictionary
	"""
	if map:
		logger.debug(map)
		return ",".join([key + ":" + value for key, value in map.items()])
	else:
		return ""


def deformat_str2dict(text):
	"""
	Convert str formated dictionary to dictionary

	Parameters
	----------
	text: str
		string formated dictionary

	Returns
	-------
	dict of str: str
		equivalent dictionary
	"""
	if text != "":
		map = {}
		for pair in text.strip().split(","):
			key, value = pair.split(":")
			map[key] = value
		return map
	else:
		return {}


def load_hugo_genes_df(data_path, hugo_genes_filename):
	"""
	Load hugo genes dataframe

	Parameters
	----------
	data_path : os.path
		path to data directory
	hugo_genes_filename : str
		hugo genes file name

	Returns
	-------
	pandas.DataFrame
		loaded hugo genes data frame
	"""
	hugo_genes_df = read_csv(join(data_path, hugo_genes_filename),
	                         sep='\t',  # field separator
	                         header=0,  # get header
	                         # index_col=[0, 1, 2],  # set as index the chr, start stop
	                         skipinitialspace=True,
	                         skip_blank_lines=True,
	                         error_bad_lines=False,
	                         warn_bad_lines=True
	                         )
	hugo_genes_df.set_index('symbol', inplace=True)
	return hugo_genes_df


def normalize_gene_info(clinvar_rec):
	"""
	Normalize ClinVar gene information

	Parameters
	----------
	clinvar_rec : dict
		ClinVar record organized in dictionary

	Returns
	-------
	list of str, int
		list of genes containing [gene symbol, gene id]
	"""
	logger.debug("Normalize gene information")
	if 'GENEINFO' in clinvar_rec.INFO:
		genes_info = []
		for gene in clinvar_rec.INFO['GENEINFO'].strip().split("|"):
			genes_info.append(gene.split(':'))
		return genes_info
	else:
		return None


def normalize_disease_id(clinvar_rec):
	"""
	Normalize disease id for input clinvar record

	Parameters
	----------
	clinvar_rec : dict
		clinvar record organized in a dictionary

	Returns
	-------
	str
		normalized disease id
	"""
	if 'CLNDISDB' in clinvar_rec.INFO:
		if not clinvar_rec.INFO['CLNDISDB']:
			disease_db_id = ','.join(clinvar_rec.INFO['CLNDISDB'])
		else:
			disease_db_id = None
	else:
		disease_db_id = None
	return disease_db_id


def convert_review_status2stars(clinvar_stars_df, star_status2int, clinvar_rev_status):
	"""
	Convert CLNREVSTAT (review status) tab from clinvar vcf file to number of review stars
	for unknown description -> star= -1
	ClinVar review status documentation: https://www.ncbi.nlm.nih.gov/clinvar/docs/review_status/

	Parameters
	----------
	clinvar_stars_df : pandas.DataFrame
		clinvar stars dataframe
	star_status2int : dict of str : int
		map clinvar star status to integer value
	clinvar_rev_status : list of str

	Returns
	-------
	int
		clinvar gold stars
	"""
	rev_status = [review_elem.replace('_', ' ') for review_elem in clinvar_rev_status]
	rev_status = ",".join(rev_status)
	if rev_status not in clinvar_stars_df.Review_status.values:  # if retrieved status not in status
		return star_status2int["unknown review status"]
	else:
		return star_status2int[
			clinvar_stars_df.loc[clinvar_stars_df["Review_status"] == rev_status][
				"Number_of_gold_stars"].iloc[0]]


def convert_status2stars(clinvar_stars_df, star_status2int, clinvar_status):
	"""
	Convert clinvar status to clinvar star using the dictionary
	to convert the clinvar star text to clinvar star integer value
	Use the output for assignment comment

	for unknown description -> star= -1

	Parameters
	----------
	clinvar_stars_df : pandas.DataFrame
		clinvar stars dataframe
	star_status2int : dict of str : int
		map clinvar star status to integer value
	clinvar_status : str
		clinvar status to be mapped

	Returns
	-------
	int
		clinvar star
	"""
	if clinvar_status not in clinvar_stars_df.Review_status.values:  # if retrieved status not in status
		return star_status2int["unknown review status"]
	else:
		return star_status2int[
			clinvar_stars_df.loc[clinvar_stars_df["Review_status"] == clinvar_status][
				"Number_of_gold_stars"].iloc[0]]


def get_clinvar_strand(hugo_genes_df, gene_symbol):
	"""
	Get strand for gene found in clinvar entry

	Parameters
	----------
	hugo_genes_df : pandas.DataFrame
		hugo genes dataframe
	gene_symbol : str
		gene symbol found in ClinVar entry

	Returns
	-------
	str
		strand of ClinVar gene
	"""
	logger.debug("Get strand for gene found in ClinVar entry")
	if gene_symbol in hugo_genes_df.index:
		return hugo_genes_df.loc[gene_symbol].strand
	else:
		return None


def deactivate_lower_strength_rules(examined_rules, ACMG_rules, ACMG_rules_comment):
	"""
	Deactivate lower strength rules to avoid double counting (PM2, PM2 supporting, BS1, BS1 supporting)
	For example if PM2 is triggered then do not trigger PM2_Supporting

	Parameters
	----------
	examined_rules: dict of str : list of int
		examined rule names mapped to list of 0/1 elements each showing activating or not from a subpopulation
	ACMG_rules: list of str
		assignment of examined rules
	ACMG_rules_comment: list of str
		comment on rule assignments

	Returns
	-------
	list of str
		updated assignment of examined rule
	list of str
		updated comments on rule assignments
	dict of list of int
		updated activating or not activating gnomAD subpopulation per rule
	"""
	logger.debug("Correct assignment of PM2 and BS1 strength rules for double counting")
	examined_rule_names = list(examined_rules.keys())
	if "PM2" in examined_rule_names and "PM2_Supporting" in examined_rule_names:
		# deactivate PM2_Supporting if PM2 is set
		pm2_idx = examined_rule_names.index("PM2")
		pm2_sup_idx = examined_rule_names.index("PM2_Supporting")
		if ACMG_rules[pm2_idx].split(":")[1].strip() == "True" and ACMG_rules[pm2_sup_idx].split(":")[
			1].strip() == "True":
			ACMG_rules[pm2_sup_idx] = "PM2_Supporting: False"
			ACMG_rules_comment[pm2_sup_idx] = ACMG_rules_comment[pm2_sup_idx] + " (deactivated as PM2 triggered)"
			examined_rules["PM2_Supporting"] = [0]
	if "BS1" in examined_rule_names and "BS1_Supporting" in examined_rule_names:
		# deactivate BS1_Supporting if BS1 is set
		bs1_idx = examined_rule_names.index("BS1")
		bs1_sup_idx = examined_rule_names.index("BS1_Supporting")
		if ACMG_rules[bs1_idx].split(":")[1].strip() == "True" and ACMG_rules[bs1_sup_idx].split(":")[
			1].strip() == "True":
			ACMG_rules[bs1_sup_idx] = "BS1_Supporting: False"
			ACMG_rules_comment[bs1_sup_idx] = ACMG_rules_comment[bs1_sup_idx] + " (deactivated as BS1 triggered)"
			examined_rules["BS1_Supporting"] = [0]

	return ACMG_rules, ACMG_rules_comment, examined_rules


def pad_seq(sequence):
	"""
	Pad sequence to multiple of 3 with N
	Credits: https://stackoverflow.com/questions/53894575/how-can-i-fix-this-error-biopythonwarning-partial-codon-lensequence-not-a

	Parameters
	----------
	sequence : Bio.Seq
		sequence object to pad with N's

	Returns
	-------
	Bio.Seq
		padded sequence object
	"""

	remainder = len(sequence) % 3
	return sequence if remainder == 0 else sequence + 'N' * (3 - remainder)


def normalize_codon_nucleotides(nucl_per_codon_pos):
	"""
	Normalize nucleotides registered for the three positions of a codon

	Parameters
	----------
	nucl_per_codon_pos: list of str
		nucleotides per codon position

	Returns
	-------
	"""
	logger.debug("Normalize nucleotides: {}, registered on one codon".format(nucl_per_codon_pos))
	seq2translate = ''
	nucl_pos = (nucl for nucl in nucl_per_codon_pos)
	while len(seq2translate) < 3:
		seq2translate += next(nucl_pos)

	logger.debug("Normalized pos: {}".format(seq2translate))
	logger.debug("Padded: {}".format(pad_seq(seq2translate)))
	return pad_seq(seq2translate)


def parse_revel_pred(variants):
	"""
	Parse REVEL score from variant row
	and return its pathogenicity prediction
	Parameters
	----------
	variants : pandas.Series
		variant row
	Returns
	-------
	float
		REVEL score predicts pathogenic variant (1), benign (0), grey zone (0.5) otherwise -1
	"""
	if isnan(variants["REVEL"]):
		return -1
	elif float(variants["REVEL"]) <= 0.15:
		return 0
	elif float(variants["REVEL"]) >= 0.7:
		return 1
	else:
		return 0.5


def parse_cadd_pred(variants):
	"""
	Parse CADD score from variant row
	and return its pathogenicity prediction

	Parameters
	----------
	variants : pandas.Series
		variant row

	Returns
	-------
	int
		CADD score predicts pathogenic mutation (1), benign mutation (0), otherwise -1
	"""
	if isnan(variants["CADD"]):
		return -1
	elif float(variants["CADD"]) > 20.0:
		return 1
	else:
		return 0


def is_variant_pathogenic(variants, predictor_name):
	"""
	Examine if variant is pathogenic using REVEL or CADD score

	Parameters
	----------
	variants: pandas.Series
		variant row
	predictor_name : str
		predictor name

	Returns
	-------
	float
		variant is pathogenic (1), benign (0), grey zone (0.5), otherwise -1
	str
		predictor name
	"""
	logger.debug("Examine predictors for variant pathogenicity")
	if predictor_name in variants:
		if predictor_name == "REVEL":
			return parse_revel_pred(variants), "REVEL"
		elif predictor_name == "CADD":
			return parse_cadd_pred(variants), "CADD"
	elif predictor_name == "REVEL" and predictor_name not in variants:
		# if the selected predictor is REVEL, but is not available try for CADD
		if "CADD" in variants:
			return parse_cadd_pred(variants), "CADD"
		else:
			return -1, None
	else:
		return -1, None


def aggregate_patho_predictions(patho_predictions):
	"""
	Aggregate pathogenicity score from different predictors
	to assess total pathogenicity
	Logic: average of all non negative pathogenicity scores,
	if positive sum => Pathogenic, else => Benign

	Parameters
	----------
	patho_predictions : dict of str to float

	Returns
	-------
	str
		after aggregation variant is voted as pathogenic, benign or unknown (if none of the predictors gave evidence)
	"""
	no_negative_predictions = []
	for predictor, prediction in patho_predictions.items():
		if prediction > -1:
			no_negative_predictions.append(prediction)
	if len(no_negative_predictions) == 0:
		return "unknown"
	elif mean(no_negative_predictions) >= 1:
		return "pathogenic"
	else:
		return "benign"


def is_variant_conserved(variants):
	"""
	Examine if variant is conserved using the PhyloP score
	and the cutoff of 1.6

	Parameters
	----------
	variants : pandas.Series
		variant row

	Returns
	-------
	int
		variant is conserved (1), is not conserved (0), value does not exist -1
	"""
	if variants["phyloP"]:
		if float(variants["phyloP"]) > 1.6:
			return 1
		else:
			return 0
	else:
		return -1


def parse_maxentscan_pred(variants):
	"""
	Parse MaxEntScan column for current variant row

	Parameters
	----------
	variants : pandas.Series
		variant row to extract MaxEntScan

	Returns
	-------
	list of float
		list of ratios that affect splicing
	"""
	if variants["MaxEntScan"] and ">" in str(variants["MaxEntScan"]):
		# if MaxEntScan not null and it is parse-able
		# calculate MaxEntScan ratio
		impacting_splicing_ratios, no_impacting_splicing_ratios = [], []
		for ref_obs_scores in str(variants["MaxEntScan"]).split(","):
			ref_score, obs_score = ref_obs_scores.split(">")
			if ref_score.isdigit() and obs_score.isdigit():
				ref_obs_ratio = abs(float(obs_score) - float(ref_score) / float(ref_score))
			else:
				ref_obs_ratio = -1.0
			if ref_obs_ratio >= 0.15:
				impacting_splicing_ratios.append(str(ref_obs_ratio))
			else:
				no_impacting_splicing_ratios.append(str(ref_obs_ratio))
	else:
		# no MaxEntScan score found
		impacting_splicing_ratios = None
	return impacting_splicing_ratios


def parse_dbscsnv_pred(variants):
	"""
	Parse dbscSNV splicing predictor column of current variant row

	Parameters
	----------
	variants : pandas.Series
		variant row to extract dbscsnv score

	Returns
	-------
	list of float
		list of score that support splicing effect
	"""
	impacting_splicing_scores = []
	if variants["dbscSNV"] and "/" in str(variants["dbscSNV"]):
		ada_score, rf_score = str(variants["dbscSNV"]).split("/")
		if ada_score.isdigit() and float(ada_score) > 0.6:
			impacting_splicing_scores.append(float(ada_score))
		if rf_score.isdigit() and float(rf_score) > 0.6:
			impacting_splicing_scores.append(float(rf_score))
	else:
		# no dbscSNV score found
		impacting_splicing_scores = None
	return impacting_splicing_scores


def aggregate_splicing_predictions(splicing_predictions):
	"""
	Aggregate splicing predictions to assess splicing impact
	Logic: all predictors, with score, should have a impacting score

	Parameters
	----------
	splicing_predictions: dict of str to list of float
		impacting splicing predictions

	Returns
	-------
	bool
		aggregation of predictor votes for splicing impact (True), otherwise (False)
	"""
	sum_impact = 0
	for predictor, impacting_splicing_scores in splicing_predictions.items():
		if impacting_splicing_scores is not None and len(impacting_splicing_scores) > 0:
			sum_impact = sum_impact + 1
		else:
			sum_impact = sum_impact - 1
	if sum_impact > 0:
		return True
	else:
		return False


def prepare_parameters(settings_file_path, root_dir):
	"""
	Prepare run parameters by parsing settings yaml file

	Parameters
	----------
	settings_file_path : str
		settings yaml file path
	root_dir : str
		analysis root directory


	Returns
	-------
	dict of str : str or int
		parsed data settings
	str
		data file path
	str
		temporary data file path
	str
		critical regions of proteins
	"""
	# read data settings yaml file
	with open(settings_file_path) as settings_file:
		data_settings = yaml.safe_load(settings_file)

	# create the central data path
	data_path = join(root_dir, "genotoscope_data")

	# create output folders
	# temporary folder
	temp_path = join(root_dir, "genotoscope_temp")

	# create version of critical regions of proteins
	critical_prot_version = data_settings["clinvar_version"] + "_" + data_settings["uniprot_version"]
	return data_settings, data_path, temp_path, critical_prot_version


def is_dna(sequence, exception_char):
	"""
	Validate if input sequence is a DNA sequence (plus the character -)

	Parameters
	----------
	sequence : str
		input sequence
	exception_char : str
		additional accepted character

	Returns
	-------
	bool
		input is DNA sequence (True), otherwise False
	"""
	for base in sequence.upper():
		if base not in IUPACData.unambiguous_dna_letters and base != exception_char:
			return False
	return True


def acmg_column2dict(acmg_column, is_rules_column):
	"""
	Convert ACMG column to dictionary format

	Parameters
	----------
	acmg_column : str
		examined rules column
	is_rules_column : bool
		input column is the rules column (True), otherwise False

	Returns
	-------
    dict of str : str
        map each rule name to examination result
	"""
	rule2assignment = {}
	is_preferred_seen = False
	if is_rules_column:
		### #### ###
		# when parsing ACMG assignments column
		# initialize preferred and alternative inheritance modes
		### #### ###
		rule2assignment['preferred_inheritance'] = 'NA'
		rule2assignment['alternative_inheritance'] = 'NA'
	for rule in acmg_column.split("||"):
		if "NGSD:" in rule or "AF:" in rule or "QUAL:" in rule:
			# exclude not ACMG rules from parsed output
			continue
		elif "=" in rule and is_rules_column:
			### ### ###
			# extract inheritance related rules
			### ### ###
			if not is_preferred_seen:
				preferred_inheritance, rules_inheritance = rule.split("=")
				rule2assignment['preferred_inheritance'] = preferred_inheritance
				for rule_inheritance in rules_inheritance.split("|"):
					rule_name, rule_assignment = split_name_assignment(rule_inheritance)
					rule2assignment[rule_name + "_preferred"] = rule_assignment
				# logger.debug(
				#	"Preferred: parsed rule\n name: {}, assignment: {}".format(rule_name+"_preferred", rule_assignment))
				is_preferred_seen = True
			else:
				alternative_inheritance, rules_inheritance = rule.split("=")
				rule2assignment['alternative_inheritance'] = alternative_inheritance
				for rule_inheritance in rules_inheritance.split("|"):
					rule_name, rule_assignment = split_name_assignment(rule_inheritance)
					rule2assignment[rule_name + "_alternative"] = rule_assignment
		# logger.debug(
		#	"Alternative: parsed rule\n name:{}, assignment: {}".format(rule_name+"_alternative", rule_assignment))
		elif "Inheritance_specific_rules:" in rule:
			# comments of ACMG rules for preferred inheritance
			if "preferred:" in rule:
				# extract preferred inheritance if two modes are available
				inheritance_info, rules_comment_preferred = rule.split("preferred:")
			else:
				# extract single inheritance mode (preferred)
				inheritance_info, rules_comment_preferred = rule.split("single_inheritance:")
			rule2assignment['inheritance_general_comment'] = inheritance_info.split("Inheritance_specific_rules: ")[1]
			rules_comment_preferred = rules_comment_preferred.split(",")[1]
			for rule_comment_preferred in rules_comment_preferred.split("|"):
				rule_name, rule_assignment = split_name_assignment(rule_comment_preferred)
				rule2assignment[rule_name + "_preferred"] = rule_assignment
		elif "alternative:" in rule:
			# comments of ACMG rules for alternative inheritance
			inheritance_info, rules_comment_alternative = rule.split("alternative:")
			rules_comment_alternative = rules_comment_alternative.split(",")[1].strip()
			for rule_comment_alternative in rules_comment_alternative.split("|"):
				rule_name, rule_assignment = split_name_assignment(rule_comment_alternative)
				rule2assignment[rule_name + "_alternative"] = rule_assignment
		elif ":" in rule:
			# assignments or comments for ACMG rules not related to inheritance
			[rule_name, rule_assignment] = rule.split(": ", 1)
			### ### ### ###
			# format rule assignment to separate total assignment #
			# from assignments per transcript                     #
			### ### ### ###

			if "|" in rule_assignment:
				if is_rules_column:
					rule_assignment = rule_assignment.replace("|", " \nTranscripts:\n", 1)
					rule_assignment = rule_assignment.replace("|","\n")
					rule_assignment = "Total:\n" + rule_assignment
				else:

					if rule_name == "PM5" or rule_name == "PS1":
						#rule_assignment = rule_assignment.replace(";","\n")
						rule_assignment = rule_assignment.replace(" Supporting pathogenic ClinVar entries: ","Supporting pathogenic ClinVar entries:\n").replace(";",'\n').replace("|E","\nE")
					else:
						rule_assignment = "Per transcript:\n" + rule_assignment.replace("|", "\n")
			rule2assignment[rule_name] = rule_assignment
	return rule2assignment


def acmg_rules2dict(rules_assignments, rules_comments):
	"""
	Convert ACMG columns (ACMG_rules & ACMG_rules_comment)
	to a merged dictionary

	Parameters
	----------
	rules_assignments : str
		rules assignments
	rules_comments : str
		rules assignments

	Returns
	-------
	dict of str : str
		merged dictionary of ACMG columns
	"""
	### ### ###
	# convert activations and comments to dictionaries
	### ### ###
	assignments = acmg_column2dict(rules_assignments, True)
	comments = acmg_column2dict(rules_comments, False)

	exists_alternative = True
	if assignments['alternative_inheritance'] == 'NA':
		exists_alternative = False

	### ### ###
	# add supporting version of
	# PM2 and BS1 for all available inheritance modes
	### ### ###
	assignments = add_inheritance_supporting_rules(assignments, True, exists_alternative)
	comments = add_inheritance_supporting_rules(comments, False, exists_alternative)

	### ### ###
	# merge the created dictionaries
	### ### ###
	all_keys = list(assignments.keys()) + list(comments.keys())
	acmg_merged = {}
	for key in all_keys:
		if key in assignments and key in comments:
			acmg_merged[key] = {'assignment': assignments[key], 'comment': comments[key]}
		elif key in assignments:
			acmg_merged[key] = assignments[key]
		else:
			acmg_merged[key] = comments[key]
	return acmg_merged


def add_inheritance_supporting_rules(acmg_dict, exists_alternative, is_rules_column):
	"""
	Add supporting rules (PM2_Supporting, BS1_Supporting) for all inheritance modes

	Parameters
	----------
	acmg_dict: dict of str : str
        map each rule name to examination result
	exists_alternative: bool
		exists alternative inheritance (True), otherwise False
	is_rules_column: bool
		input column is rules assignment (True), otherwise False

	Returns
	-------
	dict of str : str
        updated rules with PM2 and BS1 with supporting strength
	"""
	### ### ### ### ### ### ###
	# if preferred inheritance does not contain supporting rules
	# then add them as not applicable
	### ### ### ### ### ### ###
	if "PM2_Supporting_preferred" not in acmg_dict:
		if is_rules_column:
			acmg_dict["PM2_Supporting_preferred"] = 'NA'
		else:
			acmg_dict["PM2_Supporting_preferred"] = 'Not applicable'
	if "BS1_Supporting_preferred" not in acmg_dict:
		if is_rules_column:
			acmg_dict["BS1_Supporting_preferred"] = 'NA'
		else:
			acmg_dict["BS1_Supporting_preferred"] = 'Not applicable'

	### ### ### ### ### ### ###
	# if alternative exists   #
	# add supporting rules    #
	# as not applicable       #
	### ### ### ### ### ### ###
	if exists_alternative:
		if "PM2_Supporting_alternative" not in acmg_dict:
			if is_rules_column:
				acmg_dict["PM2_Supporting_alternative"] = 'NA'
			else:
				acmg_dict["PM2_Supporting_alternative"] = 'Not applicable'
		if "BS1_Supporting_alternative" not in acmg_dict:
			if is_rules_column:
				acmg_dict["BS1_Supporting_alternative"] = 'NA'
			else:
				acmg_dict["BS1_Supporting_alternative"] = 'Not applicable'
	return acmg_dict


def format_coding_splicing_column(coding_column):
	"""
	Format coding and splicing column

	Parameters
	----------
	coding_column : str
	input coding and splicing column

	Returns
	-------
	str
		formatted coding and splicing column
	"""
	return coding_column.replace(",", "\n")

def format_extra_info_columns(classified_variant_info):
	"""
	Format columns with extra info: population, effect and conservation

	Parameters
	----------
	classified_variant_info : dict of str : str
		classified variant row in dictionary format

	Returns
	-------
	dict of str : str
		classified variant row with formated the extra info column
	"""
	for extra_info in ["dbSNP","gnomAD","gnomAD_sub","CADD","REVEL","MaxEntScan","dbscSNV","PhyloP"]:
		if extra_info in classified_variant_info:
			if str(classified_variant_info[extra_info]) == 'nan':
				classified_variant_info[extra_info] = "Not available"
	return classified_variant_info

def format_clinvar_column(clinvar_column):
	"""
	Format ClinVar column of the input GSvar file

	Parameters
	----------
	clinvar_column : str
	clinvar column

	Returns
	-------
	str
		formatted ClinVar column
	"""

	if str(clinvar_column) != '':
		clinvar_column = str(clinvar_column)
		if ';' in clinvar_column:
			if clinvar_column.count(';') > 1:
				return ' , '.join(clinvar_column.split(';'))
			else:
				return clinvar_column.strip(';')
		else:
			return clinvar_column
	else:
		return 'Not available'

def format_acmg_warnings(acmg_warning_column):
	"""
	Format ACMG warning flags

	Parameters
	----------
	acmg_warning_column : str
		ACMG warning flag column

	Returns
	-------
	str
		formatted warning flags
	"""

	if str(acmg_warning_column) == '-':
		return 'NA'
	else:
		return str(acmg_warning_column).replace("Warning: ","").replace("||","\n")

def split_name_assignment(rule):
	"""
	Split rule to name and assignment

	Parameters
	----------
	rule : str
		assigned ACMG rule

	Returns
	-------
	str
		rule name
	str
		rule assignment
	"""

	rule_name, rule_assignment = rule.split(": ", 1)

	if "|" in rule_assignment:
		rule_assignment = rule_assignment.split("|")[0]
	else:
		rule_assignment = rule_assignment
	if "PVS1" in rule_assignment or "PM5" in rule_assignment:
		rule_name = rule_assignment

	return rule_name, rule_assignment


def argparse_float_range(min_value, max_value, deactivating_value):
	"""
	Function to check floating value of command line argument

	Parameters
	----------
	min_value : float
		minimum allowed value
	max_value : float
		maximum allowed value
	deactivating_value : int
		value to deactivate floating argument

	Returns
	-------
	Function
		float_range_checker() function
	"""

	def float_range_checker(arg):
		try:
			float_arg = float(arg)
		except ValueError:
			raise argparse.ArgumentTypeError("Parameter value must be floating number")
		if float_arg != deactivating_value and (float_arg < min_value or float_arg > max_value):
			raise argparse.ArgumentTypeError(
				"Parameter value must be in the range of [{},{}] or equal to {}, to deactivate option".format(min_value,
				                                                                                              max_value,
				                                                                                              deactivating_value))
		return float_arg

	return float_range_checker
