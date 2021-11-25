from genotoscope.utils import extract_all_variant_types, arrange_transcripts_per_variant_type, aggregate_examined_rules, \
	parse_coding_splicing_column, extract_codons, get_codon_index, is_var2exclude
from genotoscope.variants.VariantInfo import VariantInfo

from os.path import join
import logging
from operator import xor

from pybedtools import BedTool
import vcf
from Bio.Seq import Seq
from pyensembl import EnsemblRelease
import hgvs.parser


class RefineLossOfFunction:
	"""
	Class to examine the PVS1 (very strong rule) from ACMG guideline and refine its impact significance
	based on the paper:
	Abou Tayoun, Ahmad N., et al. "Recommendations for interpreting the
	loss of function PVS1 ACMG/AMP variant criterion." Human mutation 39.11 (2018): 1517-1524.
	DOI: 10.1002/humu.23626

	To assess NMD use rules described at subsection "Prediction algorithm of NMD-elicit mutations" in Methods of the paper:
	Hu, Zhiyuan, Christopher Yau, and Ahmed Ashour Ahmed. "A pan-cancer genome-wide analysis reveals
	tumour dependencies by induction of nonsense-mediated decay." Nature communications 8.1 (2017): 1-9.
	DOI: 10.1038/ncomms15943.
	"""

	def __init__(self, data_path, clinvar_root, clinvar_file, clinvar_stars_df,
	             min_review_stars, beds_root, critical_prot_regions_file,
	             clinical_significant_exons_file, pheno_high_af_variants, phenotype_exons_file, uniprot_domains_file,
	             hugo_genes_df):
		"""
		RefineLossOfFunction constructor

		Parameters
		----------
		data_path : str
			data path
		clinvar_root : str
			root to clinvar folder
		clinvar_file : str
			ncbi clinvar filename
		clinvar_stars_df : pandas.DataFrame
			loaded clinvar review stars dataframe
		min_review_stars : int
			minimum review stars needed to filter in ClinVar variant
		beds_root : str
			annotation bed file root folder
		critical_prot_regions_file : str
			protein critical regions annotation filename
		clinical_significant_exons_file : str
			clinical significant exons annotation filename
		pheno_high_af_variants: vcf.Reader
			phenotype pathogenic variants with high AF
		phenotype_exons_file : str
			phenotype-relevant exons filename
		uniprot_domains_file : str
			UniProt domains file name
		hugo_genes_df : pandas.DataFrame
			loaded hugo genes dataframe

		Returns
		-------
		None
		"""
		### logger ###
		self.logger = logging.getLogger("GenOtoScope_Classify.PVS1")
		self.logger.info("Initialize class to refine ACMG/AMP PVS1 rule")
		self.data_path = data_path

		### PyEnsembl ###
		# release 75 uses human reference genome GRCh37
		self.ensembl_data = EnsemblRelease(75)

		### Annotation track files: Proteins critical regions & clinical significant exons ###
		self.beds_root = beds_root
		self.critical_prot_regions = BedTool(critical_prot_regions_file).sort()
		self.clinical_significant_exons = BedTool(clinical_significant_exons_file).sort()
		self.pheno_high_af_variants = pheno_high_af_variants
		self.phenotype_exons = BedTool(phenotype_exons_file).sort()
		self.uniprot_domains = BedTool(uniprot_domains_file).sort()
		self.matched_domains_cache = []

		### ClinVar ###
		clinvar_file = join(clinvar_root, clinvar_file)
		self.load_clinvar_file(clinvar_file)
		self.clinvar_stars_df = clinvar_stars_df
		# quality star system
		self.star_status2int = {"none": 0, "one": 1, "two": 2, "three": 3, "four": 4, "unknown review status": -1}
		self.min_review_stars = min_review_stars

		### Genes to strand information ###
		self.hugo_genes_df = hugo_genes_df

		### Variant types in variant_type column ###
		self.nonsense_frameshift_types = ["stop_gained", "stop_lost", "frameshift", "inframe_insertion",
		                                  "inframe_deletion"]
		self.splice_site_types = ["splice_acceptor"]
		self.start_types = ["start_lost"]

		### Variant types set in coding_and_splicing column ###
		self.intron_variant_types = {"splice_acceptor_variant", "splice_acceptor", "splice_donor_variant",
		                             "splice_donor"}
		self.exon_variant_types = {"stop_gained", "stop_lost", "frameshift", "frameshift_variant", "inframe_insertion",
		                           "inframe_deletion",
		                           "disruptive_inframe_insertion", "disruptive_inframe_deletion",
		                           "conservative_inframe_insertion", "conservative_inframe_deletion"}
		self.hgvs_parser = hgvs.parser.Parser()

	def is_transcript_type_splice_acceptor_donor(self, transcript_types):
		"""
		Examine if transcript type is splice acceptor or splice donor

		Parameters
		----------
		transcript_types : set of str
			variant types for affected transcript

		Returns
		-------
		bool
			transcript type is splice acceptor or splice donor (True), otherwise False
		"""
		if len(transcript_types.intersection(self.intron_variant_types)) > 0 and len(
				transcript_types.intersection(self.exon_variant_types)) == 0:
			return True
		else:
			return False

	def load_clinvar_file(self, clinvar_file):
		"""
		Load clinvar file to find annotated variants
		to decide for ACMG/AMP rules
		kudos: https://github.com/jamescasbon/PyVCF/issues/201

		Parameters
		----------
		clinvar_file : str
			clinvar file name

		Returns
		-------
		None
		"""
		# print("Loading ClinVar file from {}".format(join(self.data_path, clinvar_file)))
		self.vcf_reader = vcf.Reader(filename=join(self.data_path, clinvar_file), compressed=True, encoding='UTF-8')

	def filter_transcripts_by_info(self, transcripts_info, info_name):
		"""
		Filter transcripts per info

		Parameters
		----------
		list of dict of str : str
			list of transcripts and their info to be filtered
		info_name : str
			information to filter transcripts for

		Returns
		-------
		list of dict of str: str
			filtered transcripts containing information
		"""
		self.logger.debug("Filter transcripts: {}".format(transcripts_info))
		filtered_transcripts_info = []
		for transcript_info in transcripts_info:
			transcript = self.ensembl_data.transcript_by_id(transcript_info["transcript_id"])
			if info_name == "start_codon":
				if transcript.contains_start_codon:
					filtered_transcripts_info.append(transcript_info)
			elif info_name == "exon_info":
				if transcript_info['exon'] != '':
					filtered_transcripts_info.append(transcript_info)
			elif info_name == "start_stop_codons":
				if transcript.contains_start_codon and transcript.contains_stop_codon:
					filtered_transcripts_info.append(transcript_info)
			else:
				self.logger.debug("Remove transcript with id: {}".format(transcript_info["transcript_id"]))
		return filtered_transcripts_info

	@staticmethod
	def truncated_exons2bed_interval(truncated_exons_pos, variant_info, transcript_strand):
		"""
		Convert truncated exons 2 bed interval object

		Parameters
		----------
		truncated_exons_pos : list of dict of str : str or int
			truncated exons information
		variant_info : VariantInfo
			variant information
		transcript_strand : str
			transcript strand

		Returns
		-------
		BedTool.Interval
			bedtool interval object of truncated exons
		"""
		exon_lines = []
		for truncated_exon in truncated_exons_pos:
			exon_lines.append(' '.join([variant_info.chrom, str(truncated_exon["exon_start"]), str(
				truncated_exon["exon_end"]), truncated_exon["exon_id"], '.', transcript_strand]))
		return BedTool('\n'.join(exon_lines), from_string=True)[0]

	def intersect_exon2phenotype_relevants(self, truncated_exons_pos, variant_info, transcript_strand):
		"""
		Intersect exon with phenotype relevant transcripts

		Parameters
		----------
		truncated_exons_pos : list of dict of str : str or int
			truncated exons (by variant) information
		variant_info : VariantInfo
			variant information
		transcript_strand : str
			transcript strand

		Returns
		-------
		bool
			truncated exon id is found in phenotype relevant transcript (True), otherwise (False)
		list of str
			matched exon Ensembl ids
		"""
		self.logger.debug("Examine if exon is contained in phenotype relevant transcripts")
		self.logger.debug("Affected region: {}".format(
			RefineLossOfFunction.truncated_exons2bed_interval(truncated_exons_pos, variant_info,
			                                                  transcript_strand)))
		matches = self.phenotype_exons.all_hits(
			RefineLossOfFunction.truncated_exons2bed_interval(truncated_exons_pos, variant_info,
			                                                  transcript_strand), same_strand=True)
		if len(matches) > 0:
			return True, [match.name for match in matches]
		else:
			return False, []

	def run(self, sample_name, chrom, var_start, var_end, gene_name, var_ref, var_obs, variant_type,
	        transcript_splicing_info):
		"""
		Follow figure 1 of Tayoun et al. paper and run refinement rules per variant type

		Parameters
		----------
		sample_name : str
			sample name that contains variant
		chrom : str
			chromosome harbouring variant
		var_start : str
			variant start position
		var_end : str
			variant end position
		gene_name : str
			gene name column in GSvar file
		var_ref : str
			reference bases in variant position
		var_obs : str
			observed bases in variant position
		variant_type : str
			variant type column in GSvar file
		transcript_splicing_info : str
			coding_and_splicing column in GSVar file

		Returns
		-------
		str
			assigned class based on Tayoun et al.
		str
			assignment comment
		dict of str : str
			constructed variant coding sequence per transcript
		"""
		var_chrom = chrom
		var_genomic_start = int(var_start)
		var_genomic_end = int(var_end)
		var_ref = var_ref
		var_obs = var_obs
		var_coding_seqs = None
		# read all variant types
		variant_types = extract_all_variant_types(variant_type)
		transcripts_info = parse_coding_splicing_column(transcript_splicing_info, gene_name, self.hgvs_parser)
		# create variant info object
		variant_info = VariantInfo(sample_name, gene_name, variant_type, transcripts_info, var_chrom, var_genomic_start,
		                           var_genomic_end, var_ref, var_obs)
		self.logger.debug("variant types: {}".format(variant_types))
		self.logger.debug("transcripts info: {}".format(transcripts_info))
		# filter transcripts to the ones containing start codon information
		filtered_transcripts_info = self.filter_transcripts_by_info(transcripts_info, "start_stop_codons")
		self.logger.debug("filtered transcripts: {}".format(filtered_transcripts_info))
		assigned_pvs1_per_var_type = []
		assignment_comment_per_var_type = []
		if len(filtered_transcripts_info) == 0:
			# self.logger.debug("Only one class")
			assigned_pvs1_per_var_type.append("False")
			assignment_comment_per_var_type.append("variant contains zero filtered transcripts")
		else:
			# arrange affected transcripts per variant type
			# then run for all transcripts of the same variant the PVS1 variant type specific refinement
			pvs1_variant_types = set(
				["start_lost", "stop_gained", "stop_lost", "frameshift", "inframe_insertion",
				 "inframe_deletion", "conservative_inframe_insertion", "conservative_inframe_deletion",
				 "disruptive_inframe_deletion", "disruptive_inframe_insertion", "splice_acceptor", "splice_donor"])
			transcripts_per_var_type = arrange_transcripts_per_variant_type(variant_types, filtered_transcripts_info,
			                                                                pvs1_variant_types)
			self.logger.debug("transcripts per var type: {}".format(transcripts_per_var_type))
			assigned_pvs1_per_var_type = []
			assignment_comment_per_var_type = []
			for var_type, transcripts_info in transcripts_per_var_type.items():
				current_var_type_classes, current_var_type_comments = None, None
				self.logger.debug("var_type: {}".format(var_type))
				self.logger.debug("transcripts: {}".format(transcripts_info))
				if var_type == "start_lost":
					self.logger.debug("{} to run".format(var_type))
					current_var_type_classes, current_var_type_comments, var_coding_seqs = self.refine_start_lost(
						transcripts_info,
						variant_info)
				elif var_type == "splice_acceptor" or var_type == "splice_donor":
					self.logger.debug("{} to run".format(var_type))
					current_var_type_classes, current_var_type_comments, var_coding_seqs = self.refine_splice_site(
						transcripts_info,
						variant_info)
				elif var_type == "stop_gained" or var_type == "stop_lost" or var_type == "frameshift" or "inframe" in var_type:
					self.logger.debug("{} to run".format(var_type))
					transcripts_containing_exon_info = self.filter_transcripts_by_info(transcripts_info, "exon_info")
					current_var_type_classes, current_var_type_comments, var_coding_seqs = self.refine_nonsense_frameshift(
						transcripts_containing_exon_info, variant_info)
				else:
					# PVS1 refinement is not applicable for current variant type
					for transcript_info in transcripts_info:
						assigned_pvs1_per_var_type.append("NA")
						assignment_comment_per_var_type.append(
							transcript_info[
								"transcript_id"] + ": " + "PVS1 refinement not applicable for variant type:{}".format(
								var_type))
				# self.logger.debug(
				# 	"For variant type: {}, PVS1 refinement is not applicable\n => var pos: {}".format(var_type,
				# 	                                                                                  variant_info.to_string()))

				if current_var_type_classes and current_var_type_comments:
					# current variant type classes and comments are set, then append them
					assigned_pvs1_per_var_type += current_var_type_classes
					assignment_comment_per_var_type += current_var_type_comments
				else:
					# PVS1 refinement is not applicable for current variant type
					for transcript_info in transcripts_info:
						assigned_pvs1_per_var_type.append("NA")
						assignment_comment_per_var_type.append(
							transcript_info[
								"transcript_id"] + ": " + "PVS1 refinement not applicable for variant info")
			# if after refinement there no assignment for any transcript
			# mark variant as not applicable for PVS1
			if len(assigned_pvs1_per_var_type) == 0:
				# self.logger.debug("Current variant is not applicable for PVS1")
				assigned_pvs1_per_var_type.append("NA")
				assignment_comment_per_var_type.append("PVS1 refinement not applicable for variant info")

		aggregated_pvs1, aggregated_pvs1_comment = aggregate_examined_rules(assigned_pvs1_per_var_type,
		                                                                    assignment_comment_per_var_type, "PVS1")
		return aggregated_pvs1, aggregated_pvs1_comment, var_coding_seqs

	@staticmethod
	def is_transcript_in_positive_strand(transcript):
		"""
		Check if transcript is on positive strand

		Parameters
		----------
		transcript: pyensembl.transcript
			ensemble transcript id

		Returns
		-------
		bool
			True if transcript is in positive strand, otherwise False
		"""
		strand = transcript.exons[0].to_dict()["strand"]
		if strand == "+":
			return True
		else:
			return False

	def get_transcript_exon_offsets(self, transcript, is_genomic):
		"""
		Get all exons cDNA or genomic offsets of transcript

		Parameters
		----------
		transcript : pyensembl.transcript
			ensemble transcript id
		is_genomic : bool
			create genomic position for coding exons of transcript (True), otherwise return cDNA position

		Returns
		-------
		list of list of int
			all exons start and end position
		"""
		if is_genomic:
			### ### ###
			# follow strand direction for chromosomal positions:
			# exon_i start < exon_i end for positive strand
			# exon_i start > exon_i end for negative strand
			### ### ###
			coding_exons_chrom_ranges = []
			for exon in transcript.exons:
				if RefineLossOfFunction.is_transcript_in_positive_strand(transcript):
					coding_exons_chrom_ranges.append([exon.to_dict()["start"], exon.to_dict()["end"]])
				else:
					coding_exons_chrom_ranges.append([exon.to_dict()["end"], exon.to_dict()["start"]])
			return coding_exons_chrom_ranges
		else:
			# create cDNA exon positions
			# offsets starting from 1
			exon_positions = []
			is_first_exon = True
			for exon in transcript.exons:
				if is_first_exon:
					start = 1
					is_first_exon = False
				else:
					start = end + 1
				end = start + exon.to_dict()["end"] - exon.to_dict()["start"]
				exon_positions.append([start, end])
				try:
					assert end - start == exon.to_dict()["end"] - exon.to_dict()["start"]
				except AssertionError:
					self.logger.error("Exon genomic regions and coding regions should be equal in size", exc_info=True)
			# sum = 0
			# for exon_position in exon_positions:
			# 	sum += exon_position[1] - exon_position[0] + 1
			self.logger.debug("Exon positions: {}".format(exon_positions))
			return exon_positions

	def parse_variant_intron_pos(self, var_coding):
		"""
		Parse variant offset in intron

		Parameters
		----------
		var_coding : hgvs.PosEdit
			variant in coding co-ordinates

		Returns
		-------
		list of str
			for each variant edit part (start and end):
			'+' to show that variant is just after an exon, '-' to show that the variant is before an exon
		list of int
			for each variant edit part:
			variant offset in intron
		list of int
			for each variant edit part:
			+1 if variant needs to increase its position to get to the closest exon (split_symbol = -)
			-1 if variant needs to decrease its position to get to the closest exon (split_symbol = +)
		"""
		self.logger.debug("Parse intron variant position")
		var_coding_str = str(var_coding)
		var_edit = str(var_coding.edit)
		intron_offset, direction2closest_exon = 0, 0
		split_symbols, intron_offsets, directions2exon = [], [], []
		# find the direction of the closest exon
		if '_' in var_coding_str:
			# for duplication, insertion and deletion
			# split on '_' character before finding direction
			[edit_start, edit_end] = var_coding_str.split(var_edit)[0].split('_')
			for edit_part in var_coding_str.split(var_edit)[0].split('_'):
				if "+" in edit_part:  # after exon in start of the intron
					split_symbol = "+"
					direction2closest_exon = -1
				elif "-" in edit_part:  # before exon in the end of the intron
					split_symbol = "-"
					direction2closest_exon = +1
				else:
					# print("var edit part: {} does not contain split symbol".format(edit_part))
					continue
				# get position of edit part inside the exon

				intron_offset_pos = int(edit_part.split(split_symbol)[1])
				split_symbols.append(split_symbol)
				directions2exon.append(direction2closest_exon)
				intron_offsets.append(intron_offset_pos)
		else:
			# for SNP find direction symbol
			if "+" in var_coding_str:  # after exon in start of the intron
				split_symbol = "+"
				direction2closest_exon = -1
			else:  # before exon in the end of the intron
				split_symbol = "-"
				direction2closest_exon = +1
			split_symbols.append(split_symbol)
			directions2exon.append(direction2closest_exon)
			# to get the splice site position, after splitting on - or +,
			# split on the edit to get the splice offset
			intron_offset_pos = int(var_coding_str.split(split_symbol)[-1].split(var_edit)[0])
			intron_offsets.append(intron_offset_pos)
		self.logger.debug(
			"split_symbols: {}, intron_offsets: {}, directions2exon: {}".format(split_symbols, intron_offsets,
			                                                                    directions2exon))
		try:
			assert len(split_symbols) == len(intron_offsets) == len(directions2exon)
		except AssertionError:
			self.logger.error(
				"All parts of variant edit should contain split symbol, intron offset and direction to closest exon.",
				exc_info=True)
		return split_symbols, intron_offsets, directions2exon

	def find_exon_by_ref_pos(self, transcript, ref_pos, is_genomic):
		"""
		Find exon info that contains the reference position
		the returned exon index is 0-based index

		Parameters
		----------
		transcript: pyensembl.transcript
			transcript object to extract exons from
		ref_pos: int
			reference position (1-index) to find contained exon
		is_genomic: bool
			is the position genomic (True), otherwise the position is cDNA position (False)

		Returns
		-------
		int
			exon index (0-based) containing the reference position
		int
			offset of reference position in overlapping exon
		"""
		self.logger.debug("Find exon containing reference position: {}".format(ref_pos))
		if RefineLossOfFunction.is_transcript_in_positive_strand(
				transcript):
			strand_direction = + 1
		else:
			strand_direction = - 1
		exon_positions = self.get_transcript_exon_offsets(transcript, is_genomic)
		overlap_exon_index, pos_in_overlap_exon = -1, -1

		# reference position is found using pyensembl coding sequence,
		# this coding sequence does not contain the untranslated region of the starting exons
		# thus you need to add the length of the sequence of the first exons up to start codon (ATG)
		if not is_genomic:
			# for cDNA position don't need to reverse exon start and end
			strand_direction = + 1
			len_first_exons_up_start_codon = transcript.start_codon_spliced_offsets[0]
			ref_pos = ref_pos + len_first_exons_up_start_codon

		# find reference position into exon intervals
		for exon_idx, exon_interval in enumerate(exon_positions):
			normalized_exon_interval = range(exon_interval[0], exon_interval[1] + strand_direction, strand_direction)
			if ref_pos in normalized_exon_interval:
				overlap_exon_index = exon_idx
				pos_in_overlap_exon = abs(ref_pos - exon_interval[0])

		try:
			assert overlap_exon_index >= 0 and pos_in_overlap_exon >= 0
		except AssertionError:
			self.logger.error(
				"Exon containing the reference position could not be found\n=> reference pos: {}".format(ref_pos),
				exc_info=True)
		return overlap_exon_index, pos_in_overlap_exon

	def find_exon_by_var_pos(self, transcript, transcript_info, variant_info, is_genomic):
		"""
		Find variant exon index by variant coding position, exon index is 1-based

		Parameters
		----------
		transcript : pyensembl.transcript
			transcript object contains the penultimate and final exon
		transcript_info : dict of str: str
			transcript information
		variant_info : VariantInfo
			variant basic info
		is_genomic : bool
			genomic positions will be used (True), otherwise cDNA positions will be used (False)

		Returns
		-------
		list of int
			indices of exons overlapping with variant (1-based)
		int
			variant start position as offset in first overlapping exon
		int
			variant end	position as offset in last overlapping exon
		"""
		self.logger.debug("Find exon indices containing the variant")
		if not ("exon" in transcript_info["exon"] or "intron" in transcript_info["exon"]):
			self.logger.debug("Exon information is not found in VEP => use VEP spliced offset")
			is_genomic = False

		exon_positions = self.get_transcript_exon_offsets(transcript, is_genomic)

		overlap_exon_indices = []
		var_start_offset, var_end_offset = 0, 0  # save the variant start and end offset in the overlapping exon regions
		var_coding = transcript_info["var_coding"]
		### ### ###
		# get the strand direction
		### ### ###
		if RefineLossOfFunction.is_transcript_in_positive_strand(
				transcript):
			strand_direction = + 1
		else:
			strand_direction = - 1
		if is_genomic:
			if not self.is_transcript_type_splice_acceptor_donor(transcript_info["type_variant"]):
				### ### ###
				# Exonic variant
				### ### ###
				self.logger.debug("Exonic variant type")
				var_start = variant_info.genomic_start
				var_end = variant_info.genomic_end
			else:
				self.logger.debug("Intron variant type => update as variant pos the starting position of skipped exon")
				split_symbols, intron_offsets, directions2exon = self.parse_variant_intron_pos(var_coding)

				### ### ###
				# Intronic variant
				# Parse affected exon index by VEP (1-based)
				# add splice site direction to affected exon index
				# (currently modeled that) intron variant can disrupt at most one exon, so equal variant end with its start position
				### ### ###

				if "exon" in transcript_info["exon"] and "intron" in transcript_info["exon"]:
					exon_idx = int(transcript_info["exon"].split("intron")[1].split("/")[0]) - 1
					self.logger.debug("Transcript info contains both exon and intron offset => use intron offset")
				elif "exon" in transcript_info["exon"]:
					exon_idx = int(transcript_info["exon"].split("/")[0].split("exon")[1]) - 1
				elif "intron" in transcript_info["exon"]:
					exon_idx = int(transcript_info["exon"].split("/")[0].split("intron")[1]) - 1
				else:
					try:
						assert "exon" in transcript_info["exon"] or "intron" in transcript_info["exon"]
					except AssertionError:
						self.logger.error(
							"Variant does not contain exon or intron index in VEP column\n => variant position: {}".format(
								variant_info.to_string()), exc_info=True)
				if exon_idx + directions2exon[0] >= 0:
					skipped_exon = exon_positions[exon_idx + directions2exon[0]]
				else:
					self.logger.debug("Skipping first exon")
					skipped_exon = exon_positions[0]
				var_start = skipped_exon[0]
				var_end = var_start
				self.logger.debug("Updated variant start:{}, end: {} on exon idx: {}".format(var_start, var_end,
				                                                                             exon_idx + directions2exon[
					                                                                             0]))
		else:
			# use cDNA offset for both frameshift and nonsense mutation
			var_start = int(var_coding.pos.start.base)
			var_end = int(var_coding.pos.end.base)
			if var_start < 0:
				var_start = 0
			if var_end < 0:
				if transcript_info["coding_diff_len"] > 0:
					var_end = var_start + transcript_info["coding_diff_len"]
				else:
					# deletion case
					var_end = var_start

			# VEP positions include only the coding sequence of the transcript,
			# so you need to add the length of the sequence of the first exons up to start codon (ATG)
			len_first_exons_up_start_codon = transcript.start_codon_spliced_offsets[0]
			var_start = var_start + len_first_exons_up_start_codon
			var_end = var_end + len_first_exons_up_start_codon

		### ### ###
		# find variant position into exon intervals
		### ### ###
		for exon_idx, exon_interval in enumerate(exon_positions):
			if is_genomic:
				normalized_exon_interval = range(exon_interval[0], exon_interval[1] + 2 * strand_direction,
				                                 strand_direction)
			else:
				normalized_exon_interval = range(exon_interval[0], exon_interval[1] + 1)
			self.logger.debug("Exon interval: {}".format(normalized_exon_interval))
			if var_start in normalized_exon_interval:
				overlap_exon_indices.append(exon_idx + 1)
				break
		self.logger.debug("Search var_start: {} found in exon(s): {}".format(var_start, overlap_exon_indices))

		### ### ###
		# find variant start and end offset in exons
		### ### ###
		if len(overlap_exon_indices) == 1:  # variant included in only one exon
			if self.is_transcript_type_splice_acceptor_donor(transcript_info["type_variant"]):
				# exon-skipping, get as offset the total range of skipped exon
				var_start_offset = 0
				# negative strand transcripts contain higher start position than end position
				var_end_offset = abs(exon_positions[overlap_exon_indices[0] - 1][1] - \
				                     exon_positions[overlap_exon_indices[0] - 1][0])
			else:
				if is_genomic:
					if strand_direction == 1:
						# for positive strand, compute distance from exon lower position (start)
						var_start_offset = var_start - exon_positions[overlap_exon_indices[0] - 1][0]
						var_end_offset = var_end - exon_positions[overlap_exon_indices[0] - 1][0]
					else:
						# for negative strand, compute distance from exon higher position (end)
						var_start_offset = exon_positions[overlap_exon_indices[0] - 1][0] - var_end
						var_end_offset = exon_positions[overlap_exon_indices[0] - 1][0] - var_start
				else:
					var_start_offset = var_start - exon_positions[overlap_exon_indices[0] - 1][0]
					var_end_offset = var_end - exon_positions[overlap_exon_indices[0] - 1][0]
		else:
			try:
				assert len(overlap_exon_indices) > 0
			except AssertionError:
				self.logger.error("Overlapping exon indices should be more than 0\n=> variant position: {}".format(
					variant_info.to_string()), exc_info=True)
			### ### ###
			# variant included in more than one exons,
			# so start is in the first exon, end is in the last exon
			### ### ###
			if is_genomic:
				if strand_direction == 1:
					var_start_offset = var_start - exon_positions[overlap_exon_indices[0] - 1][0]
					var_end_offset = var_end - exon_positions[overlap_exon_indices[len(overlap_exon_indices)] - 1][0]
				else:
					# for negative strand, use the highest value (exon start) to compute the offset for the start and exon position of the variant
					var_start_offset = exon_positions[overlap_exon_indices[0] - 1][0] - var_end
					var_end_offset = exon_positions[overlap_exon_indices[len(overlap_exon_indices)] - 1][0] - var_start
			else:
				var_start_offset = var_start - exon_positions[overlap_exon_indices[0] - 1][0]
				var_end_offset = var_end - exon_positions[overlap_exon_indices[len(overlap_exon_indices)] - 1][0]
		try:
			assert var_start_offset >= 0 and var_end_offset >= 0
		except AssertionError:
			self.logger.error("Variant start and end offset should be higher than 0\n=> variant position: {}".format(
				variant_info.to_string()), exc_info=True)
		self.logger.debug(
			"Overlap exon indices: {}, var start offset: {}, var end offset: {}".format(overlap_exon_indices,
			                                                                            var_start_offset,
			                                                                            var_end_offset))
		return overlap_exon_indices, var_start_offset, var_end_offset

	def assess_exon_skipping(self, transcript, transcript_info, intron_offsets, variant_info):
		"""
		Assess if exon will be skipped
		By examining if the affected intron position in on the splice sites: +/- 1,2

		Parameters
		----------
		transcript : pyensembl.transcript
			ensembl object for variant-affected transcript
		transcript_info : dict of str: str
			variant-affected transcript GSvar information
		intron_offsets : list of int
			variant offset positions in intron
		variant_info : VariantInfo
			variant basic info

		Returns
		-------
		list of int
			list of exon (1-based) indices containing the variant
		bool
			True if exon predicted to be skipped, otherwise False
		int
			if exon is skipped, variant deletion start position
		int
			if exon is skipped, variant deletion end position
		bool
			variant skips start codon exon
		bool
			variant skips stop codon exon
		bool
			variant skips coding exon
		"""
		self.logger.debug("Assess exon skipping for splice acceptor variant")
		is_exon_skipped = False
		variant_skips_start_codon_exon, variant_skips_stop_codon_exon = False, False
		variant_skips_coding_exon = False
		skipped_exon_start, skipped_exon_end = 0, 0
		self.logger.debug("intron offsets: {}".format(intron_offsets))

		for intron_offset in intron_offsets:
			if intron_offset in [1, 2]:
				# variant is disrupting donor/acceptor splicesome site
				# predict that exon will be skipped
				is_exon_skipped = True
		self.logger.debug("is_exon_skipped: {}".format(is_exon_skipped))

		if is_exon_skipped:
			self.logger.debug("exon is skipped, call find_exon_by_var_pos()")
			# if exon is skipped find its start and end
			skipped_exons, var_exon_start, var_exon_end = self.find_exon_by_var_pos(transcript,
			                                                                        transcript_info,
			                                                                        variant_info,
			                                                                        is_genomic=True)
			### ### ### ###
			# Examine if skipped exon is coding
			### ### ### ###
			# find exons containing start and stop codons
			if RefineLossOfFunction.is_transcript_in_positive_strand(transcript):
				transcript_strand = "+"
				start_codon_first_pos = transcript.start_codon_positions[0]
			else:
				transcript_strand = "-"
				start_codon_first_pos = transcript.start_codon_positions[-1]
			start_codon_exon_idx, start_codon_exon_offset = self.find_exon_by_ref_pos(transcript, start_codon_first_pos,
			                                                                          True)
			stop_codon_exon_idx, stop_codon_exon_offset = self.find_exon_by_ref_pos(transcript,
			                                                                        len(transcript.coding_sequence),
			                                                                        False)

			if skipped_exons[0] >= start_codon_exon_idx + 1:
				### ### ### ###
				# skipped exon is coding
				# find start and end position of skipped coding exons
				### ### ### ###
				variant_skips_coding_exon = True
				# get as start and end the skipped exon positions
				exon_offsets = self.get_transcript_exon_offsets(transcript, False)
				try:
					assert len(skipped_exons) == 1
				except AssertionError:
					self.logger.error(
						"Currently intron variant case is implemented only affecting one exon\n=> variant position: {}".format(
							variant_info.to_string()), exc_info=True)

				# exons offsets contain the length of the sequence of the first exons up to start codon (ATG)
				# so subtract this length to interact with pyensembl transcipt.coding_sequence
				len_first_exons_up_start_codon = transcript.start_codon_spliced_offsets[0]
				self.logger.debug("len of first exons up to start codon: {}".format(len_first_exons_up_start_codon))
				skipped_exon_offsets = exon_offsets[skipped_exons[0] - 1]
				skipped_exon_start = skipped_exon_offsets[0] - len_first_exons_up_start_codon

				if skipped_exon_start <= 1:
					# if exon start is before the start codon position,
					# make the variant start position equal to 1
					skipped_exon_start = 1
					variant_skips_start_codon_exon = True
				skipped_exon_end = skipped_exon_offsets[1] - len_first_exons_up_start_codon
				self.logger.debug("skipped exon start: {}, end:{}".format(skipped_exon_start, skipped_exon_end))

				### ### ###
				# examine if the last exon that was skipped,
				# contained the stop codon
				### ### ###
				if is_exon_skipped and skipped_exons[0] == stop_codon_exon_idx + 1:
					self.logger.debug("Search for termination codon on the skipped exon")
					# create skipped coding sequence with in-frame start
					inframe_start = (skipped_exon_start - 1) % 3
					skipped_inframe_seq = transcript.coding_sequence[
					                      skipped_exon_start - 1 + inframe_start:skipped_exon_end]
					self.logger.debug("skipped inframe coding seq: {}".format(skipped_inframe_seq))
					if self.search_termination_codon(extract_codons(skipped_inframe_seq), False):
						# for stop codon exon skipping, normalize exon end to stop codon position
						variant_skips_stop_codon_exon = True
						skipped_exon_end = len(str(transcript.coding_sequence))
		else:
			self.logger.debug("Exon is not skipped")
		return skipped_exons, is_exon_skipped, skipped_exon_start, skipped_exon_end, variant_skips_start_codon_exon, variant_skips_stop_codon_exon, variant_skips_coding_exon

	def is_genomic_pos_in_coding_exon(self, transcript, genomic_pos):
		"""
		Examine if genomic position is contained in coding exon

		Parameters
		----------
		transcript : pyensembl.transcript
			ensembl object containing the genomic position
		genomic_pos : int
			chromosomal position

		Returns
		-------
		bool
			position is contained in coding exon (True), otherwise False
		"""
		pos_in_coding_exon = False
		self.logger.debug("Examine if genomic pos: {} is contained in coding exon sequence".format(genomic_pos))
		if RefineLossOfFunction.is_transcript_in_positive_strand(transcript):
			strand_direction = +1
		else:
			strand_direction = -1
		num_coding_exons = len(transcript.coding_sequence_position_ranges)

		### ### ###
		# loop through values and check if variant overlap a coding exonic range
		### ### ###
		for coding_exon_idx, coding_interval in enumerate(transcript.coding_sequence_position_ranges):
			# set up exon coding start and end positions
			if strand_direction == 1:
				exon_coding_start = coding_interval[0]
				exon_coding_end = coding_interval[1]
			else:
				exon_coding_start = coding_interval[1]
				if coding_exon_idx + 1 == num_coding_exons:
					# update end of last exon to be the lowest chromosomal position of the stop codon
					exon_coding_end = transcript.stop_codon_positions[0]
				else:
					exon_coding_end = coding_interval[0]
			normalized_coding_interval = range(exon_coding_start, exon_coding_end + 2 * strand_direction,
			                                   strand_direction)

			self.logger.debug("normalized interval: {}".format(normalized_coding_interval))
			if genomic_pos in normalized_coding_interval:
				self.logger.debug("position in exon with offset: {}".format(coding_exon_idx))
				pos_in_coding_exon = True
				break
		return pos_in_coding_exon

	def convert_genomic2coding_pos(self, transcript, genomic_pos, variant_info):
		"""
		Convert genomic position into coding position

		Parameters
		----------
		transcript : pyensembl.transcript
			ensembl object for transcript containing genomic position
		genomic_pos : int
			genomic position
		variant_info : VariantInfo
			variant basic information

		Returns
		-------
		int
			coding sequence offset
		"""
		self.logger.debug("Convert variant chromosomal position to coding position")

		if RefineLossOfFunction.is_transcript_in_positive_strand(transcript):
			strand_direction = +1
			start_codon_first_pos = transcript.start_codon_positions[0]
		else:
			strand_direction = -1
			start_codon_first_pos = transcript.start_codon_positions[-1]

		sum_coding_len, pos_exon_offset = 0, -1
		self.logger.debug("start codon: {}".format(transcript.start_codon_positions))
		self.logger.debug("stop codon: {}".format(transcript.stop_codon_positions))
		num_coding_exons = len(transcript.coding_sequence_position_ranges)

		### ### ###
		# loop through coding sequence of exons
		# if variant chromosomal position is not included on current exon, add up of coding length
		# else compute the offset of the variant position from the current exon start (respecting strand direction)
		### ### ###
		for coding_exon_idx, coding_interval in enumerate(transcript.coding_sequence_position_ranges):
			# set up exon coding start and end positions
			if strand_direction == 1:
				exon_coding_start = coding_interval[0]
				exon_coding_end = coding_interval[1]
			else:
				exon_coding_start = coding_interval[1]
				if coding_exon_idx + 1 == num_coding_exons:
					# update end of last exon to be the lowest chromosomal position of the stop codon
					exon_coding_end = transcript.stop_codon_positions[0]
				else:
					exon_coding_end = coding_interval[0]
			normalized_coding_interval = range(exon_coding_start, exon_coding_end + strand_direction, strand_direction)
			self.logger.debug("normalized interval: {}".format(normalized_coding_interval))
			if genomic_pos in normalized_coding_interval:
				self.logger.debug("position in exon")
				# find offset of genomic position in containing exon
				if strand_direction == +1:
					pos_exon_offset = genomic_pos - exon_coding_start
				else:
					pos_exon_offset = exon_coding_start - genomic_pos
				self.logger.debug("exon offset: {}".format(pos_exon_offset))
				break
			else:
				# add coding length for each coding exon before exon containing genomic position
				sum_coding_len = sum_coding_len + abs(exon_coding_end - exon_coding_start) + 1
				self.logger.debug("sum of coding: {}".format(sum_coding_len))

		if pos_exon_offset == -1:
			self.logger.debug("Normalize position that is before the start codon")
			# check if variant contains non-translated positions before the start codon
			# update position to be the first coding position
			sum_coding_len = 0
			if strand_direction == +1:
				if genomic_pos < start_codon_first_pos:
					pos_exon_offset = 0
			else:
				if genomic_pos > start_codon_first_pos:
					pos_exon_offset = 0
		try:
			assert pos_exon_offset != -1
		except AssertionError:
			self.logger.error("Chromosome position should be on any coding exon\n => variant position: {}".format(
				variant_info.to_string()), exc_info=True)
		return sum_coding_len + pos_exon_offset

	def construct_variant_coding_seq(self, transcript, transcript_info, variant_info):
		"""
		Add variant to coding sequence to get observed coding sequence
		following hgvs recommendations: http://varnomen.hgvs.org/recommendations/general/

		Then return variant-integrated coding sequence in the forward strand
		or the reverse complement, based on the transcript directionality

		Parameters
		----------
		transcript : pyensembl.transcript
			ensembl object for variant-affected transcript
		transcript_info : dict of str: str
			variant-affected transcript GSvar information
		variant_info : VariantInfo
			variant basic info

		Returns
		-------
		str
			variant coding sequence
		int
			diff_len
		"""
		self.logger.debug(
			"Add variant to coding sequence to create observed coding sequence of transcript id: {}".format(
				transcript_info["transcript_id"]))
		var_coding = transcript_info["var_coding"]
		var_edit = str(var_coding.edit)

		is_exon_skipped = False
		diff_len = 0
		variant_skips_start_codon_exon, variant_skips_stop_codon_exon, variant_skips_coding_exon = False, False, False
		# find start end of variant positions
		if not self.is_transcript_type_splice_acceptor_donor(transcript_info["type_variant"]):
			self.logger.debug("Exonic variant")
			# subtract from variant cDNA position the length of the prime 5' end,
			# to get position in only coding parts of exons
			if int(var_coding.pos.start.base) >= 1:
				var_start = int(var_coding.pos.start.base)
			else:
				var_start = 1
			if int(var_coding.pos.end.base) >= 1:
				var_end = int(var_coding.pos.end.base)
			else:
				var_end = 1
		else:
			self.logger.debug("Intronic variant")
			# find intron variant position
			split_symbols, intron_offsets, directions2exon = self.parse_variant_intron_pos(var_coding)
			# find if variant results to exon skipping
			exons_containing_var, is_exon_skipped, var_start, var_end, variant_skips_start_codon_exon, variant_skips_stop_codon_exon, variant_skips_coding_exon = self.assess_exon_skipping(
				transcript,
				transcript_info,
				intron_offsets,
				variant_info)
			try:
				assert is_exon_skipped is True
			except AssertionError:
				self.logger.error(
					"Creating variant observed coding sequence, implemented for exon skipping case only\n=> variant position: {}".format(
						variant_info.to_string()), exc_info=True)
			self.logger.debug("Exon(s) in list: {} are skipped".format(exons_containing_var, is_exon_skipped))
			self.logger.debug("variant skips start codon: {}, stop codon: {}, skips coding exon: {}".format(
				variant_skips_start_codon_exon,
				variant_skips_stop_codon_exon, variant_skips_coding_exon))
			self.logger.debug("Skipped exon offset start: {}, end: {}".format(var_start, var_end))
		# load the coding sequence to introduce the variant
		coding_seq = str(transcript.coding_sequence)
		self.logger.debug("Coding seq: {}, len={}".format(coding_seq, len(coding_seq)))
		var_coding_seq = ""
		vep_coordinates_used, vep_contains_seq = True, False

		### ### ### ###
		# integrate variant into coding sequence, per variant type
		### ### ### ###
		if self.is_transcript_type_splice_acceptor_donor(transcript_info["type_variant"]):
			# to delete skipped exon sequence
			# use the calculated ensembl coding offsets
			vep_coordinates_used, vep_contains_seq = False, False
			if variant_skips_start_codon_exon and variant_skips_stop_codon_exon:
				### ### ### ###
				# variant skips both exon containing both start and stop codon
				### ### ### ###
				var_coding_seq = ""
			elif variant_skips_coding_exon:
				### ### ### ###
				# Simulate coding exon skipping
				# by deletion of coding exon sequence
				### ### ### ###
				if is_exon_skipped:
					# exon is skipped => delete exon from coding sequence
					self.logger.debug("Delete the skipped exon from the coding sequence")
					if var_start >= 1:
						var_coding_seq = coding_seq[0:var_start - 1] + coding_seq[var_end:len(coding_seq)]
						self.logger.debug("Deleting sequence: {}".format(coding_seq[var_start - 1:var_end]))
					else:
						try:
							assert var_start >= 1
						except AssertionError:
							self.logger.error(
								"Supplied coding region position should be >= 0\n=> variant position: {}".format(
									variant_info.to_string()), exc_info=True)

					# calculate deletion length and assert its result
					del_length = var_end - var_start + 1
					self.logger.debug(
						"len(coding)={}, len(var_coding)={}, del_length={}".format(len(coding_seq), len(var_coding_seq),
						                                                           del_length))
					try:
						assert len(coding_seq) - len(var_coding_seq) == del_length
					except AssertionError:
						self.logger.error(
							"For deletion should hold: len(reference coding) - len(variant coding) = len(deletion)\n=> variant position: {}".format(
								variant_info.to_string()), exc_info=True)
				diff_len = -1 * del_length
			else:
				var_coding_seq = coding_seq
		else:
			### ### ### ###
			# construct observed coding sequence for variant in exonic regions
			### ### ### ###
			if ">" in var_edit[
			          0:3]:
				self.logger.debug("Add SNP to coding sequence")
				# get directionally corrected reference and observed bases for SNP
				[ref_seq, obs_seq] = var_edit.split(">")
				if var_start > 1:
					var_coding_seq = coding_seq[0:var_start - 1] + obs_seq.lower() + coding_seq[
					                                                                 var_start:len(coding_seq)]
				elif var_start == 1:
					var_coding_seq = obs_seq.lower() + coding_seq[var_start:len(coding_seq)]
				else:
					try:
						assert var_start >= 0
					except AssertionError:
						self.logger.error(
							"Supplied coding region position should be >= 0\n=> variant position: {}".format(
								variant_info.to_string()), exc_info=True)
				try:
					assert len(coding_seq) == len(var_coding_seq)
				except AssertionError:
					self.logger.error(
						"For SNP the sum of length of coding exons should not change\n=> variant position: {}".format(
							variant_info.to_string()), exc_info=True)
				diff_len = 0
			elif "delins" in var_edit[0:6]:
				### ### ###
				# delins variant edits add a deletion and then an insertion
				### ### ###
				self.logger.debug("Add deletion and then insertion to coding sequence")

				### ### ###
				# introduce the deletion
				### ### ###
				self.logger.debug("Deleting: {}".format(coding_seq[var_start - 1:var_end]))
				self.logger.debug("VEP sequence to delete is: {}".format(transcript_info["var_seq"][0]))
				if var_start >= 1:
					var_coding_seq = coding_seq[0:var_start - 1] + coding_seq[var_end:len(coding_seq)]
				elif var_start == 0:
					var_coding_seq = coding_seq[var_end:len(coding_seq)]
				else:
					try:
						assert var_start >= 0
					except AssertionError:
						self.logger.error(
							"Supplied coding region position should be >= 0\n=> variant position: {}".format(
								variant_info.to_string()), exc_info=True)

				# assert deletion by resulted sequence length
				del_length = var_end - var_start + 1
				try:
					self.logger.debug("var_coding_seq: {}, len: {}".format(var_coding_seq, len(var_coding_seq)))
					assert len(coding_seq) - len(var_coding_seq) == del_length
				except AssertionError:
					self.logger.error(
						"For deletion should hold: len(reference coding) - len(variant coding) = len(deletion)\n variant position: {}".format(
							variant_info.to_string()), exc_info=True)
				### ### ###
				# introduce the insertion
				### ### ###
				obs_seq = transcript_info["var_seq"][1]
				try:
					assert len(obs_seq) >= 1
				except AssertionError:
					self.logger.error(
						"Supplied VEP does not contain insertion sequence\n=> variant position: {}".format(
							variant_info.to_string(variant_info.to_string())), exc_info=True)
				self.logger.debug("Insert sequence: {}".format(obs_seq))
				if var_start >= 1:
					var_coding_seq_ins = var_coding_seq[0:var_start - 1] + obs_seq.lower() + var_coding_seq[
					                                                                         var_start - 1:len(
						                                                                         coding_seq)]
				else:
					try:
						assert var_start >= 1
					except AssertionError:
						self.logger.error(
							"Supplied coding region position should be >= 1\n=> variant position: {}".format(
								variant_info.to_string()), exc_info=True)

				# assert insertion operation by resulted length
				ins_length = len(obs_seq)
				try:
					assert len(var_coding_seq_ins) - len(var_coding_seq) == ins_length
				except AssertionError:
					self.logger.error(
						"For insertion should hold: len(var_coding) - len(var_coding with deletion) = len(insertion)\n=> variant position: {}".format(
							variant_info.to_string()), exc_info=True)
				var_coding_seq = var_coding_seq_ins
				# assert the length difference for both operations
				diff_len = ins_length - del_length
				try:
					self.logger.debug("var_coding_seq: {}, len: {}".format(var_coding_seq, len(var_coding_seq)))
					assert len(var_coding_seq) - len(coding_seq) == diff_len
				except AssertionError:
					self.logger.error(
						"For insertion after deletion should hold: len(variant coding) - len(reference coding) = -len(deletion) + len(insertion)\n variant position: {}".format(
							variant_info.to_string()), exc_info=True)
			elif "del" in var_edit[0:3]:
				# deletion will be from the start up to end position
				self.logger.debug("Add deletion to coding sequence")
				if not self.is_transcript_type_splice_acceptor_donor(transcript_info["type_variant"]):
					# assert deleted coding sequence to equal the reference coding sequence
					# include case at which the variant start before or ends after the coding sequence
					if int(var_coding.pos.start.base) > 0 and int(var_coding.pos.end.base) > 0:
						if not (self.is_genomic_pos_in_coding_exon(transcript,
						                                           variant_info.genomic_start) and self.is_genomic_pos_in_coding_exon(
							transcript, variant_info.genomic_end)):
							# variant start or ends outside an coding exonic range
							self.logger.debug("Variant start or end not in a coding exonic range")
							self.logger.debug("var start: {}, end: {}".format(var_start, var_end))
							affected_coding_length = var_end - var_start + 1
							self.logger.debug(
								"Length of variant sequence that affects coding exonic sequence: {}".format(
									affected_coding_length))
						else:
							# variant starts and ends inside coding exonic range
							affected_coding_length = -1
				### ### ###
				# perform deletion
				### ### ###
				if var_start >= 1:
					var_coding_seq = coding_seq[0:var_start - 1] + coding_seq[var_end:len(coding_seq)]
				else:
					try:
						assert var_start >= 0
					except AssertionError:
						self.logger.error(
							"Supplied coding region position should be >= 0\n=> variant position: {}".format(
								variant_info.to_string()), exc_info=True)
				del_length = var_end - var_start + 1
				try:
					self.logger.debug("var_coding_seq: {}, len: {}".format(var_coding_seq, len(var_coding_seq)))
					assert len(coding_seq) - len(var_coding_seq) == del_length
				except AssertionError:
					self.logger.error(
						"For deletion should hold: len(reference coding) - len(variant coding) = len(deletion)\n variant position: {}".format(
							variant_info.to_string()), exc_info=True)
				diff_len = -1 * del_length
			elif "ins" in var_edit[0:3]:
				# for insertion annotation: start & end are the flanking regions
				# the insertion will be placed between the flanking regions
				self.logger.debug("Add insertion to coding sequence")
				if not (self.is_genomic_pos_in_coding_exon(transcript,
				                                           variant_info.genomic_start) and self.is_genomic_pos_in_coding_exon(
					transcript, variant_info.genomic_end)):
					# variant start or ends outside an coding exonic range
					self.logger.debug("Variant start or end not in a coding exonic range")
					affected_coding_length = var_end - var_start + 1
					self.logger.debug(
						"Length of variant sequence that affects coding exonic sequence: {}".format(
							affected_coding_length))
				else:
					# variant inside coding exonic range
					affected_coding_length = -1

				### ### ###
				# calculate sequence to be inserted
				### ### ###
				if len(transcript_info["var_seq"][0]) > 0:
					# VEP contains insertion sequence
					if "_" in transcript_info["var_seq"][0]:
						# insertion sequence is described by coding coordinates
						self.logger.info(
							"Insertion sequence is described by coding coordinates \n=> variant pos: {}".format(
								variant_info.to_string()))
						[ins_source_start, ins_source_end] = transcript_info["var_seq"][0].split("_")
						transcript_info["var_seq"][0] = coding_seq[
						                                int(ins_source_start) - 1: int(ins_source_end.strip())]
					if affected_coding_length == -1:
						# all insertion inside coding exon
						obs_seq = transcript_info["var_seq"][0]
					else:
						obs_seq = transcript_info["var_seq"][0][-affected_coding_length:]
				else:
					try:
						assert 1 == 0
					except AssertionError:
						self.logger.error(
							"vep does not contain insertion sequence\n=> variant position: {}".format(
								variant_info.to_string()), exc_info=True)
				self.logger.debug("Insert the sequence: {}".format(obs_seq))
				### ### ###
				# perform insertion
				### ### ###
				if var_start >= 1:
					if var_start == var_end:
						# start equal ends so the rest part, after the insertion, should start on end position
						if vep_coordinates_used:
							var_coding_seq = coding_seq[0:var_start] + obs_seq.lower() + coding_seq[
							                                                             var_end:len(coding_seq)]
						else:
							var_coding_seq = coding_seq[0:var_start] + obs_seq.lower() + coding_seq[
							                                                             var_end:len(coding_seq)]
					else:
						if vep_coordinates_used:
							var_coding_seq = coding_seq[0:var_start] + obs_seq.lower() + coding_seq[
							                                                             var_end - 1:len(coding_seq)]
						else:
							var_coding_seq = coding_seq[0:var_start] + obs_seq.lower() + coding_seq[
							                                                             var_end - 1:len(coding_seq)]
				else:
					try:
						assert var_start >= 1
					except AssertionError:
						self.logger.error(
							"Supplied coding region position should be >= 0\n=> variant position: {}".format(
								variant_info.to_string()), exc_info=True)
				ins_length = len(obs_seq)
				try:
					assert len(var_coding_seq) - len(coding_seq) == ins_length
				except AssertionError:
					self.logger.error(
						"For insertion should hold: len(var_coding) - len(reference_coding) = len(insertion)\n=> variant position: {}".format(
							variant_info.to_string()), exc_info=True)
				diff_len = +1 * ins_length

			elif "dup" in var_edit[0:3]:
				### ### ###
				# for duplication annotation: start & end are the duplication region
				# the duplication will be placed right-after the end position
				### ### ###
				self.logger.debug("Add duplication to coding sequence")
				if not (self.is_genomic_pos_in_coding_exon(transcript,
				                                           variant_info.genomic_start) and self.is_genomic_pos_in_coding_exon(
					transcript, variant_info.genomic_end)):
					### ### ###
					# variant starts or ends in intronic region
					### ### ###
					self.logger.debug("Variant start or end not in a coding exonic range")
					self.logger.debug("var start: {}, end: {}".format(var_start, var_end))
					affected_coding_length = var_end - var_start + 1
					self.logger.debug(
						"Length of variant sequence that affects coding exonic sequence: {}".format(
							affected_coding_length))
				else:
					### ### ###
					# variant inside coding exonic range
					### ### ###
					affected_coding_length = -1

				### ### ###
				# calculate sequence to be duplicated
				### ### ###
				if len(transcript_info["var_seq"][0]) > 0:
					if affected_coding_length == -1:
						# all duplication inside coding exon
						obs_seq = transcript_info["var_seq"][0]
					else:
						obs_seq = transcript_info["var_seq"][0][-affected_coding_length:]
				else:
					if var_start == var_end:
						obs_seq = coding_seq[var_start - 1]
					else:
						obs_seq = coding_seq[var_start - 1:var_end]
				self.logger.debug("Duplicate the sequence: {}".format(obs_seq))

				### ### ###
				# perform duplication
				### ### ###
				if var_start == var_end:
					# duplication of one nucleotide
					self.logger.debug("Duplication of one nucleotide")
					if var_start >= 1:
						var_coding_seq = coding_seq[0:var_start] + obs_seq.lower() + coding_seq[
						                                                             var_start:len(coding_seq)]
					else:
						try:
							assert var_start >= 0
						except AssertionError:
							self.logger.error(
								"Supplied coding region position should be >= 0\n=> variant position: {}".format(
									variant_info.to_string()), exc_info=True)
				else:
					# duplication of multiple nucleotides
					self.logger.debug("Duplication of multiple nucleotides")
					if var_start >= 1:
						if vep_coordinates_used:
							var_coding_seq = coding_seq[0:var_end] + obs_seq.lower() + coding_seq[
							                                                           var_end:len(coding_seq)]
						else:
							var_coding_seq = coding_seq[0:var_end + 1] + obs_seq.lower() + coding_seq[
							                                                               var_end + 1:len(coding_seq)]
					else:
						try:
							assert var_start >= 1
						except AssertionError:
							self.logger.error(
								"Supplied coding region position should be >= 0\n=> variant position: {}".format(
									variant_info.to_string()), exc_info=True)
				dupl_length = len(obs_seq)
				try:
					assert len(var_coding_seq) - len(coding_seq) == dupl_length
				except AssertionError:
					self.logger.error(
						"For duplication should hold: len(var_coding) - len(reference_coding) = len(duplication)\n variant position: {}".format(
							variant_info.to_string()), exc_info=True)
				diff_len = +1 * dupl_length

		### ### ### ###
		# after creating coding sequence with variant,
		# check that indeed is different from the reference
		### ### ### ###
		if variant_skips_coding_exon or not self.is_transcript_type_splice_acceptor_donor(
				transcript_info["type_variant"]):
			try:
				assert coding_seq.upper() != var_coding_seq.upper()
			except AssertionError:
				self.logger.error(
					"Coding sequence of reference and sample should be different\n variant position: {}".format(
						variant_info.to_string()), exc_info=True)
			self.print_ref_observed_seq(coding_seq, var_coding_seq, transcript_info, var_start, var_end,
			                            var_edit, is_exon_skipped)

			### ### ### ### ### ###
			# assert constructed variant coding sequence
			# contains start (if edit not start_lost) and stop codon
			### ### ### ### ### ###
			if variant_skips_start_codon_exon and variant_skips_stop_codon_exon:
				pass
			elif variant_skips_start_codon_exon:
				try:
					assert var_coding_seq[-3:] in ["TAG", "TAA", "TGA"]
				except AssertionError:
					self.logger.error(
						"Constructed observed coding sequence should contain stop codon\n=> variant position: {}".format(
							variant_info.to_string()), exc_info=True)
			elif variant_skips_stop_codon_exon:
				if coding_seq[0:3] == "ATG":
					try:
						assert var_coding_seq[0:3] == "ATG"
					except AssertionError:
						self.logger.debug(
							"Constructed observed coding sequence should contain start codon\n=> variant position: {}".format(
								variant_info.to_string()))
			else:
				### ### ###
				# variant does not skip exon => prepare range of variant start stop
				# to assert for existence of start and stop codons
				### ### ###
				if var_start == var_end:
					var_positions = range(var_start, var_start + 1)
				else:
					var_positions = range(var_start, var_end)

				if not ("start_lost" in transcript_info[
					"type_variant"] or "start_retained" in variant_info.variant_type):
					start_positions = set(range(1, 4))
					if not start_positions.intersection(var_positions) and coding_seq[0:3] == "ATG":
						# if variant does not change positions on the start of the coding sequence
						# and ensembl record for transcript starts by ATG
						# assert that the observed coding sequence starts with ATG
						try:
							assert var_coding_seq[0:3] == "ATG"
						except AssertionError:
							self.logger.error(
								"Constructed observed coding sequence should contain start codon\n=> variant position: {}".format(
									variant_info.to_string()), exc_info=True)

				if not ("stop_retained_variant" in transcript_info["type_variant"] or "stop_lost" in transcript_info[
					"type_variant"]):
					stop_positions = set(range(len(coding_seq) + 1 - 3, len(coding_seq) + 1))
					if not stop_positions.intersection(var_positions):
						# if variant does not change positions in the end of the coding sequence
						# assert that the observed coding sequence finishes with "TAG", "TAA" or "TGA"
						try:
							assert var_coding_seq[-3:] in ["TAG", "TAA", "TGA"]
						except AssertionError:
							self.logger.error(
								"Constructed observed coding sequence should contain stop codon\n=> variant position: {}".format(
									variant_info.to_string()), exc_info=True)

		return var_coding_seq.upper(), diff_len

	def print_ref_observed_seq(self, coding_seq, var_coding_seq, transcript_info, var_start, var_end,
	                           var_edit, is_exon_skipped):
		"""
		Print reference and observed coding sequence

		Parameters
		----------
		coding_seq : str
			reference coding sequence
		var_coding_seq : str
			observed coding sequence
		transcript_info : dict of str: str
			variant-affected transcript GSvar information
		var_start : int
			variant start position (in coding sequence)
		var_end : int
			variant end position (in coding sequence)
		var_edit : str
			variant edit information
		is_exon_skipped : bool
			variant causes an exon to be skipped (True), otherwise (False)

		Returns
		-------
		None
		"""
		self.logger.debug("Print reference and observed sequence on the variant region:")
		# set up the variant print region start and end
		if var_start >= 11:
			# for variants after the 11th genomic position, print 10 before and 10 after bases
			print_start = 25
			print_end = 25
		elif 1 < var_start < 11:
			# for variants close to start of the chromosome, print bases as much as to show variant
			# and 10 bases after
			print_start = var_start - 1
			print_end = 25
		else:
			# for variants on the first base of the chromosome, print the start of the chromosome and 10 bases after
			print_start = 0
			print_end = 25
		### ### ###
		# print variant coding sequence and coding sequence
		### ### ###
		if "del" in var_edit[0:3] or (
				self.is_transcript_type_splice_acceptor_donor(transcript_info["type_variant"]) and is_exon_skipped):
			# deletion or skipped exon case
			ref_deleted_var = coding_seq[var_start - print_start - 1:var_start - 1] + coding_seq[
			                                                                          var_start - 1:var_end].lower() + coding_seq[
			                                                                                                           var_end:var_end + print_end]
			self.logger.debug(">coding seq:\n {}".format(ref_deleted_var))
		else:
			# insertion, duplication, SNP
			self.logger.debug(">coding seq:\n {}".format(coding_seq[var_start - print_start - 1:var_end + print_end]))

		if var_start == 0:
			self.logger.debug(">var coding seq:\n {}".format(var_coding_seq[var_start:var_end + print_end]))
		else:
			self.logger.debug(
				">var coding seq:\n {}".format(var_coding_seq[var_start - print_start - 1:var_end + print_end]))

	@staticmethod
	def get_codon_frame_start_offset(exon_start_offset):
		"""
		Get the frame start offset for input exon start

		Parameters
		----------
		exon_start_offset : int
			exon start offset

		Returns
		-------
		int
			frame start offset for given exon
		"""
		# coding_sequence from pyensembl start for ATG (start codon)
		# so check how many nucleotides are left from the open reading frame (ORF) from the previous exon
		num_left_nucl = (exon_start_offset - 1) % 3
		# then return the offset of starting reading frame for the exon
		return exon_start_offset - num_left_nucl

	def assess_reading_frame_preservation(self, transcript_info):
		"""
		Check if reading frame is preserved

		Parameters
		----------
		transcript_info : dict of str: str
			transcript information

		Returns
		-------
		bool
			reading frame is preserved (True), otherwise is disrupted (False)
		"""
		# print("Check preservation of reading frame on transcript with id:{}".format(transcript_info["transcript_id"]))
		# get the difference of length of variant and reference coding sequence
		if transcript_info["coding_diff_len"] % 3 == 0:  # reading frame is preserved
			return True
		else:
			# reading frame is disrupted
			return False

	def get_codons_downstream_start(self, transcript_info, downstream_length):
		"""
		Get codons in downstream region from start codon of transcript

		Parameters
		----------
		transcript_info : dict of str: str
			transcript information
		downstream_length : int
			create codons up to given downstream length

		Returns
		-------
		list of str
			codons in downstream region of mutated start codon
		"""
		# self.logger.debug("Get codons downstream of start codon for transcript id: {}".format(transcript.id))
		return extract_codons(transcript_info["var_coding_seq"][0: downstream_length])

	def search_closest_start_codon(self, codons, variant_info):
		"""
		Search for start codon in input codons
		if exists, return the closest to the start

		Parameters
		----------
		codons : list of str
			input list of codons
		find_closest : bool
			find closest codon to start
		variant_info : VariantInfo
			variant basic info

		Returns
		-------
		int
			positive value = index of start codon in input list of codons,
			negative value = start codon is not contained in the input list of codons
		"""
		# print("Search for start codon")
		try:
			return codons.index("ATG")
		except ValueError:
			self.logger.debug(
				"Codons do not contain start codon\n=> variant position: {}".format(variant_info.to_string()))
			return -1

	def search_termination_codon(self, exon_codons, is_penultimate_exon):
		"""
		Search termination codon in exon codons

		Parameters
		----------
		exon_codons : list of str
			exon codons
		is_penultimate_exon : bool
			the input exon is the one before the last

		Returns
		-------
		bool
			termination codon exists in exon codons (True), otherwise it doesn't (False)
		int
			termination codon index (1-index), if termination codon doesn't exist equals the length of the search coding region
		"""
		if is_penultimate_exon:
			self.logger.debug("Search termination codon in penultimate exon")
		termination_codon_exists = False
		if is_penultimate_exon:
			# start from the most 3' 50 base ~= 50/3 = 17 codons
			search_start = len(exon_codons) - 17
			sequence_limit = "the most 3' 50 bases = 17 codons"
		# print("Current exon observed sequence ")
		else:
			# not penultimate exon
			search_start = 0
			sequence_limit = "all sequence"
		self.logger.debug("Current exon observed sequence ({}) is: {}".format(sequence_limit, exon_codons[
		                                                                                      search_start:len(
			                                                                                      exon_codons)]))
		# find the left most (5'-most) termination codon in sequence's codons
		term_codons_indices = [get_codon_index(exon_codons[search_start:len(exon_codons)], term_codon) for term_codon in
		                       ["TAG", "TAA", "TGA"]]
		left_most_term_index = min(term_codons_indices)

		# update termination codon existence flag and position
		if left_most_term_index == len(exon_codons[search_start:len(exon_codons)]):
			termination_codon_exists = False
			termination_codon_index = -1
		else:
			termination_codon_exists = True
			termination_codon_index = left_most_term_index
		self.logger.debug("Termination codon found: {}".format(termination_codon_exists))
		return termination_codon_exists, termination_codon_index

	def find_affected_exons_pos(self, transcript, affected_exon_idxs, start_offset, end_offset, variant_info):
		"""
		Find affected exons positions (genomic or coding)

		Parameters
		----------
		transcript : pyensembl.transcript
			transcript of which exons will be extracted
		affected_exon_idxs : list of int
			1-based indices of affected exons
		start_offset : int
			variant offset in affected exon(s)
		end_offset : int
			variant offset in affected exons(s)
		variant_info : VariantInfo
			variant basic info

		Returns
		-------
		list of dict of str: str or int
			affected exons: id, start and end (genomic positions)
		"""
		self.logger.debug("Find affected exon positions")
		self.logger.debug("Start offset: {}, end: {}".format(start_offset, end_offset))
		### ### ###
		# normalize start and end offset
		### ### ###
		end_offset = end_offset - start_offset
		start_offset = 0
		interact_exon_pos = []
		if RefineLossOfFunction.is_transcript_in_positive_strand(transcript):
			strand_direction = 1
		else:
			strand_direction = -1
		if len(affected_exon_idxs) == 1:
			### ### ###
			# one affected exon
			### ### ###
			self.logger.debug("One affected exon with index: {}".format(affected_exon_idxs))

			exon = transcript.exons[affected_exon_idxs[0] - 1]
			self.logger.debug("selected exon: {}".format(exon))
			if strand_direction == 1:
				interact_exon_pos.append({"exon_id": exon.id, "exon_start": exon.to_dict()["start"] + start_offset,
				                          "exon_end": exon.to_dict()["start"] + end_offset})
			else:
				interact_exon_pos.append({"exon_id": exon.id, "exon_start": exon.to_dict()["end"] - end_offset,
				                          "exon_end": exon.to_dict()["end"] - start_offset})
		else:
			### ### ###
			# multiple affected exons
			### ### ###
			for idx, exon in enumerate(transcript.exons):
				if idx + 1 in affected_exon_idxs:
					if affected_exon_idxs.index(idx + 1) == 0:
						# first affected exon
						if strand_direction == 1:
							interact_exon_pos.append({"exon_id": exon.id, "exon_start": exon.to_dict()["start"],
							                          "exon_end": exon.to_dict()["start"] + start_offset})
						else:
							interact_exon_pos.append(
								{"exon_id": exon.id, "exon_start": exon.to_dict()["start"],
								 "exon_end": exon.to_dict()["end"] - end_offset})
					elif affected_exon_idxs.index(idx + 1) == len(affected_exon_idxs) - 1:
						# last affected exon
						if RefineLossOfFunction.is_transcript_in_positive_strand(transcript):
							interact_exon_pos.append({"exon_id": exon.id, "exon_start": exon.to_dict()["start"],
							                          "exon_end": exon.to_dict()["start"] + end_offset})
						else:
							interact_exon_pos.append(
								{"exon_id": exon.id, "exon_start": exon.to_dict()["end"] - start_offset,
								 "exon_end": exon.to_dict()["end"]})
					else:
						# middle affected exons
						interact_exon_pos.append(
							{"exon_id": exon.id, "exon_start": exon.to_dict()["start"],
							 "exon_end": exon.to_dict()["end"]})

		# assert that created information for interacting exon have start <= end
		for interacting_exon in interact_exon_pos:
			try:
				assert interacting_exon["exon_start"] <= interacting_exon["exon_end"]
			except AssertionError:
				self.logger.error("Affected exon start is higher than end\n=> variant position: {}".format(
					variant_info.to_string()), exc_info=True)

		self.logger.debug("Genomic positions of affected exons: {}".format(interact_exon_pos))
		return interact_exon_pos

	def construct_observed_exon_seq(self, transcript, exon_idx, exons_containing_var,
	                                exon_contains_start_codon, exon_contains_stop_codon, exon_contains_coding_seq,
	                                diff_len, last_start,
	                                current_codon_length, remain_codon_length):
		"""
		Construct observed exon coding sequence
		Parameters
		transcript : pyensembl.transcript
			affected transcript object
		exon_idx : int
			affected exon index
		exons_containing_var : list of int
			list of exons containing variant
		exon_contains_start_codon : bool
			exon contains start codon (True), otherwise exon contains only UTR (False)
		exon_contains_stop_codon : bool
			exon contains stop codon (True), otherwise exon contains no stop codon (False)
		exon_contains_coding_seq : bool
			exon contains coding sequence (True), otherwise (False)
		diff_len : int
			difference of the length between reference and observed seqeunce
		last_start : int
			last coding exonic start position
		current_codon_length : int
			current exon coding length
		remain_codon_length : int
			coding length remained from previous (left) exon

		Returns
		-------
		exon_contains_start_codon : bool
			exon contains start codon (True), otherwise exon contains only UTR (False)
		exon_contains_coding_seq : bool
			exon contains coding sequence (True), otherwise (False)
		last_start : int
			updated last coding exonic start position
		current_codon_length : int
			updated current exon coding length
		remain_codon_length : int
			updated coding length remained from previous (left) exon
		"""

		self.logger.debug("Construct observed exon sequence")
		self.logger.debug("Diff len: {}".format(diff_len))
		exon = transcript.exons[exon_idx]
		exon_start, exon_end = exon.to_dict()["start"], exon.to_dict()["end"]
		# print("exon start: {}, exon end: {}".format(exon_start, exon_end))

		# set up first position of start codon based on transcript strand
		if RefineLossOfFunction.is_transcript_in_positive_strand(transcript):
			start_codon_first_pos = transcript.start_codon_positions[0]
		else:
			start_codon_first_pos = transcript.start_codon_positions[-1]

		# set up stop codon as the end of the last exon
		if exon_contains_stop_codon:
			self.logger.debug("this is last exon, update its end to be the stop codon position")
			self.logger.debug("before: exon start:{}, exon end:{}".format(exon_start, exon_end))
			# exon_end = transcript.stop_codon_positions[-1]
			if RefineLossOfFunction.is_transcript_in_positive_strand(transcript):
				exon_start = exon_start
				exon_end = transcript.stop_codon_positions[-1]
			else:
				# print("stop codon positions: {}".format(transcript.stop_codon_positions))
				exon_start = transcript.stop_codon_positions[0]
				exon_end = exon_end
		self.logger.debug("updated: exon start:{}, exon end:{}".format(exon_start, exon_end))
		coding_exonic_length = 0

		### ### ### ### ### ###
		# update current exonic coding sequence
		# and its remaining for the next exon
		### ### ### ### ### ###
		if exon_contains_coding_seq and exon_contains_start_codon:
			# first exon with coding sequence
			self.logger.debug("First exon containing coding sequence")
			self.logger.debug("last start: {}".format(last_start))
			if RefineLossOfFunction.is_transcript_in_positive_strand(transcript):
				coding_exonic_length = (exon_end - start_codon_first_pos) + 1
			else:
				coding_exonic_length = (start_codon_first_pos - exon_start) + 1
			self.logger.debug("before: adding variant edit, coding_exonic_length={}".format(coding_exonic_length))
			if exon_idx + 1 == exons_containing_var and diff_len != 0:
				self.logger.debug(
					"add variant edit len={} on the current length of exon idx:{}".format(diff_len, exon_idx + 1))
				# add variant edit on the current exonic length
				coding_exonic_length = coding_exonic_length + diff_len
			self.logger.debug(
				"after update with variant length, coding_exonic = {}, remaining = {}".format(coding_exonic_length,
				                                                                              remain_codon_length))
			self.logger.debug("remain codon_length= {}".format(remain_codon_length))
			coding_exonic_length = coding_exonic_length + remain_codon_length
			remain_codon_length = coding_exonic_length % 3
			current_codon_length = coding_exonic_length - remain_codon_length
			self.logger.debug(
				"coding exonic length={}, remain codon length={}, current codon length={}".format(coding_exonic_length,
				                                                                                  remain_codon_length,
				                                                                                  current_codon_length))
			# deactivate flag for next exon with coding sequence
			exon_contains_start_codon = False
		elif exon_contains_coding_seq:
			# exon containing coding sequence, but not the start codon
			last_start = last_start + current_codon_length
			self.logger.debug("updated last start: {}".format(last_start))
			coding_exonic_length = (exon_end - exon_start) + 1
			if coding_exonic_length < 0:
				self.logger.debug("exon end is lower than exon start, stop codon between two exons")
				coding_pos = transcript.coding_sequence_position_ranges[-1]
				self.logger.debug("coding positions: {}".format(coding_pos))
				coding_exonic_length = (coding_pos[1] - coding_pos[0]) + 1

			self.logger.debug("before: adding variant edit, coding_exonic_length={}".format(coding_exonic_length))
			if exon_idx + 1 == exons_containing_var[0] and diff_len != 0:
				self.logger.debug(
					"add variant edit len={} on the current length of exon idx:{}".format(diff_len, exon_idx + 1))
				coding_exonic_length = coding_exonic_length + diff_len

			self.logger.debug(
				"after update with variant length, coding_exonic = {}, remaining = {}".format(coding_exonic_length,
				                                                                              remain_codon_length))
			coding_exonic_length = coding_exonic_length + remain_codon_length
			remain_codon_length = coding_exonic_length % 3
			current_codon_length = coding_exonic_length - remain_codon_length
			self.logger.debug(
				"coding exonic length={}, remain codon length={}, current codon length={}".format(coding_exonic_length,
				                                                                                  remain_codon_length,
				                                                                                  current_codon_length))

		return exon_contains_coding_seq, exon_contains_start_codon, last_start, current_codon_length, remain_codon_length

	def search_stop_codon(self, transcript, transcript_info, exons_containing_var, is_exon_skipped,
	                      variant_skips_stop_codon_exon, variant_info):
		"""
		Search for stop codon on all observed exonic coding sequences

		Parameters
		----------
		transcript : pyensembl.transcript
			affected transcript object
		transcript_info : dict of str: str
			affected transcript information
		exons_containing_var : list of int
			exon indices that contain variant
		is_exon_skipped : bool
			the affected exon is skipped (True), otherwise it is not skipped
		variant_skips_stop_codon_exon : bool
			variant skips stop codon exon (True), otherwise variant does not result skipping the exon with stop codon
		variant_info : VariantInfo
			variant basic info

		Returns
		-------
		dict of int: int
			map of exon index to (left-most) termination codon position
		int
			number of exons
		"""
		self.logger.debug("Search stop codon over all exons' observed coding sequence")
		# create the variant coding sequence of the trancript

		var_coding_seq, diff_len = transcript_info["var_coding_seq"], transcript_info["coding_diff_len"]
		self.logger.debug("Difference of observed to reference, in length: {}".format(diff_len))
		current_codon_length, remain_codon_length, last_start = 0, 0, 0
		if RefineLossOfFunction.is_transcript_in_positive_strand(transcript):
			start_codon_first_pos = transcript.start_codon_positions[0]
			# find index of exon that contains stop codon
			self.logger.debug("Stop codon positions: {}".format(transcript.stop_codon_positions))
			stop_codon_first_base_exon_idx, _ = self.find_exon_by_ref_pos(transcript,
			                                                              transcript.stop_codon_positions[0], True)
			stop_codon_last_base_exon_idx, _ = self.find_exon_by_ref_pos(transcript, transcript.stop_codon_positions[2],
			                                                             True)
		else:
			start_codon_first_pos = transcript.start_codon_positions[-1]
			# find index of exon that contains stop codon
			stop_codon_first_base_exon_idx, _ = self.find_exon_by_ref_pos(transcript,
			                                                              transcript.stop_codon_positions[-1], True)
			stop_codon_last_base_exon_idx, _ = self.find_exon_by_ref_pos(transcript, transcript.stop_codon_positions[0],
			                                                             True)

		# if stop codon starts on one exon and finishes on the next one
		# keep the last exon as the position of the stop codon
		stop_codon_exon_idx = max(stop_codon_first_base_exon_idx, stop_codon_last_base_exon_idx)
		exon_contains_coding_seq, exon_contains_start_codon, exon_contains_stop_codon = False, False, False
		sum_obs_exon_coding_seq = 0
		exon2termination_codon_pos = {}
		exon_idx = 0
		num_exons = len(transcript.exon_intervals)

		while not exon_contains_stop_codon and exon_idx < num_exons:
			(exon_start, exon_end) = transcript.exon_intervals[exon_idx]
			self.logger.debug("Exon idx: {}, start={} end={}".format(exon_idx + 1, exon_start, exon_end))

			# examine if exon contains stop codon
			if exon_idx == stop_codon_exon_idx:
				exon_contains_stop_codon = True
			else:
				exon_contains_stop_codon = False

			### ### ### ### ### ###
			# skip current exon if intron variant results to skip exon
			### ### ### ### ### ###
			if exon_idx + 1 == exons_containing_var[0] and is_exon_skipped:
				self.logger.debug("Exon with idx:{} is skipped".format(exon_idx + 1))
				# we should not find that exon with start codon is skipped
				# because this case is captured by start_lost variant type
				exon_idx += 1
				continue

			### ### ### ### ### ###
			# do not construct coding sequence for exons up to the exon that contains start codon
			### ### ### ### ### ###
			if not exon_contains_coding_seq:
				if exon_start <= start_codon_first_pos <= exon_end:
					self.logger.debug("First exon with start codon has index: {}".format(exon_idx + 1))
					self.logger.debug("last start: {}".format(last_start))
					exon_contains_start_codon = True
					exon_contains_coding_seq = True

					exon_contains_coding_seq, exon_contains_start_codon, last_start, current_codon_length, remain_codon_length = self.construct_observed_exon_seq(
						transcript, exon_idx, exons_containing_var,
						exon_contains_start_codon, exon_contains_stop_codon, exon_contains_coding_seq, diff_len,
						last_start, current_codon_length, remain_codon_length)
				else:
					self.logger.debug("UTR region in exon, before start codon")
					exon_idx += 1
					continue
			else:
				# exon contains coding sequence
				# construct current exon observed coding sequence
				exon_contains_coding_seq, exon_contains_start_codon, last_start, current_codon_length, remain_codon_length = self.construct_observed_exon_seq(
					transcript, exon_idx, exons_containing_var,
					exon_contains_start_codon, exon_contains_stop_codon, exon_contains_coding_seq, diff_len,
					last_start, current_codon_length, remain_codon_length)

			if exon_contains_coding_seq:
				### ### ### ### ### ###
				# get observed current exon coding sequence and codons
				### ### ### ### ### ###
				if exon_contains_stop_codon and exon_start <= start_codon_first_pos <= exon_end:
					# add all remaining sequence if current exon contains both start and stop codons
					self.logger.debug("Exon contains both start and stop codons")
					obs_exon_coding_seq = var_coding_seq[last_start:len(var_coding_seq)]
					self.logger.debug("observed exonic codons seq= {}, length of seq= {}".format(obs_exon_coding_seq,
					                                                                             len(
						                                                                             obs_exon_coding_seq)))
				# self.logger.debug("remaining: {}".format(
				# 	var_coding_seq[last_start + current_codon_length + remain_codon_length:len(var_coding_seq)]))
				elif exon_contains_stop_codon:
					# on exon containing the stop codon, add the remaining coding sequence to the observed codon sequence
					self.logger.debug("Exon contains stop codon")
					obs_exon_coding_seq = var_coding_seq[
					                      last_start:last_start + current_codon_length + remain_codon_length]
					self.logger.debug("observed exonic codons seq= {}, length of seq= {}".format(obs_exon_coding_seq,
					                                                                             len(
						                                                                             obs_exon_coding_seq)))
					self.logger.debug("remaining_ending: {}".format(
						var_coding_seq[last_start + current_codon_length + remain_codon_length:len(var_coding_seq)]))
					self.logger.debug("transcript type: {}".format(transcript_info["type_variant"]))
				# obs_exon_coding_seq = obs_exon_coding_seq + var_coding_seq[last_start+current_codon_length+remain_codon_length:len(var_coding_seq)]
				else:
					# update observed coding sequence for coding exon
					obs_exon_coding_seq = var_coding_seq[last_start:last_start + current_codon_length]
					self.logger.debug("observed exonic codons seq= {}, length of seq= {}".format(obs_exon_coding_seq,
					                                                                             len(
						                                                                             obs_exon_coding_seq)))
				# self.logger.debug("remaining: {}".format(
				# 	var_coding_seq[last_start + current_codon_length:len(var_coding_seq)]))

				observed_codons = extract_codons(obs_exon_coding_seq)

				if exon_idx == stop_codon_exon_idx - 1:
					self.logger.debug("Process penultimate exon")
					# if exon_contains_stop_codon:
					# search termination codon on penultimate exon
					termination_codon_exists, termination_codon_index = self.search_termination_codon(observed_codons,
					                                                                                  True)
				else:
					# search termination codon on any other exon
					termination_codon_exists, termination_codon_index = self.search_termination_codon(observed_codons,
					                                                                                  False)

				if termination_codon_exists:
					# to create absolute termination codon indices
					# add the current sum of observed coding sequence
					exon2termination_codon_pos[exon_idx + 1] = sum_obs_exon_coding_seq + termination_codon_index * 3
				sum_obs_exon_coding_seq += len(obs_exon_coding_seq)
			# update exon index
			exon_idx += 1

		self.logger.debug("sum observed exon coding seq: {}".format(sum_obs_exon_coding_seq))
		self.logger.debug("length of var_coding_seq={} and sum processed coding seq={}".format(len(var_coding_seq),
		                                                                                       sum_obs_exon_coding_seq))

		if not variant_skips_stop_codon_exon:
			try:
				assert sum_obs_exon_coding_seq == len(var_coding_seq) or sum_obs_exon_coding_seq + len(
					var_coding_seq) - (last_start + current_codon_length + remain_codon_length) == len(var_coding_seq)
			except AssertionError:
				self.logger.error(
					"Sum of observed exonic codons should equal the total length of the observed coding sequence\n=> variant position: {}".format(
						variant_info.to_string()), exc_info=True)
		return exon2termination_codon_pos, num_exons

	def assess_NMD(self, transcript_info, variant_info):
		"""
		Examine if "non sense mediated mRNA decay" (NMD) will occur for current variant
		Following Tayoun et al. and Zhiyuan et al.  NMD is not predicted to occur if:
		a) premature termination codon (PTC) occurs in the last exon
		b) PTC occur in the (3') last 50 nucleotides of the penultimate exon
		c) transcript contains no introns
		d) PTC occurs inbetween the first 200 bases from start codon

		Parameters
		----------
		transcript_info : dict of str: str
			transcript information
		variant_info : VariantInfo
			variant basic info

		Returns
		-------
		bool
			NMD predicted to occur (True), not to occur (False)
		list of dict of str: str or int
			var containing exon genomic positions
		bool
			variant skips exon containing start codon
		bool
			variant skips exon containing stop codon
		bool
			variant skips coding exon
		"""

		self.logger.debug("Assess NMD for transript id: {}".format(transcript_info["transcript_id"]))
		NMD_occurs = True
		transcript = self.ensembl_data.transcript_by_id(transcript_info["transcript_id"])
		# construct exon positions using chromosome position
		exon_positions = self.get_transcript_exon_offsets(transcript, True)
		num_exon_positions = len(exon_positions)
		### ### ###
		# 4 cases: find exon containing the variant position
		# a) intronic variant that genomic start and/or end overlaps with exonic range
		# b) splice acceptor variant (that results to skipping)
		# c) exonic variant whose both start and stop overlap exonic range
		# d) exonic variant whose start or stop overlaps exonic range
		### ### ###
		if len(transcript_info["type_variant"].intersection(self.intron_variant_types)) > 0 and (
				self.is_genomic_pos_in_coding_exon(transcript,
				                                   variant_info.genomic_start) or self.is_genomic_pos_in_coding_exon(
			transcript, variant_info.genomic_end)):
			### ### ###
			# a) variant overlaps both intronic and exonic region
			### ### ###
			self.logger.debug(
				"Variant genomic position overlaps both intronic and coding exonic range => use exon offsets (VEP coordinates) to find affected exon")
			exons_containing_var, var_exon_start_offset, var_exon_end_offset = self.find_exon_by_var_pos(transcript,
			                                                                                             transcript_info,
			                                                                                             variant_info,
			                                                                                             False)
			if self.is_transcript_type_splice_acceptor_donor(transcript_info["type_variant"]):
				# if variant is splice acceptor, update exon skipping boolean flags
				var_coding = transcript_info["var_coding"]
				split_symbols, intron_offsets, directions2exon = self.parse_variant_intron_pos(var_coding)
				exon_containing_var, is_exon_skipped, var_exon_start_offset, var_exon_end_offset, variant_skips_start_codon_exon, variant_skips_stop_codon_exon, variant_skips_coding_exon = self.assess_exon_skipping(
					transcript,
					transcript_info,
					intron_offsets,
					variant_info)
		elif self.is_transcript_type_splice_acceptor_donor(transcript_info["type_variant"]):
			self.logger.debug("assess NMD for splice acceptor variant")
			### ### ###
			# b) splice acceptor
			# Currently handling only splice acceptor case that results to skipping
			# find affected exon index, start and stop and boolean flags by assess_exon_skipping()
			### ### ###
			var_coding = transcript_info["var_coding"]
			split_symbols, intron_offsets, directions2exon = self.parse_variant_intron_pos(var_coding)
			exons_containing_var, is_exon_skipped, var_exon_start_offset, var_exon_end_offset, variant_skips_start_codon_exon, variant_skips_stop_codon_exon, variant_skips_coding_exon = self.assess_exon_skipping(
				transcript,
				transcript_info,
				intron_offsets,
				variant_info)
			self.logger.debug("")
		elif not self.is_transcript_type_splice_acceptor_donor(transcript_info["type_variant"]) and (
				self.is_genomic_pos_in_coding_exon(transcript,
				                                   variant_info.genomic_start) and self.is_genomic_pos_in_coding_exon(
			transcript, variant_info.genomic_end)):
			### ### ###
			# c) exonic variant with both genomic start and stop in exonic range
			### ### ###
			self.logger.debug("Exonic type of variant and its genomic start and stop both overlap exonic range")
			exons_containing_var, var_exon_start_offset, var_exon_end_offset = self.find_exon_by_var_pos(transcript,
			                                                                                             transcript_info,
			                                                                                             variant_info,
			                                                                                             True)
		elif not self.is_transcript_type_splice_acceptor_donor(transcript_info["type_variant"]) and (
				self.is_genomic_pos_in_coding_exon(transcript,
				                                   variant_info.genomic_start) or self.is_genomic_pos_in_coding_exon(
			transcript, variant_info.genomic_end)):
			### ### ###
			# d) exonic variant with either genomic start or stop in exonic range
			### ### ###
			self.logger.debug("Exonic type of variant and its genomic start or stop overlaps exonic range")
			exons_containing_var, var_exon_start_offset, var_exon_end_offset = self.find_exon_by_var_pos(transcript,
			                                                                                             transcript_info,
			                                                                                             variant_info,
			                                                                                             False)

		if not self.is_transcript_type_splice_acceptor_donor(transcript_info["type_variant"]):
			### ### ### ### ### ### ###
			#   variant inside exon   #
			### ### ### ### ### ### ###
			# exonic variant type can't lead to exon skipping
			is_exon_skipped = False
			variant_skips_start_codon_exon, variant_skips_stop_codon_exon = False, False
			variant_skips_coding_exon = False

		### ### ### ### ### ###
		# find interacting exons genomic positions
		# if splice acceptor variant results to exon skipping
		# or the variant is exonic
		### ### ### ### ### ###
		if variant_skips_coding_exon or len(self.exon_variant_types.intersection(transcript_info["type_variant"])) > 0:
			self.logger.debug("Variant applicable for stop codon searching")
			if (variant_skips_start_codon_exon and variant_skips_stop_codon_exon) or not variant_skips_start_codon_exon:
				self.logger.debug(
					"Update exon positions for variant that skips both start and stop exons or does not skip start codon exon")
				affected_exons_pos = self.find_affected_exons_pos(transcript, exons_containing_var,
				                                                  var_exon_start_offset,
				                                                  var_exon_end_offset, variant_info)

				exon2stop_index, num_exons = self.search_stop_codon(transcript, transcript_info, exons_containing_var,
				                                                    is_exon_skipped, variant_skips_stop_codon_exon,
				                                                    variant_info)

				if num_exons == 1:
					# intronless transcript
					NMD_occurs, NMD_comment = False, "Single exon"
				else:
					exons_with_stop_codon = sorted(list(exon2stop_index.keys()))
					if len(exon2stop_index) == 0:
						# PTC not found in any exon
						# => presumably in the 3' non coding region of the last exon
						NMD_occurs, NMD_comment = False, "PTC after reference stop codon"
					elif exons_with_stop_codon[0] == num_exon_positions - 1 or exons_with_stop_codon[
						0] == num_exon_positions:
						# PTC in 3'-most 50 bases of penultimate or ultimate exon
						NMD_occurs, NMD_comment = False, "PTC in last exon or 3'-most 50 bases of penultimate exon"
					elif exon2stop_index[exons_with_stop_codon[0]] < 200:
						# PTC in the first 200 bases from start codon
						NMD_occurs, NMD_comment = False, "PTC distance from start codon < 200"
					elif exon2stop_index[exons_with_stop_codon[0]] > 200:
						# PTC after the first 200 bases from start codon
						# but before the last exon junction complex EJC
						NMD_occurs, NMD_comment = True, "PTC before last EJC"

				self.logger.debug("Transcript contains in total: {} exon(s)".format(num_exons))
				self.logger.debug("positions of stop codons: {}".format(exon2stop_index))
				self.logger.debug("NMD is predicted to occur: {}, comment: {}".format(NMD_occurs, NMD_comment))
			else:
				affected_exons_pos = []
				NMD_occurs = False
		else:
			self.logger.debug("Variant type not applicable for stop codon search")
			affected_exons_pos = []

		return NMD_occurs, affected_exons_pos, variant_skips_start_codon_exon, variant_skips_stop_codon_exon, variant_skips_coding_exon

	def intersect_affected_region2significant_exons(self, affected_exons_pos, variant_info, transcript_strand):
		"""
		Intersect affected exons with significant exons

		Parameters
		----------
		affected_exons_pos : list of dict of str: str or int
			affected exons information
		variant_info : VarianInfo object
			variant basic info
		transcript_strand : str
			affected transcript strand

		Returns
		-------
		bool
			affected exons are contained in the list of significant exons (True), otherwise (False)
		list of str
			matched clinical significant exons Ensembl ids
		"""
		# print("Examine if truncated/altered exons are included in significant exons ")
		matched_exons = self.clinical_significant_exons.all_hits(
			RefineLossOfFunction.truncated_exons2bed_interval(affected_exons_pos, variant_info,
			                                                  transcript_strand), same_strand=True)
		if len(matched_exons):
			return True, [matched_exon.name for matched_exon in matched_exons]
		else:
			return False, []

	def intersect_affected_region2uniprot_domains(self, affected_exon_pos, variant_info, transcript_strand):
		"""
		Intersect var
		Parameters
		----------
		affected_exon_pos : list of dict of str: str or int
			affected exon genomic positions
		variant_info : VariantInfo
			variant basic info
		transcript_strand : str
			affected transcript strand

		Returns
		-------
		bool
			affected region intersects UniProt domains
		list of str
			protein ids of matched domains
		"""
		### ### ### ### ###
		# First, search the cached matched uniprot domains
		### ### ### ### ###
		matches_cached = False
		if len(self.matched_domains_cache) > 0:
			cached_domains = BedTool(self.matched_domains_cache)
			overlapping_cached_domains = cached_domains.all_hits(
				RefineLossOfFunction.truncated_exons2bed_interval(affected_exon_pos, variant_info, transcript_strand),
				same_strand=True)
			if len(overlapping_cached_domains) > 0:
				matches_cached = True
				return True, [matched_domain.name for matched_domain in overlapping_cached_domains]

		### ### ### ### ###
		# if there is no match with the cached annotations,
		# search on the whole bed file with UniProt domains
		### ### ### ### ###
		if not matches_cached:
			overlapping_domains = self.uniprot_domains.all_hits(
				RefineLossOfFunction.truncated_exons2bed_interval(affected_exon_pos, variant_info, transcript_strand),
				same_strand=True)
			if len(overlapping_domains) > 0:
				# update cached matched annotations
				self.matched_domains_cache = self.matched_domains_cache + overlapping_domains
				return True, [matched_domain.name for matched_domain in overlapping_domains]
			else:
				return False, []

	def intersect_affected_region2critical_regions(self, affected_exons_pos, variant_info, transcript_strand):
		"""
		Intersect  with critical protein region

		Parameters
		----------
		affected_exons_pos : list of dict of str: str or int
			affected exon genomic positions
		variant_info : VariantInfo
			variant basic info
		transcript_strand : str
			transcript strand

		Returns
		-------
		bool
			affected exon is contained in a critical protein region
		list of str
			protein ids of matched critical regions
		"""
		# print("Examine if truncated/altered exon is included in a critical protein region")
		matched_critical_regions = self.critical_prot_regions.all_hits(
			RefineLossOfFunction.truncated_exons2bed_interval(affected_exons_pos, variant_info, transcript_strand),
			same_strand=True)
		if len(matched_critical_regions) > 0:
			return True, [matched_critical_region.name for matched_critical_region in matched_critical_regions]
		else:
			return False, []

	def calculate_prot_len_ptc(self, transcript, transcript_info):
		"""
		Calculate the protein length of the observed coding sequence based on the (premature) termination codon

		Parameters
		----------
		transcript: pyensembl.transcript
			transcript object to be used for examination
		transcript_info : dict of str: str
			transcript information

		Returns
		-------
		"""
		# after variant edit, the termination codon can be found even in the 3' UTR region
		# thus, search for the very first termination codon on the constructed observed coding sequence pluts the 3' UTR
		var_coding_3_utr = transcript_info["var_coding_seq"] + transcript.three_prime_utr_sequence
		_, premature_term_codon_index = self.search_termination_codon(extract_codons(var_coding_3_utr), False)
		if premature_term_codon_index == -1:
			return 0
		else:
			return premature_term_codon_index

	def examine_truncated_exons_significance(self, transcript, transcript_info, truncated_exons_pos, variant_info):
		"""
		Examine truncated/altered exons significance
		check if truncated region overlaps a UniProt domain
		if yes:
			* truncated region is contained in a protein critical region
		if no:
			* truncated exons are overlapping clinical significant exons
			* truncated region modifies more than 10% of the reference protein length
		Based on Abou Tayoun et al., Fig 1

		Parameters
		----------
		transcript: pyensembl.transcript
			transcript object to be used for examination
		transcript_info : dict of str: str
			transcript information
		truncated_exons_pos : list of dict of str: str or int
			truncated exon information
		variant_info : VariantInfo
			variant basic info

		Returns
		-------
		str
			refined PVS1 class
		str
			comment assigning this class
		"""
		self.logger.debug(
			"Examine exon significance for protein critical region, clinical significant exon, change on the length of the protein")
		transcript_strand = transcript.exons[0].to_dict()["strand"]
		self.logger.debug("Affected region: {}".format(
			RefineLossOfFunction.truncated_exons2bed_interval(truncated_exons_pos, variant_info,
			                                                  transcript_strand)))
		is_affected_exon_in_domain, matched_domains_prot_ids = self.intersect_affected_region2uniprot_domains(
			truncated_exons_pos, variant_info, transcript_strand)
		self.logger.debug("affected region intersects uniprot domains: {}".format(is_affected_exon_in_domain))
		is_affected_exon_in_critical_prot_region = False
		if is_affected_exon_in_domain:
			# * exon is contained in a protein critical region
			is_affected_exon_in_critical_prot_region, matched_critical_proteins = self.intersect_affected_region2critical_regions(
				truncated_exons_pos, variant_info, transcript_strand)
			self.logger.debug("affected region intersects critical protein region: {}".format(
				is_affected_exon_in_critical_prot_region))
			if is_affected_exon_in_critical_prot_region:
				self.logger.debug("Assign PVS1_Strong class")
				return "PVS1_Strong", "truncated/altered exon is contained in critical regions of protein ids=[" + ",".join(
					matched_critical_proteins) + "]"
		if not (is_affected_exon_in_domain and is_affected_exon_in_critical_prot_region):
			# * truncated exons are overlapping clinical significant exons or hl clinical significant exons
			# * examine if variant is pathogenic and phenotype relevant and it is exhibited in high AF
			# * truncated region modifies more than 10% of the reference protein length
			is_truncated_exon_in_pheno_transcripts, matched_pheno_exons = self.intersect_exon2phenotype_relevants(
				truncated_exons_pos, variant_info, transcript_strand)
			self.logger.debug("affected exon intersects phenotype relevant clinical significant exons: {}".format(
				is_truncated_exon_in_pheno_transcripts))
			is_truncated_exon_in_significant_exons, matched_significant_exons = self.intersect_affected_region2significant_exons(
				truncated_exons_pos, variant_info, transcript_strand)
			self.logger.debug(
				"affected exon intersects significant exons: {}".format(is_truncated_exon_in_significant_exons))
			if (
					is_truncated_exon_in_pheno_transcripts or is_truncated_exon_in_significant_exons) and is_truncated_exon_in_pheno_transcripts:
				# add comment about clinical significant exon or known pathogenic variant with high AF
				if is_truncated_exon_in_pheno_transcripts and not is_truncated_exon_in_significant_exons:
					clinical_significant_exons_comment = "truncated exon intersects phenotype revelant clinical significant exons [exon ids= {" + ",".join(
						matched_pheno_exons) + "}]"
				else:
					clinical_significant_exons_comment = "LoFs are not frequent in truncated/altered exons [" + ",".join(
						matched_significant_exons) + "]" + "," + "truncated exon is present in phenotype relevant transcripts [exon ids= {" + ",".join(
						matched_pheno_exons) + "}]"

				### ### ###
				# calculate the protein lenth of the observed sequence by two ways:
				# a) variant edit length & b) (premature) termination codon
				# set the min(a,b) as the observed protein length
				### ### ###
				ref_prot_len = len(str(transcript.coding_sequence)) / 3
				prot_len_edit = len(transcript_info["var_coding_seq"]) / 3
				prot_len_ptc = self.calculate_prot_len_ptc(transcript, transcript_info)
				self.logger.debug(
					"ref_prot_len: {}, prot_len_edit: {}, prot_len_ptc: {}".format(ref_prot_len, prot_len_edit,
					                                                               prot_len_ptc))

				if abs(prot_len_ptc - ref_prot_len) / ref_prot_len > 0.1:
					self.logger.debug("Affected region changes more than 10% of the protein")
					self.logger.debug("** Assign PVS1_Strong class **")
					return "PVS1_Strong", clinical_significant_exons_comment + " and truncated region is > 10% of reference protein length"
				else:
					self.logger.debug("Affected region changes less than or equal to 10% of the protein")
					self.logger.debug("** Assign PVS1_Moderate class **")
					return "PVS1_Moderate", clinical_significant_exons_comment + " and truncated region is <= 10% of reference protein length"
			else:
				self.logger.debug("Assign ACMG class type: N/A => False")
				if is_truncated_exon_in_significant_exons:
					significance_comment = "LoFs are not frequent in truncated/altered exon"
				else:
					significance_comment = "LoFs are frequent in truncated/altered exon"
				if is_truncated_exon_in_pheno_transcripts:
					pheno_relevant_comment = "exon present in phenotype-relevant transcript(s)"
				else:
					pheno_relevant_comment = "exon absent from phenotype-relevant transcript(s)"
				return "False", significance_comment + " and " + pheno_relevant_comment

	def refine_splice_site(self, transcripts_info, variant_info):
		"""
		Refine category and strength of splice site (intronic) variant
		# for each annotation of affected transcript:
		# 1) check that splice site is located on +/- 1,2
		# 2) if yes, check for frameshift and NMD
		# 3) check if affected exon is present in biological relevant transcript
		# 4) examine the importance of the affected protein region
		# a) exon is contained in a protein critical region
		# b) is a clinical significant exons"
		# c) if the truncated region modifies more than 10% of the protein length
		# then aggregate assigned classes and assignment comments

		Parameters
		----------
		transcripts_info : list of dict of str : str
			transcripts with splice site variant type
		variant_info : VariantInfo
			variant basic info

		Returns
		-------
		list of str
			assigned PVS1 per input transcript, based on Tayoun et al.
		list of str
			class assignment comment per input transcript
		dict of str : str
			map of transcript id to constructed variant coding sequence
		"""
		self.logger.debug("Refine splice site for each transcript")
		assigned_pvs1_per_transcript = []
		assignment_comment_per_transcript = []
		var_coding_seq_per_trans = {}
		for transcript_info in transcripts_info:
			self.logger.debug("=== New transcript ===")
			self.logger.debug("transcript info: {}".format(transcript_info))
			transcript = self.ensembl_data.transcript_by_id(transcript_info["transcript_id"])

			assigned_class_current_transcript = "False"
			assignment_comment_current_transcript = ""
			if "3_prime_UTR_variant" in transcript_info["type_variant"] and len(transcript_info["type_variant"]) == 1:
				# exclude PVS1 refinement for variant type not mentioned in Tayoun et al.
				assigned_class_current_transcript, assignment_comment_current_transcript = "NA", str(
					transcript.id) + ": " + "PVS1 refinement is not performed for {}".format(
					list(transcript_info["type_variant"])[0].replace("_", " "))
			else:
				# examine splice site position
				# if splice site is on +/- 1,2 positions => perform refinement
				# otherwise return unknown PVS1 class (case not discussed by Tayoun et al. work)
				var_coding = transcript_info["var_coding"]
				split_symbols, intron_offsets, directions2exon = self.parse_variant_intron_pos(var_coding)
				intron_affects_splice_sites = False
				for intron_offset in intron_offsets:
					if intron_offset in [1, 2]:
						intron_affects_splice_sites = True

				if intron_affects_splice_sites:
					self.logger.debug("Perform refinement, as variant affects splice site")
					transcript_info["var_coding_seq"], transcript_info[
						"coding_diff_len"] = self.construct_variant_coding_seq(
						transcript, transcript_info, variant_info)
					var_coding_seq_per_trans[transcript.id] = transcript_info["var_coding_seq"]
					# perform refinement if variant affects splice site
					is_reading_frame_preserved = self.assess_reading_frame_preservation(transcript_info)
					self.logger.debug("Reading frame is preserved: {}".format(is_reading_frame_preserved))
					NMD_occurs, truncated_exons_pos, var_skips_start_codon_exon, var_skips_stop_codon_exon, var_skips_coding_exon = self.assess_NMD(
						transcript_info,
						variant_info)
					self.logger.debug("Predicted to undergo NMD: {}".format(NMD_occurs))
					if not var_skips_coding_exon:
						assigned_class_current_transcript, assignment_comment_current_transcript = "False", str(
							transcript.id) + ": " + "Splice acceptor skips non-coding exon"
					elif (var_skips_start_codon_exon and var_skips_stop_codon_exon) or not var_skips_start_codon_exon:
						# if splice acceptor site does not skip the start codon exon
						# follow Fig1 of PVS1 decision tree (Tayoun et al.)
						if not is_reading_frame_preserved:
							self.logger.debug("truncated exons info: {}".format(truncated_exons_pos))
							if NMD_occurs:  # reading frame disrupted and NMD predicted to occur
								# based on Tayoun et al. check if the transcript is relevant to the phenotype (hearing loss)
								# if self.is_transcript_phenotype_relevant(transcript_info["transcript_id"]):
								affected_overlaps_pheno_rel, overlapping_pheno_exons = self.intersect_exon2phenotype_relevants(
									truncated_exons_pos, variant_info,
									transcript.exons[0].to_dict()["strand"])
								if affected_overlaps_pheno_rel:
									self.logger.debug("Transcript belongs on the known list of phenotype relevant")
									self.logger.debug("** Assign PVS1 class **")
									assigned_class_current_transcript, assignment_comment_current_transcript = "PVS1", str(
										transcript.id) + ": " + "splice_site: predicted to undergo NMD, on phenotype relevant transcript with exon ids=[" + ",".join(
										overlapping_pheno_exons) + "]"
								else:
									self.logger.debug(
										"Transcript does not belong on the known list of phenotype relevant transcripts")
									self.logger.debug("** Assign ACMG class type: N/A => False **")
									assigned_class_current_transcript, assignment_comment_current_transcript = "False", str(
										transcript.id) + ": " + "splice_site: predicted to undergo NMD, on NO phenotype relevant transcript"
							else:  # reading frame disrupted and NMD predicted not to occur
								self.logger.debug(
									"splice_site: NMD is predicted NOT to occur, so examine truncated/altered exons significance")
								assigned_class_current_transcript, truncated_region_comment = self.examine_truncated_exons_significance(
									transcript,
									transcript_info,
									truncated_exons_pos, variant_info)
								assignment_comment_current_transcript = str(
									transcript.id) + ": " + "splice_site: predicted NOT to undergo NMD" + ", " + truncated_region_comment
						else:
							# reading frame is preserved
							self.logger.debug("Reading frame is preserved")
							self.logger.debug("Examine truncated/altered exons significance")
							transcript = self.ensembl_data.transcript_by_id(transcript_info["transcript_id"])
							assigned_class_current_transcript, truncated_region_comment = self.examine_truncated_exons_significance(
								transcript,
								transcript_info,
								truncated_exons_pos, variant_info)
							assignment_comment_current_transcript = str(
								transcript.id) + ": " + "splice_site: frameshift preserved, " + truncated_region_comment
					else:
						# evaluate skipping of start codon exon
						eval_class, eval_class_comment = self.evaluate_skipping_start_codon_exon(
							transcript, transcript_info, variant_info)
						assigned_class_current_transcript, assignment_comment_current_transcript = eval_class, str(
							transcript.id) + ": " + eval_class_comment
				else:
					self.logger.debug("PVS1 refinement is implemented only for splice site on +/- 1,2 positions")
					# PVS1 refinement is implemented only for splice sites located on +/- 1,2 positions
					assigned_class_current_transcript, assignment_comment_current_transcript = "False", str(
						transcript.id) + ": " + "PVS1 refinement is performed for splice sites located in the +/- 1,2 positions"
			# append assigned classes and comments of current transcript to the total classes and comments
			self.logger.debug(
				"class: {},comment: {}".format(assigned_class_current_transcript, assigned_class_current_transcript))
			assigned_pvs1_per_transcript.append(assigned_class_current_transcript)
			assignment_comment_per_transcript.append(assignment_comment_current_transcript)
		self.logger.debug("~~~")
		return assigned_pvs1_per_transcript, assignment_comment_per_transcript, var_coding_seq_per_trans

	def refine_nonsense_frameshift(self, transcripts_info, variant_info):
		"""
		Refine category and strength of non-sense, frameshift or stop lost variant
		# for each annotation of affected transcript:
		# 1) check for NMD
		# 2) check if affected exon is present in biological relevant transcript
		# 3) examine the importance of the affected protein region
		# a) exon is contained in a protein critical region
		# b) is a clinical significant exons"
		# c) if the truncated region modifies more than 10% of the protein length
		# then aggregate assigned classes and assignment comments

		Parameters
		----------
		transcripts_info : list of dict of str : str
			transcripts matching nonsense or frameshift variant type
		variant_info : VariantInfo
			variant basic info

		Returns
		-------
		list of str
			assigned PVS1 based on Tayoun et al.
		list of str
			class assignment comment
		dict of str : str
			dictionary to map transcript id constructed variant coding sequence
		"""
		self.logger.debug("Refine nonsense, frameshift or stop lost for each transcript")
		assigned_pvs1_per_transcript = []
		assignment_comment_per_transcript = []
		var_coding_seq_per_trans = {}
		for transcript_info in transcripts_info:
			self.logger.debug("=== New transcript {} ===".format(transcript_info["transcript_id"]))
			transcript = self.ensembl_data.transcript_by_id(transcript_info["transcript_id"])
			if (self.is_genomic_pos_in_coding_exon(transcript,
			                                       variant_info.genomic_start) or self.is_genomic_pos_in_coding_exon(
				transcript, variant_info.genomic_end)):
				# variant starts or ends on coding region
				transcript_info["var_coding_seq"], transcript_info[
					"coding_diff_len"] = self.construct_variant_coding_seq(
					transcript, transcript_info, variant_info)
				var_coding_seq_per_trans[transcript.id] = transcript_info["var_coding_seq"]
				NMD_occurs, truncated_exons_pos, _, _, _ = self.assess_NMD(transcript_info, variant_info)
				self.logger.debug("Predicted to undergo NMD: {}".format(NMD_occurs))
				self.logger.debug("truncated exons info: {}".format(truncated_exons_pos))
				assigned_class_current_transcript = "False"
				assignment_comment_current_transcript = ""
				if NMD_occurs:  # based on Tayoun et al. check if the transcript is relevant to the phenotype (hearing loss)
					affected_overlaps_pheno, overlapping_pheno_exons = self.intersect_exon2phenotype_relevants(
						truncated_exons_pos, variant_info,
						transcript.exons[0].to_dict()["strand"])
					if affected_overlaps_pheno:
						self.logger.debug("Transcript belongs on the known list of phenotype relevant")
						self.logger.debug("** Assign PVS1 class **")
						assigned_class_current_transcript, assignment_comment_current_transcript = "PVS1", \
						                                                                           transcript_info[
							                                                                           "transcript_id"] + ": " + "nonsense_frameshift: predicted to undergo NMD, on phenotype relevant transcript with exon ids=[" + ",".join(
							                                                                           overlapping_pheno_exons) + "]"
					else:
						self.logger.debug(
							"Transcript does not belong on the known list of phenotype relevant transcripts")
						self.logger.debug("** Assign ACMG class type: N/A => False **")
						assigned_class_current_transcript, assignment_comment_current_transcript = "False", \
						                                                                           transcript_info[
							                                                                           "transcript_id"] + ": " + "nonsense_frameshift: predicted to undergo NMD, on NO phenotype relevant transcript"
				else:
					self.logger.debug("Examine truncated/altered exons significance")
					assigned_class_current_transcript, truncated_region_comment = self.examine_truncated_exons_significance(
						transcript, transcript_info,
						truncated_exons_pos, variant_info)
					assignment_comment_current_transcript = str(
						transcript.id) + ": " + "nonsense_frameshift_stop_lost: predicted NOT to undergo NMD" + ", " + truncated_region_comment
			else:
				# variant starts and ends on non coding region
				self.logger.debug(
					"Variant start and end positions do not affect coding region of transcript: {}".format(
						transcript.id))
				assigned_class_current_transcript, assignment_comment_current_transcript = "False", str(
					transcript.id) + ": nonsense_frameshift_stop_lost: variant start and end positions not intersecting coding region"
			assigned_pvs1_per_transcript.append(assigned_class_current_transcript)
			assignment_comment_per_transcript.append(assignment_comment_current_transcript)
		self.logger.debug("~~~")
		return assigned_pvs1_per_transcript, assignment_comment_per_transcript, var_coding_seq_per_trans

	def compare_closest_with_alternatives(self, transcript, closest_start_chrom):
		"""
		Compare closest start codon with start codons of all other transcripts

		Parameters
		----------
		transcript : pyensembl.transcript
			transcript containing closest start codon
		closest_start_chrom :
			starting/ending (positive/negative strand) chromosome position of closest start codon

		Returns
		-------
		bool
			any alternative transcript uses closest start codon (True), otherwise (False)
		"""
		transcript_id = transcript.id
		if RefineLossOfFunction.is_transcript_in_positive_strand(transcript):
			transcript_strand = "+"
			closest_start_pos = [closest_start_chrom, closest_start_chrom + 1, closest_start_chrom + 2]
		else:
			transcript_strand = "-"
			closest_start_pos = [closest_start_chrom - 2, closest_start_chrom - 1, closest_start_chrom]

		# extract start codon of other transcripts of the same gene
		alternate_transcripts_unique_start_chr_pos = []
		for transcript in transcript.gene.transcripts:
			if transcript.id != transcript_id:  # if transcript is different that the one harbouring the variant
				if transcript.contains_start_codon:  # if start codon is annotated in the transcript
					if transcript.start_codon_positions not in alternate_transcripts_unique_start_chr_pos:
						alternate_transcripts_unique_start_chr_pos.append(transcript.start_codon_positions)

		# check if closest in-frame is used by another transcript
		if alternate_transcripts_unique_start_chr_pos.count(closest_start_pos) > 0:
			return True
		else:
			return False

	def examine_alternative_start_codon(self, var_transcript):
		"""
		Check on the other transcripts of the same gene, if they use an alternative start codon

		Parameters
		----------
		var_transcript : pyensembl.transcript
			transcript object harboring the variant

		Returns
		-------
		bool
			exists alternative start codon position in any alternative transcript (True), otherwise (False)
		"""
		# self.logger.debug("Check if alternate transcripts, of the same gene containing the variant, use alternative start codon positions")
		exists_alternative_start_codon = True
		var_transcript_id = var_transcript.id
		var_start_codon_chr_pos = var_transcript.start_codon_positions
		# self.logger.debug("variant transcript start codon chr positions: {}".format(var_start_codon_chr_pos))
		alternate_transcripts_unique_start_chr_pos = []
		for transcript in var_transcript.gene.transcripts:
			if transcript.id != var_transcript_id:  # if transcript is different that the one harbouring the variant
				if transcript.contains_start_codon:  # if start codon is annotated in the transcript
					if transcript.start_codon_positions not in alternate_transcripts_unique_start_chr_pos:
						alternate_transcripts_unique_start_chr_pos.append(transcript.start_codon_positions)

		# if all alternate transcripts have the same start codon position,
		# use the same chromosomal position for start codon
		# print(alternate_transcripts_unique_start_chr_pos)
		if len(alternate_transcripts_unique_start_chr_pos) == alternate_transcripts_unique_start_chr_pos.count(
				var_start_codon_chr_pos):
			exists_alternative_start_codon = False
		else:
			exists_alternative_start_codon = True
		return exists_alternative_start_codon

	def evaluate_skipping_start_codon_exon(self, transcript, transcript_info, variant_info):
		"""
		Evaluate splice site acceptor that skips start codon exon
		# 1. Check if the closest in-frame start codon exists downstream of 200 bases from lost start codon
		# if there exists such in-frame start codon, then:
		# a) the region [lost start codon, in-frame start codon] overlaps UniProt domain
		# b) any other transcript starts with the same in-frame start codon

		Parameters
		----------
		transcript : pyensembl.transcript
			transcript with the exon skipping event
		transcript_info : dict of str: str
			transcript information
		variant_info : VariantInfo
			basic variant information

		Returns
		-------
		str
			evaluation of variant class
		str
			evaluation comment
		"""
		# print("Evaluate skipping start codon exon")
		# As Tayoun et al. paper does not refine this case, splice site that skips start codon,
		# make class to False, but get useful info for the human curator
		var_class, var_class_comment = "False", ""

		# find skipped exon (containing the start codon)
		exons_containing_var, var_exon_start_offset, var_exon_end_offset = self.find_exon_by_var_pos(transcript,
		                                                                                             transcript_info,
		                                                                                             variant_info,
		                                                                                             True)
		skipped_start_codon_exon_idx = exons_containing_var[0] - 1
		try:
			assert skipped_start_codon_exon_idx >= 0
		except AssertionError:
			self.logger.error("Skipped start codon exon index should not be negative\n=> variant position: {}".format(
				variant_info.to_string()), exc_info=True)

		# find the (original) start codon offset in the skipped exon
		exon_chrom_positions = self.get_transcript_exon_offsets(transcript, True)
		if RefineLossOfFunction.is_transcript_in_positive_strand(transcript):
			transcript_strand = "+"
			start_codon_first_pos = transcript.start_codon_positions[0]
		else:
			transcript_strand = "-"
			start_codon_first_pos = transcript.start_codon_positions[-1]
		start_codon_offset = abs(start_codon_first_pos - exon_chrom_positions[skipped_start_codon_exon_idx][0])

		# find the next start codon (ATG) position
		# get the codons in the downstream (200 bases) from the start codon
		# in the 200 bases the lost start codon do not count, thus subtract 3 more bases
		skipped_coding_length = var_exon_end_offset - start_codon_offset - 3
		downstream_length = 201 - skipped_coding_length
		downstream_length = downstream_length - downstream_length % 3
		start_downstream_codons = self.get_codons_downstream_start(transcript_info, downstream_length)

		# search for start codon in the downstream region of transcript
		closest_start_codon_index = self.search_closest_start_codon(start_downstream_codons, variant_info)
		if closest_start_codon_index < 0:
			# start codon does not exist on the next exons, after the skipped one, and in the downstream region of 200 bases
			var_class_comment = "Transcript does not contain an in-frame start codon in 200 bases downstream of lost start codon"
		else:
			### ### ###
			# check if found closest in-frame start codon
			# is used by other transcript of the same gene
			### ### ###
			exon_chrom_positions = self.get_transcript_exon_offsets(transcript, True)
			# length of coding sequence of skipped exon = var_exon_end_offset - start_codon_offset)
			# starting from the end of the skipped exon, length up to the in-frame start = (closest_start_codon_index) * 3
			closest_start_codon_offset = (var_exon_end_offset - start_codon_offset) + (
				closest_start_codon_index) * 3 + 1
			closest_start_codon_exon_idx, closest_start_codon_exon_offset = self.find_exon_by_ref_pos(transcript,
			                                                                                          closest_start_codon_offset,
			                                                                                          False)

			var_class_comment = "Transcript contains in-frame start codon in {} bases downstream of lost start codon".format(
				closest_start_codon_offset + 1)

			### ### ###
			# check if region (lost start, closest start) contains UniProt domain
			### ### ###
			# create chromosomal positions for skipped exon up to alternative start codon
			# get as start and end the skipped exon positions
			affected_chrom_pos = []
			for affected_exon_idx in range(skipped_start_codon_exon_idx, closest_start_codon_exon_idx + 1):
				if transcript_strand == "+":
					affected_exon_start = exon_chrom_positions[affected_exon_idx][0]
					affected_exon_end = exon_chrom_positions[affected_exon_idx][1]
				else:
					affected_exon_start = exon_chrom_positions[affected_exon_idx][1]
					affected_exon_end = exon_chrom_positions[affected_exon_idx][0]
				if affected_exon_idx == closest_start_codon_exon_idx:
					# replace end position of last affected exon with the found in-frame stop codon
					affected_exon_end = affected_exon_start + closest_start_codon_exon_offset + 1
				affected_chrom_pos.append(
					{"exon_id": transcript.id, "exon_start": affected_exon_start, "exon_end": affected_exon_end})

			# intersect affected region (lost start, closest in-frame) with UniProt domains
			affected_overlaps_domains, overlapping_uniprot_ids = self.intersect_affected_region2uniprot_domains(
				affected_chrom_pos, variant_info, transcript_strand)
			if affected_overlaps_domains:
				var_class_comment = var_class_comment + ", " + "region: (lost start, closest in-frame start codon) overlaps Uniprot domains on protein ids: {}".format(
					",".join(overlapping_uniprot_ids))
			else:
				var_class_comment = var_class_comment + ", " + "region: (lost start, closest in-frame start codon) does not overlap Uniprot domain"

			### ### ###
			# check if any alternative transcripts
			# uses closest start codon
			any_alternative_uses_closest_start = self.compare_closest_with_alternatives(transcript,
			                                                                            affected_chrom_pos[-1][
				                                                                            "exon_end"])
			if any_alternative_uses_closest_start:
				var_class_comment = var_class_comment + ", " + "closest in-frame start is used by another transcript of gene"
			else:
				var_class_comment = var_class_comment + ", " + "closest in-frame start is not used by any transcript of gene"
		return var_class, var_class_comment

	def evaluate_start_codon_var(self, transcript_info, variant_info):
		"""
		Evaluate start codon variant
		# for each annotation of affected transcript:
		# 1) check other transcripts of the same gene contain alternative start codon
		# 2) if no, extract the closest alternative start codon in the downstream of the start
		# 3) then examine if pathogenic variants (ClinVar) exists between alternative start codon and affected start codon

		Parameters
		----------
		transcript_info : dict of str: str
			transcript information
		variant_info : VariantInfo
			variant basic info

		Returns
		-------
		str
			assigned PVS1 rule based on Tayoun et al.
		str
			assignment comment
		str
			constructed variant coding sequence
		"""
		self.logger.debug("Evaluate start codon variant for transcript id: {}".format(transcript_info["transcript_id"]))
		assigned_class, assignment_comment = None, None

		transcript = self.ensembl_data.transcript_by_id(transcript_info["transcript_id"])

		# check if the other transcripts of the same gene contain and alternative start codon
		exists_alternative_start_codon = self.examine_alternative_start_codon(transcript)
		var_coding_seq_trans = ""
		if exists_alternative_start_codon:
			# based on Tayoun et al. then the refinement of the variant is N/A
			# print("Alternative chromosomal position for start codon exists on alternate transcript")
			# print("** Predict ACMG class type: N/A => False **")
			assigned_class = "False"
			assignment_comment = transcript_info[
				                     "transcript_id"] + ": " + "start_loss: exists alternative start codon for alternate transcript"
			return assigned_class, assignment_comment, var_coding_seq_trans
		else:
			# print("No known alternative start codon in other transcripts")
			var_coding_seq_trans, diff_len = self.construct_variant_coding_seq(transcript, transcript_info,
			                                                                   variant_info)
			transcript_info["var_coding_seq"] = var_coding_seq_trans
			transcript_info["coding_diff_len"] = diff_len
			pathogenic_clinvars = None
			# find the next start codon (ATG) position
			# get the codons in the downstream (200 bases) from the start codon
			start_downstream_codons = self.get_codons_downstream_start(transcript_info, 201)
			# print("downstream codons: {}".format(start_downstream_codons))
			# search for start codon in the downstream region of transcript
			closest_start_codon_index = self.search_closest_start_codon(start_downstream_codons, variant_info)
			if closest_start_codon_index < 0:
				# print("no ATG is found in the 201 first bases of the downstream coding sequence of the transcript")
				# print("Try to find ATG on the whole transcript length")
				all_downstream_codons = self.get_codons_downstream_start(transcript_info,
				                                                         len(transcript.coding_sequence))
				# print("all downstream codons: {}".format(all_downstream_codons))
				possible_start_codon_index = self.search_closest_start_codon(all_downstream_codons, variant_info)
				if possible_start_codon_index < 0:
					# print("there is not any alternative start codon for transcript id:{}".format(transcript.id))
					# print("** PVS1 no alternative start **")
					return "False", str(
						transcript.id) + ": " + "start_loss: there is not any alternative start codon for transcript id:{}", var_coding_seq_trans
				else:
					# print(
					# 	"ATG is found in codon={}, coding downstream position={}".format(possible_start_codon_index + 1,
					# 	                                                                 (
					# 			                                                                 possible_start_codon_index + 1) * 3))
					# examine if there is pathogenic mutation in upstream region of alternative start codon
					pathogenic_clinvars = self.extract_pathogenic_var_upstream_closest_start_codon(transcript,
					                                                                               possible_start_codon_index,
					                                                                               variant_info)
			else:
				# print("ATG is found in codon={}, coding downstream position={}".format(closest_start_codon_index + 1,
				#                                                                        (
				# 		                                                                       closest_start_codon_index + 1) * 3))
				### ### ### ### ### ###
				# examine if there is pathogenic mutation
				# in upstream region of alternative start codon
				### ### ### ### ### ###
				pathogenic_clinvars = self.extract_pathogenic_var_upstream_closest_start_codon(transcript,
				                                                                               closest_start_codon_index,
				                                                                               variant_info)

			if pathogenic_clinvars and len(pathogenic_clinvars) > 0:
				### ### ### ### ### ###
				# Assign PVS1 Moderate
				# Between affected start codon and the closest in-frame start codon,
				# there exist pathogenic ClinVars
				### ### ### ### ### ###

				return "PVS1_Moderate", str(
					transcript.id) + ": " + "start_loss: between affected start codon and the closest in-frame start codon, there exist pathogenic ClinVars", var_coding_seq_trans
			else:
				### ### ### ### ### ###
				# Assign PVS1 Supporting
				# Between affected start codon and the closest in-frame start codon,
				# there exists none pathogenic ClinVars")
				### ### ### ### ### ###
				return "PVS1_Supporting", str(
					transcript.id) + ": " + "start_loss: between affected start codon and the closest in-frame start codon, there exists none pathogenic ClinVars", var_coding_seq_trans

	def convert_exon_pos2genomic_pos(self, transcript, overlap_exon_index, overlap_exon_offset):
		"""
		Convert overlap position in exonic coordinates to genomic position

		Parameters
		----------
		transcript : pyensembl.transcript
			transcript in which overlapping position lie in
		overlap_exon_index : int
			exon index of overlapping position
		overlap_exon_offset : int
			overlapping position offset in exon

		Returns
		-------
		int
			genomic position of exon overlapping position
		"""
		# print("Convert exon position to genomic position")
		overlap_exon = transcript.exons[overlap_exon_index]

		if RefineLossOfFunction.is_transcript_in_positive_strand(transcript):
			# for positive strand, add to the exon start
			exon_start = overlap_exon.to_dict()["start"]
			overlap_genomic_pos = exon_start + overlap_exon_offset
		else:
			# for negative strand, subtract from the exon end
			exon_end = overlap_exon.to_dict()["end"]
			overlap_genomic_pos = exon_end - overlap_exon_offset
		return overlap_genomic_pos

	def set_up_start_end_pos(self, transcript, original_start_codon, closest_start_codon):
		"""
		Set up start and end positions that are between the start codon and the closest in-frame start codon

		Parameters
		----------
		transcript : pyensembl.transcript
			transcript object that positions lie in
		original_start_codon : dict of str: str or int
			original start codon information
		closest_start_codon:
			closest inframe start codon information

		Returns
		-------
		list of list of int
			start and end genomic position for all exons between the start codon and closest in-frame start codon
		"""
		start_end = []
		if original_start_codon["exon_idx"] == closest_start_codon["exon_idx"]:
			# original start and closest in-frame start codon are on the same exon
			if RefineLossOfFunction.is_transcript_in_positive_strand(transcript):
				start = original_start_codon["start"]
				end = closest_start_codon["end"]
			else:
				start = closest_start_codon["start"]
				end = original_start_codon["end"]
			start_end.append([start, end])
		elif closest_start_codon["exon_idx"] > original_start_codon["exon_idx"]:
			# closest start codon lies on exon on the downstream of the exon containing the original start codon
			if RefineLossOfFunction.is_transcript_in_positive_strand(transcript):
				transcript_strand = "+"
				start_end.append([original_start_codon["start"],
				                  transcript.exons[original_start_codon["exon_idx"]].to_dict()["end"]])
			else:
				transcript_strand = "-"
				start_end.append([transcript.exons[original_start_codon["exon_idx"]].to_dict()["start"],
				                  original_start_codon["end"]])
			for overlap_exon_idx in range(original_start_codon["exon_idx"] + 1, closest_start_codon["exon_idx"]):
				overlap_exon = transcript.exons[overlap_exon_idx]
				start_end.append([overlap_exon.to_dict()["start"], overlap_exon.to_dict()["end"]])

			if transcript_strand == "+":
				start_end.append([transcript.exons[closest_start_codon["exon_idx"]].to_dict()["start"],
				                  closest_start_codon["start"]])
			else:
				start_end.append(
					[closest_start_codon["end"], transcript.exons[closest_start_codon["exon_idx"]].to_dict()["end"]])
		else:
			try:
				assert 1 == 0
			except AssertionError:
				self.logger.error(
					"AssertError: Closest start codon lies in an exon before the original start codon\n=> original start codon: {}".format(
						original_start_codon), exc_info=True)

		for positions in start_end:
			assert positions[0] <= positions[1], "Start should be lower than end position"
		return start_end

	def extract_pathogenic_var_upstream_closest_start_codon(self, transcript, closest_start_codon_index, variant_info):
		"""
		Parameters
		----------
		transcript : pyensembl.transcript
			transcript object that variant lies in
		closest_start_codon_index : int
			closest start codon index
		variant_info : VariantInfo
			variant basic info

		Returns
		-------
		list of dict of str : int or str or list of str
			pathogenic clinvars between affected start codon and closest in-frame start codon
		"""
		# print("Extract pathogenic variants in upstream of the closest codon")
		is_genomic = False
		# find genomic position for the original start codon
		original_start_codon_exon_index, original_start_codon_exon_offset = self.find_exon_by_ref_pos(transcript, 1,
		                                                                                              is_genomic)

		# based on transcript directionality, add start,end co-ordinates
		if RefineLossOfFunction.is_transcript_in_positive_strand(transcript):
			transcript_strand = "+"
		else:
			transcript_strand = "-"

		original_start_codon_start = transcript.start_codon_positions[0]
		original_start_codon_end = transcript.start_codon_positions[-1]

		original_start_codon = {"exon_idx": original_start_codon_exon_index, "start": original_start_codon_start,
		                        "end": original_start_codon_end}
		# print("original start codon: {}".format(original_start_codon))

		# find exon position of the closest in-frame start codon
		closest_start_codon_pos = (closest_start_codon_index + 1) * 3
		closest_start_codon_exon_index, closest_start_codon_exon_offset = self.find_exon_by_ref_pos(transcript,
		                                                                                            closest_start_codon_pos,
		                                                                                            is_genomic)
		# convert exon position in genomic position for variants retrieval
		closest_start_codon_pos = self.convert_exon_pos2genomic_pos(transcript, closest_start_codon_exon_index,
		                                                            closest_start_codon_exon_offset)
		# print("closet start codon pos: {}".format(closest_start_codon_pos))

		# based on transcript directionality add start,end co-ordinates
		if transcript_strand == "+":
			closest_start_codon_start = closest_start_codon_pos
			closest_start_codon_end = closest_start_codon_pos + 2
		else:
			# closest_start_codon_start = closest_start_codon_pos
			# closest_start_codon_end = transcript.exons[closest_start_codon_exon_index].to_dict()["end"]
			closest_start_codon_start = closest_start_codon_pos - 2
			closest_start_codon_end = closest_start_codon_pos

		closest_start_codon = {"exon_idx": closest_start_codon_exon_index, "start": closest_start_codon_start,
		                       "end": closest_start_codon_end}
		# print("closest start codon: {}".format(closest_start_codon))

		# assert the start and end positions in the variant search ranges
		if transcript_strand == "+":
			try:
				assert original_start_codon["start"] < closest_start_codon["end"]
			except AssertionError:
				self.logger.error(
					"Original and alternative start on the same exon; in positive strand transcript, original's start position should be lower than closest's end position\n=> variant position: {}".format(
						variant_info.to_string()), exc_info=True)
		else:
			try:
				assert closest_start_codon["start"] < original_start_codon["end"]
			except AssertionError:
				self.logger.error(
					"Original and alternative start on the same exon; in negative strand transcript, original's start position should be higher than closest's end position\n=> variant position: {}".format(
						variant_info.to_string()), exc_info=True)
		start_codon2closest_start_ranges = self.set_up_start_end_pos(transcript, original_start_codon,
		                                                             closest_start_codon)
		self.logger.debug("start codon ranges: {}".format(start_codon2closest_start_ranges))
		clinvars = self.extract_clinvars(start_codon2closest_start_ranges, transcript.exons[0].to_dict()["strand"],
		                                 variant_info)
		if clinvars:
			filtered_clinvars = self.quality_filter_clinvars(clinvars, self.min_review_stars)
			return self.extract_pathogenic_clinvars(filtered_clinvars)
		else:
			# if there no extracted clinvars, we can't filter neither select the pathogenic ones
			# thus return None as pathogenic clinvar entries
			return None

	def extract_pathogenic_clinvars(self, clinvars):
		"""
		Extract clinvars with pathogenic or likely pathogenic significance

		Parameters
		----------
		clinvars : list of dict of str : int or str or list of str
			input clinvar entries

		Returns
		-------
		list of dict of str : int or str or list of str
			subset of clinvars with pathogenic or likely pathogenic significance
		"""
		pathogenic_subset = []
		for clinvar in clinvars:
			if len(set(clinvar["CLNSIG"]).difference(set(["Pathogenic"]))) == 0:
				pathogenic_subset.append(clinvar)
		self.logger.debug("Extracted {} pathogenic clinvars: {}".format(len(pathogenic_subset), pathogenic_subset))
		return pathogenic_subset

	def quality_filter_clinvars(self, clinvars, min_review_stars):
		"""
		Filter ClinVar entries by quality

		Parameters
		----------
		clinvars : list of dict of str: int or str or list of str
			input clinvar entries
		min_review_stars : int
			minimum number of review stars

		Returns
		-------
		list of dict of str: int or str or list of str
			filtered clinvar entries
		"""
		filtered_clinvars = []
		for clinvar in clinvars:
			if clinvar['CLNREVSTAT'] >= min_review_stars:
				filtered_clinvars.append(clinvar)
		self.logger.debug("Filtered in {} ClinVar entries".format(len(filtered_clinvars)))
		return filtered_clinvars

	def refine_start_lost(self, transcripts_info, variant_info):
		"""
		Refine category and strength of start loss variant

		Parameters
		----------
		transcripts_info : list of dict of str: str
			transcripts matching start lost variant type
		variant_info : VariantInfo
			variant basic info

		Returns
		-------
		list of str
			assigned PVS1 per input transcript
		list of str
			assignment comment per input transcript
		dict of str : str
			transcript id to constructed variant coding sequence
		"""
		# self.logger.debug("Refine start loss for each transcript")
		assigned_pvs1_per_transcript = []
		assignment_comment_per_transcript = []
		var_coding_seq_per_trans = {}
		for transcript_info in transcripts_info:
			# create the variant coding sequence of the transcript
			assigned_pvs1, assignment_comment, var_coding_seq = self.evaluate_start_codon_var(transcript_info,
			                                                                                  variant_info)
			# todo: check what you add
			self.logger.debug(assigned_pvs1)
			assigned_pvs1_per_transcript.append(assigned_pvs1)
			assignment_comment_per_transcript.append(assignment_comment)
			var_coding_seq_per_trans[transcript_info["transcript_id"]] = var_coding_seq
		# print("~~~")
		return assigned_pvs1_per_transcript, assignment_comment_per_transcript, var_coding_seq_per_trans

	def get_clinvar_strand(self, gene_symbol):
		"""
		Get strand for gene found in clinvar entry

		Parameters
		----------
		gene_symbol : str
			gene symbol found in ClinVar entry

		Returns
		-------
		str
			strand of ClinVar gene
		"""
		# print("Get strand for gene found in ClinVar entry")
		if gene_symbol in self.hugo_genes_df.index:
			return self.hugo_genes_df.loc[gene_symbol].strand
		else:
			return None

	def convert_review_status2stars(self, clinvar_rev_status):
		"""
		Convert CLNREVSTAT (review status) tab from clinvar vcf file to number of review stars
		for unknown description -> star= -1
		ClinVar review status documentation: https://www.ncbi.nlm.nih.gov/clinvar/docs/review_status/

		Parameters
		----------
		clinvar_rev_status : list of str

		Returns
		int
			clinvar gold stars
		"""
		rev_status = [review_elem.replace('_', ' ') for review_elem in clinvar_rev_status]
		rev_status = ",".join(rev_status)
		if rev_status not in self.clinvar_stars_df.Review_status.values:  # if retrieved status not in status
			return self.star_status2int["unknown review status"]
		else:
			return self.star_status2int[
				self.clinvar_stars_df.loc[self.clinvar_stars_df["Review_status"] == rev_status][
					"Number_of_gold_stars"].iloc[0]]

	def normalize_disease_id(self, clinvar_rec):
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
		if "CLNDISDB" in clinvar_rec.INFO:
			if not clinvar_rec.INFO["CLNDISDB"]:
				disease_db_id = ",".join(clinvar_rec.INFO["CLNDISDB"])
			else:
				disease_db_id = None
		else:
			disease_db_id = None
		return disease_db_id

	def normalize_gene_info(self, clinvar_rec):
		"""
		Normalize ClinVar gene information

		Parameters
		----------
		clinvar_rec : vcf.model._Record
			ClinVar record organized in dictionary

		Return :
		list of str, int
			list of genes containing [gene symbol, gene id]
		"""
		if "GENEINFO" in clinvar_rec.INFO:
			genes_info = []
			for gene in clinvar_rec.INFO["GENEINFO"].strip().split("|"):
				genes_info.append(gene.split(":"))
			return genes_info
		else:
			return None

	def extract_clinvars(self, variant_search_ranges, transcript_strand, variant_info):
		"""
		Extract ClinVar entries found in the positions in the range of the annotation

		Parameters
		----------
		variant_search_ranges : list of list of int
			variant search start stop ranges
		transcript_strand : str
			transcript strand
		variant_info : VariantInfo
			variant basic info

		Returns
		-------
		list of dict of str : int or str or list of str
			extracted clivnar entries
		"""
		self.logger.debug("Extract ClinVar variants for input annotation")
		matched_clinvars = []
		for start_stop in variant_search_ranges:
			chr, start, end = variant_info.chrom.split("chr")[1], start_stop[0] - 1, start_stop[1]

			for clinvar_rec in list(self.vcf_reader.fetch(chr, start, end)):
				matched_clinvar = {'id': clinvar_rec.ID, 'pos': int(clinvar_rec.POS), 'ref': str(clinvar_rec.REF),
				                   'alt': ','.join([str(alt) for alt in clinvar_rec.ALT]),
				                   'CLNDISDB': self.normalize_disease_id(clinvar_rec),
				                   'CLNREVSTAT': self.convert_review_status2stars(clinvar_rec.INFO['CLNREVSTAT']),
				                   'CLNSIG': clinvar_rec.INFO['CLNSIG'],
				                   'gene_info': self.normalize_gene_info(clinvar_rec)}
				matched_clinvars.append(matched_clinvar)

		if len(matched_clinvars) > 0:
			self.logger.debug("Found {} matching ClinVar entries".format(len(matched_clinvars)))
			extracted_clinvars = []
			for matched_clinvar in matched_clinvars:
				genes_info = matched_clinvar['gene_info']
				if genes_info:  # get the symbol for the first registered gene in ClinVar
					gene_symbol, gene_id = genes_info[0]
				else:
					gene_symbol = None
				# save a clinvar entry if
				# a) the strand is the same as for the UniProt annotation
				# b) the significance field is filled
				if self.get_clinvar_strand(gene_symbol) == transcript_strand and 'CLNSIG' in clinvar_rec.INFO:
					if 'None' not in matched_clinvar['alt']:
						# do not extract a ClinVar record that contains the None nucleotide as alternate
						extracted_clinvars.append(matched_clinvar)
			if len(extracted_clinvars) == 0:
				# if there is no extracted clinvar entries, then return None
				extracted_clinvars = None
		else:
			extracted_clinvars = None
		return extracted_clinvars
