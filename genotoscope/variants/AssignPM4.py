from genotoscope.utils import extract_all_variant_types, arrange_transcripts_per_variant_type, parse_coding_splicing_column, \
	aggregate_examined_rules, extract_codons, get_codon_index
from genotoscope.variants.VariantInfo import VariantInfo

from pybedtools import BedTool
from pyensembl import EnsemblRelease
import hgvs.parser

from os.path import join
import logging


class AssignPM4:
	"""
	Class to examine evidence for PM4 ACMG rule
	PM4, Pathogenic - moderate: "Protein length changes as a result of in-frame deletions/insertions in a nonrepeat region or stop-loss variants"

	Full paper:
	Richards, Sue, et al. "Standards and guidelines for the interpretation of sequence variants: a joint consensus
	recommendation of the American College of Medical Genetics and Genomics and the Association for Molecular Pathology."
	Genetics in medicine 17.5 (2015): 405-423.
	DOI: 10.1038/gim.2015.30
	"""

	def __init__(self, data_path, beds_root, repeats_file):
		"""
		Assign PM4 constructor

		Parameters
		----------
		data_path : os.path object
			data path
		beds_root : os.path object
			beds root path
		repeats_file : BedTool
			uniprot repeats filename

		Returns
		-------
		None
		"""
		self.logger = logging.getLogger("GenOtoScope_Classify.PM4")
		self.logger.info("Initialize class to examine PM4 rule")
		self.data_path = data_path
		# release 75 uses human reference genome GRCh37
		self.ensembl_data = EnsemblRelease(75)

		### Annotation tracks: UniProt repeat regions ###
		self.beds_root = beds_root
		self.repeats = BedTool(repeats_file).sort()

		### Cache of matched annotation ###
		self.matched_annotations_cache = []
		self.hgvs_parser = hgvs.parser.Parser()

	def filter_transcripts_by_info(self, variant_info, info_name):
		"""
		Filter transcripts by their information

		Parameters
		----------
		variant_info : variantInfo
			current variant info object
		info_name : str
			information to filter transcripts for

		Returns
		-------
		list of dict of str: str
			filtered transcripts containing information
		"""
		# self.logger.debug("Filter out transcripts by {}".format(info_name.replace("_"," ")))
		filtered_transcripts_info = []
		for transcript_info in variant_info.transcripts_info:
			transcript = self.ensembl_data.transcript_by_id(transcript_info["transcript_id"])
			if info_name == "stop_codon":
				if transcript.contains_stop_codon:
					filtered_transcripts_info.append(transcript_info)
		return filtered_transcripts_info

	def run(self, sample_name, chrom, var_start, var_end, gene_name, var_ref, var_obs, variant_type,
	        transcript_splicing_info, var_coding_seqs):
		"""
		Process affected-variant transcript information and assess PM4

		Parameters
		----------
		sample_name : str
			sample name containing variant
		chrom : str
			variant-affected chromosome
		var_start : int
			variant start position (genomic)
		var_end : int
			variant end position (genomic)
		gene_name : str
			variant-affected gene
		var_ref : str
			variant reference bases
		var_obs : str
			variant observd bases
		variant_type : str
			variant_type column in GSvar file
		transcript_splicing_info : str
			coding_and_splicing column in GSvar file
		var_coding_seqs : dict of str : str
			constructed variant coding sequence per transcript

		Returns
		-------
		str
			assigned PM4 rule based on Richards, Sue, et al.
		str
			assignment comment
		"""
		self.logger.debug("Examine PM4 for variant with coding info: {}".format(transcript_splicing_info))

		### ### ### ###
		# Preprocess  #
		### ### ### ###
		# read all variant types
		variant_types = extract_all_variant_types(variant_type)
		# parse transcripts information
		transcripts_info = parse_coding_splicing_column(transcript_splicing_info, gene_name, self.hgvs_parser)
		# save variant basic info
		variant_info = VariantInfo(sample_name, gene_name, variant_type, transcripts_info, chrom, var_start, var_end,
		                           var_ref,
		                           var_obs)
		### ### ###
		# filter out transcripts not containing start codon information
		### ### ###
		if variant_type == "stop_lost":
			filtered_transcripts_info = self.filter_transcripts_by_info(variant_info, "stop_codon")
		else:
			filtered_transcripts_info = transcripts_info

		if len(filtered_transcripts_info) == 0:
			PM4_class = "PM4: False"
			assignment_comment_PM4 = "PM4: all affected transcripts did not contain stop codon"
		else:
			pm4_variant_types = set(["inframe_insertion", "inframe_deletion", "conservative_inframe_insertion",
			                         "conservative_inframe_deletion", "disruptive_inframe_deletion",
			                         "disruptive_inframe_insertion", "stop_lost"])
			transcripts_per_var_type = arrange_transcripts_per_variant_type(variant_types, filtered_transcripts_info,
			                                                                pm4_variant_types)
			PM4_class, assignment_comment_PM4 = "PM4: False", "PM4: not applicable"
			for var_type, transcripts_info in transcripts_per_var_type.items():
				if "inframe" in var_type:
					# assess PM4 rule
					PM4_class, assignment_comment_PM4 = self.assess_PM4_inframe_indel(transcripts_info, variant_info,
					                                                                  var_coding_seqs)

				elif var_type == "stop_lost":
					# stop lost variant type
					PM4_class, assignment_comment_PM4 = self.assess_PM4_stop_lost(transcripts_info, variant_info,
					                                                              var_coding_seqs)

		return PM4_class, assignment_comment_PM4

	def calculate_observed_protein_length(self, variant_info, new_stop_distance=0):
		if variant_info.variant_type == "frameshift_deletion":
			return len(variant_info.ref_base) / 3
		elif variant_info.variant_type == "frameshift_insertion":
			return len(variant_info.obs_base) / 3
		else:
			# stop lost
			return new_stop_distance / 3

	def find_inframe_stop(self, codons):
		"""
		Find in-frame stop codon in the 3 prime UTR region

		Parameters
		----------
		codons : list of str
			list of codons in 3' UTR region

		Returns
		-------
		closest_stop_codon : str
			closest in-frame stop codon (TAA or TAG or TGA)
		closet_stop_codon_idx : int
			0-based index of closest in-frame stop codon
		"""
		self.logger.debug("Find in-frame stop codon in 3' UTR region")
		stop_codon_indices = {"TAA": -1, "TAG": -1, "TGA": -1}
		closest_stop_codon_idx, closest_stop_codon = len(codons), None
		for stop_codon in stop_codon_indices:
			if stop_codon in codons:
				stop_codon_indices[stop_codon] = codons.index(stop_codon)
				if stop_codon_indices[stop_codon] < closest_stop_codon_idx:
					closest_stop_codon_idx = stop_codon_indices[stop_codon]
					closest_stop_codon = stop_codon
		return closest_stop_codon_idx, closest_stop_codon

	def assess_PM4_stop_lost(self, transcripts_info, variant_info, var_coding_seqs):
		"""
		Assess PM4 rule for stop lost variants

		Parameters
		----------
		transcripts_info: list of dict of str : str
			variant-affected transcripts to assess PM4 rule for
		variant_info : VariantInfo object
			variant basic info
		var_coding_seqs : dict of str : str
			variant coding sequence per transcript

		Returns
		-------
		str
			PM4 rule assignments
		str
			comment for PM4 rule assignment
		"""
		assigned_PM4_per_transcript = []
		assignment_comment_per_transcript = []
		for transcript_info in transcripts_info:
			# 1) assign PM4 for each affected transcript
			# 2) process stop codon info
			# a) find the next in-frame stop codon in UTR
			# b) examine if the region changes the protein length by more than 10%

			# 2. a) find next in-frame stop codon
			transcript = self.ensembl_data.transcript_by_id(transcript_info["transcript_id"])
			closest_stop_codon_idx, closest_stop_codon = self.find_inframe_stop(transcript.three_prime_utr_sequence)

			if closest_stop_codon:
				# if in-frame stop codon exists in 3' UTR region
				# calculate the protein change ratio
				ref_prot_len = len(transcript.protein_sequence)
				closest_stop_codon_dist = (closest_stop_codon_idx + 1) * 3
				obs_prot_len = self.calculate_observed_protein_length(variant_info,
				                                                      new_stop_distance=closest_stop_codon_dist)
				prot_change_ratio = abs(obs_prot_len - ref_prot_len) / ref_prot_len
				if prot_change_ratio > 0.1:
					# in-frame stop lost increases the protein length more than 10%
					# then, PM4 is triggered
					self.logger.debug("in-frame stop lost increases protein > 10%")
					current_transcript_class = "True"
					current_transcript_comment = str(
						transcript.id) + ": (lost stop, in-frame stop) region results to greater than 10% protein length change"
				else:
					self.logger.debug("in-frame stop lost increases protein length <= 10%")
					# in-frame stop lost increases the protein length less than or equal to 10%
					current_transcript_class = "False"
					current_transcript_comment = str(
						transcript.id) + ": (lost stop, in-frame stop) region results to less than or equal to 10% protein length change"
			else:
				self.logger.debug("3' UTR region does not contain stop codon")
				# 3' UTR region does not contain stop codon
				current_transcript_class = "False"
				current_transcript_comment = str(
					transcript.id) + ": in-frame stop codon does not exist in 3' UTR"
			# append assignment of current transcript
			assigned_PM4_per_transcript.append(current_transcript_class)
			assignment_comment_per_transcript.append(current_transcript_comment)
		return aggregate_examined_rules(assigned_PM4_per_transcript, assignment_comment_per_transcript, "PM4")

	def calculate_prot_len_ptc(self, transcript, var_coding_seq):
		"""
		Calculate the protein length of the observed coding sequence based on the (premature) termination codon

		Parameters
		----------
		transcript: pyensembl.transcript
			transcript object to be used for examination
		var_coding_seq : str
			constructed variant coding sequence

		Returns
		-------
		int
			premature termination codon index
		"""
		self.logger.debug("Calculate protein length by finding termination codon in coding region + 3'UTR region")
		# after variant edit, the termination codon can be found even in the 3' UTR region
		# thus, search for the very first termination codon on the constructed observed coding sequence pluts the 3' UTR
		var_coding_3_utr = var_coding_seq + transcript.three_prime_utr_sequence
		_, premature_term_codon_index = self.search_termination_codon(extract_codons(var_coding_3_utr), False)
		if premature_term_codon_index == -1:
			return 0
		else:
			return premature_term_codon_index

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

	def assess_PM4_inframe_indel(self, transcripts_info, variant_info, var_coding_seqs):
		"""
		Assess PM4 for a inframe insertion/deletion variant type

		Parameters
		----------
		transcripts_info : list of dict of str : str
			variant-affected transcripts to assess PM4 rule for
		variant_info : VariantInfo object
			variant basic info
		var_coding_seqs : dict of str : str
			constructed variant coding sequence per transcript

		Returns
		-------
		str
			PM4 rule assignments
		str
			comment for PM4 rule assignment
		"""
		assigned_PM4_per_transcript = []
		assignment_comment_per_transcript = []

		for transcript_info in transcripts_info:
			# 1) assign PM4 for each affected transcript
			# 2) check the variant edit length on protein product
			# if it is higher than 10%, then:
			# 3) check if variant does not intersect a repeat region
			# a) create BedTool interval object for selected variant
			# b) check if there is a hit between the self.repeats and the variant interval object
			transcript = self.ensembl_data.transcript_by_id(transcript_info["transcript_id"])
			transcript_strand = transcript.exons[0].to_dict()["strand"]

			### ### ###
			# 2. a) calculate observed protein length
			### ### ###
			if transcript_info["transcript_id"] in var_coding_seqs:
				var_coding_seq = var_coding_seqs[transcript_info["transcript_id"]]
			else:
				var_coding_seq = None
			if var_coding_seq:
				self.logger.debug("Variant coding sequence exists")
				obs_prot_len = self.calculate_prot_len_ptc(transcript, var_coding_seq)
			else:
				obs_prot_len = self.calculate_observed_protein_length(variant_info)
			ref_prot_len = len(transcript.protein_sequence)
			prot_change_ratio = abs(obs_prot_len - ref_prot_len) / ref_prot_len
			self.logger.debug("Protein change ratio: {}".format(prot_change_ratio))

			if transcript.protein_sequence and prot_change_ratio > 0.1:
				self.logger.debug("Variant changes protein length in more than 10%")
				# reference protein length is changed more than 10% by variant
				# check if variant intersects repeats
				variant_interval = \
					BedTool(variant_info.create_bed_line(transcript_strand, "6columns"), from_string=True)[0]

				### ### ### ### ###
				# First, search the cached matched repeat annotations
				### ### ### ### ###
				matches_cached = False
				if len(self.matched_annotations_cache) > 0:
					matched_annotations = BedTool(self.matched_annotations_cache)
					overlapping_matched_annotations = matched_annotations.all_hits(variant_interval,
					                                                               same_strand=True)
					if len(overlapping_matched_annotations) > 0:
						current_transcript_class = "False"
						current_transcript_comment = str(
							transcript.id) + ": variant in repeat region, repeats in protein ids=[" + ",".join(
							[overlapping_matched_annotation.name for overlapping_matched_annotation in
							 overlapping_matched_annotations]) + "]"
						matches_cached = True

				### ### ### ### ###
				# if there is no match with the cached annotations,
				# search on the whole bed file
				### ### ### ### ###
				if not matches_cached:
					annotation_hits = self.repeats.all_hits(variant_interval, same_strand=True)
					if len(annotation_hits) > 0:
						current_transcript_class = "False"
						current_transcript_comment = str(
							transcript.id) + ": variant in repeat region, repeats in protein ids=[" + ",".join(
							[annotation_hit.name for annotation_hit in annotation_hits]) + "]"
						# update cached matched annotations
						self.matched_annotations_cache = self.matched_annotations_cache + annotation_hits
					else:
						current_transcript_class = "True"
						current_transcript_comment = str(
							transcript.id) + ": variant not intersecting repeats and results to greater than 10% protein length change"
			else:
				# reference protein length is changed less than or equal to 10% by variant
				self.logger.debug("Variant changes protein length in less than or equal to 10%")
				current_transcript_class = "False"
				current_transcript_comment = str(
					transcript.id) + ": variant results to less than or equal to 10% protein length change"
			# save PM4 assignment for all affected transcript having the same transcript strand
			# strand2PM4[transcript_strand] = [current_transcript_class, current_transcript_comment]
			# append assignment of current transcript
			assigned_PM4_per_transcript.append(current_transcript_class)
			assignment_comment_per_transcript.append(current_transcript_comment)
		return aggregate_examined_rules(assigned_PM4_per_transcript, assignment_comment_per_transcript, "PM4")
