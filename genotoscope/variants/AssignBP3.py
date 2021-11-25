from genotoscope.utils import extract_all_variant_types, arrange_transcripts_per_variant_type, parse_coding_splicing_column, \
	aggregate_examined_rules
from genotoscope.variants.VariantInfo import VariantInfo
from os.path import join
import logging

from pybedtools import BedTool
from pyensembl import EnsemblRelease
import hgvs.parser


class AssignBP3:
	"""
	Class to examine evidence for BP3 ACMG rule
	BP3, Benign - supporting: "In-frame deletions/insertions in a repetitive region without a known function"

	Full paper:
	Richards, Sue, et al. "Standards and guidelines for the interpretation of sequence variants: a joint consensus
	recommendation of the American College of Medical Genetics and Genomics and the Association for Molecular Pathology."
	Genetics in medicine 17.5 (2015): 405-423.
	DOI: 10.1038/gim.2015.30
	"""

	def __init__(self, data_path, beds_root, repeats_no_domains_file):
		"""
		AssignBP3 constructor

		Parameters
		----------
		data_path : str
			data path
		ensembl_data : str
			pyensembl data for human reference genome
		beds_root : str
			annotation bed file root folder
		repeats_no_domains_file : str
			repeats without known function filename

		Returns
		-------
		None
		"""
		self.logger = logging.getLogger("GenOtoScope_Classify.BP3")
		self.logger.info("Initialize class to examine BP3 rule")
		self.data_path = data_path
		# release 75 uses human reference genome GRCh37
		self.ensembl_data = EnsemblRelease(75)

		### Annotation tracks: UniProt domains and repeat regions ###
		self.beds_root = beds_root
		self.repeats_unknown_function = BedTool(repeats_no_domains_file).sort()

		### Cache of matched annotation ###
		self.matched_annotations_cache = []

		###         HGVS parser         ###
		self.hgvs_parser = hgvs.parser.Parser()

	def run(self, sample_name, chrom, var_start, var_end, gene_name, var_ref, var_obs, variant_type,
	        transcript_splicing_info):
		"""
		Process affected-variant transcript information and assess BP3

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

		Returns
		-------
		str
			assigned class based on Richards, Sue, et al.
		str
			assignment comment
		"""
		self.logger.debug("Examine BP3 rule for variant with coding info:{}".format(transcript_splicing_info))
		# read all variant types
		variant_types = extract_all_variant_types(variant_type)
		transcripts_info = parse_coding_splicing_column(transcript_splicing_info, gene_name, self.hgvs_parser)
		# save all variant info
		variant_info = VariantInfo(sample_name, gene_name, variant_type, transcripts_info, chrom, var_start, var_end,
		                           var_ref,
		                           var_obs)

		# arrange affected transcripts per variant type
		# then run BP3 rule examination for all transcripts
		bp3_variant_types = set(
			["inframe_insertion", "inframe_deletion", "conservative_inframe_insertion", "conservative_inframe_deletion",
			 "disruptive_inframe_insertion", "disruptive_inframe_deletion", "stop_gained"])
		transcripts_per_var_type = arrange_transcripts_per_variant_type(variant_types, transcripts_info,
		                                                                bp3_variant_types)
		BP3_class, assignment_comment_BP3 = "BP3: NA", "BP3: not applicable"
		for var_type, transcripts_info in transcripts_per_var_type.items():
			if "inframe" in var_type:
				# assess BP3 rule
				BP3_class, assignment_comment_BP3 = self.assess_BP3(transcripts_info, variant_info)
		return BP3_class, assignment_comment_BP3

	def assess_BP3(self, transcripts_info, variant_info):
		"""
		Assess BP3 rule for input filtered variant-affected transcripts
		BP3: "In-frame deletions/insertions in a repetitive region without a known function"
		Credits:
		http://daler.github.io/pybedtools/3-brief-examples.html#example-2-intersections-for-a-3-way-venn-diagram
		https://daler.github.io/pybedtools/autodocs/pybedtools.bedtool.BedTool.any_hits.html#pybedtools.bedtool.BedTool.any_hits

		Parameters
		----------
		transcripts_info : list of dict of str : str
			variant-affected transcripts to assess BP3 rule for
		variant_info : VariantInfo object
			variant basic info

		Returns
		-------
		str
			BP3 rule assignments
		str
			comment for BP3 rule assignment
		"""
		assigned_BP3_per_transcript = []
		assignment_comment_per_transcript = []
		strand2BP3 = {}
		for transcript_info in transcripts_info:
			# 1) assign BP3 per strand of affected transcript
			# (rule is not related to each transcript info, thus strand info separates BP3 outcome only)
			# 2) check if variant intersects a repeat region and does not intersect a domain region
			# a) create BedTool interval object for selected variant
			# b) check if there is a hit between the self.repeats_unknown_function and the variant interval object
			transcript = self.ensembl_data.transcript_by_id(transcript_info["transcript_id"])

			# get affected transcript strand
			transcript_strand = transcript.exons[0].to_dict()["strand"]
			if transcript_strand not in strand2BP3:
				# assess BP3 rule for unseen strand
				# check if variant intersects repeat but not domain regions
				variant_interval = \
					BedTool(variant_info.create_bed_line(transcript_strand, "6columns"), from_string=True)[0]

				### ### ### ### ###
				# First, search the cached matched annotations
				### ### ### ### ###
				matches_cached = False
				if len(self.matched_annotations_cache) > 0:
					matched_annotations = BedTool(self.matched_annotations_cache)
					overlapping_matched_annotations = matched_annotations.all_hits(variant_interval,
					                                                               same_strand=True)
					if len(overlapping_matched_annotations) > 0:
						current_transcript_class = "True"
						current_transcript_comment = str(
							transcript.id) + ": variant in repeat region but not in any domain, repeats in protein ids=[" + ",".join(
							[overlapping_matched_annotation.name for overlapping_matched_annotation in
							 overlapping_matched_annotations]) + "]"
						matches_cached = True

				### ### ### ### ###
				# if there is no match with the cached annotations,
				# search on the whole bed file
				### ### ### ### ###
				if not matches_cached:
					annotation_hits = self.repeats_unknown_function.all_hits(variant_interval, same_strand=True)
					if len(annotation_hits) > 0:
						current_transcript_class = "True"
						current_transcript_comment = str(
							transcript.id) + ": variant in repeat region but not in any domain, repeats in protein ids=[" + ",".join(
							[annotation_hit.name for annotation_hit in annotation_hits]) + "]"
						# update cached matched annotations
						self.matched_annotations_cache = self.matched_annotations_cache + annotation_hits
					else:
						current_transcript_class = "False"
						current_transcript_comment = str(
							transcript.id) + ": variant not intersecting repeats with unknown function"

				# save BP3 assignment for all affected transcript having the same transcript strand
				strand2BP3[transcript_strand] = [transcript.id, current_transcript_class, current_transcript_comment]

			### ### ###
			# append assignment of current transcript
			### ### ###
			if transcript.id == strand2BP3[transcript_strand][0]:
				# the BP3 was evaluated for this transcript
				assigned_BP3_per_transcript.append(strand2BP3[transcript_strand][1])
				assignment_comment_per_transcript.append(str(transcript.id) + strand2BP3[transcript_strand][2])
			else:
				# the BP3 was evaluated for another transcript of same gene
				assigned_BP3_per_transcript.append(strand2BP3[transcript_strand][1])
				assignment_comment_per_transcript.append(
					str(transcript.id) + ": BP3 was evaluated for transcript id: " + str(
						strand2BP3[transcript_strand][0]))
		return aggregate_examined_rules(assigned_BP3_per_transcript, assignment_comment_per_transcript, "BP3")
