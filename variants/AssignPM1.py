from genotoscope.utils import extract_all_variant_types, arrange_transcripts_per_variant_type, parse_coding_splicing_column, \
	aggregate_examined_rules
from genotoscope.variants.VariantInfo import VariantInfo
import logging
from os.path import join

from pyensembl import EnsemblRelease
from pybedtools import BedTool
import hgvs.parser


class AssignPM1:
	"""
	Class to examine evidence for PM1 ACMG rule
	PM1, Pathogenic - moderate: "Located in a mutational hot spot and/or critical and well-established functional domain (e.g., active site of an enzyme) without beningn variation"

	Full papers:
	1. General rule:
	Richards, Sue, et al. "Standards and guidelines for the interpretation of sequence variants: a joint consensus
	recommendation of the American College of Medical Genetics and Genomics and the Association for Molecular Pathology."
	Genetics in medicine 17.5 (2015): 405-423.
	DOI: 10.1038/gim.2015.30

	2. Hearing loss specific regions:
	Oza, Andrea M., et al. "Expert specification of the ACMG/AMP variant interpretation guidelines for genetic hearing loss."
	Human mutation 39.11 (2018): 1593-1613.
	DOI: 10.1002/humu.23630
	"""

	def __init__(self, beds_root, critical_prot_regions_no_benign_file, pm1_regions_file):
		"""
		Assign PM1 constructor

		Parameters
		----------
		beds_root : str
			bed files root path
		critical_prot_regions_no_benign_file : str
			filename of critical protein regions without benign mutations
		pm1_regions_file : str
			filename of hearing loss specific protein regions

		Returns
		-------
		None
		"""
		self.logger = logging.getLogger("GenOtoScope_Classify.PM1")
		self.logger.info("Initialize class to examine PM1 rule")

		### Annotation tracks: UniProt domains and repeat regions ###
		self.beds_root = beds_root
		self.critical_regions_no_benign = BedTool(critical_prot_regions_no_benign_file).sort()
		self.pm1_regions = BedTool(pm1_regions_file).sort()

		# release 75 uses human reference genome GRCh37
		self.ensembl_data = EnsemblRelease(75)

		###         HGVS parser         ###
		self.hgvs_parser = hgvs.parser.Parser()

	def run(self, sample_name, chrom, var_start, var_end, gene_name, var_ref, var_obs, variant_type,
	        transcript_splicing_info):
		"""
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
			assigned class based on Richards, Sue, et al. and Oza et al.
		str
			assignment comment
		"""
		# read all variant types
		variant_types = extract_all_variant_types(variant_type)
		transcripts_info = parse_coding_splicing_column(transcript_splicing_info, gene_name, self.hgvs_parser)
		# save all variant info
		variant_info = VariantInfo(sample_name, gene_name, variant_type, transcripts_info, chrom, var_start, var_end,
		                           var_ref,
		                           var_obs)
		self.logger.debug(
			"Examine PM1 for all transcripts affected by variant in: {}".format(variant_info.to_string()))
		# set up variant types applicable for PM1
		pm1_variant_types = set(["missense"])
		transcripts_per_var_type = arrange_transcripts_per_variant_type(variant_types, transcripts_info,
		                                                                pm1_variant_types)
		PM1_class, assignment_comment_PM1 = "PM1: NA", "PM1: not applicable"
		for var_type, transcripts_info in transcripts_per_var_type.items():
			if var_type == "missense":
				# assess BP3 rule
				PM1_class, assignment_comment_PM1 = self.assess_PM1(transcripts_info, variant_info)
		return PM1_class, assignment_comment_PM1

	def assess_PM1(self, transcripts_info, variant_info):
		"""
		Assess PM1 rule for input filtered affected variant-affected transcripts
		PM1: "Located in a mutational hot spot and/or critical and well-established functional domain (e.g., active site of an enzyme) without benign variation"

		Parameters
		----------
		transcripts_info: list of dict of str : str
			variant-affected transcripts to assess PM1 rule for
		variant_info : VariantInfo object
			variant basic info

		Returns
		-------
		str
			PM1 rule assignment per transcript
		str
			comment for PM1 rule assignment per transcript
		"""
		assigned_PM1_per_transcript = []
		assignment_comment_per_transcript = []
		strand2PM1 = {}
		for transcript_info in transcripts_info:
			self.logger.debug("Assess PM1 for transcript id: {}".format(transcript_info["transcript_id"]))
			# 1) assign PM1 per strand of affected transcript
			# (rule is not related to each transcript info, thus strand info separates PM1 outcome only)
			# 2) check if variant intersects a protein critical region without benign variants
			# a) create BedTool interval object for selected variant
			# b) examine for hit between the self.critical_prot_regions_no_benign and the variant interval object
			# 3) check if variant intersects a hearing loss specific protein region
			# a) examine for hit between the self.pm1_regions and the variant interval object
			# b) variant overlaps collagen Gly-X-Y motifs, check if disrupts a Glycine site (Oza et al. DOI: 10.1002/humu.23630)
			transcript = self.ensembl_data.transcript_by_id(transcript_info["transcript_id"])

			# get affected transcript strand
			transcript_strand = transcript.exons[0].to_dict()["strand"]
			if transcript_strand not in strand2PM1:
				# assess PM1 rule for unseen strand
				# create bed variant interval object
				variant_interval = \
					BedTool(variant_info.create_bed_line(transcript_strand, "6columns"), from_string=True)[0]
				pm1_hit_comments = []
				# 2) examine overlap with critical regions without benign
				hits2critical_regions = self.critical_regions_no_benign.all_hits(variant_interval, same_strand=True)

				### ### ###
				# 3) examine overlap with hearing loss specific protein domains
				### ### ###
				hits2pm1_regions = self.pm1_regions.all_hits(variant_interval, same_strand=True)
				var_overlaps_pm1_regions = False
				if len(hits2pm1_regions) > 0:
					# 3b) if variant overlaps collagen regions then check if affects Glysine residue on this region
					hit_in_KCNQ4, hit_disrupts_GlyXY = False, False
					for hit in hits2pm1_regions:
						# if hit on collagen Gly-X-Y domains in COL11A2, COL4A3, COL4A4, COL4A5 genes
						if hit.name in ["P53420", "A0A2R8Y2F0", "A0A2R8Y6E5", "Q01955", "A0A0C4DFS1", "H0YIS1",
						                "Q4VXY6", "P29400", "H0Y9R8"]:
							self.logger.debug("Variant in collagen region")
							if transcript_info["var_protein"]:
								if transcript_info["var_protein"][0:3] == "Gly":
									hit_disrupts_GlyXY = True
									pm1_hit_comments.append("variant disrupts Glysine residue in " + hit.name)
						else:
							# variant overlaps pore-forming intramembrane region of KCNQ4 gene
							hit_in_KCNQ4 = True
							pm1_hit_comments.append("variant overlaps region " + hit.name + " in KCNQ4 gene")
					if hit_in_KCNQ4 or hit_disrupts_GlyXY:
						var_overlaps_pm1_regions = True
					else:
						var_overlaps_pm1_regions = False

				### ### ###
				# assign PM1, based on found overlaps
				### ### ###
				if var_overlaps_pm1_regions or len(hits2critical_regions):
					current_transcript_class = "True"
					if not var_overlaps_pm1_regions:
						# add comment for overlap of critical regions without benign mutation
						pm1_hit_comments.append("variant overlaps critical regions in protein ids= [" + ",".join(
							[hit.name for hit in hits2critical_regions]) + "] without benign mutation")
					current_transcript_comment = ": " + " ".join(pm1_hit_comments)
				else:
					current_transcript_class = "False"
					current_transcript_comment = ": variant does not intersect critical regions without benign, neither hearing loss specified protein domains"
				# save PM1 assignment for all affected transcript having the same transcript strand
				strand2PM1[transcript_strand] = [transcript.id, current_transcript_class, current_transcript_comment]

			### ### ###
			# append assignment of current transcript
			### ### ###
			if transcript.id == strand2PM1[transcript_strand][0]:
				# PM1 was evaluated for this transcript
				assigned_PM1_per_transcript.append(strand2PM1[transcript_strand][1])
				assignment_comment_per_transcript.append(str(transcript.id) + strand2PM1[transcript_strand][2])
			else:
				# PM1 was evaluated for another transcript of the same gene
				assigned_PM1_per_transcript.append(strand2PM1[transcript_strand][1])
				assignment_comment_per_transcript.append(
					str(transcript.id) + ": PM1 was evaluated for transcript id: " + str(
						strand2PM1[transcript_strand][0]))
		return aggregate_examined_rules(assigned_PM1_per_transcript, assignment_comment_per_transcript, "PM1")
