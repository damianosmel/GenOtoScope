from genotoscope.variants.RefineLossOfFunction import RefineLossOfFunction
from genotoscope.variants.VariantInfo import VariantInfo
from genotoscope.utils import extract_all_variant_types, arrange_transcripts_per_variant_type, \
	aggregate_examined_rules, update_assignment, parse_coding_splicing_column, compose_clinvar_change_signature, \
	normalize_gene_info, convert_review_status2stars, get_clinvar_strand, normalize_codon_nucleotides

from os.path import join
from math import ceil
import logging

from Bio.Seq import Seq
from Bio.Data import IUPACData
import vcf
from pyensembl import EnsemblRelease
import hgvs.parser


class ProcessClinVar:
	"""
	Class to process ACMG/AMP ClinVar evidence

	PS1, Pathogenic - strong: Same amino acid change as previously established pathogenic variant regardless of nucleotide change
	PM5, Pathogenic - moderate: Novel missense change at an amino acid residue where a different missense change determined to be pathogenic has been seen before
	(For hearing loss specification only:)
	PM5, Pathogenic - strong: Novel missense change at same amino acid residue, as two different pathogenic missense variants

	Full papers:
	Richards, Sue, et al. "Standards and guidelines for the interpretation of sequence variants: a joint consensus
	recommendation of the American College of Medical Genetics and Genomics and the Association for Molecular Pathology."
	Genetics in medicine 17.5 (2015): 405-423.
	DOI: 10.1038/gim.2015.30

	Oza, Andrea M., et al. "Expert specification of the ACMG/AMP variant interpretation guidelines for genetic hearing loss."
	Human mutation 39.11 (2018): 1593-1613.
	DOI: 10.1002/humu.23630
	"""

	def __init__(self, data_path, clinvar_root, clinvar_file, clinvar_stars_df, min_review_stars,
	             hugo_genes_df):
		"""
		ProcessClinVar constructor

		Parameters
		----------
		data_path : str
			absolute path where genomic variants lie
		clinvar_root : str
			clinvar root path
		clinvar_file : str
			clinvar file
		clinvar_stars_df : pandas.DataFrame
			clinvar review stars dataframe
		min_review_stars : int
			minimum review stars to filter in a ClinVar entry
		ensembl_data : pyensembl.EnsemblRelease
			ensembl data
		hugo_genes_df : pandas.DataFrame
			hugo genes dataframe

		Returns
		-------
		None
		"""
		### logger ###
		self.logger = logging.getLogger("GenOtoScope_Classify.PS1_PM5")
		self.logger.info("Initialize class to examine PS1 and PM5 rules")
		self.data_path = data_path
		self.clinvar_root = clinvar_root

		### PyEnsembl ###
		# release 75 uses human reference genome GRCh37
		self.ensembl_data = EnsemblRelease(75)

		### ClinVar ###
		self.load_clinvar_file(clinvar_file)
		# quality star system
		self.star_status2int = {"none": 0, "one": 1, "two": 2, "three": 3, "four": 4, "unknown review status": -1}
		self.clinvar_stars_df = clinvar_stars_df

		self.min_review_stars = min_review_stars
		self.clinvars_cache = {}

		### Genes to strand information ###
		self.hugo_genes_df = hugo_genes_df

		self.hgvs_parser = hgvs.parser.Parser()

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
		self.logger.debug("Loading ClinVar file from {}".format(join(self.clinvar_root, clinvar_file)))
		self.vcf_reader = vcf.Reader(filename=join(self.clinvar_root, clinvar_file), compressed=True, encoding='utf-8')

	def correct_double_counting(self, ps1_class, pm5_class, pm5_comment):
		"""
		Correct double counting by deactivate PM5_Moderate if PS1 is triggered as well

		Parameters
		----------
		ps1_class : list of str
			assigned PS1 class per examined transcript
		pm5_class: list of str
			assigned PM5 class per examined transcript
		pm5_comment: list of str
			list of comments per PM5 assignment

		Returns
		-------
		list of str
			corrected PM5 assignments
		list of str
			updated comment for PM5 assignments
		"""
		self.logger.debug("Correct double counting if PS1 is triggered")
		if "True" in ps1_class and "PM5_Moderate" in pm5_class:
			### ### ###
			# deactivate PM5_Moderate if PS1 is triggered
			### ### ###
			self.logger.debug("Deactivate PM5_Moderate, because PS1 is triggered")
			pm5_class_corrected, pm5_comment_corrected = [], []
			for idx, assigned_class_transcript in enumerate(pm5_class):
				if assigned_class_transcript == "PM5_Moderate":
					pm5_class_corrected.append("False")
					pm5_comment_corrected.append(pm5_comment[idx] + " (PM5_Moderate deactivated as PS1 triggered)")
				else:
					pm5_class_corrected.append(assigned_class_transcript)
					pm5_comment_corrected.append(pm5_comment[idx])
			return pm5_class_corrected, pm5_comment_corrected
		else:
			return pm5_class, pm5_comment

	def run(self, sample_name, chrom, var_start, var_end, gene_name, var_ref, var_obs, variant_type,
	        transcript_splicing_info):
		"""
		Preprocess variant information to extract variant-affected codon genomic positions
		then assess PS1 and PM5 for the found codon

		Parameters
		----------
		sample_name : str
			sample name
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
			coding and transcript column in GSVar file

		Returns
		-------
		str
			assigned PS1 and PM5 rules based on Richards, Sue, et al.
		str
			assignment comment
		"""
		self.logger.debug("Examine PS1 and PM5 rules for variant transcript info: {}".format(transcript_splicing_info))

		var_chrom = chrom
		var_genomic_start = int(var_start)
		var_genomic_end = int(var_end)
		gene_name = gene_name
		var_ref = var_ref
		var_obs = var_obs

		# read all variant types
		variant_types = extract_all_variant_types(variant_type)

		transcripts_info = parse_coding_splicing_column(transcript_splicing_info, gene_name, self.hgvs_parser)
		# create variant info object
		variant_info = VariantInfo(sample_name, gene_name, variant_type, transcripts_info, var_chrom, var_genomic_start,
		                           var_genomic_end, var_ref, var_obs)
		# filter out transcripts not containing start codon information
		filtered_transcripts_info = self.filter_transcripts_by_info(variant_info, "start_codon")

		# arrange affected transcripts per variant type
		# then run PS1 and PM5 rules examination for all transcripts
		ps1_pm5_variant_types = set(["missense"])
		transcripts_per_var_type = arrange_transcripts_per_variant_type(variant_types, filtered_transcripts_info,
		                                                                ps1_pm5_variant_types)
		aggregated_PS1_PM5_rules, aggregated_PS1_PM5_comment = "PS1: NA||PM5: NA", "PS1 not applicable || PM5 not applicable for this variant type"

		for var_type, transcripts_info in transcripts_per_var_type.items():
			current_var_type_classes, current_var_type_comments = None, None
			if var_type == "missense":
				# collect information for variant's codon
				transcripts_var_codon_info = self.extract_var_codon_info(transcripts_info, variant_info)
				clinvars_per_transcript = self.extract_clinvar_records(transcripts_var_codon_info, variant_info)
				# assess PS1 and PM5
				PS1_class, assignment_comment_PS1 = self.assess_PS1(clinvars_per_transcript,
				                                                    transcripts_var_codon_info)
				PM5_class, assignment_comment_PM5 = self.assess_PM5(clinvars_per_transcript,
				                                                    transcripts_var_codon_info)
				'''
				PM5_class_corrected, assignment_comment_PM5_corrected = self.correct_double_counting(PS1_class,
				                                                                                     PM5_class,
				                                                                                     assignment_comment_PM5)
				'''
				aggregated_PS1_PM5_rules = update_assignment(PS1_class, PM5_class)
				aggregated_PS1_PM5_comment = update_assignment(assignment_comment_PS1, assignment_comment_PM5)

		return aggregated_PS1_PM5_rules, aggregated_PS1_PM5_comment

	def normalize_codon_exonic_pos(self, transcript, variant_info, codon_genomic_positions, transcript_strand):
		"""
		Correct codon position after investigating intersection with an intron

		Parameters
		----------
		transcript : pyensembl.transcript
			transcript ensembl object
		variant_info : VariantInfo
			basic variant information
		codon_genomic_positions : list of int
			initial codon genomic positions
		transcript_strand : str
			strand of transcript

		Returns
		-------
		codon_pos_corrected : list of list of int
			corrected codon positions
		codon_intersects_intron_at : int
			position that codon intersects with an intron (0-index)
		"""
		self.logger.debug("Normalize codon in exonic positions")
		### ### ###
		# create coding sequence ranges respecting transcript's strand direction
		### ### ###
		coding_seq_positions = []
		if transcript_strand == "+":
			transcript_strand = +1
			codon_start, codon_end = codon_genomic_positions[0], codon_genomic_positions[2]
		else:
			transcript_strand = -1
			codon_start, codon_end = codon_genomic_positions[2], codon_genomic_positions[0]
		codon_middle = codon_genomic_positions[1]

		# for coding_range in transcript.coding_sequence_position_ranges:
		for exon in transcript.exons:
			if transcript_strand == +1:
				coding_seq_positions.append([exon.to_dict()["start"], exon.to_dict()["end"]])
			else:
				coding_seq_positions.append([exon.to_dict()["end"], exon.to_dict()["start"]])
		num_coding_sequences = len(coding_seq_positions)

		### ### ### ### ### ### ### ### ### ### ### ### ###
		# find at which exon (or between which two exons) #
		# the codon lies in                               #
		### ### ### ### ### ### ### ### ### ### ### ### ###
		codon_pos_corrected = []
		codon_intersects_intron_at = -1
		for cod_seq_idx, coding_seq_range in enumerate(coding_seq_positions):
			coding_seq_start, coding_seq_end = coding_seq_range[0], coding_seq_range[1]
			normalized_exon_interval = range(coding_seq_start, coding_seq_end + transcript_strand, transcript_strand)

			# self.logger.debug("Coding start: {}, end: {}, range: {}".format(coding_seq_start, coding_seq_end,
			#                                                                 normalized_exon_interval))
			# if codon_start >= coding_seq_start and codon_end <= coding_seq_end:
			if codon_start in normalized_exon_interval and codon_end in normalized_exon_interval:
				self.logger.debug("codon in exon")
				# codon inside coding seq =>
				# ( --- coding_seq --- )
				#      [ -codon- ]
				codon_pos_corrected = codon_genomic_positions
				codon_intersects_intron_at = -1
				break
			elif codon_start in normalized_exon_interval:
				self.logger.debug("start in exon")
				### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
				# codon may intersect with current coding seq range in the coding seq range's right side
				# ( --- current_coding_seq --- ) ( -intron- ) ...
				#                             [ -codon- ]
				### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
				if codon_start == coding_seq_end:
					self.logger.debug("Codon start is the last base of exon")
					if transcript_strand == +1:
						codon_pos_corrected.append([codon_start])
						codon_pos_corrected.append([coding_seq_positions[cod_seq_idx + 1][0],
						                            coding_seq_positions[cod_seq_idx + 1][0] + 1])
					else:
						# codon start is the one before the last base of the exon
						codon_pos_corrected.append(
							[coding_seq_positions[cod_seq_idx + 1][0] - 1, coding_seq_positions[cod_seq_idx + 1][0]])
						codon_pos_corrected.append([codon_start])
					# intersects on the 2nd position (0-index = 1)
					codon_intersects_intron_at = 1
				elif codon_start == coding_seq_end - 1 * transcript_strand:
					self.logger.debug("Codon start is the penultimate base of exon")
					if transcript_strand == +1:
						codon_pos_corrected.append([codon_start, codon_middle])
						codon_pos_corrected.append([coding_seq_positions[cod_seq_idx + 1][0]])
					else:
						codon_pos_corrected.append([coding_seq_positions[cod_seq_idx + 1][0]])
						codon_pos_corrected.append([codon_middle, codon_start])
					# intersects on the 3rd position (0-index = 2)
					codon_intersects_intron_at = 2
				break
			elif codon_middle in normalized_exon_interval and codon_end in normalized_exon_interval:
				self.logger.debug("middle,end in exon")
				### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
				# 	codon middle and end position intersects current exon
				# (--- previous_coding_seq ---) ( -intron- )   (--- current_coding_seq ---)
				#                       start               middle,end
				#                       [-------------------codon-----]
				### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
				if transcript_strand == +1:
					codon_pos_corrected.append([coding_seq_positions[cod_seq_idx - 1][1]])
					codon_pos_corrected.append([codon_middle, codon_end])
				else:
					codon_pos_corrected.append([codon_end, codon_middle])
					codon_pos_corrected.append([coding_seq_positions[cod_seq_idx - 1][1]])
				# genomic position of codon intersects 1st position (0-index = 0)
				codon_intersects_intron_at = 0
				break
			elif codon_end in normalized_exon_interval:
				self.logger.debug("end in exon")
				### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
				# codon end position intersects exon
				# (--- previous_coding_seq ---) ( -intron- )   (--- current_coding_seq ---)
				#                  start,middle                 end
				#                    [--------------------codon----]
				### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
				if transcript_strand == +1:
					codon_pos_corrected.append(
						[coding_seq_positions[cod_seq_idx - 1][1] - 1, coding_seq_positions[cod_seq_idx - 1][1]])
					codon_pos_corrected.append([codon_end])
				else:
					codon_pos_corrected.append([codon_end])
					codon_pos_corrected.append(
						[coding_seq_positions[cod_seq_idx - 1][1], coding_seq_positions[cod_seq_idx - 1][1] + 1])
				# genomic position of codon intersects 1st position (0-index = 0)
				codon_intersects_intron_at = 0
				break
		self.logger.debug("normalized positions: {}".format(codon_pos_corrected))
		self.logger.debug("codon intersects intron at: {}".format(codon_intersects_intron_at))
		if codon_intersects_intron_at == -1:
			try:
				assert len(codon_pos_corrected) == 3
			except AssertionError:
				self.logger.error(
					"AssertionError: Codon does not intersects intron, thus corrected codon positions should be 3 \n=> variant position: {}".format(
						variant_info.to_string()), exc_info=True)
			# assert corrected positions are sorted
			try:
				assert sorted(codon_pos_corrected) == codon_pos_corrected
			except AssertionError:
				self.logger.error(
					"AssertionError: Corrected position of codon are not sorted: {}\n=> variant position: {}".format(
						codon_pos_corrected, variant_info.to_string()), exc_info=True)
		else:
			try:
				assert len(codon_pos_corrected) == 2
			except AssertionError:
				self.logger.error(
					"AssertionError: Codon intersect intron at {}, thus codon positions should two lists \n=> variant position: {}".format(
						codon_intersects_intron_at, variant_info.to_string()), exc_info=True)

			# check that codon positions are increasing
			start_pos = 0
			for codon_positions in codon_pos_corrected:
				for pos in codon_positions:
					try:
						assert pos > start_pos
					except AssertionError:
						self.logger.error(
							"AssertionError: Codon positions are not increasing: {}\n=> variant position: {}".format(
								codon_pos_corrected, variant_info.to_string()), exc_info=True)
					start_pos = pos
		return codon_pos_corrected, codon_intersects_intron_at

	def filter_transcripts_by_info(self, variant_info, info_name):
		"""
		Filter transcripts per info

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
		filtered_transcripts_info = []
		for transcript_info in variant_info.transcripts_info:
			transcript = self.ensembl_data.transcript_by_id(transcript_info["transcript_id"])
			if info_name == "start_codon":
				if transcript.contains_start_codon:
					filtered_transcripts_info.append(transcript_info)
				elif variant_info.chrom == "chrMT":
					# affected variant on mitochondrial gene
					# => coding sequence is the whole transcript sequence
					filtered_transcripts_info.append(transcript_info)
		# else:
		# print("exclude transcript id:{}, as it does not contain start codon".format(
		# 	transcript_info["transcript_id"]))
		return filtered_transcripts_info

	def extract_var_codon_info(self, transcripts_info, variant_info):
		"""
		Extract variant-located codon information

		Parameters
		----------
		transcripts_info: list of dict of str : str
			Selected transcripts to extract var codon information for
		variant_info : VariantInfo
			current variant info object

		Returns
		-------
		dict of str: dict of str : str
			first-level dictionary key: transcript id, value: second-level dictionary
			second-level dictionary: information for codon, where variant is located
		"""
		self.logger.debug("Extract codon information from variant-affected genomic position")
		transcripts_var_codon_info = {}
		for transcript_info in transcripts_info:
			self.logger.debug("=== New transcript id = {} ===".format(transcript_info["transcript_id"]))
			transcript = self.ensembl_data.transcript_by_id(transcript_info["transcript_id"])

			# find at which codon is the variant located
			# based on the variant location in the coding cDNA
			var_start = int(transcript_info['var_coding'].pos.start.base)
			codon_index = var_start // 3
			var_codon_pos = var_start % 3

			try:
				var_start == codon_index * 3 + var_codon_pos
			except AssertionError:
				self.logger.error(
					"AssertionError: codon_index * 3 + var_codon_pos should be equal to codon index\n=> variant position: {}".format(
						variant_info.to_string()), exc_info=True)
			if var_codon_pos == 0:
				var_codon_pos = 3
			self.logger.debug(
				"Var start pos: {} is at codon index: {} and variant is the codon position of {}".format(var_start,
				                                                                                         codon_index,
				                                                                                         var_codon_pos))
			var_codon_genomic_pos, var_codon_coding_pos = [], []
			self.logger.debug("var_codon_pos = {}".format(var_codon_pos))

			# understand variant edit type
			# create genomic position and amino acid per codon
			if var_codon_pos == 3:  # variant is at the third (last) position of a codon
				try:
					var_start >= 3
				except AssertionError:
					self.logger.error(
						"Variant can't be at the last position of codon and not be at least at the coding position 3\n=> variant position: {}".format(
							variant_info.to_string()), exc_info=True)
				var_codon_genomic_pos = [variant_info.genomic_start - 2, variant_info.genomic_start - 1,
				                         variant_info.genomic_start]
				var_codon_coding_pos = [var_start - 2, var_start - 1, var_start]
			elif var_codon_pos == 1:  # variant is at the first position of a codon
				var_codon_genomic_pos = [variant_info.genomic_start, variant_info.genomic_start + 1,
				                         variant_info.genomic_start + 2]
				var_codon_coding_pos = [var_start, var_start + 1, var_start + 2]
			elif var_codon_pos == 2:  # variant is at the second position of a codon
				try:
					var_start >= 2
				except AssertionError:
					self.logger.error(
						"Variant can't be at the second position of a codon\n=> variant position: {}".format(
							variant_info.to_string()), exc_info=True)
				var_codon_genomic_pos = [variant_info.genomic_start - 1, variant_info.genomic_start,
				                         variant_info.genomic_start + 1]
				var_codon_coding_pos = [var_start - 1, var_start, var_start + 1]
			self.logger.debug("var_codon_genomic_pos: {}".format(var_codon_genomic_pos))

			# get the transcript strand
			if RefineLossOfFunction.is_transcript_in_positive_strand(transcript):
				var_strand = "+"
			else:
				var_strand = "-"

			# correct variant codon genomic positions
			var_codon_genomic_pos_corrected, codon_intersects_intron_at = self.normalize_codon_exonic_pos(transcript,
			                                                                                              variant_info,
			                                                                                              var_codon_genomic_pos,
			                                                                                              var_strand)

			### ### ### ### ### ### ### ### ###
			# create reference codon sequence #
			### ### ### ### ### ### ### ### ###
			try:
				# use the coding sequence attribute
				codon_seq_ref = transcript.coding_sequence[var_codon_coding_pos[0] - 1:var_codon_coding_pos[-1]]
			except ValueError:
				# solve value error by using the sequence attribute instead
				# Pyensembl returned ValueError for current transcript with id (seen for mitochondrial genes)
				codon_seq_ref = transcript.sequence[var_codon_coding_pos[0] - 1: var_codon_coding_pos[-1]]
				self.logger.debug("Use sequence attribute instead of coding_sequence for transcript id: {}".format(
					transcript_info["transcript_id"]))
			### ### ### ### ### ### ### ### ###
			# create observed codon sequence  #
			### ### ### ### ### ### ### ### ###
			if var_codon_pos == 3:  # variant on the last position of the codon
				codon_seq_obs = [char for char in codon_seq_ref[0:2]] + [variant_info.obs_base]
			elif var_codon_pos == 1:  # variant on the first position of the codon
				codon_seq_obs = [variant_info.obs_base] + [char for char in codon_seq_ref[1:3]]
			elif var_codon_pos == 2:  # variant on the middle (second) position of the codon
				codon_seq_obs = [codon_seq_ref[0], variant_info.obs_base, codon_seq_ref[2]]
			codon_seq_obs = "".join(codon_seq_obs)

			# assert that the codon size is multiple of 3
			try:
				len(codon_seq_ref) % 3 == 0 and len(codon_seq_obs) % 3 == 0
			except AssertionError:
				self.logger.error(
					"The codon sequence for reference or observed sequence is not multiple of 3\n=> variant position: {}".format(
						variant_info.to_string()), exc_info=True)
			# create amino of protein for affected codon
			prot_var_start = ceil(var_start / 3)
			codon_amino_ref, codon_amino_obs = [], []
			codon_amino_ref.append(str(Seq(codon_seq_ref).translate()))
			codon_amino_obs.append(str(Seq(codon_seq_obs).translate()))
			codon_amino_ref.append(''.join(ProcessClinVar.convert_1to3_aa(codon_amino_ref[0])))
			codon_amino_obs.append(''.join(ProcessClinVar.convert_1to3_aa(codon_amino_obs[0])))

			transcripts_var_codon_info[transcript_info["transcript_id"]] = {
				"var_start": var_start,
				"genomic_pos": var_codon_genomic_pos_corrected,
				"coding_pos": var_codon_coding_pos,
				"intersects_intron_at": codon_intersects_intron_at,
				"strand": var_strand,
				"seq_ref": codon_seq_ref,
				"seq_obs": codon_seq_obs,
				"prot_start": prot_var_start,
				"amino_ref": codon_amino_ref,
				"amino_obs": codon_amino_obs
			}
			self.logger.debug("Var codon info per transcript: {}".format(transcripts_var_codon_info))
		return transcripts_var_codon_info

	def match_pos2cache(self, chr, pos, strand):
		"""
		Match position to ClinVar cache

		Parameters
		----------
		chr: str
			variant chromosome
		pos : int
			variant position
		strand : int
			variant strand

		Returns
		-------
		bool
			variant position found in cache (True), otherwise False
		list of dict of str: str
			cached ClinVars (or None)
		"""
		self.logger.debug("Search pos: chr {}, {}, {} in cache".format(chr, pos, strand))
		cache_hit, cached_clinvars = False, None
		if chr in self.clinvars_cache:
			if pos in self.clinvars_cache[chr]:
				if strand in self.clinvars_cache[chr][pos]:
					cache_hit = True
					cached_clinvars = self.clinvars_cache[chr][pos][strand]
		return cache_hit, cached_clinvars

	def cache_clinvars(self, chr, pos, strand, extracted_clinvars):
		"""
		Save extracted clinvars to cache

		Parameters
		----------
		chr : str
			variant chromosome
		pos : int
			variant position
		strand : str
			variant strand
		extracted_clinvars : list of dict of str: str
			extracted clinvars

		Returns
		-------
		None
		"""
		self.logger.debug("Cache clinvars for pos: chr{},{},{}".format(chr, pos, strand))
		if chr in self.clinvars_cache:
			if pos in self.clinvars_cache[chr]:
				if strand in self.clinvars_cache[chr][pos]:
					self.clinvars_cache[chr][pos][strand] = self.clinvars_cache[chr][pos][strand] + extracted_clinvars
				else:
					self.clinvars_cache[chr][pos][strand] = extracted_clinvars
			else:
				self.clinvars_cache[chr][pos] = {strand: extracted_clinvars}
		else:
			self.clinvars_cache[chr] = {pos: {strand: extracted_clinvars}}

	def add_info2extracted_clinvars(self, clinvar_records):
		"""
		Add information for extracted clinvar records
		Each extracted clinvar record is a dictionary of
		ID -> clinvar id
		CLNDISDB -> clinvar disease db
		CLNREVSTAT -> converted gold stars from clinvar review status (https://www.ncbi.nlm.nih.gov/clinvar/docs/review_status/)
		CLNSIG -> clinvar clinical significance

		Parameters
		----------
		clinvar_records : list of vcf.model._Record
			input clinvar records
		Returns
		-------
		list of dict
			clinvar records with needed information from vcf
		"""
		self.logger.debug("Parse needed information for extracted ClinVar records")
		clinvars_info, uniq_ids = [], []
		for clinvar_rec in clinvar_records:
			if "CLNSIG" not in list(clinvar_rec.INFO.keys()):
				# if significance is not specified, pass
				continue
			# parse disease db
			if "CLNDISDB" in clinvar_rec.INFO:
				if clinvar_rec.INFO["CLNDISDB"][0]:
					rec_dis_db = ','.join(clinvar_rec.INFO['CLNDISDB'])
				else:
					rec_dis_db = "not_specified"
			else:
				rec_dis_db = "not_specified"
			clinvar_dict = {'id': clinvar_rec.ID, 'pos': int(clinvar_rec.POS), 'ref': str(clinvar_rec.REF),
			                'alt': ','.join([str(alt) for alt in clinvar_rec.ALT]),
			                'CLNDISDB': rec_dis_db,
			                'CLNREVSTAT': convert_review_status2stars(self.clinvar_stars_df, self.star_status2int,
			                                                          clinvar_rec.INFO['CLNREVSTAT']),
			                'CLNSIG': clinvar_rec.INFO['CLNSIG'], 'gene_info': normalize_gene_info(clinvar_rec)}
			if 'None' not in clinvar_dict['alt'] and clinvar_dict['id'] not in uniq_ids:
				# extract all unique clinvars that do not contain the None nucl as alternate
				clinvars_info.append(clinvar_dict)
				uniq_ids.append(clinvar_dict['id'])
		return clinvars_info

	def extract_clinvar_records(self, transcripts_var_codon_info, variant_info):
		"""
		Extract ClinVar records matching codon genomic positions per transcript
		Each extracted clinvar record is a dictionary of
		ID -> clinvar id
		CLNDISDB -> clinvar disease db
		CLNREVSTAT -> converted gold stars from clinvar review status (https://www.ncbi.nlm.nih.gov/clinvar/docs/review_status/)
		CLNSIG -> clinvar clinical significance

		Parameters
		----------
		transcripts_var_codon_info : dict of str : dict of str : str
			two level dictionary, map each transcript id to codon information
		variant_info : VariantInfo
			current variant info object

		Returns
		-------
		dict of str : list of dict of str: str
			dictionary mapping transcript id to all clinvar record found for the specific codon of this transcript
		"""
		self.logger.debug("Extract ClinVar records that are found at the same codon as input variant")
		codon_clinvars = {}
		for transcript_id, transcript_codon_info in transcripts_var_codon_info.items():
			self.logger.debug(" === Transcript {} === ".format(transcript_id))
			chr = variant_info.chrom.split("chr")[1]
			transcript = self.ensembl_data.transcript_by_id(transcript_id)
			transcript_strand = transcript.exons[0].to_dict()["strand"]
			### ### ### ###
			# for each genomic position,
			# first search the cached clinvars
			# if no hit, then perform vcf fetch
			# vcf fetch uses 0-index and open upper limit (stop position)
			### ### ### ###
			if transcript_codon_info["intersects_intron_at"] == -1:
				self.logger.debug("Codon does not intersect intron")
				# codon does not intersect intron
				start, end = transcript_codon_info["genomic_pos"][0] - 1, transcript_codon_info["genomic_pos"][2]
				matched_clinvars = []
				for pos in range(start, end):
					case_hit, cached_clinvars = self.match_pos2cache(chr, pos, transcript_strand)
					if not case_hit:
						extracted_clinvars = self.add_info2extracted_clinvars(
							list(self.vcf_reader.fetch(chr, pos, pos + 1)))
						self.logger.debug("Extracted clinvars: {}".format(extracted_clinvars))
						if len(extracted_clinvars) > 0:
							# cache: save strand and quality filter ClinVars
							# uniq_clinvars = self.uniq_id_filter(extracted_clinvars)
							strand_filtered_clinvars = self.strand_filter(extracted_clinvars, transcript_strand)
							quality_filtered_clinvars = self.quality_filter(strand_filtered_clinvars,
							                                                self.min_review_stars)
							unique_clinvars = self.uniq_new_filter(quality_filtered_clinvars)
							self.cache_clinvars(chr, pos, transcript_strand, unique_clinvars)
							matched_clinvars = matched_clinvars + unique_clinvars
						else:
							# cache: no ClinVar matching this chr,pos,strand
							self.cache_clinvars(chr, pos, transcript_strand, [])
							matched_clinvars = matched_clinvars + []
					else:
						# use cached ClinVar entries
						try:
							assert len(cached_clinvars) >= 0
						except AssertionError:
							self.logger.error(
								"Number of extracted clinvars should be >= 0\n=> variant position: {}".format(
									variant_info.to_string()), exc_info=True)
						matched_clinvars = matched_clinvars + cached_clinvars
			else:
				### ### ###
				# aggregate clinvar entries for the 'separated' genomic positions of the codon
				### ### ###
				self.logger.debug("Variant on two different exons")
				try:
					len(transcript_codon_info["genomic_pos"]) == 2
				except AssertionError:
					self.logger.error(
						"Corrected codon genomic position should be contained into two lists\n=> variant position: {}".format(
							variant_info.to_string()), exc_info=True)
				matched_clinvars = []
				for codon_genomic_range in transcript_codon_info["genomic_pos"]:
					if len(codon_genomic_range) == 2:  # two positions are contained in the current codon partition
						start, end = codon_genomic_range[0] - 1, codon_genomic_range[1]
					else:  # one position is contained in the current codon partition
						start, end = codon_genomic_range[0] - 1, codon_genomic_range[0]
					for pos in range(start, end):
						case_hit, cached_clinvars = self.match_pos2cache(chr, pos, transcript_strand)
						if not case_hit:
							extracted_clinvars = self.add_info2extracted_clinvars(
								list(self.vcf_reader.fetch(chr, pos, pos + 1)))
							self.logger.debug("Extracted clinvars: {}".format(extracted_clinvars))
							if len(extracted_clinvars) > 0:
								# cache: save strand and quality filter ClinVars
								strand_filtered_clinvars = self.strand_filter(extracted_clinvars, transcript_strand)
								quality_filtered_clinvars = self.quality_filter(strand_filtered_clinvars,
								                                                self.min_review_stars)
								unique_clinvars = self.uniq_new_filter(quality_filtered_clinvars)
								self.cache_clinvars(chr, pos, transcript_strand, unique_clinvars)
								matched_clinvars = matched_clinvars + unique_clinvars
							else:
								# cache: no ClinVar matching this chr,pos,strand
								self.cache_clinvars(chr, pos, transcript_strand, [])
								matched_clinvars = matched_clinvars + []
						else:
							# use cached ClinVars
							try:
								assert len(cached_clinvars) >= 0
							except AssertionError:
								self.logger.error(
									"Number of extracted clinvars should be >= 0\n=> variant position: {}".format(
										variant_info.to_string()), exc_info=True)
							matched_clinvars = matched_clinvars + cached_clinvars
			if len(matched_clinvars) > 0:
				codon_clinvars[transcript_id] = self.uniq_new_filter(matched_clinvars)
		return codon_clinvars

	def strand_filter(self, matched_clinvars, transcript_strand):
		"""
		Filter matched ClinVars to ensure the same strand with affected transcript

		Parameters
		----------
		matched_clinvars: list of dict of str: str
			matched clinvar records (for chrom,pos,strand)
		transcript_strand : str
			affected transcript strand

		Returns
		-------
		list of dict of str: str
			filtered clinvars by strand information
		"""
		self.logger.debug("\nFilter by strand the {} matching clinvars".format(len(matched_clinvars)))
		filtered_clinvars = []
		for candidate_clinvar in matched_clinvars:
			if candidate_clinvar["gene_info"]:
				# get the symbol for the first registered gene in ClinVar
				gene_symbol, gene_id = candidate_clinvar["gene_info"][0]
			else:
				gene_symbol = None
			self.logger.debug("clinvar's gene: {}".format(gene_symbol))

			### ### ###
			# save a clinvar entry if
			# a) the strand is the same as for the PyEnsembl annotation
			# b) the significance field is filled
			### ### ###
			if get_clinvar_strand(self.hugo_genes_df,
			                      gene_symbol) == transcript_strand and 'CLNSIG' in candidate_clinvar:
				if 'None' not in candidate_clinvar['alt']:
					# do not extract a clinvar record that contains the None nucl as alternate
					filtered_clinvars.append(candidate_clinvar)
		self.logger.debug("Strand-filtered clinvars: {}".format(filtered_clinvars))
		return filtered_clinvars

	def uniq_new_filter(self, new_clinvars):
		"""
		Filter clinvar records to keep the ones not already found in matched list

		Parameters
		----------
		new_clinvars : list of dict of str : str
			newly matched clinvars

		Returns
		-------
		list of dict of str : str
			unique clinvars
		"""
		self.logger.debug("\nFilter all new clinvars to keep the unique ones only")
		ids, uniq_clinvars = [], []
		for clinvar in new_clinvars:
			if clinvar["id"] not in ids:
				uniq_clinvars.append(clinvar)
				ids.append(clinvar["id"])
		self.logger.debug("Unique clinvars: {}".format(uniq_clinvars))
		return uniq_clinvars

	def quality_filter(self, matched_clinvars, min_quality_stars):
		"""
		Filter clinvar records and keep the ones that have at least the minimum quality of stars

		Parameters
		----------
		matched_clinvars : list of dict of str: str
			matched clinvar records (for chrom,pos,strand)
		min_quality_stars : int
			minimum number of quality stars needed to filter in a clinical assertion

		Returns
		-------
		list of dict of str: str
			filtered clinvars by quality per transcript id
		"""
		self.logger.debug("\nFilter ClinVar entries by quality")
		filtered_clinvars = []
		if len(matched_clinvars) > 0:
			for clinvar in matched_clinvars:
				if clinvar['CLNREVSTAT'] >= min_quality_stars:
					filtered_clinvars.append(clinvar)
		self.logger.debug("Review-stars filtered clinvars: {}".format(filtered_clinvars))
		return filtered_clinvars

	def compact_clinvar_entries(self, clinvar_records):
		"""
		Compact ClinVar entries

		Parameters
		----------
		clinvar_records : list of dict
			list of clinvar records

		Returns
		-------
		str
			compacted (all available) clinvar records for all transcript
		"""
		self.logger.debug("Compact clinvar entries")
		clinvar_attributes_ordered = ['id', 'pos',
		                              'ref', 'alt',
		                              'CLNDISDB', 'CLNREVSTAT', 'CLNSIG']

		compacted_clinvars = []
		for clinvar in clinvar_records:
			clinvar_str = []
			for attr in clinvar_attributes_ordered:
				if attr != "CLNSIG":
					clinvar_str.append(attr + ":" + str(clinvar[attr]))
				else:
					clinvar_str.append(
						attr + ":" + '::'.join([str(signif).replace("_", ' ') for signif in clinvar[attr]]))
			# save compacted clinvar and append to all compacts
			compacted_clinvar = ",".join(clinvar_str)
			self.logger.debug("compacted clinvar: {}".format(compacted_clinvar))
			compacted_clinvars.append(compacted_clinvar)
		self.logger.debug("compacted clinvars: {}".format(compacted_clinvars))
		return ";".join(compacted_clinvars)

	@staticmethod
	def convert_1to3_aa(amino_acids):
		"""
		Convert 1 letter amino acid genotoscope to 3 letter equivalent
		For termination codon the return protein change will contain the 3-letter protein letter 'Ter'
		as discussed: https://varnomen.hgvs.org/recommendations/protein/variant/frameshift/

		Parameters
		----------
		amino_acids : list of str
			list of 1-letter amino acids

		Returns
		-------
		list of str
			list of 3-letters amino acids
		"""
		aa_3code = []
		for aa in amino_acids:
			try:
				letter3 = IUPACData.protein_letters_1to3[aa]
			except KeyError:
				letter3 = 'Ter'
			aa_3code.append(letter3)
		return aa_3code

	def construct_clinvar_prot_change(self, clinvar_rec, var_codon_info):
		"""
		Construct protein change for clinvar record

		For termination codon the return protein change will contain the 3-letter protein letter 'Ter'
		as discussed: https://varnomen.hgvs.org/recommendations/protein/variant/frameshift/

		Parameters
		----------
		clinvar_rec: dict of str : str
			clinvar record information
		var_codon_info: dict of str: str
			variant codon information

		Returns
		-------
		int
			protein change position
		str
			clinvar record's protein change
		"""
		self.logger.debug("Construct protein change for clinvar record: {}".format(clinvar_rec))
		# get position of clinvar nucleotide
		var_pos_idx = var_codon_info["genomic_pos"].index(clinvar_rec["pos"])
		self.logger.debug("ClinVar record is found in codon position: {}".format(var_pos_idx))
		### ### ###
		# construct the reference coding sequence for clinvar record
		### ### ###
		clinvar_ref_seq, clinvar_alt_seq = [], []
		for idx, nucl in enumerate(var_codon_info["seq_ref"]):
			if idx == var_pos_idx:
				clinvar_ref_seq.append(clinvar_rec["ref"])
			else:
				clinvar_ref_seq.append(nucl)
		self.logger.debug("clinvar ref: {}".format(clinvar_ref_seq))
		### ### ###
		# and then construct the alternate sequence
		### ### ###
		for idx, nucl in enumerate(var_codon_info["seq_ref"]):
			if idx == var_pos_idx:
				clinvar_alt_seq.append(clinvar_rec["alt"])
			else:
				clinvar_alt_seq.append(nucl)
		self.logger.debug("clinvar alt: {}".format(clinvar_alt_seq))

		### ### ###
		# translate these ref and alt sequence to protein edit
		### ### ###
		clinvar_ref_translated = str(Seq(normalize_codon_nucleotides(clinvar_ref_seq)).translate())
		clinvar_ref_aa = ''.join(ProcessClinVar.convert_1to3_aa(clinvar_ref_translated))
		clinvar_alt_translated = str(Seq(normalize_codon_nucleotides(clinvar_alt_seq)).translate())
		clinvar_alt_aa = ''.join(ProcessClinVar.convert_1to3_aa(clinvar_alt_translated))
		self.logger.debug(
			"ref_translated: {}, alt_translated: {}".format(clinvar_ref_translated, clinvar_alt_translated))
		return var_codon_info['prot_start'], clinvar_ref_aa + str(var_codon_info['prot_start']) + clinvar_alt_aa

	def assess_PS1(self, codon_clinvars, transcripts_var_codon_info):
		"""
		Assess PS1 rule, that is:
		PS1 Strong: Same amino acid change as previously established pathogenic variant regardless of nucleotide change

		Parameters
		----------
		codon_clinvars : dict of str : dict of str : str
			two-level dictionary mapping transcript id to clinvar entries of variant codon on that transcript

		transcripts_var_codon_info : dict of str : dict of str : str
			two level dictionary, map each transcript id to codon information

		Returns
		-------
		str
			PS1 rule assignments
		str
			comment for PS1 rule assignment
		"""
		self.logger.debug("Assess PS1 rule")
		### ### ###
		# for each transcript, examine each filtered clinvar record for PS1 rule
		# then aggregate the assigned classes and comments
		### ### ###
		assigned_class_per_trans, assignment_comment_per_trans = [], []
		for transcript_id, var_codon_info in transcripts_var_codon_info.items():
			current_transcript_class, current_transcript_comment = None, None
			is_PS1_triggered = False
			supporting_clinvar_data = {}
			self.logger.debug(" === Transcript id = {} === ".format(transcript_id))
			if transcript_id in codon_clinvars:
				var_protein_change = var_codon_info['amino_ref'][1] + str(var_codon_info['prot_start']) + \
				                     var_codon_info['amino_obs'][1]
				self.logger.debug("variant protein change: {}".format(var_protein_change))
				self.logger.debug("affected codon positions: {}".format(var_codon_info["genomic_pos"]))
				'''
				process each clinvar record to get its protein change
				if the protein change is the same as the variant annotation and has pathogenic significance,
				then PS1 is triggered
				'''
				self.logger.debug("Filtered clinvars: {}".format(codon_clinvars[transcript_id]))
				for clinvar_rec in codon_clinvars[transcript_id]:
					is_clinvar_pathogenic = False
					# if clinvar genomic position is contained in the variant codon (genomic position)
					# => create clinvar protein change and compare it to the variant prot change
					if clinvar_rec["pos"] in var_codon_info["genomic_pos"]:
						clinvar_prot_change_pos, clinvar_prot_change = self.construct_clinvar_prot_change(clinvar_rec,
						                                                                                  var_codon_info)
						if clinvar_prot_change == var_protein_change:
							self.logger.debug("Clinvar record has same protein change, now check for significance")
							if len(set(clinvar_rec["CLNSIG"]).difference(
									set(["Pathogenic", "Likely_pathogenic", "Pathogenic/Likely_pathogenic"]))) == 0:
								# if only pathogenic significance is recorded for this clinvar
								# flag the clinvar
								is_clinvar_pathogenic = True

							if is_clinvar_pathogenic:
								# single existence of pathogenic CliVar with same protein change, fires PS1
								self.logger.debug("Pathogenic ClinVar record with same protein change found")
								is_PS1_triggered = True
								self.logger.debug("** PS1 is triggered **")
								# keep all supporting clinvar entries
								if transcript_id not in supporting_clinvar_data:  # supporting data for PS1 rule has not yet processed for transcript id
									supporting_clinvar_data[transcript_id] = [clinvar_rec]
								else:
									supporting_clinvar_data[transcript_id].append(clinvar_rec)
				if is_PS1_triggered:
					current_transcript_class = "True"
					current_transcript_comment = transcript_id + ": " + "Supporting pathogenic ClinVar entries: " + self.compact_clinvar_entries(
						supporting_clinvar_data[transcript_id])
				else:
					current_transcript_class = "False"
					current_transcript_comment = transcript_id + ": " + "No evidence to support PS1 rule"
			else:
				# no ClinVar entries for current transcript
				current_transcript_class = "False"
				current_transcript_comment = transcript_id + ": " + "No matching ClinVar to support PS1 rule"
			# append assignment of current transcript
			assigned_class_per_trans.append(current_transcript_class)
			assignment_comment_per_trans.append(current_transcript_comment)
		return aggregate_examined_rules(assigned_class_per_trans, assignment_comment_per_trans, "PS1")

	def assess_PM5(self, codon_clinvars, transcripts_var_codon_info):
		"""
		Assess PM5 rule, that is:
		PM5, Moderate: Novel missense change at an amino acid residue where a different missense change determined
		to be pathogenic has been seen before

		Based on Oza et al. (DOI: 10.1002/humu.23630)
		PM5, Strong: Located at an amino acid residue with known pathogenic variation
		(at least 2 other variants at the same site meet pathogenic criteria based on independent data)

		Parameters
		----------
		codon_clinvars : dict of str : dict of str : str
			two-level dictionary mapping transcript id to clinvar entries of variant codon on that transcript

		transcripts_var_codon_info : dict of str : dict of str : str
			two level dictionary, map each transcript id to codon information
		Returns
		-------
		bool
			is PM5 rule triggered
		str
			comment for PM5 rule assignment
		"""
		self.logger.debug("Assess PM5 rule")
		'''
		for each transcript, examine each filtered clinvar record for PM5 rule
		then aggregate the assigned class and assignment comments for all transcripts
		'''
		assigned_class_per_trans, assignment_comment_per_trans = [], []
		for transcript_id, var_codon_info in transcripts_var_codon_info.items():
			self.logger.debug(" === Transcript id = {} === ".format(transcript_id))
			current_transcript_class, current_transcript_comment = None, None
			is_pm5_triggered, is_pm5_strong_triggered = False, False
			supporting_clinvar_data = {}
			if transcript_id in codon_clinvars:
				var_protein_change_pos = var_codon_info['prot_start']
				self.logger.debug("variant protein change_pos: {}".format(var_protein_change_pos))
				### ### ###
				# for each clinvar record, construct the protein change and extract the change position
				# if protein change position is the same as the variant annotation and has pathogenic significance,
				# then PM5 is triggered
				### ### ###
				for clinvar_rec in codon_clinvars[transcript_id]:
					is_clinvar_pathogenic = False
					# if clinvar genomic position is contained in the variant codon (genomic position)
					# => create clivnar protein change and compare it to the variant prot change
					if clinvar_rec["pos"] in var_codon_info["genomic_pos"]:
						clinvar_prot_change_pos, clinvar_prot_change = self.construct_clinvar_prot_change(clinvar_rec,
						                                                                                  var_codon_info)

						if clinvar_prot_change_pos == var_protein_change_pos:
							self.logger.debug(
								"Clinvar record has same position for protein change, now check for significance")
							if len(set(clinvar_rec["CLNSIG"]).difference(
									set(['Pathogenic', 'Likely_pathogenic', 'Pathogenic/Likely_pathogenic']))) == 0:
								# if pathogenic significance is recorded for this clinvar
								# flag the entry
								is_clinvar_pathogenic = True

					if is_clinvar_pathogenic:
						# keep all supporting clinvar entries
						if transcript_id not in supporting_clinvar_data:
							# supporting data for PM5 rule has not yet processed for transcript id
							supporting_clinvar_data[transcript_id] = [clinvar_rec]
						else:
							supporting_clinvar_data[transcript_id].append(clinvar_rec)
				### ### ###
				# examine PM5 strength
				### ### ###
				self.logger.debug("Supporting Clinvar data: {}".format(supporting_clinvar_data))
				if len(supporting_clinvar_data.keys()):
					uniq_pathogenic_dna_signatures = set([compose_clinvar_change_signature(supporting_clinvar["ref"],
					                                                                       supporting_clinvar["pos"],
					                                                                       supporting_clinvar["alt"])
					                                      for
					                                      supporting_clinvar in supporting_clinvar_data[transcript_id]])
				else:
					self.logger.debug("No supporting Clinvars for PM5")
					uniq_pathogenic_dna_signatures = None

				if not uniq_pathogenic_dna_signatures:
					# no unique pathogenic dna signatures => no supporting Clinvars for PM5
					# => PM5 is not triggered
					current_transcript_class = "False"
					current_transcript_comment = transcript_id + ": " + "No evidence to support PM5 rule"
				elif len(uniq_pathogenic_dna_signatures) >= 2:
					current_transcript_class = "PM5_Strong"
					current_transcript_comment = transcript_id + ": " + "Supporting pathogenic ClinVar entries: " + self.compact_clinvar_entries(
						supporting_clinvar_data[transcript_id])
				elif len(uniq_pathogenic_dna_signatures) == 1:
					current_transcript_class = "PM5_Moderate"
					current_transcript_comment = transcript_id + ": " + "Supporting pathogenic ClinVar entries: " + self.compact_clinvar_entries(
						supporting_clinvar_data[transcript_id])

			else:
				# no ClinVar entries to support rule
				current_transcript_class = "False"
				current_transcript_comment = transcript_id + ": " + "No matching ClinVar to support PM5 rule"
			# append assignment of current transcript
			assigned_class_per_trans.append(current_transcript_class)
			assignment_comment_per_trans.append(current_transcript_comment)
		return aggregate_examined_rules(assigned_class_per_trans, assignment_comment_per_trans, "PM5")
