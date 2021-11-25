from genotoscope.utils import create_dir
from os.path import join
import numpy as np
import vcf
from pyensembl import EnsemblRelease


class FindClinicalSignificantExons:
	"""
	Class to identify clinical significant exons using gnomAD data sets
	Following 2.1.4 section of the research work:
	Xiang, Jiale, et al. "AutoPVS1: An automatic classification tool for PVS1 interpretation of null variants."
	Human Mutation 41.9 (2020): 1488-1498.
	"""

	def __init__(self, data_path, gnomAD_version, gnomAD_exomes_file, output_path):
		"""
		FindClinicalSignificantExons constructor

		Parameters
		----------
		data_path : str
			absolute path where gnomAD files lie
		gnomAD_version : str
			gnomAD version
		gnomAD_exomes_file : str
			gnomAD exomes absolute file name
		output_path : str
			output path, where resulted bed files will be saved

		Returns
		-------
		None
		"""
		# input and output data paths
		self.data_path = data_path
		output_root = join(self.data_path, output_path)
		self.output_path = join(output_root, gnomAD_version)
		create_dir(self.output_path)

		### gnomAD ###
		self.gnomAD_exomes_reader = self.load_vcf_file(gnomAD_exomes_file)

		### PyEnsembl ###
		# release 75 uses human reference genome GRCh37
		self.ensembl_data = EnsemblRelease(75)

	def load_vcf_file(self, file_path):
		"""
		Load data from VCF file

		Parameters
		----------
		file_path : str
			vcf file path
		Returns
		-------
		vcf.Reader object
			loaded vcf
		"""
		print("Load vcf file from: {}".format(file_path))
		# return vcf.Reader(filename=file_path, encoding='UTF-8')
		return vcf.VCFReader(open(file_path, 'rb'))

	def run(self, gnomAD_pLoF_file):
		"""
		Run function to find clinical significant exons

		Parameters
		----------
		gnomAD_pLoF_file : str
			gnomAD all pLoF variants absolute file name

		Returns
		-------
		None
		"""
		# aggregate PLoF variants per exon
		exons2pLoF_info = self.aggregate_variants_per_exon(gnomAD_pLoF_file)
		'''
		for exon_id, exon2var in exons2pLoF_info.items():
			print("exon id= {}".format(exon_id))
			print("exon chr={}, start={}, end={}".format(exon2var['chr'], exon2var['start'], exon2var['end']))
			print("var_info: \ntranscripts: {}\npLoF_info:{}".format(exon2var['transcripts'], exon2var['pLoF']))
		'''

		# identify clinical significant exons
		# following 2.1.4 from AutoPVS1 research work
		self.find_clinical_significant_exons(exons2pLoF_info)
		print("*** ***")

	@staticmethod
	def convert_pLoF2fields(pLoF_line):
		"""
		Convert pLoF variant line into dictionary with fields
		using as header:
		chrom	pos	ref	alt	most_severe_consequence	gene_ids	gene_symbols	transcript_ids

		tested pLoF file:
		https://storage.googleapis.com/gnomad-public/papers/2019-flagship-lof/v1.0/gnomad.v2.1.1.all_lofs.txt.bgz

		Parameters
		----------
		pLoF_line: str
			pLoF variant line

		Returns
		-------
		dict of str : int, str or list of str
			pLoF variant fields
		"""

		fields = pLoF_line.strip().split('\t')
		return {'chr': fields[0], 'pos': int(fields[1]), 'ref': fields[2], 'alt': fields[3].split(","),
		        'most_severe_consequence': fields[4], 'gene_ids': fields[5].split(","),
		        'gene_symbols': fields[6].split(","), "transcript_ids": fields[7].split(",")}

	def is_var_in_exon(self, var_fields, exon):
		"""
		Examine if variant is contained in exonic range (-/+2 (for splice acceptor variants))

		Parameters
		----------
		var_fields : dict of str : int, str or list of str
			variant fields
		exon : pyensembl.exon
			ensembl object of input exon

		Returns
		-------
		bool
			variant is contained in exonic range (True), otherwise False
		"""
		# pass exon strand direction into an integer
		if exon.to_dict()["strand"] == "+":
			strand_direction = +1
		else:
			strand_direction = -1

		if var_fields['most_severe_consequence'] == "splice_acceptor_variant":
			# add the -2 and +2 positions to capture splice acceptor variants
			exon_genom_range = range(exon.to_dict()["start"] - 2, exon.to_dict()["end"] + 2 + 1)
		else:
			# normal exon range
			exon_genom_range = range(exon.to_dict()["start"], exon.to_dict()["end"] + 1,
			                         strand_direction)
		if var_fields['pos'] in exon_genom_range:
			return True
		else:
			return False

	def find_var_exon_ids(self, var_fields):
		"""
		Find exon ids where variant lies in

		Parameters
		----------
		var_fields : dict of str : int, str or list of str
			variant fields

		Returns
		-------
		dict of str : dict of str : str or set of str
			map of variant including exons accompanying with their transcript id
		"""
		print("Find exon ids overlapping with variant: {}".format(var_fields))
		### ### ###
		# get all transcript ids associated with the variant
		### ### ###
		exons2transcripts = {}
		transcripts = []
		for transcript_id in var_fields['transcript_ids']:
			transcripts.append(self.ensembl_data.transcript_by_id(transcript_id))

		### ### ###
		# find exon id that variant lies in
		# save the mapping to its transcript id
		### ### ###
		for transcript in transcripts:
			print("Investigating transcript id: {}".format(transcript.id))
			for exon in transcript.exons:
				if self.is_var_in_exon(var_fields, exon):
					print("PLoF in exon id: {}".format(exon.id))
					if var_fields['most_severe_consequence'] == "splice_acceptor_variant":
						# extend exon start and stop to contain splice acceptor sites (+/- 1,2)
						overlapping_exon = {'chr': var_fields['chr'], 'start': exon.to_dict()["start"] - 2,
						                    'end': exon.to_dict()["end"] + 2,
						                    'strand': exon.to_dict()["strand"],
						                    'transcripts': {transcript.id}}
					else:
						# all other supported variants are inside the exon
						overlapping_exon = {'chr': var_fields['chr'], 'start': exon.to_dict()["start"],
						                    'end': exon.to_dict()["end"], 'strand': exon.to_dict()["strand"],
						                    'transcripts': {transcript.id}}

					if exon.id not in exons2transcripts:
						exons2transcripts[exon.id] = overlapping_exon
					else:
						# make sure that the version of the exon
						# that contains the splice acceptor sites will be saved finally
						same_exon_length = exons2transcripts[exon.id]['end'] - exons2transcripts[exon.id]['start'] + 1
						current_length = overlapping_exon['end'] - overlapping_exon['start'] + 1
						if current_length > same_exon_length:
							exons2transcripts[exon.id] = overlapping_exon
						# save transcript that contains the overlapping exon
						exons2transcripts[exon.id]['transcripts'].add(transcript.id)
		return exons2transcripts

	def aggregate_variants_per_exon(self, gnomAD_pLoF_file):
		"""
		Aggregate pLoF variants per exon

		Parameters
		----------
		gnomAD_pLoF_file : str
			gnomAD pLoF variants file path

		Returns
		-------
		dict of str : str or int or set of str
			exons mapped to variants
		"""
		print("Aggregate pLoF variants per exon")
		exons2var_info = {}
		num_pLoFs = 0
		with open(gnomAD_pLoF_file, 'r') as pLoF_handle:
			is_header = True
			# process each pLoF variant to get variant fields
			# find the exon ids where the variant lies in
			# aggregate exon ids
			for pLoF_line in pLoF_handle.readlines():
				if is_header:  # pass the first line
					is_header = False
					continue
				else:
					num_pLoFs += 1
					pLoF_fields = FindClinicalSignificantExons.convert_pLoF2fields(pLoF_line)
					exons2transcripts = self.find_var_exon_ids(pLoF_fields)
					for exon_id, exon2transcript_info in exons2transcripts.items():
						if exon_id not in exons2var_info:
							exons2var_info[exon_id] = {'chr': exon2transcript_info['chr'],
							                           'start': exon2transcript_info['start'],
							                           'end': exon2transcript_info['end'],
							                           'strand': exon2transcript_info['strand'],
							                           'transcripts': exon2transcript_info['transcripts'],
							                           'pLoF': [pLoF_fields]}
						else:
							exons2var_info[exon_id]['transcripts'] = exons2var_info[exon_id]['transcripts'].union(
								exon2transcript_info['transcripts'])
							exons2var_info[exon_id]['pLoF'].append(pLoF_fields)
		print("Num of all investigated pLoFs: {}".format(num_pLoFs))
		return exons2var_info

	def extract_matching_gnomAD_exome_rows(self, target_var_fields):
		"""
		Extract gnomAD exome variants matching target variant.
		The matching is on the position and the REF and ALT sequences

		Parameters
		----------
		target_var_fields : dict of str : str or int
			target variant's fields to be used to find matching gnomAD exome variants

		Returns
		-------
		list of dict of str : str or int
			matched gnomAD variants to target variant
		"""
		# print("Find gnomAD variants matching the pLoFs aggregated for exon")
		matched_gnomADs_equal_seq = []

		# find matching gnomAD variants by position
		# pyVCF uses 0-index for coordinate fetching
		matched_gnomAD_rows = list(
			self.gnomAD_exomes_reader.fetch(target_var_fields['chr'], target_var_fields['pos'] - 1,
			                                target_var_fields['pos']))

		# refine matched gnomAD variants by filtering in the ones with same ref and alt seq
		for gnomAD_row in matched_gnomAD_rows:
			gnomAD_row_alt_char = [str(c) for c in gnomAD_row.ALT]  # convert the alt elements to str
			if gnomAD_row.REF == target_var_fields['ref'] and set(gnomAD_row_alt_char) == set(target_var_fields['alt']):
				matched_gnomADs_equal_seq.append(gnomAD_row)
		return matched_gnomADs_equal_seq

	def quality_filter_variants(self, variants):
		"""
		Filter variants by gnomAD vcf quality values:
		* FILTER field -> PASS

		Parameters
		----------
		variants : list of dict of str : str or int
			input variants to be filtered

		Returns
		-------
		list of dict of str : str or int
			filtered in variants
		"""
		# print("Filter variants by quality")
		filtered_variants = []
		for var in variants:
			if len(var.FILTER) == 0:
				filtered_variants.append(var)
		return filtered_variants

	def sum_AF_per_subpopulation(self, variants):
		"""
		Sum AF per subpopulation for each quality filtered-in pLoF variants
		All variants are aggregated by exon id

		Parameters
		----------
		variants : list of dict of str: str or int
			input variants to sum AF for

		Returns
		-------
		dict of str: float
			sum of AF, per subpopulation, for all pLoFs of an exon
		"""
		print("Summing AF, per subpopulation, of each quality filtered-in pLoF variants")
		AF_subpopulations = {'AF_afr': 0.0, 'AF_asj': 0.0,
		                     'AF_amr': 0.0,
		                     'AF_sas': 0.0, 'AF_eas_jpn': 0.0, 'AF_eas_kor': 0.0, 'AF_eas_oea': 0.0, 'AF_eas': 0.0,
		                     'AF_nfe_nwe': 0.0, 'AF_nfe_seu': 0.0, 'AF_nfe_bgr': 0.0, 'AF_nfe_onf': 0.0,
		                     'AF_nfe_est': 0.0, 'AF_nfe': 0.0, 'AF_fin': 0.0,
		                     'AF_oth': 0.0}

		for var in variants:
			for population in AF_subpopulations.keys():
				if population in var.INFO:
					AF_subpopulations[population] = AF_subpopulations[population] + var.INFO[population][0]
		print("Sum AF subpopulations of all variants matching PLoF regions of exon: {}".format(AF_subpopulations))
		return AF_subpopulations

	def filter_sum_AF(self, AF_subpopulations):
		"""
		Filter variants by allele frequency
		if sum of AF for exons variants <= 0.001 in any population (Africa, America, Asia, Europe, Other) found in gnomAD vcf
		then filter is passed

		Parameters
		----------
		AF_subpopulations : dict of str: float
			sum of AF per subpopulation, for all variants of an exon

		Returns
		-------
		boolean
			sum AF of exon variants filter is passed (True), otherwise False
		"""
		print("Compute filter for sum AF of the exon variants")
		is_filter_passed = False
		passed_populations = []
		# for each subpopulation, check if the sum AF of exon variants is above the pre-defined cut-off
		# if yes, turn on filter passed flag
		for population, sum_AF in AF_subpopulations.items():
			if sum_AF <= 0.001:
				is_filter_passed = True
				passed_populations.append(population)
		return is_filter_passed

	def find_clinical_significant_exons(self, exons2vars_info):
		"""
		Find clinical significant exons

		Parameters
		----------
		exons2vars_info : dict of str : str or int or set of str
			aggregated pLoF variants per exon

		Returns
		-------
		None
		"""
		### ### ###
		# for each pLoF variant in an exon
		#  extract the gnomAD exome variants matching pLoF position and seq info
		#  filter the matching gnomAD variants by quality

		# then aggregate all quality filtered in gnomAD variants and:
		#  sum AF, per sub-population, for all aggregated variants of exon
		#  if sum AF is above pre-defined threshold
		#  => save exon in clinical significant exon bed file
		### ### ###
		exons_filename = "clinical_significant_exons.bed"
		exons_header = "#chr\tstart\tend\texon\ttranscripts\tstrand"
		clinical_significant_exons = []
		num_filtered_AF_gnomAD_per_exon = np.array([])
		print("For each exon: ")
		print("1) match gnomAD variants to pLoF variants aggregated for the exon")
		print("2) quality filter the matched gnomAD variants")
		print("3) If filter by sum AF subpopulations is PASS => save exon as (potentially) clinical significant")
		for exon_id, exon_var_info in exons2vars_info.items():
			print("Exon id: {}".format(exon_id))
			pLoF_variants = exon_var_info['pLoF']
			exon_qual_filtered_gnomAD_vars = []
			# extract all gnomAD variants matching pLoFs for the current exon
			# then quality filtered the exon gnomAD variants
			matched_gnomAD_vars_exon = []
			for pLoF_fields in pLoF_variants:
				matched_gnomAD_vars_exon += self.extract_matching_gnomAD_exome_rows(pLoF_fields)
			print("Number of matched gnomAD variants= {}".format(len(matched_gnomAD_vars_exon)))
			qual_filtered_gnomAD_vars_exon = self.quality_filter_variants(matched_gnomAD_vars_exon)
			print("Number of quality-filtered gnomAD variants= {}".format(len(qual_filtered_gnomAD_vars_exon)))

			# from all quality filtered variants of exons, sum the AF
			sum_AF_subpopulations = self.sum_AF_per_subpopulation(qual_filtered_gnomAD_vars_exon)
			is_sum_AF_filter_passed = self.filter_sum_AF(sum_AF_subpopulations)

			# if sum AF filter is passed, save exon as bed row
			if is_sum_AF_filter_passed:
				print("Filter is passed, save exon as clinical significant")
				exon_row = "\t".join(
					['chr' + exon_var_info['chr'], str(exon_var_info['start']), str(exon_var_info['end']), exon_id,
					 ','.join(exon_var_info['transcripts']), exon_var_info['strand']])
				clinical_significant_exons.append(exon_row)
				# keep the number of quality and AF filter gnomAD variants per exon
				num_filtered_AF_gnomAD_per_exon = np.append(num_filtered_AF_gnomAD_per_exon,
				                                            len(qual_filtered_gnomAD_vars_exon))
			print(" --- ")
		# finally save all clinical significant exon rows in a bed file
		print("Identified {} possibly clinical significant exons".format(len(clinical_significant_exons)))
		# print basic statistics for number of gnomAD variants per clinical significant exons
		print(
			"Basic statistics for filtered gnomAD variants (by quality and sum AF subpopulation) matching pLoFs per clinical significant exons:\nmin: {}, mean: {}, median: {}, max: {}".format(
				np.min(num_filtered_AF_gnomAD_per_exon), np.mean(num_filtered_AF_gnomAD_per_exon),
				np.median(num_filtered_AF_gnomAD_per_exon), np.max(num_filtered_AF_gnomAD_per_exon)))
		self.write_bed(exons_header, clinical_significant_exons, exons_filename)

	def write_bed(self, header_bed, bed_rows, bed_filename):
		"""
		Write bed file

		Parameters
		----------
		header_bed : str
			header for bed file
		bed_rows : list of str
			bed rows
		bed_filename : str
			bed file name

		Returns
		-------
		None
		"""
		# print("Writing bed file, containing {} rows".format(len(bed_rows)))
		with open(join(self.output_path, bed_filename), 'w') as bed_handle:
			bed_handle.write(header_bed + "\n")
			for bed_row in bed_rows:
				bed_handle.write(bed_row + "\n")
