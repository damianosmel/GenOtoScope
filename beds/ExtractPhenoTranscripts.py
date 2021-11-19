from genotoscope.utils import create_dir
from os.path import join
from pandas import read_excel

from pyensembl import EnsemblRelease
from pybedtools import BedTool


class ExtractPhenoTranscripts:
	"""
	Read hearing loss relevant transcripts and clinical significant exons from the followings two papers:

	1. DiStefano, Marina T., et al. "Curating clinically relevant transcripts for the interpretation of sequence variants."
	The Journal of Molecular Diagnostics 20.6 (2018): 789-801.
	DOI: 10.1016/j.jmoldx.2018.06.005.
	2. DiStefano, Marina T., et al. "ClinGen expert clinical validity curation of 164 hearing loss gene–disease pairs."
	Genetics in Medicine 21.10 (2019): 2239-2247.
	DOI: doi: 10.1038/s41436-019-0487-0.
	"""

	def __init__(self, data_path, output_root, output_dir):
		"""

		Parameters
		----------
		data_path : str
			hearing loss specific info data path
		output_root : str
			output root, where the produced bed files will be saved
		output_dir : str
			output directory where the bed file will be placed

		Returns
		-------
		None
		"""
		# input and output directories
		self.data_path = data_path
		self.output_dir = join(output_root, output_dir)
		create_dir(self.output_dir)

		### PyEnsembl ###
		# release 75 uses human reference genome GRCh37
		self.ensembl_data = EnsemblRelease(75)

	def get_gene_chr(self, gene_name):
		"""
		Get chromosome number from gene name

		Parameters
		----------
		gene_name: str
			gene name

		Returns
		-------
		str
			chromosome at which gene lies in
		"""
		known_chromosomes = [str(i) for i in range(1, 23)]
		known_chromosomes += ['X', 'Y', 'MT']
		if len(gene_name.strip().split(" ")) > 1:
			[gene, gene_alternative] = gene_name.strip().split(" ")
		else:
			gene = gene_name.strip()
			gene_alternative = None
		try:
			# print("try to get locus for gene name: {}".format(gene))
			gene_locus = self.ensembl_data.loci_of_gene_names(gene)
		except ValueError:
			if gene_alternative:  # if gene alternative name exists
				try:
					gene_alternative = gene_alternative.replace("(", "")
					gene_alternative = gene_alternative.replace(")", "")
					# print("gene locus could not be found, try with alternative name: {}".format(gene_alternative))
					gene_locus = self.ensembl_data.loci_of_gene_names(gene_alternative)
				except ValueError:
					print("Could not find locus for input gene name (first and alternative): {}".format(gene_name))
					gene_locus = None
			else:
				# alternative name does not exist, and gene locus was not found for the first gene name
				print("Could not find locus for input gene name: {}".format(gene_name))
				gene_locus = None

		# convert gene locus to chromosome number
		if gene_locus:
			first_locus = gene_locus[0]
			chr = first_locus.to_dict()['contig']
			print("Gene is at chromosome: {}".format(chr))
			assert chr in known_chromosomes, "AssertionError: retrieved chromosome number should be one of [1:22] + X,Y"
			return chr
		else:
			return None

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
		print("Writing bed file, containing {} rows".format(len(bed_rows)))
		with open(join(self.output_dir, bed_filename), 'w') as bed_handle:
			bed_handle.write(header_bed + "\n")
			for bed_row in bed_rows:
				bed_handle.write(bed_row + "\n")

	@staticmethod
	def parse_exons_index(exons_column):
		"""
		Parse phenotype relevant exon indices

		Parameters
		----------
		exons_column : str
			clinical significant exon indices column

		Returns
		-------
		dict of str : list of int
			dictionary of clinical significant exons (to be included or excluded)
		"""
		print("Parse exons column: {}".format(exons_column))
		clinical_exons = {"include": [], "exclude": []}
		if exons_column[0].isdigit():
			### ### ###
			# exons to include
			### ### ###
			if "-" in exons_column:
				# at least two ranges of exon indices to be included
				for exons_range in exons_column.split(","):
					if "-" in exons_range:
						exon_start, exon_end = exons_range.strip().split("-")
						clinical_exons["include"] = clinical_exons["include"] + [exon_idx for exon_idx in
						                                                         range(int(exon_start) - 1,
						                                                               int(exon_end))]
					else:
						clinical_exons["include"] = clinical_exons["include"] + [int(exons_range.strip()) - 1]
			else:
				# only one exon to include
				clinical_exons["include"] = clinical_exons["include"] + [int(exons_column) - 1]
		elif exons_column[0] == "-":
			if "(" in exons_column and ")" in exons_column:
				# at least two ranges of exon indices to be excluded
				exons_column = exons_column[2:-1]
				for exons_range in exons_column.split(","):
					if "-" in exons_range:
						exon_start, exon_end = exons_range.strip().split("-")
						clinical_exons["exclude"] = clinical_exons["exclude"] + [exon_idx for exon_idx in
						                                                         range(int(exon_start) - 1,
						                                                               int(exon_end))]
					else:
						clinical_exons["exclude"] = clinical_exons["exclude"] + [int(exons_range.strip()) - 1]
			else:
				# only one exon to be excluded
				clinical_exons["exclude"] = clinical_exons["exclude"] + [-1 * (int(exons_column) + 1)]
		else:
			# all exons to be included
			# equivalent exclude a non existent exon index
			clinical_exons["exclude"] = [-1]
		print("Clinical exons to process: {}".format(clinical_exons))
		### ### ###
		# assert that exons are only included or only excluded
		### ### ###
		try:
			assert bool(len(clinical_exons["include"]) > 0) ^ bool(len(clinical_exons["exclude"]) > 0)
		except AssertionError:
			print(
				"AssertionError: Please make clinical_exons columns to described only inclusion of exons or only exclusion of them")

		return clinical_exons

	def extract_pheno_relevant_exons(self, pheno_relevant_transcripts_file):
		"""
		Create phenotype relevant exons using transcript and gene information

		Parameters
		----------
		pheno_relevant_transcripts_file : str
			full path to phenotpy relevant transcripts list

		Returns
		-------
		list of str
			phenotype relevant exons, one per list element
		"""
		print("Extract information for phenotype relevant exons")
		pheno_transcripts_df = read_excel(join(self.data_path, pheno_relevant_transcripts_file), header=0, index_col=0)
		# for each gene get the chromosome number
		# and then save all exons of the transcripts to an list
		pheno_relevant_exons = []
		genes_transcripts_df = pheno_transcripts_df.iloc[:, 0:2]
		gene_transcripts_exons_df = pheno_transcripts_df.iloc[:, 0:4]
		extract_stats = {"uniq_genes": set(), "uniq_trans": set(), "uniq_exons": set(), "grch38_only": set(),
		                 "no_distinct_exons": set(), "already_added": set()}

		for gene_name, row in gene_transcripts_exons_df.iterrows():
			if row["comment"] == "only_grch38":
				# exclude current transcript, as is found only on GRCh38
				extract_stats["grch38_only"].add(row["ensembl GRCh37.p13"].strip())
			elif row["comment"] == "transcript_already_included":
				# exclude current transcript, as it is already added
				extract_stats["already_added"].add(row["ensembl GRCh37.p13"].strip())
			elif row["clinical_exons"] == "no_distinct_coding_exons":
				# exclude current transcript, as it does not contain distinct exons
				extract_stats["no_distinct_exons"].add(row["ensembl GRCh37.p13"].strip())
			else:
				extract_stats["uniq_genes"].add(gene_name)
				# get transcript and chromosome info
				if len(row["ensembl GRCh37.p13"].split(".")[0]) > 1:
					transcript_id = row["ensembl GRCh37.p13"].split(".")[0]
				else:
					transcript_id = row["ensembl GRCh37.p13"]
				extract_stats["uniq_trans"].add(transcript_id)
				transcript_chrom = self.get_gene_chr(gene_name)
				transcript = self.ensembl_data.transcript_by_id(transcript_id)

				# get clinical exons
				clinical_exons = ExtractPhenoTranscripts.parse_exons_index(str(row["clinical_exons"]))
				# understand if you need to only include exons or only exlude
				if len(clinical_exons["include"]) > 0:
					include_or_exclude = 0
				elif len(clinical_exons["exclude"]) > 0:
					include_or_exclude = 1
				for exon_idx, exon in enumerate(transcript.exons):
					if include_or_exclude == 0:
						if exon_idx in clinical_exons["include"]:
							extract_stats["uniq_exons"].add(exon.id)
							pheno_relevant_exons.append(
								"\t".join(
									["chr" + transcript_chrom, str(exon.to_dict()["start"]), str(exon.to_dict()["end"]),
									 exon.id, transcript_id, exon.to_dict()["strand"]]))
					else:
						if exon_idx not in clinical_exons["exclude"]:
							extract_stats["uniq_exons"].add(exon.id)
							pheno_relevant_exons.append(
								"\t".join(
									["chr" + transcript_chrom, str(exon.to_dict()["start"]), str(exon.to_dict()["end"]),
									 exon.id, transcript_id, exon.to_dict()["strand"]]))
			print("--")
		print("### ### ###")
		print("Extraction statistics:")
		for k, v in extract_stats.items():
			print(k + " contains: {} elements, which are:".format(len(v)))
			print(v)
		print("### ### ###")
		return pheno_relevant_exons

	def run(self, pheno_relevant_transcripts_file, header_bed, out_filename):
		"""
		Extract phenotype relevant transcripts and write their clinical significant exons on bed file

		Parameters
		----------
		pheno_relevant_transcripts_file : str
			phenotype relevant transcripts xslx file
		header_bed : str
			bed header for extracted phenotype relevant transcripts
		out_filename: str
			bed annotation filename

		Returns
		-------
		None
		"""
		print("Extract phenotype relevant transcripts and write their clinical significant exons on bed file")
		# extract phenotype relevant transcripts and their clinical significant exons
		pheno_relevant_exons = self.extract_pheno_relevant_exons(pheno_relevant_transcripts_file)
		# write them on a bed file
		self.write_bed(header_bed, pheno_relevant_exons, out_filename)
		# sort bed file
		print("Sort resulted bed and save sorted version")
		pheno_relevant_exons = BedTool(join(self.output_dir, out_filename)).sort()
		pheno_relevant_exons.saveas(join(self.output_dir, out_filename.split(".bed")[0] + "_sorted" + ".bed"))
