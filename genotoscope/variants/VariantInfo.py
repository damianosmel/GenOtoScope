from Bio.Seq import Seq
from datetime import datetime
from os.path import join

class VariantInfo:
	"""
	Class to keep variant information
	"""

	def __init__(self, sample_name, gene_name, variant_type, transcripts_info, var_chrom, var_genomic_start,
	             var_genomic_end,
	             var_ref, var_obs):
		"""
		Parameters
		----------
		sample_name : str
			sample name containing the specific variant
		gene_name : str
			name of gene's harbouring variant
		variant_type : str
			variant type column in GSvar file
		transcripts_info : list of dict of str : str
			variant information for each affected transcript
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
		"""
		self.sample_name = sample_name
		self.gene_name = gene_name
		self.variant_type = variant_type
		self.transcripts_info = transcripts_info
		self.chrom = var_chrom
		self.genomic_start = var_genomic_start
		self.genomic_end = var_genomic_end
		self.ref_base = var_ref
		self.obs_base = var_obs

	def to_string(self):
		return ",".join(
			["sample_name: " + self.sample_name, "chrom: " + self.chrom, "start: " + str(self.genomic_start),
			 "end: " + str(self.genomic_end)])

	def create_bed_line(self, transcript_strand, format):
		"""
		Create bed line to represent variant info
		Credits: http://daler.github.io/pybedtools/intervals.html

		Parameters
		----------
		transcript_strand : str
			strand of variant-affected transcript

		Returns
		-------
		str
			variant bed line
		"""
		if format == "6columns":
			return ' '.join(
				[self.chrom, str(self.genomic_start), str(self.genomic_end), self.gene_name, '.', transcript_strand])

	def is_same_seq_edit(self, var2compare):
		"""
		Check if self variant and variant to compare have the same sequence edit

		Parameters
		----------
		var2compare:

		Returns
		-------
		bool
			sequence edits are equal (True), otherwise False
		"""
		equal_ref_base, equal_obs_base = False, False
		if str(self.ref_base) == "-" and str(var2compare.ref_base) == "-":
			equal_ref_base = True
		elif str(self.ref_base) == str(var2compare.ref_base) or str(self.ref_base) == str(
				Seq(var2compare.ref_base).reverse_complement()):
			# check for sequence equality on positive strand or negative strand
			equal_ref_base = True

		if str(self.obs_base) == "-" and str(var2compare.obs_base) == "-":
			equal_obs_base = True
		elif str(self.obs_base) == str(var2compare.obs_base) or str(self.obs_base) == str(
				Seq(var2compare.obs_base).reverse_complement()):
			# check for sequence equality on positive strand or negative strand
			equal_obs_base = True
		return equal_ref_base and equal_obs_base

	def is_same_var(self, var2compare):
		"""
		Check if variant is the same
		to the variant characterized by the input arguments

		Parameters
		----------
		var2compare: VariantInfo
			object of variant to compare to the first

		Returns
		-------
		bool
			variant is the same as the input one (True), otherwise False
		"""
		if self.chrom == var2compare.chrom and int(self.genomic_start) == int(var2compare.genomic_start) and int(
				self.genomic_end) == int(var2compare.genomic_end) and self.is_same_seq_edit(var2compare):
			return True
		else:
			return False

	def save_vcf(self,source,reference,output_dir):
		"""
		Save variant in VCF format
		VCF version: 4.1

		Parameters
		----------
		source : str
			source metadata field
		reference :
			reference metadata field
		output_dir :
			output directory where the vcf file will be saved

		Returns
		-------
		None
		"""
		# convert the spaces with tabs and create the vcf name from the variant description
		vcf_name = "_".join([self.chrom,self.genomic_start,self.genomic_end,self.ref_base,self.obs_base]) + ".vcf"
		datetime.today().strftime("%d-%m-%Y")
		header_rows = ["##fileformat=VCFv4.1","##fileDate="+datetime.today().strftime("%d-%m-%Y"),"##Source="+source,"##reference="+reference,"##ID=<Description=\"ClinVar Variation ID\">","#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"]
		# replace stop position with clinvar id (currently use '.')
		# add "." (unknown) for QUAL, FILTER and INFO fields
		description_tabs = [self.chrom.split("chr")[-1],self.genomic_start,".",self.ref_base,self.obs_base] + [".",".","."]
		# write variant into the vcf file
		with open(join(output_dir,vcf_name),'w') as vcf_handle:
			for header_row in header_rows:
				vcf_handle.write(header_row+"\n")
			vcf_handle.write("\t".join(description_tabs)+"\n")