from genotoscope.utils import create_dir,load_hugo_genes_df, normalize_gene_info,normalize_disease_id, convert_review_status2stars,get_clinvar_strand
from os.path import join, isfile
from pandas import read_csv
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import vcf
from pybedtools import BedTool
from bioservices import HGNC
from pyensembl import EnsemblRelease
from math import ceil
from scipy.stats.mstats import mquantiles
import pickle

plt.style.use('seaborn-paper')


class FindCriticalProteinRegions:
	"""
	Use ClinVar vcf file and UniProt domain annotation track to find critical regions for protein function
	Following section 2.1.3 of research work:
	Xiang, Jiale, et al. "AutoPVS1: An automatic classification tool for PVS1 interpretation of null variants."
	Human Mutation 41.9 (2020): 1488-1498.

	Also follow VarSome definition of critical regions of protein:
	https://varsome.com/about/resources/acmg-implementation/#pm1
	publication:
	Kopanos, Christos, et al. "VarSome: the human genomic variant search engine." Bioinformatics 35.11 (2019): 1978.
	"""

	def __init__(self, data_path, clinvar_path, clinvar_file, clinvar_review_stars_file, clinvar_version, uniprot_domains_path,
	             uniprot_version,hugo_genes_file, output_path):
		"""
		FindCriticalProteinRegions contructor

		Parameters
		----------
		data_path : str
			absolute path where clinvar and uniprot annotation files lie
		clinvar_path : str
			Clinvar root path
		clinvar_file: str
			clinvar vcf file
		clinvar_review_stars_file : str
			clinvar review stars file
		clinvar_version : str
			clinvar version
		uniprot_domains_path : str
			absolute path of uniprot domains file
		uniprot_version : str
			uniprot version
		hugo_genes_file : str
			hugo genes file name
		output_path : str
			output path, where resulted bed files will be saved

		Returns
		-------
		None
		"""

		# input and output data paths
		self.data_path = data_path

		output_root = join(self.data_path, output_path)
		self.output_path = join(output_root, clinvar_version + "_" + uniprot_version)
		create_dir(self.output_path)

		### ClinVar ###
		self.clinvar_root = join(data_path,clinvar_path)
		self.clinvar_path = join(self.clinvar_root, clinvar_version)
		self.load_clinvar_file(clinvar_file)
		# quality star system
		self.star_status2int = {"none": 0, "one": 1, "two": 2, "three": 3, "four": 4, "unknown review status": -1}
		self.review2color = {-1: "black", 0: "grey", 1: "orange", 2: "goldenrod", 3: "gold", 4: "yellow"}
		self.read_clinvar_review_stars(clinvar_review_stars_file)

		### UniProt domains ###
		self.uniprot_domains = self.load_uniprot_annotation_file(uniprot_domains_path)

		### PyEnsembl ###
		# release 75 uses human reference genome GRCh37
		self.ensembl_data = EnsemblRelease(75)

		### Hugo genes ###
		self.hugo_genes_df = load_hugo_genes_df(self.data_path, hugo_genes_file)

	def load_uniprot_annotation_file(self, uniprot_annotations_path):
		"""
		Load uniprot annotation into BedTool object

		Parameters
		----------
		uniprot_annotations_path : str
			UniProt annonations path

		Returns
		-------
		uniprot_annot: BedTool
			BedTool object containing uniprot annotations
		"""
		print("Loading UniProt annotation file from {}".format(uniprot_annotations_path))
		return BedTool(uniprot_annotations_path)

	def load_clinvar_file(self, clinvar_file):
		"""
		Load clinvar file to find annotated variants
		to decide for ACMG/AMP rules

		Parameters
		----------
		clinvar_file : str
			clinvar file name

		Returns
		-------
		None
		"""
		print("Loading ClinVar file from {}".format(join(self.clinvar_path, clinvar_file)))
		self.vcf_reader = vcf.Reader(filename=join(self.clinvar_path, clinvar_file), encoding='UTF-8')

	def read_clinvar_review_stars(self, clinvar_review_stars_file):
		"""
		Load clinvar review stars into dataframe

		Parameters
		----------
		clinvar_stars_file : str
		clinvar stars file name

		Returns
		-------
		None
		"""
		self.clinvar_stars_df = read_csv(join(self.clinvar_root, clinvar_review_stars_file), sep="\t", header=0)
		self.clinvar_stars_df = self.clinvar_stars_df.rename(lambda x: x.strip().replace(' ', '_'), axis='columns')
		self.clinvar_stars_df['Review_status'] = self.clinvar_stars_df['Review_status'].str.strip()
		self.clinvar_stars_df['Number_of_gold_stars'] = self.clinvar_stars_df['Number_of_gold_stars'].str.strip()

	def run(self, min_review_stars):
		"""
		Run function to find critical protein regions

		Parameters
		----------
		min_review_stars : int
			filter in clinvar entries with at least minimum number of review stars

		Returns
		-------
		None
		"""
		print("Find critical protein regions")

		### ### ###
		# extract clinvars intersecting with domains
		### ### ###
		domain_clinvars = self.aggregate_clinvars_per_domain(min_review_stars)
		self.calculate_dom_clinvar_stats(domain_clinvars)

		'''
		### create the significance summary of the ClinVars intersecting domains ###
		significance_summary = self.summarize_significance(domain_clinvars)
		# for significance, counts_per_review in significance_summary.items():
		# 	print("significance: {}".format(significance))
		# 	print("counts: {}".format(counts_per_review))

		### plot the significance summary ###
		use_review_status = True
		self.plot_significance(significance_summary, use_review_status)
		use_review_status = False
		self.plot_significance(significance_summary, use_review_status)
		'''
		### ### ###
		# identify the important truncated/altered regions
		### ### ###
		print("")
		criterion, min_ratio_cutoff = "no_benign_atleast5pathogenic", None
		self.identify_important_truncated_regions(domain_clinvars, criterion, min_ratio_cutoff)

		print("")
		criterion, min_ratio_cutoff = "pathogenic_probability", 0.51
		self.identify_important_truncated_regions(domain_clinvars, criterion, min_ratio_cutoff)

		### ### ###
		# identify well-known domains without benign variant
		### ### ###
		print("")
		criterion, min_ratio_cutoff = "pathogenic_probability", 0.51
		self.identify_domains_without_benign(domain_clinvars, criterion, min_ratio_cutoff)

		plot_ratio = True
		self.calculate_domains_pathogenic_distribution(domain_clinvars, plot_ratio)

		print("*** ***")

	def calculate_dom_clinvar_stats(self, domain_clinvars, print_stats=True):
		"""
		Calculate basic count statistics for UniProt domains and ClinVar entries

		Parameters
		----------
		domain_clinvars : list of dict
			filtered clinvars per domain
		print_stats : bool
			print calculated statistics

		Returns
		-------
		None
		"""
		total_domains, total_clinvars = 0, 0
		clinvars_per_dom = np.array([])
		for dom_clinvar in domain_clinvars:
			total_domains += 1
			total_clinvars += len(dom_clinvar.values())
			clinvars_per_dom = np.append(clinvars_per_dom, len(dom_clinvar.values()))
		if print_stats:
			print("\n# ### ### ### ### ### ### ### ### ### #")
			print(" UniProt domains & ClinVar entries ")
			print("        ###  Statistics  ###         ")
			print("# total UniProt domains= {}".format(total_domains))
			print("# total ClinVar entries= {}".format(total_clinvars))
			print("# ClinVars per domain=> min={}, mean={}, max={}".format(np.min(clinvars_per_dom),
			                                                               np.mean(clinvars_per_dom),
			                                                               np.max(clinvars_per_dom)))
			print("# ### ### ### ### ### ### ### ### ### #\n")

	def parse_bed_row(self, bed_row):
		"""
		Parse bed row and convert it to a custom annotation dictionary

		Parameters
		----------
		bed_row : str
			bed annotation row

		Returns
		-------
		dict of str : int or str
			dictionary to saved selected fields from input annotation row
		"""
		fields = bed_row.strip().split('\t')
		if len(fields[-1].split(";")) > 1:
			description = fields[-1].split(";")[1].lstrip().replace(" ", "_")
		else:
			description = fields[-1]

		### ### ###
		# accept annotation with correct chromosome name
		### ### ###
		if fields[0].split("chr")[1].isdigit():
			annotation = {'chr': fields[0], 'start': int(fields[1]), 'end': int(fields[2]), 'strand': fields[5].strip(),
		              'uniprot_id': fields[3],
		              'uniprot_description': description}
			return annotation
		else:
			return None

	def extract_clinvars(self, annotation):
		"""
		Extract ClinVar entries found in the positions in the range of the annotation

		Parameters
		----------
		annotation : dict of str: int or str
			uniprot annotation to restrict ClinVars extraction

		Returns
		-------
		list of dict of str: int or str or list of str
			extracted clivnar entries
		"""
		print("Extract ClinVar variants for input annotation")
		if annotation:
			chr, start, end = annotation["chr"].split("chr")[1], annotation["start"] - 1, annotation["end"]
			print("chr={},start={},end={}".format(chr,start,end))
			matched_clinvars = list(self.vcf_reader.fetch(chr, start, end))

			if len(matched_clinvars) > 0:
				print("Found {} matching ClinVar entries".format(len(matched_clinvars)))
				extracted_clinvars = []
				for clinvar_rec in matched_clinvars:
					genes_info = normalize_gene_info(clinvar_rec)
					if genes_info:  # get the symbol for the first registered gene in ClinVar
						gene_symbol, gene_id = genes_info[0]
					else:
						gene_symbol = None
					# save a clinvar entry if
					# a) the strand is the same as for the UniProt annotation
					# b) the significance field is filled
					if get_clinvar_strand(self.hugo_genes_df,gene_symbol) == annotation["strand"] and 'CLNSIG' in clinvar_rec.INFO:
						clinvar_entry = {'ID': clinvar_rec.ID, 'POS': int(clinvar_rec.POS), 'REF': str(clinvar_rec.REF),
						                 'ALT': ','.join([str(alt) for alt in clinvar_rec.ALT]),
						                 'CLNDISDB': normalize_disease_id(clinvar_rec),
						                 'CLNREVSTAT': convert_review_status2stars(self.clinvar_stars_df,self.star_status2int, clinvar_rec.INFO['CLNREVSTAT']),
						                 'CLNSIG': clinvar_rec.INFO['CLNSIG'],
						                 'ENTREZ_GENE_ID': gene_id}
						if 'None' not in clinvar_entry['ALT']:
							# do not extract a ClinVar record that contains the None nucleotide as alternate
							extracted_clinvars.append(clinvar_entry)
				if len(extracted_clinvars) == 0:
					# if there is no extracted clinvar entries, then return None
					extracted_clinvars = None
			else:
				extracted_clinvars = None
		else:
			# domain annotation was not available
			extracted_clinvars = None
		return extracted_clinvars

	@staticmethod
	def label_bar(ax, rects):
		"""
		Attach a text label above each bar in *rects*, displaying its height.
		Credits: https://matplotlib.org/3.1.1/gallery/lines_bars_and_markers/barchart.html#sphx-glr-gallery-lines-bars-and-markers-barchart-py

		"""
		for rect in rects:
			height = rect.get_height()
			if height != 0:
				ax.annotate('{}'.format(height), xy=(rect.get_x(), height + rect.get_width() / 2), xytext=(0, 3),
				            # 3 points vertical offset
				            textcoords="offset points", ha='center', va='bottom')

	@staticmethod
	def label_barh(ax, rects):
		"""
		Label horizontial bars

		Parameters
		----------
		ax : matplotlib.axes
			axes object
		rects :
		horizontial bars to be annotated

		Credits : https://matplotlib.org/3.3.1/gallery/statistics/barchart_demo.html#sphx-glr-gallery-statistics-barchart-demo-py
		"""
		# Lastly, write in the ranking inside each bar to aid in interpretation
		for rect in rects:
			# Rectangle widths are already integer-valued but are floating
			# type, so it helps to remove the trailing decimal point and 0 by
			# converting width to int type
			width = int(rect.get_width())
			height = int(rect.get_height())

			# The bars aren't wide enough to print the ranking inside
			if width < 40:
				# Shift the text to the right side of the right edge
				xloc = 5
				# Black against white background
				clr = 'black'
				align = 'left'
			else:
				# Shift the text to the left side of the right edge
				xloc = -5
				# White on magenta
				clr = 'white'
				align = 'right'

			# Center the text vertically in the bar
			yloc = rect.get_y() + height / 2
			if height != 0:
				label = ax.annotate(
					height, xy=(width, yloc), xytext=(xloc, 0),
					textcoords="offset points",
					horizontalalignment=align, verticalalignment='center',
					color=clr, weight='bold', clip_on=True)

	def plot_significance(self, significance_summary, show_review_status):
		"""
		Plot significance of ClinVar entries (also separated by review status)
		Credits: # https://matplotlib.org/3.1.1/gallery/lines_bars_and_markers/bar_stacked.html#sphx-glr-gallery-lines-bars-and-markers-bar-stacked-py
		Parameters
		----------
		dict of str: dict of int: int
			significance summary aggregated for review quality
		show_review_status: bool
			separate each significance bar per review status (True), otherwise show aggregated significance
		Returns
		-------
		None
		"""
		fig, ax = plt.subplots()
		if show_review_status:
			plot_file_name = "significance_summary_review_stars.png"
			x_labels_loc = np.arange(len(significance_summary))
			width = 0.3
			all_significances = list(significance_summary.keys())
			# aggregate per review
			review2significance = {}
			for sig, counts_per_review in significance_summary.items():
				for review, counts in counts_per_review.items():
					if review not in review2significance:
						review2significance[review] = {}
						if sig not in review2significance[review]:
							review2significance[review][sig] = counts
					else:
						if sig not in review2significance[review]:
							review2significance[review][sig] = counts

			# normalize review2sig dictionary to save all possible significances
			for review_status, significances in review2significance.items():
				for sig in all_significances:
					if sig not in list(significances):
						review2significance[review_status][sig] = 0
			boxes = []
			review_idx = 0
			for review, significance in review2significance.items():
				significances_per_review = ax.barh(x_labels_loc + review_idx * width, significance.values(), width,
				                                   align='center', color=mcolors.CSS4_COLORS[self.review2color[review]],
				                                   edgecolor='k', alpha=0.8)
				boxes.append(significances_per_review)
				review_idx += 1

			plt.legend([box[0] for box in boxes], [review for review in list(review2significance.keys())],
			           title="Review status stars", loc='upper right')
			plt.ylabel('Significance', fontsize=14)
			plt.xlabel('Num. of ClinVar entries $\cap$ UniProt domains', fontsize=14)
			plt.yticks(range(len(all_significances)),
			           [sig.replace("/", "/\n").replace("_", " ") for sig in all_significances])
			sig = [list(counts_per_review_status.values()) for significance, counts_per_review_status in
			       significance_summary.items()]
			plt.xlim(0, max([sum(s) for s in sig]) + 10)

			ax.xaxis.grid(True, linestyle='--', which='major',
			              color='grey', alpha=.25)
		else:
			plot_file_name = "significance_summary.png"
			summary_all_reviews = {sig: sum(counts_per_review.values()) for (sig, counts_per_review) in
			                       significance_summary.items()}

			recs = plt.bar(range(len(summary_all_reviews)), list(summary_all_reviews.values()), align='center',
			               color='tab:green', edgecolor='k', alpha=0.8)
			plt.xlabel('Significance', fontsize=14)
			plt.ylabel('Num. of ClinVar entries $\cap$ UniProt domains', fontsize=14)
			plt.xticks(range(len(summary_all_reviews)),
			           [sig.replace("/", "/\n").replace("_", " ") for sig in list(summary_all_reviews.keys())],
			           rotation=45)
			plt.ylim(0, max(summary_all_reviews.values()) + 10)
			ax.yaxis.grid(True, linestyle='--', which='major',
			              color='grey', alpha=.25)
		print("Plotting significance summary at {}".format(join(self.output_path, plot_file_name)))
		fig.savefig(join(self.output_path, plot_file_name), bbox_inches='tight', dpi=600)

	@staticmethod
	def calculate_pathogenic_probability(clinvars):
		"""
		Calculate pathogenic over benign ratio
		p_{pathogenic} = (count(pathogenic or likely pathogenic) + delta) / Sum_{x=pathogenic,benign,VUS} (count(x) + delta)
		Credits: https://en.wikipedia.org/wiki/Additive_smoothing

		Parameters
		----------
		clinvars : list of dict
			clinvar entries each organized in dict

		Returns
		-------
		float
			pathogenic probability
		"""
		num_B_LB, num_P_LP, num_VUS = FindCriticalProteinRegions.calculate_pathogenic_class_counts(clinvars)
		delta = 10 ** (-6)

		# apply smoothing to compute the ratio
		return (num_P_LP + delta) / sum([num_B_LB + delta, num_P_LP + delta, num_VUS + delta])

	@staticmethod
	def calculate_pathogenic_class_counts(clinvars):
		"""
		Calculate pathogenic classes counts (benign/likely benign, pathogenic/likely pathogenic, VUS)

		Parameters
		----------
		clinvars : list of dict
			clinvar entries each organized in dict

		Returns
		-------
		int
			number of benign/likely benign ClinVar records
		int
			number of pathogenic/likely pathogenic ClinVar records
		int
			number of variants of uncertain significance (VUS) ClinVar records
		"""
		num_B_LB, num_P_LP, num_VUS = 0, 0, 0

		for clinvar in clinvars:
			current_significance = clinvar['CLNSIG'][0]
			if "benign" in current_significance.lower() or "likely_benign" in current_significance.lower():
				num_B_LB += 1
			elif current_significance.lower() == "pathogenic" or current_significance.lower() == "likely_pathogenic" or current_significance.lower() == "pathogenic/likely_pathogenic":
				num_P_LP += 1
			elif "uncertain_significance" in current_significance.lower():
				num_VUS += 1
		return num_B_LB, num_P_LP, num_VUS

	def calculate_domains_pathogenic_distribution(self, domain_clinvars, plot_ratio):
		"""
		Calculate domains pathogenic/likely pathogenic over benign/likely benign ratio

		Parameters
		----------
		domain_clinvars : list of dict
			bed annotation and clinvar entries per domain
		plot_ratio : bool
			Plot computed ratio and save it on file (True), otherwise do not plot ratio (False)

		Returns
		-------
		None
		"""
		domains_pathogenic_distribution = []
		# for each domain, check its ClinVar entries to calculate its pathogenic probability
		for dom_clinvar in domain_clinvars:
			domains_pathogenic_distribution.append(
				FindCriticalProteinRegions.calculate_pathogenic_probability(dom_clinvar['clinvar']))

		if plot_ratio:
			fig = plt.figure()
			n, bins, patches = plt.hist(domains_pathogenic_distribution, histtype='bar', color='g',
			                            align='left', edgecolor='k', alpha=0.8)

			plt.xlabel("P/LP empirical probability distribution", fontsize=14)
			plt.ylabel("Num. of UniProt domains", fontsize=14)
			plt.xticks([bin_value for bin_value in bins])
			plt.xlim(min(bins), max(bins))
			ratio_plot_name = "pathogenic_prob_histogram.png"
			fig.savefig(join(self.output_path, ratio_plot_name), bbox_inches='tight', dpi=600)
			print("Plotted Pathogenic/likely pathogenic probability empirical distribution at {}".format(
				join(self.output_path, ratio_plot_name)))

	@staticmethod
	def is_region_important(clinvars, criterion, min_cutoff):
		"""
		Investigate if domain region is important by criteria discussed in
		2.1.3 section of AutoPVS1 research work

		Parameters
		----------
		clinvars : list of dict
			clinvar entries each organized in dict
		criterion : str
			criterion to decide for region importance
		min_cutoff: float
			minimum cutoff to decide for a region, if criterion is using cutoff value

		Returns
		-------
		float
			region score
		bool
			True: region is important, otherwise False
		"""
		is_important_region = False
		if criterion == "no_benign_atleast5pathogenic":
			num_B_LB, num_P_LP, num_VUS = FindCriticalProteinRegions.calculate_pathogenic_class_counts(clinvars)
			region_score = 0.0
			if num_B_LB == 0 and num_P_LP >= 5:
				region_score = float(num_P_LP)
				is_important_region = True
		elif criterion == "pathogenic_probability":
			region_score = FindCriticalProteinRegions.calculate_pathogenic_probability(clinvars)
			if region_score >= min_cutoff:
				is_important_region = True
		return region_score, is_important_region

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
		with open(join(self.output_path, bed_filename), 'w') as bed_handle:
			bed_handle.write(header_bed + "\n")
			for bed_row in bed_rows:
				bed_handle.write(bed_row + "\n")

	def identify_important_truncated_regions(self, domain_clinvars, criterion, min_ratio_cutoff):
		"""
		Identify important truncated/altered protein regions
		then write them on bed file using 6-column format,

		Following:
		* Section 2.1.3 from AutoPVS1 research work
		* bedtools guidelines:
		https://felixfan.github.io/bedtools/

		Parameters
		----------
		domain_clinvars: list of dict
			bed annotation and clinvar entries per domain
		criterion : str
			criterion to decide for region importance
		min_ratio_cutoff: float
			minimum ratio cutoff to decide for a region, if criterion is using the ratio of B_LB over P_LP

		Returns
		-------
		None
		"""
		print("Identify important truncated/altered protein regions using criterion: {}".format(criterion))
		important_region_rows = []
		regions_filename = "critical_regions_" + criterion + ".bed"
		regions_header = "#chr\tstart\tend\tuniprot_id\tscore\tstrand"

		# for each domain investigate if it is of critical importance using the supplied criterion
		for dom_clinvar in domain_clinvars:
			clinvars = dom_clinvar['clinvar']
			region_score, is_region_important = FindCriticalProteinRegions.is_region_important(clinvars, criterion,
			                                                                                   min_ratio_cutoff)
			if is_region_important:
				important_region_rows.append('\t'.join(
					[dom_clinvar['bed']['chr'], str(dom_clinvar['bed']['start']), str(dom_clinvar['bed']['end']),
					 dom_clinvar['bed']['uniprot_id'],
					 "{:10.4f}".format(region_score), dom_clinvar['bed']['strand']]))

		print("Identified {} possibly critical protein regions".format(len(important_region_rows)))
		self.write_bed(regions_header, important_region_rows, regions_filename)

	def identify_domains_without_benign(self, domain_clinvars, criterion, min_ratio_cutoff):
		"""
		Identify well-known functional domains without benign variants
		to use for ACMG PM1 rule, following VarSome definition found at:
		https://varsome.com/about/resources/acmg-implementation/#pm1

		domain_clinvars : list of dict
			bed annotation and clinvar entries per domain
		criterion : str
			criterion to decide for region importance
		min_ratio_cutoff: float
			minimum ratio cutoff to decide for a region, if criterion is using the ratio of B_LB over P_LP

		Returns
		-------
		None
		"""
		print("Identify well-known functional domains without benign mutations")
		domains_no_benign = []
		regions_filename = "critical_regions_" + criterion + "_no_benign.bed"
		regions_header = "# chr\tstart\tend\tuniprot_id\tscore\tstrand"

		# for each domain investigate if its pathogenic probability >= 0.5 and there no benign mutations
		for dom_clinvar in domain_clinvars:
			clinvars = dom_clinvar['clinvar']
			num_B_LB, num_P_LP, num_VUS = FindCriticalProteinRegions.calculate_pathogenic_class_counts(clinvars)
			region_score, is_region_important = FindCriticalProteinRegions.is_region_important(clinvars, criterion,
			                                                                                   min_ratio_cutoff)
			if is_region_important and num_B_LB == 0:
				domains_no_benign.append('\t'.join(
					[dom_clinvar['bed']['chr'], str(dom_clinvar['bed']['start']), str(dom_clinvar['bed']['end']),
					 dom_clinvar['bed']['uniprot_id'], "{:10.4f}".format(region_score), dom_clinvar['bed']['strand']]))
		print("Identified {} well-known functional domains without a benign mutation".format(len(domains_no_benign)))
		self.write_bed(regions_header, domains_no_benign, regions_filename)

	def summarize_significance(self, domain_clinvars):
		"""
		Summarize signifance per review status

		Parameters
		----------
		domain_clinvars: list of dict
			bed annotation and clinvar entries per domain

		Returns
		-------
		dict of str: dict of int: int
			significance summary aggregated for review quality
		"""
		print("Summarize clinvars per significance")
		significance_summary = {}
		for dom_clinvar in domain_clinvars:
			clinvars = dom_clinvar['clinvar']
			# summarize over significance
			for clinvar in clinvars:
				current_significance = clinvar['CLNSIG'][0]
				if current_significance not in significance_summary:
					significance_summary[current_significance] = {}
				current_review_status = clinvar['CLNREVSTAT']
				if current_review_status not in significance_summary[current_significance]:
					significance_summary[current_significance][current_review_status] = 1
				else:
					significance_summary[current_significance][current_review_status] += 1
		# print("Significance summary: {}".format(significance_summary))
		return significance_summary

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
		print("Filter ClinVar by review stars")
		filtered_clinvars = []
		for clinvar in clinvars:
			if clinvar['CLNREVSTAT'] >= min_review_stars:
				filtered_clinvars.append(clinvar)
		print("Filtered in {} ClinVar entries".format(len(filtered_clinvars)))
		return filtered_clinvars

	def aggregate_clinvars_per_domain(self, min_review_stars):
		"""
		Aggregate clinvar entries per domain

		Parameters
		----------
		min_review_stars : int
			number of minimum review stars

		Returns
		-------
		list of dict
			filtered clinvars per domain
		"""
		print("Aggregate ClinVar entries per domain")
		### ### ###
		# For each UniProt domain
		# a) extract clinvar entries found on the same genomic region and strand
		# b) filter these extracted clinvar entries by clinvar status stars
		### ### ###
		domain_clinvars = []
		for row in self.uniprot_domains:
			domain_annotation = self.parse_bed_row(str(row))
			# a) extract clinvar entries for this domain range
			clinvars = self.extract_clinvars(domain_annotation)
			if clinvars:
				print("dom: {} {} {}".format(domain_annotation['chr'], domain_annotation['start'], domain_annotation['end']))
				# b) filter extracted entries by review quality
				filtered_clinvars = self.quality_filter_clinvars(clinvars, min_review_stars)
				if len(filtered_clinvars) > 0:
					domain_clinvars.append({'bed': domain_annotation, 'clinvar': filtered_clinvars})
		print("---")
		return domain_clinvars
