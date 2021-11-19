from pybedtools import BedTool
from pyensembl import EnsemblRelease

from os.path import join
from collections import Counter
import matplotlib.pyplot as plt
import numpy as np

from pandas import read_excel


# ### ### ### ### ### ### ### ### ### ### ### #
# test created annotations by intersecting    #
# critical regions to significant exons       #
# ### ### ### ### ### ### ### ### ### ### ### #


def intersect_exons_with_protein_regions(critical_regions_file, significant_exons_file, ensembl_data):
	"""
	Intersect exons with protein regions to get number of intersecting exons per region
	for example if there are 10 regions without intersecting exons, 1 region with one intersecting exon,
	1 region with 2 intersecting exons, 1 region with 3 intersecting exons
	=> return Counter({0: 10, 1:2, 3:1})

	Parameters
	----------
	critical_regions_file : str
		critical regions filename
	significant_exons_file : str
		significant exons filename
	ensembl_data : pyensembl.EnsemblRelease
		ensemble data object

	Returns
	-------
	Counter
		regions group by number of interesting domains
	"""
	critical_regions = BedTool(critical_regions_file)
	significant_exons = BedTool(significant_exons_file)

	# for each exon examine overlap with each region
	# to find intersecting exons
	intersect_regions2exons = {}
	for exon in significant_exons:
		for region in critical_regions:
			if region[0].split('chr')[1] == exon[0]:  # same chr
				uniprot_id = region[4]
				if int(exon[1]) >= int(region[1]) and int(exon[2]) <= int(region[2]):
					# examine if exon and critical region are on the same strand
					region_strand = region[3]
					exon_strand = ensembl_data.exon_by_id(exon[3]).to_dict()["strand"]
					if region_strand == exon_strand:
						# print("match")
						if uniprot_id not in intersect_regions2exons:
							intersect_regions2exons[uniprot_id] = [exon]
						else:
							intersect_regions2exons[uniprot_id].append(exon)

	print("### Intersecting regions with exons ###")
	num_intersecting_exons_per_region = []
	for region_id, exons in intersect_regions2exons.items():
		print('critical region of uniprot id: {} contains intersecting exons: {}'.format(region_id, len(exons)))
		num_intersecting_exons_per_region.append(len(exons))
	# add the regions with no intersecting exons
	num_critical_regions_all = critical_regions.count()
	assert num_critical_regions_all >= len(list(
		intersect_regions2exons.keys())), "AssertionError: Regions with intersecting exons can't be more than all available regions"
	num_regions_no_intersecting_exons = num_critical_regions_all - len(list(intersect_regions2exons.keys()))
	regions_no_intersecting_exons_percentage = (num_regions_no_intersecting_exons / num_critical_regions_all) * 100
	print("Regions without intersecting exons: {} out of all regions: {} (percentage= {:.4f}%)".format(
		num_regions_no_intersecting_exons, num_critical_regions_all, regions_no_intersecting_exons_percentage))

	num_intersecting_exons_per_region += [0] * num_regions_no_intersecting_exons
	print("### ### ### ### ###")
	return Counter(num_intersecting_exons_per_region)


def plot_intersecting_bar(regions_intersecting_exons, regions_name, exons_name, data_path):
	"""
	Plot cumulative bar plot of number of intersecting exons with protein regions

	Parameters
	----------
	regions_intersecting_exons:
		Counter object containing num of intersecting exons per region
	regions_name : str
		regions type information
	exons_name : str
		exons type information
	data_path : str
		data root path

	Returns
	-------
	None
	"""
	print("Plot cumulative bar plot of number of intersecting exons with protein regions")
	fig = plt.figure()
	# sort counter by number of intersecting exons
	num_intersecting_exons_sorted, num_regions_sorted = [], []
	for num_intersecting_exons, num_regions in sorted(regions_intersecting_exons.items()):
		num_intersecting_exons_sorted.append(num_intersecting_exons)
		num_regions_sorted.append(num_regions)
	# calculate the cumulative density of the number of regions
	cumulative_density_num_regions = np.cumsum(num_regions_sorted) / np.sum(num_regions_sorted)
	plt.bar(num_intersecting_exons_sorted, cumulative_density_num_regions, align='center', color='tab:green',
	        edgecolor='k', alpha=0.8)

	if regions_name == "important_trunc_regions_no_benign_atleast5pathogenic":
		regions_title = "protein critical regions (no B & >= 5 P/LP)"
	elif regions_name == "important_trunc_regions_pathogenic_probability":
		regions_title = "protein critical regions (P_pathogenic >= 0.5)"
	exons_title = exons_name.replace("_", " ")
	title = "Cumulative histogram of intersecting " + exons_title + " per " + regions_title
	plt.title(title)
	plt.xlabel("Num. intersecting exons per region", fontsize=14)
	plt.ylabel("Percentage of regions", fontsize=14)
	if regions_name == "important_trunc_regions_pathogenic_probability":
		plt.xticks(num_intersecting_exons_sorted, num_intersecting_exons_sorted, rotation=90)
	# plt.show()
	bar_name = exons_name + "_in_" + regions_name + "_cumulative_density" + ".png"
	fig.savefig(join(data_path, bar_name), bbox_inches='tight', dpi=600)


### data file paths ###
# input
data_path = "/home/damian/Documents/L3S/projects/hearing4all/human_genetics/GenoPipeAutomate/shared_files/beds4autoPVS1"
# autoPVS1
# critical_regions_file = join(data_path, 'important_trunc_regions_no_benign_atleast5pathogenic.bed')
# P_pathogenic >= 0.5
critical_regions_file = join(data_path, 'critical_regions_pathogenic_probability.bed')
significant_exons_PASS_file = join(data_path, 'clinical_significant_exons_PASS.bed')
significant_exons_file = join(data_path,'clinical_significant_exons.bed')
pheno_relevant_transcripts_file = join("/home/damian/Documents/L3S/projects/hearing4all/human_genetics/data",
                                       'gene_transcripts.xlsx')
# output
output_path = "/home/damian/Documents/L3S/projects/hearing4all/human_genetics/data/beds4autoPVS1"
# ensembl object
ensembl_data = EnsemblRelease(75)

### Intersect called critical regions of proteins with called clinical significant exons ###
# num_intersecting_regions_per_exon = intersect_exons_with_protein_regions(critical_regions_file,significant_exons_file,ensembl_data)
# plot_intersecting_bar(num_intersecting_regions_per_exon, "important_trunc_regions_pathogenic_probability", "clinical_significant_exons",output_path)

### Intersect critical regions of proteins with phenotype relevant transcripts ###

# todo: remove when no needed anymore
"""
# 1) extract phenotype relevant exons from transcripts information
pheno_relevant_exons = extract_pheno_relevant_exons(pheno_relevant_transcripts_file, ensembl_data)
exons_header = "#chr\tstart\tend\texon\ttranscripts\tstrand"
pheno_relevant_exons_bed_name = "pheno_relevant_exons.bed"
write_bed(exons_header, pheno_relevant_exons, pheno_relevant_exons_bed_name, output_path)
"""

'''
# 2) intersect phenotype relevant exons with critical regions (no benign at least 5 pathogenic)
pheno_relevant_exons_bed_name = "pheno_relevant_exons.bed"
pheno_relevant_exons_file = join(output_path, pheno_relevant_exons_bed_name)
num_intersecting_regions_per_exon = intersect_exons_with_protein_regions(critical_regions_file,
                                                                         pheno_relevant_exons_file, ensembl_data)
plot_intersecting_bar(num_intersecting_regions_per_exon, "important_trunc_regions_no_benign_atleast5pathogenic",
                      "phenotype_relevant_exons", output_path)
'''
'''
### After change on quality filter, examine if clinical significant exons have changed ###
significant_exons_DP_QUAL = BedTool(significant_exons_file)
significant_exons_PASS = BedTool(significant_exons_PASS_file)
print("num of significant exons: {}".format(significant_exons_DP_QUAL.count()))
print("num of signigicant exons PASS: {}".format(significant_exons_PASS.count()))
#only_old_type = significant_exons.subtract(significant_exons_PASS)
exons_DP_QUAL_only = significant_exons_DP_QUAL - significant_exons_PASS
print("num of significant exons for old type of quality filter: {}".format(exons_DP_QUAL_only.count()))
'''