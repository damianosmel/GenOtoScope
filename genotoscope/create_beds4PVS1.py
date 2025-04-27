import os, sys

sys.path.insert(0, os.path.abspath(".."))

from genotoscope.beds.FindCriticalProteinRegions import FindCriticalProteinRegions
from genotoscope.beds.FindClinicalSignificantExons import FindClinicalSignificantExons
from os.path import join

# ### ### ### ### ### ### ### ### ### ### ### ### #
#  Create annotation tracks for autoPVS1 rules    #
# ### ### ### ### ### ### ### ### ### ### ### ### #

### ### ###
# input and output paths
### ### ###

### ### ###
#  input  #
### ### ###
### ### ###
# ClinVar.vcf + UniProt domain annotations => critical regions of protein
### ### ###

# local machine #
data_path = "/home/damian/Documents/L3S/projects/hearing4all/human_genetics/genotoscope_data"
# server #
# data_path = "/home/melidis/hearing_genomics/data"

### ClinVar ###
clinvar_path = "ClinVar"
clinvar_version = "clinvar_apr_2025"
clinvar_file = "clinvar_20250409.vcf.gz"
clinvar_stars_file = "clinvar_review_stars.tsv"

### UniProt ###
uniprot_dir = join(data_path, "uniprot")
uniprot_version = "uniprot_feb_2025"
uniprot_dir = join(uniprot_dir,uniprot_version)
uniprot_domains_path = join(uniprot_dir, "UP000005640_9606_domain.bed")

### Hugo genes ###
hugo_genes_directory = join("misc", "hgnc_genes_info")
hugo_genes_directory = join(hugo_genes_directory, "apr_2025")
hugo_genes_file = "hgnc_complete_set_strand.txt"
hugo_genes_path = join(hugo_genes_directory,hugo_genes_file)

### ### ###
# Parameters
### ### ###
min_review_stars = 2  # ClinVar minimum review stars

### ### ###
# output  #
### ### ###

output_path = 'annotation_beds'
### ### ###
# Find critical protein region
### ### ###
critical_protein_regions = FindCriticalProteinRegions(data_path, clinvar_path, clinvar_file, clinvar_stars_file,
                                                     clinvar_version,uniprot_domains_path, uniprot_version, hugo_genes_path,
                                                     output_path)
critical_protein_regions.run(min_review_stars)

'''
###
# gnomAD.exomes.vcf + all pLoF variants => clinical significant exons
###

gnomAD_root = join(data_path, "gnomAD")
gnomAD_version = "gnomAD_v2.1"
gnomAD_path = join(gnomAD_root, gnomAD_version)
### GRCh37 all at server ###
gnomAD_exomes_file = join(gnomAD_path, "gnomad.exomes.r2.1.1.sites.vcf.gz")
gnomAD_pLoF_file = join(gnomAD_path, "gnomad.v2.1.1.all_lofs.txt")

### test at local machine ###
# gnomAD_exomes_file = join(gnomAD_path, "gnomad.exomes.r2.1.1.sites.1.vcf.gz")
# gnomAD_pLoF_file = join(gnomAD_path, "gnomad.v2.1.1.all_lofs.txt")
# gnomAD_exomes_file = join(gnomAD_path, "gnomad.exomes.r2.1.1.sites.1.10000.vcf.gz")
# gnomAD_pLoF_file = join(gnomAD_path, "pLoF_30.txt")
# gnomAD_exomes_file = join(gnomAD_path, "gnomad.exomes.r2.1.1.sites.1.1000.vcf.gz")


FindClinicalSignificantExons = FindClinicalSignificantExons(data_path, gnomAD_version, gnomAD_exomes_file, output_path)
FindClinicalSignificantExons.run(gnomAD_pLoF_file)
'''