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
# local machine #
# data_path = "/home/damian/Documents/L3S/projects/hearing4all/human_genetics/data"
# server #
data_path = "/home/melidis/hearing_genomics/data"
output_path = "beds4autoPVS1"

'''
### ### ###
# ClinVar.vcf + UniProt domain annotations => critical regions of protein
### ### ###

### ClinVar ###
clinvar_path = "ClinVar"
clinvar_version = "clinvar_mar_2021"
clinvar_file = "clinvar_20210323.vcf.gz"
clinvar_stars_file = "clinvar_review_stars.tsv"

### UniProt ###
uniprot_dir = join(data_path, "uniprot_annotations")
uniprot_version = "uniprot_feb_2021"
uniprot_dir = join(uniprot_dir,uniprot_version)
uniprot_domains_path = join(uniprot_dir, "UP000005640_9606_domain_hg19.bed")

### Hugo genes ###
hugo_genes_file = "hgnc_complete_set_strand.txt"

### ### ###
# Parameters
### ### ###
min_review_stars = 2  # ClinVar minimum review stars

### ### ###
# Find critical protein region
### ### ###
CriticalProteinsRegions = FindCriticalProteinRegions(data_path, clinvar_path, clinvar_file, clinvar_stars_file,
                                                     clinvar_version,uniprot_domains_path, uniprot_version, hugo_genes_file,
                                                     output_path)
CriticalProteinsRegions.run(min_review_stars)
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
