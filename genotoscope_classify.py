import os, sys
import logging
import argparse
from genotoscope.variants.ClassifyVariants import ClassifyVariants
from genotoscope.utils import prepare_parameters, argparse_float_range

sys.path.insert(0, os.path.abspath(".."))

# ### ### ### ### ### ### ### ### ### ### ### ### #
#  GSvar + workflow --> variants classification   #
# ### ### ### ### ### ### ### ### ### ### ### ### #


# input = data_path where you have placed:
# 1. folder with gsvar files
# 2. hearing loss relevant genes inheritance from DOI:  10.1038/s41436-019-0487-0
# 3. clinvar.vcf.gz as downloaded from https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh37/
# (please download also the clinvar.vcf.gz.tbi that contains the tabix indexes)
# 4. clinvar_stars.tsv as downloaded from https://www.ncbi.nlm.nih.gov/clinvar/docs/review_status/
# 5. UniProt domains and repeats
# 7. GnomAD pLoF
# 8. hearing loss relevant well characterized functional regions of proteins (for PM1) from DOI: 10.1002/humu.23630
# 8. hearing loss relevant transcripts from DOI: https://pubmed.ncbi.nlm.nih.gov/30096381/
# 9. complete HGNC genes text file augmented with strand information. Original downloaded from https://www.genenames.org/download/statistics-and-files/

# ### ### ### ### ### ### ### ### ### ### ### ### #
#      Set up GenOtoScope classification run      #
# ### ### ### ### ### ### ### ### ### ### ### ### #

parser = argparse.ArgumentParser(description="Classify variants using GenOtoScope.")
parser.add_argument("--analysis_root",
                    help="Absolute path to directory under which GenOtoScope data is placed")
parser.add_argument("--input_variants", help="Absolute path to patients variants or single file path with variants")
parser.add_argument("--output_path", help="Absolute path to output directory")
parser.add_argument("--settings_file", help="Settings file (.yaml) path")
parser.add_argument("--filter_HL_genes", type=int,
                    help="Apply classification on variants affecting ONLY known HL genes (1) or all presented genes (0)",
                    choices=range(0, 2))
parser.add_argument("--min_pathogenicity", type=argparse_float_range(0.0, 1.0, -1),
                    help="Minimum pathogenicity probability, in the range of [0.0,1.0], for variants classified with class 3 (ACMG), along with all variants classified with 4/5 class, to be extracted in a separate 'summary' file. Choose -1 to deactivate extraction step.")
parser.add_argument("--logging", choices=["INFO", "DEBUG"], help="Logging level")
parser.add_argument("--threads", type=int, help="Number of threads")
args = parser.parse_args()

# read up input arguments
genotoscope_root = args.analysis_root
variants_path = args.input_variants
output_path = args.output_path
yaml_file_path = args.settings_file
filter_hl_genes = args.filter_HL_genes
min_pathogenicity = args.min_pathogenicity
logging_level = args.logging
num_threads = args.threads

### ### #### #### ### ###
# check input arguments #
### ### #### #### ### ###
try:
	assert os.path.isdir(genotoscope_root)
except AssertionError:
	logging.error("Please provide an existing GenOtoScope root directory.", exc_info=True)
	sys.exit(os.EX_NOINPUT)

try:
	assert os.path.isdir(variants_path) or os.path.isfile(variants_path)
except AssertionError:
	logging.error("Please provide an existing directory or existing single file as input variants path.", exc_info=True)
	sys.exit(os.EX_NOINPUT)

try:
	assert os.path.isdir(os.path.dirname(output_path))
except AssertionError:
	logging.error("Please provide an output folder, whose parent folder exists.", exc_info=True)
	sys.exit(os.EX_NOINPUT)

try:
	assert os.path.isfile(yaml_file_path) and os.path.splitext(yaml_file_path)[1] == ".yaml"
except AssertionError:
	logging.error("Please provide an existing YAML file", exc_info=True)
	sys.exit(os.EX_NOINPUT)

if filter_hl_genes == 0:
	filter_hl_genes = False
else:
	filter_hl_genes = True
# prepare parameters for run
data_settings, data_path, temp_path, critical_prot_version = prepare_parameters(yaml_file_path, genotoscope_root)

### #### #### #### #### #### ###
#   Run GenOtoScope Classify   #
### #### #### #### #### #### ###
logging.info("==== ====")
GenOtoScopeClassify = ClassifyVariants(data_path, variants_path, data_settings["genes_inheritance_file"],
                                       data_settings["thresholds_inheritance_file"],
                                       data_settings["clinvar_root"], data_settings["clinvar_version"],
                                       data_settings["clinvar_file"],
                                       data_settings["clinvar_stars_file"], data_settings["min_review_stars"],
                                       data_settings["beds_root"], critical_prot_version,
                                       data_settings["critical_prot_regions_file"],
                                       data_settings["critical_regions_no_benign_file"],
                                       data_settings["gnomAD_version"],
                                       data_settings["clinical_exons_file"], data_settings["uniprot_root"],
                                       data_settings["uniprot_version"], data_settings["uniprot_domains_file"],
                                       data_settings["uniprot_repeat_file"], data_settings["hearing_loss_annot_root"],
                                       data_settings["hl_relevant_exons_file"],
                                       data_settings["hl_pm1_regions_file"],
                                       data_settings["hugo_genes_file"], data_settings["pheno_high_af_variants_file"],
                                       temp_path,
                                       output_path, filter_hl_genes, min_pathogenicity, num_threads, logging,
                                       logging_level)
if os.path.isdir(variants_path):
	# run genotoscope for each GSvar annotated variants file (one file per patient)
	GenOtoScopeClassify.run_multiple_samples()
elif os.path.isfile(variants_path):
	# run genotoscope for a single GSvar annotated variant file
	GenOtoScopeClassify.run_single_sample()
logging.info("=== * ===")
logging.info("== *** ==")
