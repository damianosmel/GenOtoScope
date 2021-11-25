import argparse
import logging
import os,sys
from os.path import join

from genotoscope.variants.AnnotateVariants import AnnotateVariants

sys.path.insert(0, os.path.abspath(".."))
# ### ### ### ### ### ### ### ### ### ### ### ### #
#  VCF + megSAP --> annotated GSvar files         #
# ### ### ### ### ### ### ### ### ### ### ### ### #

# input = data_path where you have placed:
# folder with gsvar files

# docker parameters
# docker_root = name of mounted root of the docker container
# docker_root_host = absolute host path to the docker container
# docker_image_id = image id of the configured docker container



# ### ### ### ### ### ### ### ### ### ### ### ### #
#   Set up GenOtoScope variant annotation run     #
# ### ### ### ### ### ### ### ### ### ### ### ### #

parser = argparse.ArgumentParser(description="Annotate variants using megSAP.")
parser.add_argument("--input_variants", help="Absolute path to patients variants or single file path with variants")
parser.add_argument("--output_path", help="Absolute path to output directory (it will be created, if not existing)")
parser.add_argument("--docker_root", help="Name of mounted docker root directory")
parser.add_argument("--docker_root_host", help="Absolute path to docker directory in host")
parser.add_argument("--docker_image_id", help="Configured docker image ID")
parser.add_argument("--logging", choices=["INFO"], help="Logging level")
args = parser.parse_args()

# read up input arguments
variants_path = args.input_variants
output_path = args.output_path
docker_root_path = args.docker_root
docker_root_host_path = args.docker_root_host
docker_image_id = args.docker_image_id
logging_level = args.logging

### ### #### #### ### ###
# check input arguments #
### ### #### #### ### ###

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
	assert os.path.isdir(os.path.dirname(docker_root_host_path))
except AssertionError:
	logging.error("Please provide an docker root folder that exists in the host filesystem.",exc_info=True)
	sys.exit(os.EX_NOINPUT)

### ### ###
# configure data path and logging folder
### ### ###
# create the central data path
#data_path = join(root_dir, "data")
# logging folder
#logging_path = join(analysis_root, "genotoscope_annotate_logs")
temp_path = join(output_path,"genotoscope_annotate_temp")

### #### #### #### #### #### ###
#   Run GenOtoScope Annotate   #
### #### #### #### #### #### ###
logging.info("==== ====")
GenOtoScopeAnnotate = AnnotateVariants(variants_path,output_path,temp_path,docker_root_path,docker_root_host_path,docker_image_id,logging,logging_level)

if os.path.isdir(variants_path):
	# annotate each VCF file, one file per patient, contained in the input path
	GenOtoScopeAnnotate.annotate_multiple_samples()
elif os.path.isfile(variants_path):
	# annotate a single VCF file
	GenOtoScopeAnnotate.annotate_single_sample()
logging.info("=== * ===")
logging.info("== *** ==")