import time
from rq import get_current_job
import subprocess
from os.path import join, split, abspath, isfile
from os import pardir
from genotoscope.utils import split_name_from_extension
from genotoscope.website import website


def example(args_dict):
	print("Task started")
	job = get_current_job()

	for i in range(args_dict["seconds"]):
		print(args_dict["message"])
		time.sleep(1)
		job.meta['progress'] = 100.0 * i / args_dict["seconds"]
		job.save_meta()
	return 'Printed message {} for {} seconds'.format(args_dict["message"], args_dict["seconds"])


def annotate_and_classify_variant(args_dict):
	website.logger.info("=== ==== ===")
	website.logger.info("Annotate&Classify: Start - Annotate and classify variant")
	website.logger.info(
		"Annotate&Classify: input variant file: {}".format(join(args_dict["uploads"], args_dict["vcf_file"])))

	### ### ### ### ### ### ### ### ### ### ###
	#                  megSAP                 #
	### ### ### ### ### ### ### ### ### ### ###

	### ### ### ### ### ### ### ### ### ### ###
	# 1) VeP annotate vcf file
	# save on annotated vcf in tmp directory
	### ### ### ### ### ### ### ### ### ### ###

	# construct the parts of the subprocess command
	[vcf_filename, _] = split_name_from_extension(args_dict["vcf_file"])
	vep_annotate_command = 'php /megSAP/src/NGS/an_vep.php -in ' + args_dict['docker_root'] + '/host_data/uploads/' + \
	                       args_dict["vcf_file"] + ' -out ' + args_dict[
		                       'docker_root'] + '/host_data/tmp/' + vcf_filename + '_annot.vcf' + ' -no_splice >& ' + \
	                       args_dict['docker_root'] + '/host_data/logs/' + vcf_filename + '_annot.log'

	# run the subprocess command
	an_vep_return_code = subprocess.call(
		['docker', 'run', '--volume', args_dict['docker_root_host'] + ':' + args_dict['docker_root'],
		 '--volume', args_dict['web_data_host'] + ':' + args_dict['docker_root'] + '/host_data',
		 args_dict['docker_image_id'],
		 '-c', vep_annotate_command])
	website.logger.info("Annotate&Classify: an_vep return code: {}".format(an_vep_return_code))

	### ### ### ### ### ### ### ### ### ### ###
	# 2) convert to gsvar, downloads
	# save gsvar file in downloads directory
	### ### ### ### ### ### ### ### ### ### ###

	if isfile(join(args_dict['tmp_folder'], vcf_filename + '_annot.vcf')):
		# construct the parts of the subprocess command
		vcf2gsvar_command = 'php /megSAP/src/NGS/vcf2gsvar.php -in ' + args_dict[
			'docker_root'] + '/host_data/tmp/' + vcf_filename + '_annot.vcf' + ' -out ' + args_dict[
			                    'docker_root'] + '/host_data/tmp/' + vcf_filename + '_annot.GSvar' + ' >& /megSAP/host_data/logs/' + vcf_filename + '_gsvar.log'

		# run the subprocess command
		vcf2gsvar_return_code = subprocess.call(
			['docker', 'run', '--volume', args_dict['docker_root_host'] + ':' + args_dict['docker_root'],
			 '--volume', args_dict['web_data_host'] + ':' + args_dict['docker_root'] + '/host_data',
			 args_dict['docker_image_id'],
			 '-c', vcf2gsvar_command])
		website.logger.info("Annotate&Classify: vcf2gsvar return code: {}".format(vcf2gsvar_return_code))

		### ### ### ### ### ### ### ### ### ### ###
		# 3) Run GenOtoScope
		# save classified variant in downloads directory
		### ### ### ### ### ### ### ### ### ### ###
		website_dir, _ = split(abspath(__file__))
		genotoscope_scripts_dir = abspath(join(website_dir, pardir))
		genotoscope_classify_script_path = join(genotoscope_scripts_dir, "genotoscope_classify.py")
		annot_gsvar_file_name = vcf_filename + '_annot.GSvar'
		classification_process = subprocess.run(["python",
		                                         genotoscope_classify_script_path,
		                                         "--analysis_root",
		                                         args_dict["genotoscope_script_analysis_root"],
		                                         "--input_variants",
		                                         join(args_dict["tmp_folder"], annot_gsvar_file_name),
		                                         "--output_path", args_dict['genotoscope_web_out'], "--settings_file",
		                                         args_dict['genotoscope_settings'],"--filter_HL_genes", "0", "--min_pathogenicity", "-1",
		                                         "--logging", "INFO", "--threads",
		                                         "1"])
		website.logger.info(
			"Annotate&Classify: Classification process returncode: {}".format(classification_process.returncode))
		classified_gsvar_file_name = annot_gsvar_file_name.split(".GSvar")[0] + "_genotoscope" + ".GSvar"
		move_result = subprocess.run(
			["mv", join(args_dict['genotoscope_web_out'], classified_gsvar_file_name), args_dict['downloads']])
		website.logger.info(
			"Annotate&Classify: Moving classified gsvar to downloads returncode: {}".format(move_result.returncode))
	else:
		website.logger.info(
			"Annotate&Classify: Error: an_vep could not produce annotated vcf file => terminate annotation and classification")
	website.logger.info("Annotate&Classify: End - Annotate and classify variant")
	website.logger.info("=== ==== ===")
	return ''
