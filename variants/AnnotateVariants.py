from genotoscope.utils import create_dir, split_name_from_extension, sec2hour_min_sec
from os import scandir
from os.path import join, basename, dirname, normpath, isfile
import logging
from datetime import datetime
from timeit import default_timer as timer
import subprocess


class AnnotateVariants:
	"""
	Class to annotate VCF file(s) using megSAP application. After annotation, run again megSAP application to convert to GSvar file(s)
	megSAP: https://github.com/imgag/megSAP
	"""

	def __init__(self, variants_path, output_path, temp_path, docker_root, docker_root_host, docker_image_id,
	             logging_root, logging_level):
		"""
		AnnotateVariants constructor

		Parameters
		----------

		variants_path : str
			input VCFs path
		output_path : str
			output absolute path
		temp_path : str
			absolute path of temporary files

		docker_root : str
			absolute path to (mounted) docker root
		docker_root_host : str
			absolute path to docker root in host
		docker_image_id : str
			configured image id
		logging_root : logging
			logging root initialized by genotoscope_annotate script
		logging_level : str
			logging level
		Returns
		-------
		None
		"""
		# self.analysis_root = analysis_root
		self.variants_path = variants_path
		self.out_path = output_path
		self.temp_path = temp_path
		self.docker_root = docker_root
		self.docker_root_host = docker_root_host
		self.docker_image_id = docker_image_id

		# create output and logging folders
		create_dir(output_path)
		create_dir(temp_path)

		### Logging ###
		self.init_logger(logging_root, "GenOtoScope_Annotate", logging_level)

	def init_logger(self, logging_root, logger_name, logging_level):
		"""
		Initiliaze logger
		Credits: Credits: https://stackoverflow.com/questions/20240464/python-logging-file-is-not-working-when-using-logging-basicconfig

		Parameters
		----------
		logging_root : logging
			logging root object
		logger_name : str
			logger name
		logging_level : str
			logging level

		Returns
		-------
		None
		"""
		self.logger = logging_root.getLogger(logger_name)
		self.logger.setLevel(logging_level)

		# create file handler
		fh = logging_root.FileHandler(join(self.out_path,
		                                   'genotoscope_annotate_' + datetime.today().strftime("%d_%m_%Y") + '.log'))
		fh.setLevel(logging_level)
		# create console handler
		ch = logging.StreamHandler()
		ch.setLevel(logging.ERROR)

		# create format
		formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
		fh.setFormatter(formatter)
		ch.setFormatter(formatter)
		# add handler to the logger
		self.logger.addHandler(fh)
		self.logger.addHandler(ch)

	def annotate_convert2gsvar(self, input_vcf_path):
		"""
		Run megSAP docker image to
		a) annotate vcf file
		b) convert it to GSvar file

		Parameters
		----------
		input_vcf_path : str
			path to input vcf file

		Returns
		-------
		None
		"""
		self.logger.info("Run megSAP to annotate vcf file and convert it to GSvar")
		### ### ### ### ### ### ### ### ### ### ###
		#                  megSAP                 #
		### ### ### ### ### ### ### ### ### ### ###

		### ### ### ### ### ### ### ### ### ### ###
		# 1) VeP annotate vcf file
		# save on annotated vcf and log in tmp directory
		### ### ### ### ### ### ### ### ### ### ###

		# construct the parts of the subprocess command
		input_vcf_parent = dirname(normpath(input_vcf_path))
		vcf_file_name = basename(input_vcf_path)
		[vcf_name, _] = split_name_from_extension(input_vcf_path)
		vep_annotate_command = 'php /megSAP/src/NGS/an_vep.php -in ' + self.docker_root + '/host_data_in/' + \
		                       vcf_file_name + ' -out ' + self.docker_root + '/host_data_out/genotoscope_annotate_temp/' + vcf_name + '_annot.vcf' + ' -no_splice >& ' + \
		                       self.docker_root + '/host_data_out/genotoscope_annotate_temp/' + vcf_name + '_annot.log'

		self.logger.info("Execute command: {}".format(
			' '.join(['docker', 'run', '--volume', self.docker_root_host + ':' + self.docker_root,
			          '--volume', input_vcf_parent + ':' + self.docker_root + '/host_data_in', '--volume',
			          self.out_path + ':' + self.docker_root + "/host_data_out",
			          self.docker_image_id,
			          '-c', vep_annotate_command])))
		# run the subprocess command
		an_vep_return_code = subprocess.call(
			['docker', 'run', '--volume', self.docker_root_host + ':' + self.docker_root,
			 '--volume', input_vcf_parent + ':' + self.docker_root + '/host_data_in', '--volume',
			 self.out_path + ':' + self.docker_root + "/host_data_out",
			 self.docker_image_id,
			 '-c', vep_annotate_command])
		self.logger.info("an_vep return code: {}".format(an_vep_return_code))

		if isfile(join(self.temp_path, vcf_name + '_annot.vcf')):
			### ### ### ### ### ### ### ### ### ### ###
			# 2) convert to annotated vcf to gsvar file
			# save gsvar file in output directory
			### ### ### ### ### ### ### ### ### ### ###

			# construct the parts of the subprocess command
			vcf2gsvar_command = 'php /megSAP/src/NGS/vcf2gsvar.php -in ' + self.docker_root + '/host_data_out/genotoscope_annotate_temp/' + vcf_name + '_annot.vcf' + ' -out ' + self.docker_root + '/host_data_out/' + vcf_name + '_annot.GSvar' + ' >& /megSAP/host_data_out/genotoscope_annotate_temp/' + vcf_name + '_gsvar.log'
			# run the subprocess command
			self.logger.info("Execute command: {}".format(
				' '.join(['docker', 'run', '--volume', self.docker_root_host + ':' + self.docker_root,
				          '--volume', self.out_path + ':' + self.docker_root + '/host_data_out',
				          self.docker_image_id,
				          '-c', vcf2gsvar_command])))
			vcf2gsvar_return_code = subprocess.call(
				['docker', 'run', '--volume', self.docker_root_host + ':' + self.docker_root,
				 '--volume', self.out_path + ':' + self.docker_root + '/host_data_out',
				 self.docker_image_id,
				 '-c', vcf2gsvar_command])

			self.logger.info("vcf2gsvar return code: {}".format(vcf2gsvar_return_code))
			self.logger.info("Saving annotated variants GSvar %s", join(self.out_path, vcf_name + '_annot.GSvar'))
		else:
			self.logger.info("Error: an_vep did not produce annotated variants vcf file => terminate annotation")

	def annotate_single_sample(self):
		"""
		Annotate single vcf file

		Parameters
		----------
		None

		Returns
		-------
		None
		"""
		self.logger.info("Process single file by annotating VCF and converting it to GSvar file")

		self.logger.info("==== ==== Start ==== ====")
		file_name, file_extension = split_name_from_extension(self.variants_path)
		if file_extension == ".vcf":
			self.logger.info("######################### ~ #########################")
			self.logger.info("### Start:  %s ###", basename(self.variants_path))
			time_start = timer()
			self.annotate_convert2gsvar(self.variants_path)
			time_end = timer()
			self.logger.info(
				"Elapsed time to process GSvar: {}".format(sec2hour_min_sec(time_end - time_start)))
			self.logger.info("### End:  %s  ###", basename(self.variants_path))
		self.logger.info("==== ====  End  ==== ====")

	def annotate_multiple_samples(self):
		"""
		Run annotation and convertion to gsvar for each sample/VCF file

		Returns
		-------
		None
		"""
		self.logger.info("Process multiple VCF files by annotating and converting to GSvar files")
		excluded_file_names = []

		self.logger.info("==== ==== Start ==== ====")
		with scandir(self.variants_path) as vcf_dir:
			for content in vcf_dir:
				file_name, file_extension = split_name_from_extension(content.name)
				if file_extension == ".vcf" and content.is_file():
					self.logger.info("######################### ~ #########################")
					self.logger.info("### Start:  %s ###", content.name)
					time_start = timer()
					self.annotate_convert2gsvar(join(self.variants_path, content.name))
					time_end = timer()
					self.logger.info(
						"Elapsed time to process GSvar: {}".format(sec2hour_min_sec(time_end - time_start)))
					self.logger.info("### End:  %s  ###", content.name)
				else:
					excluded_file_names.append(content.name)
		self.logger.info("==== ====  End  ==== ====")

		if len(excluded_file_names) > 0:
			self.logger.info("=== Excluded ===")
			self.logger.info("# Excluding following contents as they don't have vcf extension:")
			self.logger.info("\n".join(excluded_file_names))
