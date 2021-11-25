import os
from os.path import join
import configparser


class Config(object):
	### ### ### ### ### ### ###
	# parse input config file #
	### ### ### ### ### ### ###
	config = configparser.ConfigParser()
	config_file_dir, _ = os.path.split(os.path.abspath(__file__))
	# local ubuntu box
	# config.read(join(config_file_dir,'website_settings.ini'))
	# server
	config.read(join(config_file_dir, 'website_settings_server.ini'))

	### ### ### ##### ### ### ###
	# set up website parameters #
	### ### ### ##### ### ### ###
	SECRET_KEY = os.environ.get('SECRET_KEY') or config['APP_INFO']['secret_key']
	REDIS_URL = os.environ.get('REDIS_URL') or config['APP_INFO']['redis_url']
	# print(REDIS_URL)
	DOWNLOADS_URL = config['APP_INFO']['downloads_url']
	UPLOADS_URL = config['APP_INFO']['uploads_url']
	TEMP_FOLDER_URL = config['APP_INFO']['temp_folder_url']
	WEB_DATA_HOST = config['APP_INFO']['web_data_host']
	LOGS_DIR_ROOT = config['APP_INFO']['logs_root_dir']

	# docker variables
	DOCKER_ROOT_HOST = config['APP_INFO']['docker_root_host']
	DOCKER_ROOT = config['APP_INFO']['docker_root']
	DOCKER_IMAGE_ID = config['APP_INFO']['docker_image_id']

	# genotoscope variables
	GENOTOSCOPE_WEB_OUT = config['APP_INFO']['genotoscope_web_out']
	GENOTOSCOPE_SETTINGS = config['APP_INFO']['genotoscope_settings']
	GENOTOSCOPE_ANALYSIS_ROOT = config['APP_INFO']['genotoscope_analysis_root']

	UPLOAD_EXTENSIONS = ['vcf']
