from flask import Flask,request
from flask_bootstrap import Bootstrap
from genotoscope.config_web import Config
from logging.handlers import RotatingFileHandler
import logging
import os
from redis import Redis
import rq


website = Flask(__name__)
website.config.from_object(Config)
bootstrap = Bootstrap(website)
redis_connection = Redis.from_url(website.config['REDIS_URL'])
task_queue = rq.Queue('genotoscope-tasks', connection=redis_connection)
# website data folder variables
downloads_path = website.config['DOWNLOADS_URL']
uploads_path = website.config['UPLOADS_URL']
tmp_path = website.config['TEMP_FOLDER_URL']
web_data_host = website.config["WEB_DATA_HOST"]
upload_extensions = website.config['UPLOAD_EXTENSIONS']

# docker variables
docker_root_host = website.config["DOCKER_ROOT_HOST"]
docker_root = website.config["DOCKER_ROOT"]
docker_image_id = website.config["DOCKER_IMAGE_ID"]

# genotoscope variables
genotoscope_web_out_path = website.config['GENOTOSCOPE_WEB_OUT']
genotoscope_settings_file_path = website.config['GENOTOSCOPE_SETTINGS']
genotoscope_analysis_root = website.config['GENOTOSCOPE_ANALYSIS_ROOT']


if not website.debug:
	data_path = website.config['LOGS_DIR_ROOT']
	logging_path = os.path.join(data_path,'logs')
	if not os.path.exists(logging_path):
		os.mkdir(logging_path)
	file_handler = RotatingFileHandler(os.path.join(logging_path,'genotoscope_web.log'), maxBytes=10240,
                                       backupCount=10)
	file_handler.setFormatter(logging.Formatter(
        '%(asctime)s %(levelname)s: %(message)s [in %(pathname)s:%(lineno)d]'))
	file_handler.setLevel(logging.INFO)
	website.logger.addHandler(file_handler)

	website.logger.setLevel(logging.INFO)

	website.logger.info('#### #### #### #### #### #### #### #### ####')
	website.logger.info('GenOtoScope Website Startup')

from genotoscope.website import routes, errors, tasks
