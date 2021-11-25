from genotoscope.website import website
from flask import render_template, flash, redirect, request, Response, url_for
from werkzeug.utils import secure_filename

from genotoscope.website.forms import VariantForm

from genotoscope.website import task_queue, redis_connection, downloads_path, uploads_path, tmp_path, web_data_host, \
	docker_root_host, docker_root, docker_image_id, genotoscope_web_out_path, \
	genotoscope_settings_file_path, genotoscope_analysis_root
from genotoscope.utils import acmg_rules2dict, format_coding_splicing_column, format_clinvar_column, \
	split_name_from_extension, format_acmg_warnings, format_extra_info_columns

from rq.job import Job
from pandas import read_csv
from os.path import join, exists
import time


@website.route("/", methods=["GET"])
@website.route("/home", methods=["GET"])
def home():
	variant_form = VariantForm(posted_request=None, time_stamp=None)
	return render_template("home.html", title="Home", form=variant_form)


@website.route("/", methods=["POST"])
@website.route("/home", methods=["POST"])
def vcf_upload():
	### ### ### ### ### ### ### ### ### ###
	# create an instance of the vcf form
	# and validate updated input file
	### ### ### ### ### ### ### ### ### ###
	website.logger.info("#### #### #### #### #### ####")
	website.logger.info("VCF_UPLOAD: User uploaded file")
	current_timestamp = int(time.time())
	variant_form = VariantForm(posted_request=request, time_stamp=current_timestamp)
	if variant_form.validate_on_submit():
		website.logger.info("VCF_UPLOAD: vcf file is validated")
		# get the secure name of the uploaded vcf file
		vcf_file_name = secure_filename(variant_form.posted_request.files['variant_vcf_file'].filename)
		vcf_unique_filename = split_name_from_extension(vcf_file_name)[0] + "_time" + str(current_timestamp) + ".vcf"
		### ### ### ### ### ### ### ### ### ### ###
		# start a background process
		# to annotate and classify the variant in vcf
		### ### ### ### ### ### ### ### ### ### ###
		task_arg = {'vcf_file': vcf_unique_filename, 'uploads': uploads_path, 'downloads': downloads_path,
		            'tmp_folder': tmp_path, 'web_data_host': web_data_host,
		            'docker_root_host': docker_root_host, 'docker_root': docker_root,
		            'docker_image_id': docker_image_id,
		            'genotoscope_web_out': genotoscope_web_out_path,
		            'genotoscope_settings': genotoscope_settings_file_path,
		            'genotoscope_script_analysis_root': genotoscope_analysis_root}
		job = task_queue.enqueue("genotoscope.website.tasks." + "annotate_and_classify_variant", task_arg)

		# after processing, get the name of the resulted file
		classified_vcf_file = split_name_from_extension(vcf_unique_filename)[0] + "_annot" + "_genotoscope" + ".GSvar"

		# return a processing page where progress is posted
		return render_template('processing.html', title='Processing', job_id=job.id, vcf_file=classified_vcf_file)
	return render_template("home.html", title="Home", form=variant_form)


@website.route('/progress/<job_id>')
def progress(job_id):
	def get_status():
		job = Job.fetch(job_id, connection=redis_connection)
		return job.get_status()

	return Response('{"status":"' + get_status() + '"}', mimetype='application/json')


@website.route('/processing/<job_id>')
def processing(job_id, vcf_file):
	return render_template('processing.html', title='Processing', job_id=job_id, vcf_file=vcf_file)


@website.route('/classification/<gsvar_file>')
def classification(gsvar_file, genotoscope_results=None):
	### ### ###
	# todo: vcf_file is the input to the URL, then add the suffix and the GSvar and present to the user
	### ### ###
	website.logger.info("DISPLAY_CLASSIFICATION: render html for {}".format(join(downloads_path, gsvar_file)))
	if exists(join(downloads_path, gsvar_file)):
		# parse variant classification output from the gsvar file
		gsvar_df = read_csv(join(downloads_path, gsvar_file),
		                    sep='\t',  # field separator
		                    comment='#',  # comment
		                    header=0,  # get header
		                    # index_col=[0, 1, 2],  # set as index the chr, start stop
		                    skipinitialspace=True,
		                    skip_blank_lines=True,
		                    error_bad_lines=False,
		                    warn_bad_lines=True
		                    )
		if gsvar_df.shape[0] == 1:
			### ### ### ### ### ### ### ###
			# if the resulted classification file
			# contains a single variant row
			# => display results
			### ### ### ### ### ### ### ###
			# get the classification results of the single variant
			genotoscope_results = gsvar_df.iloc[0].to_dict()

			### ### ###
			# Format coding and splicing columns
			### ### ###
			# format coding and splicing column to be shown better on the webpage
			genotoscope_results['coding_and_splicing'] = format_coding_splicing_column(
				genotoscope_results['coding_and_splicing'])
			genotoscope_results['coding_and_splicing_refseq'] = format_coding_splicing_column(
				genotoscope_results['coding_and_splicing_refseq'])
			# save zygosity field to easier-to-reference key
			genotoscope_results['zygosity'] = gsvar_df.iloc[0, 5]
			# format ClinVar column
			genotoscope_results['ClinVar'] = format_clinvar_column(genotoscope_results['ClinVar'])
			# format extra info columns
			genotoscope_results = format_extra_info_columns(genotoscope_results)
			# format warnings from ACMG criteria examination
			genotoscope_results['acmg_warnings'] = format_acmg_warnings(genotoscope_results['ACMG_rule_flags'])

			# merge the dictionary of the assignment and the comments of the evidence-based criteria
			genotoscope_results = {**genotoscope_results, **acmg_rules2dict(genotoscope_results['ACMG_rules'],
			                                                                genotoscope_results['ACMG_rules_comment'])}
			website.logger.info("#### #### #### #### #### ####")
			return render_template("classification.html", title="ACMG Classification",
			                       genotoscope_results=genotoscope_results)
		else:
			### ### ### ### ### ### ### ### ###
			# resulted classification file
			# does not contain any variant row
			# => display error page
			### ### ### ### ### ### ### ### ###
			website.logger.info("#### #### #### #### #### ####")
			user_input_vcf_file_name = gsvar_file.split("_time")[0] + ".vcf"
			return render_template("500.html", user_input_file=user_input_vcf_file_name)
	else:
		website.logger.info("#### #### #### #### #### ####")
		user_input_vcf_file_name = gsvar_file.split("_time")[0] + ".vcf"
		return render_template("500.html", user_input_file=user_input_vcf_file_name)


@website.route("/about")
def about():
	return render_template("about.html", title="About")


@website.route("/documentation")
def documentation():
	return render_template("documentation.html", title="Documentation")
