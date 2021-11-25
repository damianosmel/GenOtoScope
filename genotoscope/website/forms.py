from genotoscope.website import upload_extensions
from genotoscope.website import website

from flask_wtf import FlaskForm
from werkzeug.utils import secure_filename
from wtforms import SubmitField
from flask_wtf.file import FileField, FileRequired
from wtforms.validators import ValidationError

from genotoscope.website import uploads_path
import logging
from os.path import splitext, join
from codecs import decode


class VariantForm(FlaskForm):
	variant_vcf_file = FileField('VCF file:', validators=[FileRequired()])
	submit = SubmitField('Annotate & Classify!')

	def __init__(self, time_stamp, *args, **kwargs):
		posted_request = kwargs.pop('posted_request')  # cache the user object you pass in
		self.time_stamp = time_stamp
		if posted_request:
			self.posted_request = posted_request
		super(VariantForm, self).__init__(*args, **kwargs)  # and carry on to init the form

	def validate_variant_vcf_file(self, variant_vcf_file):
		"""
		Validate vcf file field on VariantForm

		Parameters
		----------
		variant_vcf_file : FileField
			form's file field

		Returns
		-------
		None
		"""
		posted_file = self.posted_request.files['variant_vcf_file']
		website.logger.info("=== ==== ===")
		website.logger.info("Validate_VCF: Start - Validate input file")
		variant_file_name = secure_filename(variant_vcf_file.data.filename)

		if variant_file_name != '':
			# if variant file is not empty
			input_name, input_extension = splitext(variant_file_name)
			if input_extension[1:] not in upload_extensions:
				# check if input extension is VCF
				website.logger.info("Validate_VCF: Error - User uploaded a file with no vcf extension")
				website.logger.info("Validate_VCF: End - Validate input file")
				website.logger.info("=== ==== ===")
				raise ValidationError("Please upload a file in VCF format only.")
			else:
				# count the number of variants
				vcf_content = posted_file.read()
				num_variants_minimum = 0
				for line in decode(vcf_content).split("\n"):
					if len(line) > 0 and line[0] != "#":
						num_variants_minimum += 1
					if num_variants_minimum > 1:
						break

				if num_variants_minimum == 0:
					# check if the vcf contains only the header
					website.logger.info("Validate_VCF: Error - User uploaded an vcf file containing no variant row")
					website.logger.info("Validate_VCF: End - Validate input file")
					website.logger.info("=== ==== ===")
					raise ValidationError("Please upload a non-empty vcf file, containing a single variant.")
				elif num_variants_minimum > 1:
					website.logger.info(
						"Validate_VCF: Error - User uploaded an vcf file containing more than one variant row")
					website.logger.info("Validate_VCF: End - Validate input file")
					website.logger.info("=== ==== ===")
					raise ValidationError("Please upload a vcf file containing a single variant.")
				else:
					# save content of valid vcf file in uploads path
					variant_unique_file = input_name + "_time" + str(self.time_stamp) + ".vcf"
					with open(join(uploads_path, variant_unique_file), 'w') as out_file:
						out_file.write(decode(vcf_content))
					website.logger.info("Validate_VCF: Success - User uploaded a vcf file containing one variant row")
					website.logger.info("Validate_VCF: End - Validate input file")
					website.logger.info("=== ==== ===")
