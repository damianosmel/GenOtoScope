from flask import render_template
from genotoscope.website import website

@website.errorhandler(404)
def not_found_error(error):
	return render_template('404.html'),404

@website.errorhandler(500)
def internal_error(error):
	return render_template('500.html'),500