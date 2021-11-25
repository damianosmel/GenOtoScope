from os.path import join

"""
Format UniProt bed files to contain the same number of columns,
in order bedtools utilities, such intersect(), run without errors
"""

def restrict_bed_columns(input_dir,filename_in, max_columns_num):
	"""
	Resrict columns of input bed file

	Parameters
	----------
	input_dir : input directory, where the bed file lies in
	filename_in : file name
	max_columns_num : maximum number of columns

	Returns
	-------
	string
		name of output bed file with restricted number of columns
	"""
	filename_out = filename_in.split(".bed")[0] + "_" + str(max_columns_num) +"columns" + ".bed"
	print("Restrict input: {} to {} max number of columns, outputing: {}".format(filename_in,max_columns_num,filename_out))
	with open(join(input_dir,filename_in),"r") as file_in, open(join(input_dir,filename_out),"w") as file_out:
		for line_index,line in enumerate(file_in):
			line_splits = line.strip("\n").split("\t")
			if len(line_splits) > max_columns_num:
				# restrict the number of the columns to maximum number
				print("restrict line with index: {}".format(line_index))
				line_splits = line_splits[0:max_columns_num]
			elif len(line_splits) < max_columns_num:
				print("append '.' to line with index: {}".format(line_index))
				# append with '.' to create a line with max number of columns
				while len(line_splits) < max_columns_num:
					line_splits.append(".")
			file_out.write("\t".join(line_splits)+"\n")
	return filename_out

### #### ###
#  Input   #
### #### ###
input_dir = "/home/damian/Documents/L3S/projects/hearing4all/human_genetics/genotoscope_data/uniprot/uniprot_feb_2021"
domains_file = "UP000005640_9606_domain_hg19.bed"
repeats_file = "UP000005640_9606_repeat_hg19.bed"

### #### ###
#  Output  #
### #### ###
domains_15_col_file = restrict_bed_columns(input_dir,domains_file,15)
repeats_16_col_file = restrict_bed_columns(input_dir,repeats_file,16)