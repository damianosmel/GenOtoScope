from os.path import join
from pandas import read_csv
from pyensembl import EnsemblRelease


def extract_strand(hugo_gene_row, ensembl_db):
	"""
	Extract strand for gene row

	Parameters
	----------
	hugo_gene_row : pandas.Series
		hugo gene row
	ensembl_db: pyensembl.EnsemblRelease
		ensembl data object

	Returns
	-------
	str
		gene strand
	"""
	### ### ### ###
	# get ensembl id of gene row
	# use the ensembl id to get its trancripts
	# and finally extract the strand of the gene's transcript
	### ### ### ###
	ensembl_gene_id = hugo_gene_row["ensembl_gene_id"]
	try:
		transcript_ids = ensembl_db.transcript_ids_of_gene_id(ensembl_gene_id)
	except ValueError as ve:
		transcript_ids = []
		print("Pyensembl does not contain transcript ids for gene id: {}".format(ensembl_gene_id))
	finally:
		if len(transcript_ids) > 0:
			transcript = ensembl_data.transcript_by_id(transcript_ids[0])
			return transcript.exons[0].to_dict()["strand"]
		else:
			return "unknown"


def extract_gene_strands(data_path, hugo_genes_filename, ensembl_db, update_file=False):
	"""
	Extract strand information for each gene row of hugo genes file

	Parameters
	----------
	data_path : str
		data root path
	hugo_genes_filename : str
		hugo genes file name
	ensembl_db : pyensembl
		ensembl data
	update_file : bool
		Update hugo gene file with strand information (True), otherwise (False)

	Returns
	-------
	None
	"""
	### ### ### ###
	# apply function to get strand for a gene symbol through ensembl
	# if specified, save this hugo_genes_df['strand']
	### ### ### ###

	print("Extract strand information for each gene")
	hugo_genes_df = read_csv(join(data_path, hugo_genes_filename),
	                         sep='\t',  # field separator
	                         header=0,  # get header
	                         # index_col=[0, 1, 2],  # set as index the chr, start stop
	                         skipinitialspace=True,
	                         skip_blank_lines=True,
	                         error_bad_lines=False,
	                         warn_bad_lines=True
	                         )

	print("Hugo genes df head:\n {}".format(hugo_genes_df.head(5)))
	print("Hugo genes df columns:\n {}".format(hugo_genes_df.columns.values))

	hugo_genes_df["strand"] = hugo_genes_df.apply(extract_strand, ensembl_db=ensembl_db, axis=1)
	print("### ### Extracted gene strand statistics ### ###")
	print("Number of all genes: {}".format(hugo_genes_df.shape[0]))
	print("Assigned strand information:\n{}".format(hugo_genes_df["strand"].value_counts()))
	print("### ### ### ###")
	if update_file:
		[filename, file_extension] = hugo_genes_filename.split(".")
		strand_filename = filename + "_strand" + "." + file_extension
		print("Updated Hugo genes file with strand saved in {}".format(join(data_path, strand_filename)))
		hugo_genes_df.to_csv(join(data_path, strand_filename), sep="\t", header=True, index=False)


### data ###
path_root = "/home/damian/Documents/L3S/projects/hearing4all/human_genetics"
data_path = join(path_root, "data")
hugo_genes_file = "hgnc_complete_set.txt"

### PyEnsembl ###
# release 75 uses human reference genome GRCh37
ensembl_data = EnsemblRelease(75)

extract_gene_strands(data_path, hugo_genes_file, ensembl_data, update_file=True)
