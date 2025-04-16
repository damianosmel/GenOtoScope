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
	# use the ensembl id to get its transcripts
	# and finally extract the strand of the gene's transcript
	### ### ### ###
	ensembl_gene_id = hugo_gene_row["ensembl_gene_id"]
	gene_strand = 'unknown'
	try:
		gene = ensembl_db.gene_by_id(ensembl_gene_id)
		gene_strand = gene.strand
	except ValueError as ve:
		print(f'Pyensembl does not contain gene with input gene id: {ve}')
	return gene_strand


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
path_root = "/home/damian/Documents/L3S/projects/hearing4all/human_genetics/genotoscope_data/misc/hgnc_genes_info"
data_path = join(path_root, "apr_2025")
hugo_genes_file = "hgnc_complete_set.txt"

### PyEnsembl ###
# Ensembl release 108 (Oct. 2022) uses the human genome reference GRCh38
ensembl_data = EnsemblRelease(108)

extract_gene_strands(data_path, hugo_genes_file, ensembl_data, update_file=True)