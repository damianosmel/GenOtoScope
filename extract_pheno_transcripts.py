from genotoscope.beds.ExtractPhenoTranscripts import ExtractPhenoTranscripts

### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# Extract clinical significant exons                      #
# from phenotype relevant transcripts                     #
# HL transcripts from: DOI: 10.1016/j.jmoldx.2018.06.005. #
### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

### ### ###
# input and output paths
### ### ###
# local machine
# input
data_path = "/home/damian/Documents/L3S/projects/hearing4all/human_genetics/data/misc"
hl_transcripts_file = "hl_relevant_transcripts.xlsx"
output_path = "/home/damian/Documents/L3S/projects/hearing4all/human_genetics/data/annotation_beds"
output_dir = "hearing_loss"

### ### ###
# Extract HL exons into bed file
### ### ###
# set up parameters for extraction
bed_6columns_header = "#chr\tstart\tend\texon\ttranscripts\tstrand"
extracted_exons_file = "pheno_relevant_exons.bed"
hl_transcripts_extractor = ExtractPhenoTranscripts(data_path,output_path,output_dir)
hl_transcripts_extractor.run(hl_transcripts_file,bed_6columns_header,extracted_exons_file)