#!/bin/env nextflow 

// Enable DSL-2 syntax
nextflow.enable.dsl=2

// Define the process
process PREPROCESS {	
	cpus "${params.cpus}"
	publishDir "${params.output_folder}", mode: 'copy'
	conda "${projectDir}/envs/environment.yml"
	memory "${params.memory}"
	beforeScript "${params.before_script}"

	// See: https://www.nextflow.io/docs/latest/process.html#inputs
	// each input needs to be placed on a new line
	input:
	path qupath_data

	// See: https://www.nextflow.io/docs/latest/process.html#outputs
	// each new output needs to be placed on a new line
	output:
	path ("report.html")
	path ("${params.batch_name}_cell_type_labels.csv")
	path ("${params.batch_name}_images.csv")
	path ("${params.batch_name}_preprocessed_input_data.csv")
	path ("${params.batch_name}_decoder.json")
	
	// this is an example of some code to run in the code block 
	shell:
	if (params.additional_meta_data == "") { flag_a = "" } else { flag_a = "-a '${params.additional_meta_data}'" }
	if (params.cell_types_to_remove == "") { flag_l = "" } else { flag_l = "-l '${params.cell_types_to_remove}'" }
	if (params.change_to == "") { flag_t = "" } else { flag_t = "-t '${params.change_to}'" }
	if (params.unwanted_markers == "") { flag_m = "" } else { flag_m = "-m '${params.unwanted_markers}'" }
	if (params.unwanted_compartments == "") { flag_c = "" } else { flag_c = "-c '${params.unwanted_compartments}'" }
	if (params.unwanted_statistics == "") { flag_s = "" } else { flag_s = "-s '${params.unwanted_statistics}'" }
	'''
	python3 "!{projectDir}/scripts/preprocess-training-data.py" \\
		-d "$(realpath !{qupath_data})" \\
		-o "$(realpath .)" \\
		-n "!{params.batch_name}" \\
		!{flag_a} !{flag_l} !{flag_t} !{flag_m} !{flag_m} !{flag_c} !{flag_s} \\
		> report.qmd

	# add output locations (script has no knowledge of publishDir)
	echo "" >> report.qmd
	echo '# Output Paths' >> report.qmd
	echo "" >> report.qmd
	echo "**Cell type labels:**" >> report.qmd
	echo "" >> report.qmd
	echo "\\`\\`\\`" >> report.qmd
	echo "!{params.output_folder}/!{params.batch_name}_cell_type_labels.csv" >> report.qmd
	echo "\\`\\`\\`" >> report.qmd
	echo "" >> report.qmd
	echo "**Image list:**" >> report.qmd
	echo "" >> report.qmd
	echo "\\`\\`\\`" >> report.qmd
	echo "!{params.output_folder}/!{params.batch_name}_images.csv" >> report.qmd
	echo "\\`\\`\\`" >> report.qmd
	echo "" >> report.qmd
	echo "**Decoder:**" >> report.qmd
	echo "" >> report.qmd
	echo "\\`\\`\\`" >> report.qmd
	echo "!{params.output_folder}/!{params.batch_name}_decoder.json" >> report.qmd
	echo "\\`\\`\\`" >> report.qmd
	echo "" >> report.qmd
	echo "**Preprocessed data:**" >> report.qmd
	echo "" >> report.qmd
	echo "\\`\\`\\`" >> report.qmd
	echo "!{params.output_folder}/!{params.batch_name}_preprocessed_input_data.csv" >> report.qmd
	echo "\\`\\`\\`" >> report.qmd
	echo "" >> report.qmd

	quarto render report.qmd --to html
	'''
}