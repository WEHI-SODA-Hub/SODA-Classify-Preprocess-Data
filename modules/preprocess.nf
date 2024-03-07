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
	val batch_name
	path qupath_data
	val additional_meta_data
	val cell_types_to_remove
	val change_to
	val unwanted_markers
	val unwanted_compartments
	val unwanted_statistics
	path preprocess_script

	// See: https://www.nextflow.io/docs/latest/process.html#outputs
	// each new output needs to be placed on a new line
	output:
	path ("report.html")
	path ("${batch_name}_*_labels.csv")
	path ("${batch_name}_images.csv")
	path ("${batch_name}_preprocessed_input_data.csv")
	path ("${batch_name}_decoder.json")
	
	// this is an example of some code to run in the code block 
	shell:
	flag_a = additional_meta_data == "" ? "" : "-a '${additional_meta_data}'"
	flag_l = cell_types_to_remove == "" ? "" : "-l '${cell_types_to_remove}'"
	flag_t = change_to == "" ? "" : "-t '${change_to}'"
	flag_m = unwanted_markers == "" ? "" : "-m '${unwanted_markers}'"
	flag_c = unwanted_compartments == "" ? "" : "-c '${unwanted_compartments}'"
	flag_s = unwanted_statistics == "" ? "" : "-s '${unwanted_statistics}'"
	'''
	python3 "!{preprocess_script}" !{params.target} \\
		-d "$(realpath !{qupath_data})" \\
		-o "$(realpath .)" \\
		-n "!{batch_name}" \\
		!{flag_a} !{flag_l} !{flag_t} !{flag_m} !{flag_m} !{flag_c} !{flag_s} \\
		> report.qmd

	# add output locations (script has no knowledge of publishDir)
	echo "" >> report.qmd
	echo '# Output Paths' >> report.qmd
	echo "" >> report.qmd
	echo "**Cell type labels:**" >> report.qmd
	echo "" >> report.qmd
	echo "\\`\\`\\`" >> report.qmd
	echo "!{params.output_folder}/!{batch_name}_cell_type_labels.csv" >> report.qmd
	echo "\\`\\`\\`" >> report.qmd
	echo "" >> report.qmd
	echo "**Image list:**" >> report.qmd
	echo "" >> report.qmd
	echo "\\`\\`\\`" >> report.qmd
	echo "!{params.output_folder}/!{batch_name}_images.csv" >> report.qmd
	echo "\\`\\`\\`" >> report.qmd
	echo "" >> report.qmd
	echo "**Decoder:**" >> report.qmd
	echo "" >> report.qmd
	echo "\\`\\`\\`" >> report.qmd
	echo "!{params.output_folder}/!{batch_name}_decoder.json" >> report.qmd
	echo "\\`\\`\\`" >> report.qmd
	echo "" >> report.qmd
	echo "**Preprocessed data:**" >> report.qmd
	echo "" >> report.qmd
	echo "\\`\\`\\`" >> report.qmd
	echo "!{params.output_folder}/!{batch_name}_preprocessed_input_data.csv" >> report.qmd
	echo "\\`\\`\\`" >> report.qmd
	echo "" >> report.qmd

	quarto render report.qmd --to html
	'''
}