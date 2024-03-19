#!/bin/env nextflow 

// Enable DSL-2 syntax
nextflow.enable.dsl=2

process PREPROCESS {	
	cpus "${params.cpus}"
	publishDir "${params.output_folder}", mode: 'copy'
	conda "${projectDir}/envs/environment.yml"
	memory "${params.memory}"
	beforeScript "${params.before_script}"

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

	output:
	path ("${batch_name}_report.html")
	path ("${batch_name}_cell_type_labels.csv")
	path ("${batch_name}_images.csv")
	path ("${batch_name}_preprocessed_input_data.csv")
	path ("${batch_name}_decoder.json", optional: true)
	
	shell:
	flag_a = additional_meta_data == "" ? "" : "-a '${additional_meta_data}'"
	flag_l = cell_types_to_remove == "" ? "" : "-l '${cell_types_to_remove}'"
	flag_t = change_to == "" ? "" : "-t '${change_to}'"
	flag_m = unwanted_markers == "" ? "" : "-m '${unwanted_markers}'"
	flag_c = unwanted_compartments == "" ? "" : "-c '${unwanted_compartments}'"
	flag_s = unwanted_statistics == "" ? "" : "-s '${unwanted_statistics}'"
	'''
	REPORT_QMD=!{batch_name}_report.qmd
	python3 "!{preprocess_script}" \\
		-d "$(realpath !{qupath_data})" \\
		-o "$(realpath .)" \\
		-n "!{batch_name}" \\
		!{flag_a} !{flag_l} !{flag_t} !{flag_m} !{flag_m} !{flag_c} !{flag_s} \\
		> "$REPORT_QMD"

	# add output locations (script has no knowledge of publishDir)
	echo "" >> "$REPORT_QMD"
	echo '# Output Paths' >> "$REPORT_QMD"
	echo "" >> "$REPORT_QMD"
	echo "**Cell type labels:**" >> "$REPORT_QMD"
	echo "" >> "$REPORT_QMD"
	echo "\\`\\`\\`" >> "$REPORT_QMD"
	echo "!{params.output_folder}/!{batch_name}_cell_type_labels.csv" >> "$REPORT_QMD"
	echo "\\`\\`\\`" >> "$REPORT_QMD"
	echo "" >> "$REPORT_QMD"
	echo "**Image list:**" >> "$REPORT_QMD"
	echo "" >> "$REPORT_QMD"
	echo "\\`\\`\\`" >> "$REPORT_QMD"
	echo "!{params.output_folder}/!{batch_name}_images.csv" >> "$REPORT_QMD"
	echo "\\`\\`\\`" >> "$REPORT_QMD"
	echo "" >> "$REPORT_QMD"
	echo "**Decoder:**" >> "$REPORT_QMD"
	echo "" >> "$REPORT_QMD"
	echo "\\`\\`\\`" >> "$REPORT_QMD"
	[ -f "!{batch_name}_decoder.json" ] && echo "!{params.output_folder}/!{batch_name}_decoder.json" >> "$REPORT_QMD" || echo "None" >> "$REPORT_QMD"
	echo "\\`\\`\\`" >> "$REPORT_QMD"
	echo "" >> "$REPORT_QMD"
	echo "**Preprocessed data:**" >> "$REPORT_QMD"
	echo "" >> "$REPORT_QMD"
	echo "\\`\\`\\`" >> "$REPORT_QMD"
	echo "!{params.output_folder}/!{batch_name}_preprocessed_input_data.csv" >> "$REPORT_QMD"
	echo "\\`\\`\\`" >> "$REPORT_QMD"
	echo "" >> "$REPORT_QMD"

	quarto render "$REPORT_QMD" --to html
	'''
}