#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Import subworkflows to be run in the workflow
include { PREPROCESS } from './modules/preprocess'

/// Print a header
log.info """\

=======================================================================================
MIBI-preprocess - nf 
=======================================================================================

Research Computing Platform, WEHI

Find documentation and more info @ https://github.com/BioimageAnalysisCoreWEHI/MIBI-preprocess-data

Cite this pipeline @ INSERT DOI

Log issues @ https://github.com/BioimageAnalysisCoreWEHI/MIBI-preprocess-data/issues

=======================================================================================
Workflow run parameters 
=======================================================================================
batch_name            : ${params.batch_name}
output_folder         : ${params.output_folder}
input_data            : ${params.input_data}
additional_meta_data  : ${params.additional_meta_data}
cell_types_to_remove  : ${params.cell_types_to_remove}
change_to             : ${params.change_to}
unwanted_markers      : ${params.unwanted_markers}
unwanted_compartments : ${params.unwanted_compartments}
unwanted_statistics   : ${params.unwanted_statistics}
preprocess_script     : ${params.preprocess_script}
workDir               : ${workflow.workDir}
RAM                   : ${params.memory}
=======================================================================================

"""

// Help function
def helpMessage() {
    log.info"""
  Usage:  nextflow run main.nf

  Required Arguments:

  --batch_name BATCH_NAME
        Batch name used to label output files.
  --output_folder OUTPUT_FOLDER
        Where preprocessed files will be stored. The folder will be created if it doesn't already exist.
  --input_data QUPATH_DATA
        The raw data exported from QuPath to be preprocessed.

  Optional Arguments:

  --additional_meta_data ADDITIONAL_METADATA
        A comma-delimited list of additional metadata columns you wish to keep.
  --cell_types_to_remove CELL_TYPES_TO_REMOVE
        A comma-delimited list of cell types identified that you wish to remove. E.g., "B cells,CD4 T cells".
  --change_to CHANGE_UNWANTED_CELLTYPES_TO
        The label assigned to celltypes that you have flagged for removal. Default: Other.
  --unwanted_markers UNWANTED_MARKERS
        A comma-delimited list of markers you want to remove from the phenotyping.
  --unwanted_compartments UNWANTED_COMPARTMENTS
        A comma-delimited list of compartments you want to remove from the phenotyping.
  --unwanted_statistics UNWANTED_STATISTICS
        A comma-delimited list of statistics you want to remove from the phenotyping.
  --memory
        The RAM to allocate for the preprocessing. Include units e.g. "2 GB"
  --before_script
        Command or script to run before the process runs e.g. to make conda available.
	
""".stripIndent()
}

workflow {

	// Show help message if --help is run or if any required params are not  provided at runtime
	if ( params.help || params.batch_name == "" || params.output_folder == "" || params.input_data == ""){   
		helpMessage()
		exit 1
	} else {
		// Run preprocessing process
		(qmd, cell_type_labels, images, results, decoder) = PREPROCESS(
                  params.batch_name,
                  Channel.fromPath(params.input_data),
                  params.additional_meta_data,
                  params.cell_types_to_remove,
                  params.change_to,
                  params.unwanted_markers,
                  params.unwanted_compartments,
                  params.unwanted_statistics,
                  Channel.fromPath(params.preprocess_script)
            )
	}
}

workflow.onComplete {
summary = """
=======================================================================================
Workflow execution summary
=======================================================================================

Duration    : ${workflow.duration}
Success     : ${workflow.success}
workDir     : ${workflow.workDir}
Exit status : ${workflow.exitStatus}
outDir      : ${params.output_folder}

=======================================================================================
  """
println summary

}
