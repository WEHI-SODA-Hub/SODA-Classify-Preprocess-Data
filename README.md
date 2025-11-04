# MIBI Preprocess Data Pipeline

This Nextflow pipeline is a sub-pipeline in the MIBI suite. It is used to preprocess output QuPath
data in preparation for XGBoost model training or application. This README contains WEHI-specific
as well as general usage instructions.

## Introduction

This pipeline is a single-process pipeline, and is mostly used as an interface to Nextflow tower. It
takes output from QuPath and applies the following transformations:

1. Try to fix punctuation issues in headers from newer QuPath data
   * e.g., on export QuPath will represent `MHC I (HLA-DR): Membrane: Percentile: 98.0` as `MHC.I..HLA.DR...Membrane..Percentile..98.0`.
1. Removes redundant text from the class column like "Editied: ", and "Immune cells: "
2. Converts pixel length measurements to micrometre (if relevant).
4. Removes "Target: " from column names and replaces underscores with spaces
5. Attempts to merge duplicate columns (after transformation 3.)
6. Replaces Cytoplasm measurements with membrane measurements
7. Replace Nucleus measurements with Cell measurements in the case where there are Nucleus measurement values missing
8. Remove user-defined unwanted cell types, markers, compartments, or statistics

At the end of pipeline, a report is generated which supplies information about the input data, which
you can use to decide which cell types, markers, compartments, and/or statistics to discard in the
next run.

## Usage

Parameters:

```
  Usage:  nextflow run main.nf

  Required Arguments:

  --batch_name BATCH_NAME
        Batch name used to label output files.
  --target TARGET
        Whether to preprocess the data for the "main-cell-type", 
        "fm-markers-only", or the "fm-with-celltype" pipeline.
  --output_folder OUTPUT_FOLDER
        Where preprocessed files will be stored. The folder will be created if 
        it doesn't already exist.
  --input_data QUPATH_DATA
        The raw data exported from QuPath to be preprocessed.

  Optional Arguments:

  --additional_meta_data ADDITIONAL_METADATA
        A comma-delimited list of additional metadata columns you wish to keep.
  --cell_types_to_remove CELL_TYPES_TO_REMOVE
        A comma-delimited list of cell types identified that you wish to remove. 
        E.g., "B cells,CD4 T cells".
  --change_to CHANGE_UNWANTED_CELLTYPES_TO
        The label assigned to celltypes that you have flagged for removal. 
        Default: Other.
  --unwanted_markers UNWANTED_MARKERS
        A comma-delimited list of markers you want to remove from the 
        phenotyping.
  --unwanted_compartments UNWANTED_COMPARTMENTS
        A comma-delimited list of compartments you want to remove from the 
        phenotyping.
  --unwanted_statistics UNWANTED_STATISTICS
        A comma-delimited list of statistics you want to remove from the 
        phenotyping.
  --memory
        The RAM to allocate for the preprocessing. Include units e.g. "2 GB"
  --before_script
        Command or script to run before the process runs e.g. to make conda 
        available.
```

When you're preprocessing "fresh" data, at first, you only need to supply `--batch_name`, `--output_folder`, and `--input_data`.
After inspecting the results, you can supply the optional arguments as needed.

If you feel comfortable with the command line, you can run the preprocessing Python script directly.

```
$ conda env create -f envs/environment.yml
$ conda activate xgboost-cell-classification
$ python python3 bin/mibi-preprocess.py --help
usage: MIBI-preprocess-data [-h] -n BATCH_NAME [-o OUTPUT_FOLDER] -d QUPATH_DATA [-a ADDITIONAL_METADATA_TO_KEEP]
                            [-l UNWANTED_CELLTYPES] [-t CHANGE_UNWANTED_CELLTYPES_TO] [-m UNWANTED_MARKERS]
                            [-c UNWANTED_COMPARTMENTS] [-s UNWANTED_STATISTICS]
                            {cell-type,fm-markers-only,fm-with-celltype}

This script is for preprocessing annotated data which has been exported from QuPath. The data will be exported for XGBoost
training or any supervised machine learning method of choice.

positional arguments:
  {main-cell-type,fm-markers-only,fm-with-celltype}
                        Whether to preprocess the data for the "cell type" classification pipeline, or the "functional marker"
                        classification pipeline

optional arguments:
  -h, --help            show this help message and exit
  -n BATCH_NAME, --batch-name BATCH_NAME
                        Batch name used to label output files.
  -o OUTPUT_FOLDER, --output-folder OUTPUT_FOLDER
                        Where preprocessed files will be stored. The folder will be created if it doesn't already exist.
  -d QUPATH_DATA, --qupath-data QUPATH_DATA
                        The raw data exported from QuPath to be preprocessed.
  -a ADDITIONAL_METADATA_TO_KEEP, --additional-metadata-to-keep ADDITIONAL_METADATA_TO_KEEP
                        A comma-delimited list of additional metadata columns you wish to keep.
  -l UNWANTED_CELLTYPES, --unwanted-celltypes UNWANTED_CELLTYPES
                        A comma-delimited list of cell types identified that you wish to remove. E.g., "B cells,CD4 T cells".
  -t CHANGE_UNWANTED_CELLTYPES_TO, --change-unwanted-celltypes-to CHANGE_UNWANTED_CELLTYPES_TO
                        The label assigned to celltypes that you have flagged for removal. Default: Other.
  -m UNWANTED_MARKERS, --unwanted-markers UNWANTED_MARKERS
                        A comma-delimited list of markers you want to remove from the phenotyping.
  -c UNWANTED_COMPARTMENTS, --unwanted-compartments UNWANTED_COMPARTMENTS
                        A comma-delimited list of compartments you want to remove from the phenotyping.
  -s UNWANTED_STATISTICS, --unwanted-statistics UNWANTED_STATISTICS
                        A comma-delimited list of statistics you want to remove from the phenotyping.
```

### Example Usage

Download and unzip example input from [here](https://github.com/BioimageAnalysisCoreWEHI/MIBI/blob/main/annotationsTest.csv.zip)

```
nextflow run main.nf \
    --batch_name test \
    --target main-cell-type \
    --output_folder /tmp/mibi-test-run-output \
    --input_data annotationsTest.csv \
    --cell_types_to_remove Unknown \
    --change_to Other \
    --unwanted_markers dsDNA,Beta-Tubulin,CD39,CD49a,Tantalum

# Using default values, equivalent to

nextflow run main.nf \
    --batch_name test \
    --output_folder /tmp/mibi-test-run-output \
    --input_data annotationsTest.csv \
    --unwanted_markers dsDNA,Beta-Tubulin,CD39,CD49a,Tantalum
```

Outputs will be located in `/tmp/mibi-test-run-output`.

## Pipeline Output

The pipeline will produce an HTML report which you can view in your browser. This report povides detailed information for you to
inspect. This report is named `report.html` and is saved in the directory as specified by `--output_folder`.

The pipeline will also produce 4 files:

* `<batch name>_preprocessed_input_data.csv` - the fully pre-processed results.
* `<batch name>_cell_type_labels.csv` - the integer cell type label corresponding to each row of the preprocessed data CSV.
* `<batch name>_decoder.json` - the dictionary to convert cell type integer values back to the text labels. This won't be produced when preprocessing training data i.e., when the `CLASS` column is empty.
* `<batch name>_images.csv` - the image file and centroid coordinates associated with each cell obvservation.
* `<batch name>_binarized_labels.csv` - Classifications converted to 1/0 based on +/- (only if using `fm-*` targets).

If using the Python script directly, the same 4 output files are produced, but the report is printed to the terminal in markdown
format. This can be rendered by quarto, if so desired (but not strictly necessary).

## Credits

The core functionality of the MIBI pipeline was developed by Kenta Yotoke (@yokotenka) under the supervision of Claire Marceaux 
(@ClaireMarceaux). The pipeline was adapted to Nextflow by Edward Yang (@edoyango) and maintained by Michael Mckay (@mikemcka) and Michael Milton (@multimeric).

## Citation

TBC
