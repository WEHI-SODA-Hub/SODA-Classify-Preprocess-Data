#!/usr/bin/env python3

import pandas as pd
import numpy as np
import tabulate

import re
import os
import json
import datetime

class mibi_reporter:
    
    def print_report(self):
        print(self.report_template_str.format(
            report_date = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            batch_name = self.batch_name,
            output_folder = self.output_folder,
            expression_mat_path = self.expression_mat_path,
            cell_types_to_remove_table = self.cell_types_to_remove_table,
            change_to = self.change_to,
            additional_meta_data_to_keep_table = self.additional_meta_data_to_keep_table,
            unwanted_markers_table = self.unwanted_markers_table,
            unwanted_compartments_table = self.unwanted_compartments_table,
            unwanted_statistics_table = self.unwanted_statistics_table,
            found_cell_types_table = self.found_cell_types_table,
            cell_types_table = self.cell_types_table,
            encoding_table = self.encoding_table,
            cell_type_count_table = self.cell_type_count_table,
            markers_table = self.markers_table,
            null_columns_table = self.null_columns_table,
            warnings = self.warnings
        ))
    def __init__(self):
        self.report_template_str = """
---
title: MIBI data preprocessing summary report
date: {report_date}
format:
  html:
    number-sections: true
    embed-resources: true
    theme: cosmo
---

# Input Summary

## Batch name:

{batch_name}

## Output folder:

```
{output_folder}
```

## QuPath data to be preprocessed:
        
```
{expression_mat_path}
```

## Unwanted cell types:

{cell_types_to_remove_table}

## Changing unwanted cell types to:

{change_to}

## Additional metadata to keep:

{additional_meta_data_to_keep_table}

## Unwanted markers:

{unwanted_markers_table}

## Unwanted compartments:

{unwanted_compartments_table}

## Unwanted statistics:

{unwanted_statistics_table}

# Results

---

## Cell types found:

{found_cell_types_table}

## Cell types after removing user-defined cells:

{cell_types_table}

## Encoding:

{encoding_table}

## Count of each cell type:

{cell_type_count_table}
        
## Markers found:

{markers_table}

## User-defined markers to remove:

{unwanted_markers_table}

## Columns with NA values:

{null_columns_table}

**If there are columns with NA values:**

This will be an issue with the measurement names across different images, and cohorts. 
If the problem is due to different measurement names across different images, this can be fixed by changing the names for the columns
in the images where this is a problem. 

Some things to check:

* Do all of my images have the same channel names on QuPath?
* Did I change the channel names before or after the segmentation? 
    * If after, the measurements would have been created using the previous channel names.
    * You can either change the names of the columns (best option) or if the channel names were
    completely different and you don't know which corresponds to which, you should re run the segmentation
    using the new channel names. 

# Warnings

{warnings}
"""

def list_2_md_table(input_list, columns=3) -> str:
    """
    Convert a 1D list to a markdown table with specified no. of columns.
    """

    # could've used if input_list, but input_list might be a pandas series
    if len(input_list) > 0:
        list2 = [
            input_list[i : i + columns] for i in range(0, len(input_list), columns)
        ]
        return tabulate.tabulate(list2)
    else:
        return str(None)


def setup(output_folder, expression_mat_path):
    os.makedirs(output_folder, exist_ok=True)
    expression_df = pd.read_csv(expression_mat_path)
    return expression_df

def generate_warnings(expression_df) -> str:
    """
    Checks input dataframe for known issues:
        * Duplicate column names after removing underscores
    """
    warning_str = ''

    # duplicate column names
    column_names_orig = expression_df.columns
    column_names = column_names_orig.str.replace("Target:", "")
    column_names = column_names.str.replace('_', ' ')
    duplicated_column_names = column_names[column_names.duplicated()]

    if len(duplicated_column_names) > 0:
        warning_str += """## Duplicate Columns\n
Column names were transformed by removing \"Target:\" and replacing underscores with spaces. \
After this transformation, duplicate columns were merged. Merging involves averaging the values in the duplicated columns. \
Post-transformation duplicate columns found:

"""
        # create dataframe matching new duplicate names to original names
        duplicated_columns = {}
        for name in duplicated_column_names:
            idx = np.where(column_names == name)
            duplicated_columns[name] = column_names_orig[idx].to_list()
        duplicated_columns_df = pd.DataFrame(duplicated_columns).T

        # label the dataframe and save table
        duplicated_columns_df.index.name = "Duplicate Column Name"
        colnames = [f"Original Column Name {int(i)+1}" for i in duplicated_columns_df.columns]
        duplicated_columns_df.columns = colnames
        warning_str += duplicated_columns_df.to_markdown(tablefmt="simple")

    return warning_str

def preprocess_celltypecolumn(
    expression_df, cell_types_to_remove=["Unknown"], change_to="Other"
):
    """
    Preprocess the cell type column
    cell types which you want to remove. By remove, the cell type will be changed to what ever you set to the variable change_to
    """

    # Check that all the cell types are there
    # remove the Edited prefix which may have occured from the qupath script
    expression_df.loc[:, "Class"] = expression_df.loc[:, "Class"].str.replace(
        "Edited: ", ""
    )
    expression_df.loc[:, "Name"] = expression_df.loc[:, "Name"].str.replace(
        "Edited: ", ""
    )

    expression_df.loc[:, "Class"] = expression_df.loc[:, "Class"].str.replace(
        "Immune cells: ", ""
    )
    expression_df.loc[:, "Name"] = expression_df.loc[:, "Name"].str.replace(
        "Immune cells: ", ""
    )

    found_cell_types = sorted(expression_df.loc[:, "Class"].unique())

    expression_df.loc[:, "Class"] = expression_df.loc[:, "Class"].replace(
        cell_types_to_remove, change_to
    )
    expression_df.loc[:, "Name"] = expression_df.loc[:, "Name"].replace(
        cell_types_to_remove, change_to
    )

    cell_types = expression_df.loc[:, "Class"].unique()
    cell_types = sorted(cell_types)

    return found_cell_types, cell_types


def create_encoder_decoder(cell_types, output_folder, batch_name):
    # encoder for converting your labels
    encoder = {cell_types[i]: i for i in range(len(cell_types))}

    # decoder for decoding the results of the model. Save somewhere safe.
    decoder = {i: cell_types[i] for i in range(len(cell_types))}

    with open(
        os.path.join(output_folder, f"{batch_name}_decoder.json"), "w"
    ) as json_file:
        json.dump(decoder, json_file, indent=4)

    return encoder, decoder


def save_encoded_labels(expression_df, encoder, output_folder, batch_name):
    """
    Save the labels as a separate csv file. The labels will be encoded with the above encoding
    """

    filename = os.path.join(output_folder, "{}_cell_type_labels.csv".format(batch_name))
    labels = expression_df.loc[:, ["Name"]]
    labels = labels.replace({"Name": encoder})
    labels.to_csv(filename, index=False)


def convert_pixels_to_micrometre(expression_df, pixel_size=0.3906):
    """
    If for some reason, the centroid measurements are done in pixels and not µm, this will convert the pixel values to microns.

    The pixel_size variable is the microns/pixel. This information should be available somewhere idk.
    """

    pixel_size = 0.3906  # fixed size (for now)

    for dim in ["X", "Y"]:
        try:
            null_arr = expression_df.loc[:, "Centroid {} µm".format(dim)].isnull()
            if null_arr.any() != False:
                expression_df.loc[null_arr.values, "Centroid {} µm".format(dim)] = (
                    expression_df.loc[null_arr.values, "Centroid {} px".format(dim)]
                    * pixel_size
                )
                expression_df.drop(["Centroid {} px".format(dim)], axis=1)
        except:
            expression_df.loc[:, "Centroid {} µm".format(dim)] = (
                expression_df.loc[:, "Centroid {} px".format(dim)] * pixel_size
            )
            expression_df = expression_df.drop(["Centroid {} px".format(dim)], axis=1)

    return expression_df


def save_image_coordinate_columns(
    expression_df, additional_meta_data, output_folder, batch_name
):
    """
    Save the image and coordinate columns. This is for when we want to import the results back into qupath
    """

    image_coord_cols = [
        "Image",
        "Centroid X µm",
        "Centroid Y µm",
    ] + additional_meta_data
    image_coord_df = expression_df.loc[:, image_coord_cols]
    image_coord_file_name = os.path.join(
        output_folder, "{}_images.csv".format(batch_name)
    )
    image_coord_df.to_csv(image_coord_file_name, index=False)


def remove_prefixes_underscores(expression_df):
    """
    Remove unnecessary prefixes and underscores.
    """
    expression_df.columns = expression_df.columns.str.replace("Target:", "")
    expression_df.columns = expression_df.columns.str.replace("_", " ")

    return expression_df

def remove_duplicate_columns(expression_df):
    """
    Removes duplicate columns and columns with the same name, but less data
    """
    # find duplicate column names
    duplicated_cols = expression_df.columns[expression_df.columns.duplicated()]
    # copy duplicated column data to new dataframe
    expression_df_duplicated = expression_df[duplicated_cols]
    # attempt to merge duplicated columns
    ## Net effect:
    ## * rows with same value are preserved
    ## * rows with one empty (NaN) and another with a value takes the value
    ## * rows with both empty stay empty
    ## * average is taken for rows with both values present
    expression_df_duplicated = expression_df_duplicated.groupby(
        expression_df_duplicated.columns, axis=1
    ).mean()
    expression_df = expression_df.drop(columns=duplicated_cols)
    expression_df = expression_df.join(expression_df_duplicated)

    return expression_df


def collect_markers(expression_df):
    """
    Collects all of the markers in this cohort
    """
    # markers to include
    markers = [
        col.replace(": Cell: Mean", "")
        for col in expression_df.columns
        if "Cell: Mean" in col
    ]
    return markers


def drop_markers(expression_df, markers, excluded_markers):
    """
    In this step, markers which do not help in determining the cell type should be removed. For example, dsDNA will not help in determining
    cell types.

    Any markers where the staining did not work should also be removed.
    Keep only the columns with the markers you want to keep.
    """
    # return [marker for marker in markers if marker not in excluded_markers]
    markers_ = [s + ": " for s in markers if s not in excluded_markers]
    measurement_columns = [
        col for col in expression_df.columns if any(map(col.__contains__, markers_))
    ]
    return expression_df.loc[:, measurement_columns]


def replace_cytoplasm_with_membrane(expression_df):
    """
    Due to the segmentation, some cells will not have a cytoplasm compartment. That is because the nuclei boundary and the cell boundary
    are the same pixels. This usually occurs in densely packed tumours where the nuclei and cell boundary merge.

    Because of this, some cells will have missing values in the cell cytoplasm measurements. We will therefore, instead of imputing the
    missing values with a 0, we will use the membrane measurement. This is a more representative way to impute the missing measurements.
    """

    for col in expression_df.columns:
        null_arr = expression_df.loc[:, col].isnull()
        if null_arr.values.any():
            if "Cytoplasm" in col:
                new_col = col.replace("Cytoplasm", "Membrane", 1)
                expression_df.loc[null_arr.values, col] = expression_df.loc[
                    null_arr.values, new_col
                ]

    return expression_df


def use_cell_measurement(expression_df):
    """
    Due to the segmentation, some cells will not have a nucleus area too small.

    Because of this, some cells will have missing values in the nucleus measurements. We will therefore, instead of imputing the
    missing values with a 0, we will use the cell measurement. This is a more representative way to impute the missing measurements.
    """

    for col in expression_df.columns:
        null_arr = expression_df.loc[:, col].isnull()
        if null_arr.values.any():
            if "Nucleus" in col:
                new_col = col.replace("Nucleus", "Cell", 1)
                expression_df.loc[null_arr.values, col] = expression_df.loc[
                    null_arr.values, new_col
                ]

    return expression_df


def remove_unwanted_compartments(expression_df, unwanted_compartments):
    """
    If you believe that some compartments should not be considered during the phenotyping, remove them here
    """

    compartments_cols_to_remove = [
        col
        for col in expression_df.columns
        if any(map(col.__contains__, unwanted_compartments))
    ]
    expression_df = expression_df.drop(columns=compartments_cols_to_remove)

    return expression_df


def remove_statistics(expression_df, unwanted_stats):
    """
    If you believe that some statistics should not be considered during the phenotyping, remove them here
    """
    statistics_cols_to_remove = [
        col
        for col in expression_df.columns
        if any(map(col.__contains__, unwanted_stats))
    ]
    expression_df = expression_df.drop(columns=statistics_cols_to_remove)

    return expression_df


def save_preprocessed_data(expression_df, output_folder, batch_name):
    """
    Save the input preprocessed data
    """
    output_expression_df_path = os.path.join(
        output_folder, "{}_preprocessed_input_data.csv".format(batch_name)
    )
    expression_df.to_csv(output_expression_df_path, index=False)


def preprocess_training_data(
    batch_name,
    output_folder,
    expression_mat_path,
    cell_types_to_remove,
    change_to,
    additional_meta_data_to_keep,
    unwanted_markers,
    unwanted_compartments,
    unwanted_statistics
) -> mibi_reporter:
    
    output_mibi_reporter = mibi_reporter()
    output_mibi_reporter.batch_name = batch_name
    output_mibi_reporter.output_folder = output_folder
    output_mibi_reporter.expression_mat_path = expression_mat_path
    output_mibi_reporter.cell_types_to_remove_table = list_2_md_table(cell_types_to_remove)
    output_mibi_reporter.change_to = change_to
    output_mibi_reporter.additional_meta_data_to_keep_table = list_2_md_table(additional_meta_data_to_keep)
    output_mibi_reporter.unwanted_markers_table = list_2_md_table(unwanted_markers)
    output_mibi_reporter.unwanted_compartments_table = list_2_md_table(unwanted_compartments)
    output_mibi_reporter.unwanted_statistics_table = list_2_md_table(unwanted_statistics)

    expression_df = setup(output_folder, expression_mat_path)

    output_mibi_reporter.warnings = generate_warnings(expression_df)

    cell_types, found_cell_types = preprocess_celltypecolumn(
        expression_df, cell_types_to_remove, change_to
    )

    output_mibi_reporter.found_cell_types_table = list_2_md_table(found_cell_types)
    output_mibi_reporter.cell_types_table = list_2_md_table(cell_types)

    encoder, decoder = create_encoder_decoder(cell_types, output_folder, batch_name)

    output_mibi_reporter.encoding_table = tabulate.tabulate([[k] for k in encoder.keys()], showindex="always")

    save_encoded_labels(expression_df, encoder, output_folder, batch_name)

    expression_df = convert_pixels_to_micrometre(expression_df)

    save_image_coordinate_columns(
        expression_df, additional_meta_data_to_keep, output_folder, batch_name
    )

    output_mibi_reporter.cell_type_count_table = expression_df.loc[:, "Class"].value_counts().to_markdown(tablefmt="simple")

    expression_df = remove_prefixes_underscores(expression_df)

    expression_df = remove_duplicate_columns(expression_df)

    markers = collect_markers(expression_df)

    output_mibi_reporter.markers_table = list_2_md_table(markers)
    output_mibi_reporter.unwanted_markers_table = list_2_md_table(unwanted_markers)

    expression_df = drop_markers(expression_df, markers, unwanted_markers)

    expression_df = replace_cytoplasm_with_membrane(expression_df)

    expression_df = use_cell_measurement(expression_df)

    expression_df = remove_unwanted_compartments(expression_df, unwanted_compartments)

    expression_df = remove_statistics(expression_df, unwanted_statistics)

    output_mibi_reporter.null_columns_table = list_2_md_table(expression_df.columns[expression_df.isna().any()].values, 2)

    save_preprocessed_data(expression_df, output_folder, batch_name)

    return output_mibi_reporter


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        prog="MIBI-preprocess-data",
        description="""This script is for preprocessing annotated data which has been exported from QuPath.
The data will be exported for XGBoost training or any supervised machine learning method of choice.""",
    )

    # declare input arguments.
    parser.add_argument(
        "-n",
        "--batch-name",
        required=True,
        help="Batch name used to label output files.",
    )
    parser.add_argument(
        "-o",
        "--output-folder",
        default="output",
        help="Where preprocessed files will be stored. The folder will be created if it doesn't already exist.",
    )
    parser.add_argument(
        "-d",
        "--qupath-data",
        required=True,
        help="The raw data exported from QuPath to be preprocessed.",
    )
    parser.add_argument(
        "-a",
        "--additional-metadata-to-keep",
        help="A comma-delimited list of additional metadata columns you wish to keep.",
    )
    parser.add_argument(
        "-l",
        "--unwanted-celltypes",
        help='A comma-delimited list of cell types identified that you wish to remove. E.g., "B cells,CD4 T cells".',
    )
    parser.add_argument(
        "-t",
        "--change-unwanted-celltypes-to",
        default="Other",
        help="The label assigned to celltypes that you have flagged for removal. Default: Other.",
    )
    parser.add_argument(
        "-m",
        "--unwanted-markers",
        help="A comma-delimited list of markers you want to remove from the phenotyping.",
    )
    parser.add_argument(
        "-c",
        "--unwanted-compartments",
        help="A comma-delimited list of compartments you want to remove from the phenotyping.",
    )
    parser.add_argument(
        "-s",
        "--unwanted-statistics",
        default=" Nucleus: Mean, Nucleus: Median, Nucleus: Min, Nucleus: Max, Nucleus: Std.Dev, Nucleus: Percentile: 91.0, Nucleus: Percentile: 92.0, Nucleus: Percentile: 93.0, Nucleus: Percentile: 94.0, Nucleus: Percentile: 96.0,Nucleus: Percentile: 97.0, Nucleus: Percentile: 98.0, Nucleus: Percentile: 99.0, Nucleus: Percentile: 99.5, Nucleus: Percentile: 99.9, Nucleus: Percentile: 95.0, Nucleus: Percentile: 90.0, Nucleus: Percentile: 80.0, Nucleus: Percentile: 70.0",
        help="A comma-delimited list of statistics you want to remove from the phenotyping.",
    )

    args = parser.parse_args()

    # assign args to variables
    batch_name = args.batch_name
    output_folder = args.output_folder
    expression_mat_path = args.qupath_data
    change_to = args.change_unwanted_celltypes_to

    # might be None.
    try:
        cell_types_to_remove = args.unwanted_celltypes.split(",")
    except:
        cell_types_to_remove = []

    # might be None
    try:
        additional_meta_data_to_keep = args.additional_metadata_to_keep.split(",")
    except:
        additional_meta_data_to_keep = []

    # might be None
    try:
        unwanted_markers = args.unwanted_markers.split(",")
    except:
        unwanted_markers = []

    # might be None
    try:
        unwanted_compartments = args.unwanted_compartments.split(",")
    except:
        unwanted_compartments = []

    # might be None
    try:
        unwanted_statistics = args.unwanted_statistics.split(",")
    except:
        unwanted_statistics = []

    output_mibi_reporter = preprocess_training_data(
        batch_name,
        output_folder,
        expression_mat_path,
        cell_types_to_remove,
        change_to,
        additional_meta_data_to_keep,
        unwanted_markers,
        unwanted_compartments,
        unwanted_statistics
    )

    output_mibi_reporter.print_report()