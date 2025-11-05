#!/usr/bin/env python3

import pandas as pd
import numpy as np
import tabulate

import os
import json
import datetime


# class used to store report string and relevant data to be presented (as a string)
# the member attributes are populated in the functions below
class mibi_reporter:
    def print_report(self):
        print(
            self.report_template_str.format(
                report_date=datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                batch_name=self.batch_name,
                output_folder=self.output_folder,
                expression_mat_path=self.expression_mat_path,
                cell_types_to_remove_table=self.cell_types_to_remove_table,
                change_to=self.change_to,
                additional_meta_data_to_keep_table=self.additional_meta_data_to_keep_table,
                unwanted_markers_table=self.unwanted_markers_table,
                unwanted_compartments_table=self.unwanted_compartments_table,
                unwanted_statistics_table=self.unwanted_statistics_table,
                found_cell_types_table=self.found_cell_types_table,
                cell_types_table=self.cell_types_table,
                encoding_table=self.encoding_table,
                cell_type_count_table=self.cell_type_count_table,
                markers_table=self.markers_table,
                after_drop_markers_table=self.after_drop_markers_table,
                null_columns_table=self.null_columns_table,
                warnings=self.warnings,
            )
        )

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
    page-layout: full
    toc: true
    toc-location: left
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

## Markers after removing user-defined markers:

{after_drop_markers_table}

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

    try:
        # len(input_list) > 0 in case pandas series
        if len(input_list) > 0:
            list2 = [
                input_list[i : i + columns] for i in range(0, len(input_list), columns)
            ]
            return tabulate.tabulate(list2)
        else:
            return str(None)
    except TypeError: # catch input_list = None
        return str(None)


def setup(output_folder, expression_mat_path):
    """
    Create output folder and read in the input CSV or Parquet file
    """
    os.makedirs(output_folder, exist_ok=True)
    
    # Determine file format from extension
    if expression_mat_path.endswith('.parquet'):
        expression_df = pd.read_parquet(expression_mat_path)
    elif expression_mat_path.endswith('.csv'):
        try:
            expression_df = pd.read_csv(expression_mat_path)
        except UnicodeDecodeError:
            expression_df = pd.read_csv(expression_mat_path, encoding="cp1252")
    else:
        raise ValueError(f"Unsupported file format. Expected .csv or .parquet, got: {expression_mat_path}")
    
    return expression_df


def remove_dots(expression_df) -> pd.DataFrame:
    """
    Tries to translate newer QuPath data to older, more sensible format.
    """

    cols = expression_df.columns.copy()
    # let's start easy...
    cols = cols.str.replace("Âµm", "µm")
    cols = cols.str.replace("µm.2", "µm^2", regex=False)

    # these are "known" specific replacements
    specific_matches = (
        ("MHC.I..", "MHC I ("),
        ("MHC.II..", "MHC II ("),
        ("MHC_I_.", "MHC_I_("),
        ("MHC_II_.", "MHC_II_("),
        ("Target.", "Target:"),
        ("Beta.Tubulin", "Beta-Tubulin"),
        ("IFN.y", "IFN-y"),
        ("HLA.DR", "HLA-DR")
    )
    for m, r in specific_matches:
        cols = cols.str.replace(m, r, regex=False)

    # once all the known specific replacements are performed, we can be a little more presumptuous...
    cols = cols.str.replace("...", "): ", regex=False)
    cols = cols.str.replace("..", ": ", regex=False)
    # replaces periods with spaces when the period isn't between two numbers
    # first protect Std.Dev.
    cols = cols.str.replace("Std.Dev.", "STD_DEV_PLACEHOLDER", regex=False)
    cols = cols.str.replace("(?<!\d)\.(?!\d)", " ", regex=True)
    # return Std.Dev.
    cols = cols.str.replace("STD_DEV_PLACEHOLDER", "Std.Dev.", regex=False)

    expression_df.columns = cols

    return expression_df

def generate_warnings(expression_df) -> str:
    """
    Checks input dataframe for known issues:
        * Duplicate column names after removing underscores
    """
    warning_str = ""

    # duplicate column names
    column_names_orig = expression_df.columns
    column_names = column_names_orig.str.replace("Target:", "")
    column_names = column_names.str.replace("_", " ")
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
        colnames = [
            f"Original Column Name {int(i)+1}" for i in duplicated_columns_df.columns
        ]
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

    if expression_df["Class"].notnull().any():
        # Check that all the cell types are there
        # remove the Edited prefix which may have occured from the qupath script
        expression_df.loc[:, "Class"] = expression_df.loc[:, "Class"].str.replace(
            "Edited: ", ""
        )
        try:
            expression_df.loc[:, "Name"] = expression_df.loc[:, "Name"].str.replace(
                "Edited: ", ""
            )
        except KeyError:
            pass

        expression_df.loc[:, "Class"] = expression_df.loc[:, "Class"].str.replace(
            "Immune cells: ", ""
        )
        try:
            expression_df.loc[:, "Name"] = expression_df.loc[:, "Name"].str.replace(
                "Immune cells: ", ""
            )
        except KeyError:
            pass

        found_cell_types = sorted(expression_df.loc[:, "Class"].unique())

        expression_df.loc[:, "Class"] = expression_df.loc[:, "Class"].replace(
            cell_types_to_remove, change_to
        )
        try:
            expression_df.loc[:, "Name"] = expression_df.loc[:, "Name"].replace(
                cell_types_to_remove, change_to
            )
        except KeyError:
            pass

        cell_types = expression_df.loc[:, "Class"].unique()
        cell_types = sorted(cell_types)
    else:
        found_cell_types = None
        cell_types = None

    return found_cell_types, cell_types


def create_encoder_decoder(cell_types, output_folder, batch_name):
    """
    Creates a dictionary that maps between integer labels and string labels for the cell types
    """
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
        # case 1: Centroid X/Y µm" column exists but so does "Centroid X/Y px". Try to merge the two
        um_col = f"Centroid {dim} µm"
        px_col = f"Centroid {dim} px"
        cols = expression_df.columns
        if um_col in cols and px_col in cols:
            # branch 1: both um and px measurements are present.
            # fill empty um measurements with available px measurements, then drop px measurements
            expression_df = expression_df[um_col].fillna(
                expression_df[px_col] * pixel_size
            )
            expression_df.drop([px_col], axis=1)
        elif um_col not in cols and px_col in cols:
            # branch 2: only px measurements are present.
            # create new um column using px measurements. Then drop px measurements
            expression_df.loc[:, um_col] = (
                expression_df.loc[:, px_col] * pixel_size
            )
            expression_df = expression_df.drop([px_col], axis=1)
        elif um_col in cols and px_col not in cols:
            # branch 3: only um measurements are present.
            # do nothing
            pass
        else:
            # branch 4: neither um or px measurements are present.
            # throw error.
            raise RuntimeError(
                "X/Y centroid measurements (in either pixels or µm) are missing!"
            )

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


def save_preprocessed_data(expression_df, output_folder, batch_name, output_format='csv'):
    """
    Save the input preprocessed data in CSV or Parquet format
    """
    if output_format == 'parquet':
        output_expression_df_path = os.path.join(
            output_folder, "{}_preprocessed_input_data.parquet".format(batch_name)
        )
        expression_df.to_parquet(output_expression_df_path, index=False)
    else:  # csv
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
    unwanted_statistics,
    output_format='csv',
) -> mibi_reporter:
    output_mibi_reporter = mibi_reporter()
    output_mibi_reporter.batch_name = batch_name
    output_mibi_reporter.output_folder = output_folder
    output_mibi_reporter.expression_mat_path = expression_mat_path
    output_mibi_reporter.cell_types_to_remove_table = list_2_md_table(
        cell_types_to_remove
    )
    output_mibi_reporter.change_to = change_to
    output_mibi_reporter.additional_meta_data_to_keep_table = list_2_md_table(
        additional_meta_data_to_keep
    )
    output_mibi_reporter.unwanted_markers_table = list_2_md_table(unwanted_markers)
    output_mibi_reporter.unwanted_compartments_table = list_2_md_table(
        unwanted_compartments
    )
    output_mibi_reporter.unwanted_statistics_table = list_2_md_table(
        unwanted_statistics
    )

    expression_df = setup(output_folder, expression_mat_path)

    expression_df = remove_dots(expression_df)

    output_mibi_reporter.warnings = generate_warnings(expression_df)

    found_cell_types, cell_types = preprocess_celltypecolumn(
        expression_df, cell_types_to_remove, change_to
    )

    output_mibi_reporter.found_cell_types_table = list_2_md_table(found_cell_types)
    output_mibi_reporter.cell_types_table = list_2_md_table(cell_types)

    if cell_types:
        encoder, decoder = create_encoder_decoder(cell_types, output_folder, batch_name)

        output_mibi_reporter.encoding_table = tabulate.tabulate(
            [[k] for k in encoder.keys()], showindex="always"
        )

    else:

        encoder = None
        decoder = None
        output_mibi_reporter.encoding_table = list_2_md_table(None)

    save_encoded_labels(expression_df, encoder, output_folder, batch_name)

    expression_df = convert_pixels_to_micrometre(expression_df)

    save_image_coordinate_columns(
        expression_df, additional_meta_data_to_keep, output_folder, batch_name
    )

    output_mibi_reporter.cell_type_count_table = (
        expression_df.loc[:, "Class"].value_counts().to_markdown(tablefmt="simple")
    )

    expression_df = remove_prefixes_underscores(expression_df)

    expression_df = remove_duplicate_columns(expression_df)

    markers = collect_markers(expression_df)

    output_mibi_reporter.markers_table = list_2_md_table(markers)

    expression_df = drop_markers(expression_df, markers, unwanted_markers)

    after_drop_markers = collect_markers(expression_df)

    output_mibi_reporter.after_drop_markers_table = list_2_md_table(after_drop_markers)

    expression_df = replace_cytoplasm_with_membrane(expression_df)

    expression_df = use_cell_measurement(expression_df)

    expression_df = remove_unwanted_compartments(expression_df, unwanted_compartments)

    expression_df = remove_statistics(expression_df, unwanted_statistics)

    output_mibi_reporter.null_columns_table = list_2_md_table(
        expression_df.columns[expression_df.isna().any()].values, 2
    )

    save_preprocessed_data(expression_df, output_folder, batch_name, output_format)

    return output_mibi_reporter
