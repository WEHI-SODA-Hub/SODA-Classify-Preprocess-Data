#!/usr/bin/env python3

import sklearn.preprocessing
import pandas as pd
import tabulate

import os
import json
import datetime

from preprocess_cell_type_classification import (
    list_2_md_table,
    remove_dots,
    generate_warnings,
    preprocess_celltypecolumn,
    convert_pixels_to_micrometre,
    save_image_coordinate_columns,
    remove_prefixes_underscores,
    remove_duplicate_columns,
    collect_markers,
    replace_cytoplasm_with_membrane,
    use_cell_measurement,
    remove_unwanted_compartments,
    remove_statistics,
    save_preprocessed_data,
)


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
                cell_type_count_table=self.cell_type_count_table,
                encoding_table=self.encoding_table,
                classification_count_table=self.classification_count_table,
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

## Count of each cell type:

{cell_type_count_table}

## Encoding:

{encoding_table}

## Count of each classification:

{classification_count_table}
        
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


def setup(output_folder, expression_mat_path):
    """
    Create output folder and read in the input CSV or Parquet file
    """
    os.makedirs(output_folder, exist_ok=True)
    
    # Determine file format from extension
    if expression_mat_path.endswith('.parquet'):
        expression_df = pd.read_parquet(expression_mat_path)
    elif expression_mat_path.endswith('.csv'):
        expression_df = pd.read_csv(expression_mat_path, index_col=0)
    else:
        raise ValueError(f"Unsupported file format. Expected .csv or .parquet, got: {expression_mat_path}")
    
    return expression_df


def binarize_and_save_fm(expression_df, output_folder, batch_name) -> tuple:
    """
    Converts +/- to 1/0 in the Classification column (representing functional marker of interest).
    """

    binarized = expression_df["Classification"].copy()

    labels = sorted(binarized.unique())
    decoder = {1: labels[0], 0: labels[1]}
    encoder = {labels[0]: 1, labels[1]: 0}

    filename = os.path.join(output_folder, f"{batch_name}_decoder.json")
    with open(filename, "w") as json_file:
        json.dump(decoder, json_file, indent=4)

    # replace "+" with 1, otherwise, 0.
    binarized = binarized.map(lambda x: int("+" in x))

    filename = os.path.join(output_folder, f"{batch_name}_binarized_labels.csv")
    binarized.to_csv(filename, index=False)

    return encoder, decoder


def one_hot_encode_cell_types(expression_df, cell_types) -> pd.DataFrame:
    """
    Convert column of cell types to one-hot-encoded array

    For example, column of [Epithelial Cells, Dendritic Cells, Other] becomes:

    | Dendritic Cells, Epithelial Cells, Other |
    |               0,                 1,    0 |
    |               1,                 0,    0 |
    |               0,                 0,    1 |
    """

    lb = sklearn.preprocessing.LabelBinarizer()
    cell_types_s = expression_df["Class"]
    expression_df = pd.concat(
        (
            expression_df,
            pd.DataFrame(
                lb.fit_transform(cell_types_s),
                columns=cell_types,
                index=expression_df.index,
            ),
        ),
        axis=1,
    )
    expression_df = expression_df.drop("Class", axis=1)

    return expression_df


def drop_markers(expression_df, markers, excluded_markers, cell_types):
    """
    In this step, markers which do not help in determining the cell type should be removed. For example, dsDNA will not help in determining
    cell types.

    Any markers where the staining did not work should also be removed.
    Keep only the columns with the markers you want to keep.
    """
    # return [marker for marker in markers if marker not in excluded_markers]
    markers_ = [s + ": " for s in markers if s not in excluded_markers]

    # identify which columns contain measurements of the markers we want to keep
    measurement_columns = [
        col for col in expression_df.columns if any(map(col.__contains__, markers_))
    ]

    # return dataframe with measurement columns, but retain also retain one-hot-enconded cell types
    # if cell-types aren't one-hot-encoded, then return just the measurements
    try:
        return expression_df.loc[:, measurement_columns + cell_types]
    except KeyError:
        return expression_df.loc[:, measurement_columns]


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
    with_celltype = True,
    output_format='csv'
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

    encoder, decoder = binarize_and_save_fm(expression_df, output_folder, batch_name)

    output_mibi_reporter.encoding_table = tabulate.tabulate(encoder.items())

    output_mibi_reporter.classification_count_table = (
        expression_df["Classification"].value_counts().to_markdown(tablefmt="simple")
    )

    expression_df = convert_pixels_to_micrometre(expression_df)

    save_image_coordinate_columns(
        expression_df, additional_meta_data_to_keep, output_folder, batch_name
    )

    output_mibi_reporter.cell_type_count_table = (
        expression_df["Class"].value_counts().to_markdown(tablefmt="simple")
    )

    if with_celltype: 
        expression_df = one_hot_encode_cell_types(expression_df, cell_types)
    else:
        expression_df.drop("Class", axis=1, inplace=True)

    expression_df = remove_prefixes_underscores(expression_df)

    expression_df = remove_duplicate_columns(expression_df)

    markers = collect_markers(expression_df)

    output_mibi_reporter.markers_table = list_2_md_table(markers)

    expression_df = drop_markers(expression_df, markers, unwanted_markers, cell_types)

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
