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
    parser.add_argument(
        "target",
        help='Whether to preprocess the data for the "cell type" classification pipeline, or the "functional marker" classification pipeline',
        choices=["cell-type", "functional-marker"],
    )

    args = parser.parse_args()

    # assign args to variables
    batch_name = args.batch_name
    output_folder = args.output_folder
    expression_mat_path = args.qupath_data
    change_to = args.change_unwanted_celltypes_to

    # might be None.
    try:
        cell_types_to_remove = [s.strip() for s in args.unwanted_celltypes.split(",")]
    except:
        cell_types_to_remove = []

    # might be None
    try:
        additional_meta_data_to_keep = [
            s.strip() for s in args.additional_metadata_to_keep.split(",")
        ]
    except:
        additional_meta_data_to_keep = []

    # might be None
    try:
        unwanted_markers = [s.strip() for s in args.unwanted_markers.split(",")]
    except:
        unwanted_markers = []

    # might be None
    try:
        unwanted_compartments = [
            s.strip() for s in args.unwanted_compartments.split(",")
        ]
    except:
        unwanted_compartments = []

    # might be None
    try:
        unwanted_statistics = [s.strip() for s in args.unwanted_statistics.split(",")]
    except:
        unwanted_statistics = []

    if args.target == "cell-type":
        from preprocess_cell_type_classification import preprocess_training_data
    elif args.target == "functional-marker":
        from preprocess_functional_marker_classification import preprocess_training_data

    output_mibi_reporter = preprocess_training_data(
        batch_name,
        output_folder,
        expression_mat_path,
        cell_types_to_remove,
        change_to,
        additional_meta_data_to_keep,
        unwanted_markers,
        unwanted_compartments,
        unwanted_statistics,
    )

    output_mibi_reporter.print_report()
