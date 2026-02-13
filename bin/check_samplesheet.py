#!/usr/bin/env python

"""Provide a command line tool to validate and transform tabular samplesheets."""

import argparse
import csv
import logging
import sys
from collections import Counter
from pathlib import Path

logger = logging.getLogger()


class RowChecker:
    """
    Define a service that can validate and transform each given row.

    Attributes:
        modified (list): A list of dicts, where each dict corresponds to a previously
            validated and transformed row. The order of rows is maintained.

    """

    VALID_TYPES = (
        "native",
        "ivt",
    )

    def __init__(
        self,
        sample_col="sample",
        fastq_col="fastq",
        type_col="type",
        replicate_col="replicate",
        fast5_dir_col="fast5_dir",
        **kwargs,
    ):
        """
        Initialize the row checker with the expected column names.

        Args:
            sample_col (str): The name of the column that contains the sample name
                (default "sample").
            fastq_col (str): The name of the column that contains the FastQ file path
                (default "fastq").
            type_col (str): The name of the column that contains the sample type
                (default "type").
            replicate_col (str): The name of the column that contains the replicate
                identifier (default "replicate").
            fast5_dir_col (str): The name of the column that contains the FAST5
                directory path (default "fast5_dir").

        """
        super().__init__(**kwargs)
        self._sample_col = sample_col
        self._fastq_col = fastq_col
        self._type_col = type_col
        self._replicate_col = replicate_col
        self._fast5_dir_col = fast5_dir_col
        self._seen = set()
        self.modified = []

    def validate_and_transform(self, row):
        """
        Perform all validations on the given row and insert the read pairing status.

        Args:
            row (dict): A mapping from column headers (keys) to elements of that row
                (values).

        """
        self._validate_sample(row)
        self._validate_fastq(row)
        self._validate_type(row)
        self._validate_replicate(row)
        self._validate_fast5_dir(row)
        self._seen.add((row[self._sample_col]))
        self.modified.append(row)

    def _validate_sample(self, row):
        """Assert that the sample name exists and is unique."""
        if len(row[self._sample_col]) <= 0:
            raise AssertionError("Sample input is required.")
        # Enforce that the sample name is unique
        if row[self._sample_col] in self._seen:
            raise AssertionError(f"The sample name {row[self._sample_col]} is duplicated.")

    def _validate_fastq(self, row):
        """Assert that the FastQ entry is non-empty and has the right format."""
        if len(row[self._fastq_col]) <= 0:
            raise AssertionError("The FastQ file is required.")
        self._validate_fastq_format(row[self._fastq_col])

    def _validate_fastq_format(self, filename):
        """Assert that a given filename has a valid FastQ extension."""
        if not any(filename.endswith(extension) for extension in [".fq", ".fq.gz", ".fastq", ".fastq.gz"]):
            raise AssertionError(
                f"The FastQ file has an unrecognized extension: {filename}\n"
                f"It should be '.fq', '.fq.gz', '.fastq' or '.fastq.gz'"
            )

    def _validate_type(self, row):
        """Assert that the sample type is either 'native' or 'ivt'."""
        if row[self._type_col] not in self.VALID_TYPES:
            raise AssertionError(
                f"Sample type must be 'native' or 'ivt', not '{row[self._type_col]}'"
            )

    def _validate_replicate(self, row):
        """Assert that the replicate identifier exists."""
        if len(row[self._replicate_col]) <= 0:
            raise AssertionError("Replicate identifier is required.")

    def _validate_fast5_dir(self, row):
        """Assert that the FAST5 directory path exists."""
        if len(row[self._fast5_dir_col]) <= 0:
            raise AssertionError("FAST5 directory path is required.")


def read_head(handle, num_lines=10):
    """Read the specified number of lines from the current position in the file."""
    lines = []
    for idx, line in enumerate(handle):
        if idx == num_lines:
            break
        lines.append(line)
    return "".join(lines)


def sniff_format(handle):
    """
    Detect the tabular format.

    Args:
        handle (text file): A handle to a `text file`_ object. The read position is
        expected to be at the beginning (index 0).

    Returns:
        csv.Dialect: The detected tabular format.

    .. _text file:
        https://docs.python.org/3/glossary.html#term-text-file

    """
    peek = read_head(handle)
    handle.seek(0)
    sniffer = csv.Sniffer()
    dialect = sniffer.sniff(peek)
    return dialect


def check_samplesheet(file_in, file_out):
    """
    Check that the tabular samplesheet has the structure expected by the pipeline.

    Validate the general shape of the table, expected columns, and each row. Also add
    an additional column which records whether one or two FASTQ reads were found.

    Args:
        file_in (pathlib.Path): The given tabular samplesheet. The format can be either
            CSV, TSV, or any other format automatically recognized by ``csv.Sniffer``.
        file_out (pathlib.Path): Where the validated and transformed samplesheet should
            be created; always in CSV format.

    Example:
        This function checks that the samplesheet follows the following structure,
        see also the `pipeline documentation
        <https://github.com/bhargava-morampalli/rnamodbench/blob/dev/docs/usage.md#samplesheet-input>`_::

            sample,fastq,type,replicate,fast5_dir
            SAMPLE_1,SAMPLE_1.fastq.gz,native,rep1,/path/to/fast5
            SAMPLE_2,SAMPLE_2.fastq.gz,ivt,rep1,/path/to/fast5

    """
    required_columns = {"sample", "fastq", "type", "replicate", "fast5_dir"}
    # See https://docs.python.org/3.9/library/csv.html#id3 to read up on `newline=""`.
    with file_in.open(newline="") as in_handle:
        reader = csv.DictReader(in_handle, dialect=sniff_format(in_handle))
        # Validate the existence of the expected header columns.
        if not required_columns.issubset(reader.fieldnames):
            req_cols = ", ".join(required_columns)
            logger.critical(f"The sample sheet **must** contain these column headers: {req_cols}.")
            sys.exit(1)
        # Validate each row.
        checker = RowChecker()
        for i, row in enumerate(reader):
            try:
                checker.validate_and_transform(row)
            except AssertionError as error:
                logger.critical(f"{str(error)} On line {i + 2}.")
                sys.exit(1)
    header = list(reader.fieldnames)
    # See https://docs.python.org/3.9/library/csv.html#id3 to read up on `newline=""`.
    with file_out.open(mode="w", newline="") as out_handle:
        writer = csv.DictWriter(out_handle, header, delimiter=",")
        writer.writeheader()
        for row in checker.modified:
            writer.writerow(row)


def parse_args(argv=None):
    """Define and immediately parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Validate and transform a tabular samplesheet.",
        epilog="Example: python check_samplesheet.py samplesheet.csv samplesheet.valid.csv",
    )
    parser.add_argument(
        "file_in",
        metavar="FILE_IN",
        type=Path,
        help="Tabular input samplesheet in CSV or TSV format.",
    )
    parser.add_argument(
        "file_out",
        metavar="FILE_OUT",
        type=Path,
        help="Transformed output samplesheet in CSV format.",
    )
    parser.add_argument(
        "-l",
        "--log-level",
        help="The desired log level (default WARNING).",
        choices=("CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG"),
        default="WARNING",
    )
    return parser.parse_args(argv)


def main(argv=None):
    """Coordinate argument parsing and program execution."""
    args = parse_args(argv)
    logging.basicConfig(level=args.log_level, format="[%(levelname)s] %(message)s")
    if not args.file_in.is_file():
        logger.error(f"The given input file {args.file_in} was not found!")
        sys.exit(2)
    args.file_out.parent.mkdir(parents=True, exist_ok=True)
    check_samplesheet(args.file_in, args.file_out)


if __name__ == "__main__":
    sys.exit(main())
