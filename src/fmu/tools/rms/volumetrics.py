"""Module for handling volumetrics text files from RMS"""

import sys
import logging
import argparse
import signal
from pathlib import Path

import pandas as pd

logger = logging.getLogger(__name__)


def rmsvolumetrics_txt2df(
    txtfile,
    columnrenamer=None,
    phase=None,
    outfile=None,
    regionrenamer=None,
    zonerenamer=None,
):
    # pylint: disable=too-many-arguments
    """Parse the volumetrics txt file from RMS as Pandas dataframe

    Columns will be renamed according to FMU standard,
    https://wiki.equinor.com/wiki/index.php/FMU_standards

    Args:
        txtfile (str): path to file emitted by RMS Volumetrics job.
            Can also be a Path object.
        columnrenamer (dict): dictionary for renaming column. Will be merged
            with a default renaming dictionary (anything specified here will
            override any defaults)
        phase (str): stating typically 'GAS', 'OIL' or 'TOTAL', signifying
            what kind of data is in the file. Will be appended to column names,
            and is guessed from filename if not provided.
        outfile (str): filename to write CSV data to.
            If directory does not exist, it will be made.
        regionrenamer: a function that when applied on strings, return a
            new string. If used, will be applied to every region value,
            using pandas.Series.apply()
        zonerenamer: ditto for the zone column

    Return:
        pandas.DataFrame

    The renamer functions could be defined like this::

        def myregionrenamer(s):
            return s.replace('Equilibrium_region_', '')

    or the same using a lambda expression.
    """
    # First find out which row the data starts at:
    headerline = 0  # 0 is the first line
    for line in Path(txtfile).read_text().splitlines():
        if "Zone" in line or "Region" in line or "Facies" in line:
            break
        headerline = headerline + 1
    vol_df = pd.read_csv(txtfile, sep=r"\s\s+", skiprows=headerline, engine="python")

    # Enforce FMU standard:
    # https://wiki.statoil.no/wiki/index.php/FMU_standards
    # on column names

    # The Real column from RMS is not real.. Ignore it.
    if "Real" in vol_df.columns:
        vol_df.drop("Real", axis=1, inplace=True)

    if not phase:
        phase = guess_phase(txtfile)

    columns = {
        "Zone": "ZONE",
        "Region index": "REGION",
        "Facies": "FACIES",
        "License boundaries": "LICENSE",
        "Bulk": "BULK_" + phase,
        "Net": "NET_" + phase,
        "Hcpv": "HCPV_" + phase,
        "Pore": "PORV_" + phase,
        "Stoiip": "STOIIP_" + phase,
        "Giip": "GIIP_" + phase,
        "Assoc.Liquid": "ASSOCIATEDOIL_" + phase,
        "Assoc.Gas": "ASSOCIATEDGAS_" + phase,
    }
    if columnrenamer:
        # Overwrite with user supplied column conversion
        columns.update(columnrenamer)
    vol_df.rename(columns, axis=1, inplace=True)

    # Work on the data itself:
    if regionrenamer:
        vol_df["REGION"] = vol_df["REGION"].apply(regionrenamer)
    if zonerenamer:
        vol_df["ZONE"] = vol_df["ZONE"].apply(zonerenamer)

    # Remove the Totals rows in case they are present.
    #  (todo: do this for all columns that are not not of numeric type)
    checkfortotals = ["ZONE", "REGION", "LICENSE", "FACIES"]

    totalsrows = pd.Series([False] * len(vol_df))
    for col in checkfortotals:
        if col in vol_df.columns:
            totalsrows = totalsrows | (vol_df[col] == "Totals")
    vol_df = vol_df[~totalsrows].reset_index(drop=True)

    if outfile:
        Path(outfile).parent.mkdir(exist_ok=True, parents=True)
        vol_df.to_csv(outfile, index=False)

    return vol_df


def guess_phase(text):
    """From a text-file, guess which phase the text file
    concerns, oil, gas or the "total" phase.

    Args:
        text (str): Multiline

    Returns:
        str: "OIL", "GAS" or "TOTAL"

    Raises:
        ValueError if guessing fails.
    """
    if "oil" in str(text).lower():
        return "OIL"
    if "gas" in str(text).lower():
        return "GAS"
    if "total" in str(text).lower():
        return "TOTAL"
    raise ValueError("Not able to guess phase")


def get_parser():
    """Set up parser for command line utility"""
    parser = argparse.ArgumentParser(
        description="Convert single RMS volumetrics files to CSV"
    )

    parser.add_argument("rmsvolumetricsfile", type=str)
    parser.add_argument(
        "--phase",
        type=str,
        help=(
            "Set typically to 'GAS', 'OIL' or 'TOTAL'. "
            "Will be appended to column headers. "
            "If not given, it will be guessed from the filename"
        ),
        default=None,
    )
    parser.add_argument(
        "-o",
        "--output",
        type=str,
        help="Name of output CSV file. Use '-' for stdout (default)",
        default="-",
    )
    parser.add_argument("-v", "--verbose", action="store_true", help="Be verbose")
    return parser


def rmsvolumetrics2csv_main():
    """Endpoint for command line utility for converting one file at a time"""
    parser = get_parser()
    args = parser.parse_args()

    logging.basicConfig()

    if args.verbose:
        logger.setLevel(logging.INFO)

    vol_df = rmsvolumetrics_txt2df(args.rmsvolumetricsfile, phase=args.phase)

    if args.output == "-":
        # Ignore pipe errors when writing to stdout.
        signal.signal(signal.SIGPIPE, signal.SIG_DFL)
        vol_df.to_csv(sys.stdout, index=False)
    else:
        Path(args.output).parent.mkdir(exist_ok=True, parents=True)
        vol_df.to_csv(args.output, index=False)
        print("Wrote to " + args.output)
