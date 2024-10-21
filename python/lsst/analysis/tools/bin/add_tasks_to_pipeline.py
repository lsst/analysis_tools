# This file is part of analysis_tools.
#
# Developed for the LSST Data Management System.
# This product includes software developed by the LSST Project
# (https://www.lsst.org).
# See the COPYRIGHT file at the top-level directory of this distribution
# for details of code ownership.
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

from __future__ import annotations

import logging
import os
import time
from argparse import SUPPRESS, ArgumentParser

from ..utils import add_tasks_to_pipeline


def build_argparser():
    """Build an argument parser for this script."""
    parser = ArgumentParser(
        description="""Add task(s) to a reference pipeline. This script takes a
        reference pipeline definition file in YAML format and adds tasks from a
        list of input pipeline files.
        """,
        epilog="More information is available at https://pipelines.lsst.io.",
        add_help=False,
        argument_default=SUPPRESS,
    )
    parser.add_argument(
        "-r",
        "--reference-pipeline",
        type=str,
        help="Location of a reference pipeline definition YAML file.",
        required=True,
        metavar="FILE",
    )
    parser.add_argument(
        "-i",
        "--input-pipelines",
        type=str,
        help="Location(s) of input pipeline definition YAML file(s). Tasks from input_pipelines will be added"
        " to reference_pipeline.",
        required=True,
        metavar="FILE",
        nargs="+",
    )
    parser.add_argument(
        "-s",
        "--subset_name",
        type=str,
        help="""All tasks from input_pipelines will be added to this subset. If the
        subset does not exist it will be created.""",
        metavar="FILE",
    )
    parser.add_argument(
        "-d",
        "--new_subset_description",
        type=str,
        help="The description for the new subset.",
        metavar="FILE",
    )
    parser.add_argument(
        "-f",
        "--filename",
        help="Path to save a modified pipeline definition YAML file.",
        metavar="FILE",
    )
    parser.add_argument(
        "--overwrite",
        help="Overwrite the output saved pipeline definition file if it already exists.",
        action="store_true",
    )
    parser.add_argument(
        "--instrument",
        type=str,
        help="Add instrument overrides. Must be a fully qualified class name.",
        metavar="instrument",
    )
    parser.add_argument(
        "-h",
        "--help",
        action="help",
        help="Show this help message and exit.",
    )
    return parser


def main():
    """Use this as the main entry point when calling from the command line."""
    # Set up logging.
    tz = time.strftime("%z")
    logging.basicConfig(
        format="%(levelname)s %(asctime)s.%(msecs)03d" + tz + " - %(message)s",
        datefmt="%Y-%m-%dT%H:%M:%S",
    )
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.DEBUG)

    args = vars(build_argparser().parse_args())
    filename = args.pop("filename", None)
    overwrite = args.pop("overwrite", None)
    pipeline = add_tasks_to_pipeline(**args)

    if filename:
        if os.path.exists(filename) and not overwrite:
            raise RuntimeError(f"File {filename} already exists; use --overwrite to write anyway.")
        pipeline.write_to_uri(filename)
        logger.info(
            "Modified pipeline definition YAML file saved at %s.",
            os.path.realpath(filename),
        )
    else:
        print("\n", str(pipeline), sep="")
