#!/usr/bin/env python
import argparse
import logging

import common.linearity as lin
import common.theoretical as theo
import common.multiple as many
import common.graph_data as graph_data
import common.graph_t as graph_t
from common.log import setup_logging


logger = logging.getLogger("script")

actions = {
    "linearity": {
        "message": "Perform linearity check for the given data.",
        "main": lin.main,
        "output": lin.default_output_file,
        "arguments": {
            "x": {
                "help": "Path to x data file"
            },
            "y": {
                "help": "Path to y data file"
            },
        },
        "options": {
            ("--start", "-s"): {
                "help": "The start index for the analysis",
                "type": int,
                "default": 0
            },
            ("--end", "-e"): {
                "help": "The end index for the analysis (non-inclusive)",
                "type": int,
                "default": None
            },
            ("--shift", "-i"): {
                "help": "The shift between y and x. For a positive value x is ahead",
                "type": int,
                "default": 0
            },
            ("--multiple", "-m"): {
                "help": "Run multiple fits of size specified by --fit-size.",
                "type": bool,
                "default": False
            },
            ("--plot", "-p"): {
                "help": "Plot resulting regression.",
                "type": bool,
                "default": False
            },
            ("--fit-size", "-f"): {
                "help": "When performing multiple fits, what size to use.",
                "type": int,
                "default": 100
            },
            ("--fit-increment", "-j"): {
                "help": "When performing multiple fits, how many steps to jump.",
                "type": int,
                "default": 10
            },
            ("--fit-count", "-c"): {
                "help": "When performing multiple fits, top how many results to display.",
                "type": int,
                "default": 10
            },
        }
    },
    "analyze": {
        "message": "Calculate the causality using both the theoretical and bootstrap method.",
        "main": theo.main,
        "output": theo.default_output_file,
        "arguments": {
            "x": {
                "help": "Path to x data file"
            },
            "y": {
                "help": "Path to y data file"
            },
        },
        "options": {
            ("--start", "-s"): {
                "help": "The start index for the analysis",
                "type": int,
                "default": 0
            },
            ("--end", "-e"): {
                "help": "The end index for the analysis (non-inclusive)",
                "type": int,
                "default": None
            },
            ("--shift", "-i"): {
                "help": "The shift between y and x. For a positive value x is ahead",
                "type": int,
                "default": 0
            },
            ("--alpha", "-a"): {
                "help": "Significance level (alpha) used in the analysis.",
                "type": float,
                "default": 0.05
            },
            ("--bootstrap-iter", "-b"): {
                "help": "Number of bootstraping rounds to perform.",
                "type": int,
                "default": 1000
            },
            ("--bootstrap-conf-val", "-c"): {
                "help": "Confidence value for the bootstraping method.",
                "type": float,
                "default": 1.96
            },
        }
    },
    "analyze-many": {
        "message": "Calculate the causality for many subwindows of the provided data.",
        "main": many.main,
        "output": many.default_output_file,
        "arguments": {
            "x": {
                "help": "Path to x data file"
            },
            "y": {
                "help": "Path to y data file"
            },
        },
        "options": {
            ("--alpha", "-a"): {
                "help": "Significance level (alpha) used in the analysis.",
                "type": float,
                "default": 0.05
            },
            ("--bootstrap-iter", "-b"): {
                "help": "Number of bootstraping rounds to perform.",
                "type": int,
                "default": 1000
            },
            ("--bootstrap-conf-val", "-c"): {
                "help": "Confidence value for the bootstraping method.",
                "type": float,
                "default": 1.96
            },
            ("--window-size", "-w"): {
                "help": "Size of the windows in steps to perform the analysis on.",
                "type": int,
                "default": None
            },
            ("--window-start", "-ws"): {
                "help": "First starting position in steps when perfoming windowed analysies",
                "type": int,
                "default": 0
            },
            ("--window-end", "-we"): {
                "help": "Last starting position in steps when perfoming windowed analysies (non-inclusive)",
                "type": int,
                "default": None
            },
            ("--window-step", "-wt"): {
                "help": "The amount of steps to move the window between analysies",
                "type": int,
                "default": 5
            },
            ("--shift-start", "-ss"): {
                "help": "First shift value to analyze",
                "type": int,
                "default": 0
            },
            ("--shift-end", "-se"): {
                "help": "Last shift value to analyze (inclusive)",
                "type": int,
                "default": 0
            },
            ("--shift-step", "-st"): {
                "help": "The number of steps to skip between shift analysies",
                "type": int,
                "default": 1
            },
        }
    },
    "graph-data": {
        "message": "Graph the raw data.",
        "main": graph_data.main,
        "output": graph_data.default_output_file,
        "arguments": {
            "x": {
                "help": "Path to x data file",
            },
            "y": {
                "help": "Path to y data file",
                "default": None,
                "nargs": '?'
            },
            "z": {
                "help": "Path to z data file",
                "default": None,
                "nargs": '?'
            },
            "w": {
                "help": "Path to w data file",
                "default": None,
                "nargs": '?'
            },
        },
        "options": {
            ("--start", "-s"): {
                "help": "The start index for the analysis",
                "type": int,
                "default": 0
            },
            ("--end", "-e"): {
                "help": "The end index for the analysis (non-inclusive)",
                "type": int,
                "default": None
            },
            ("--reverse", "-r"): {
                "help": "Reverse data",
                "type": bool,
                "default": False
            },
            ("--t-conv", "-tc"): {
                "help": "Conversion factor for time",
                "type": float,
                "default": 1
            },
            ("--t-unit", "-tu"): {
                "help": "Unit of time",
                "type": str,
                "default": "kyr"
            },
            ("--t-label", "-tl"): {
                "help": "Label time",
                "type": str,
                "default": "Time"
            },
            ("--x-conv", "-xc"): {
                "help": "Conversion factor for x data",
                "type": float,
                "default": 1
            },
            ("--x-unit", "-xu"): {
                "help": "Unit of x data",
                "type": str,
                "default": "unit"
            },
            ("--x-label", "-xl"): {
                "help": "Label of x data",
                "type": str,
                "default": ""
            },
            ("--x-domain", "-xd"): {
                "help": "Group for sharing y axis range",
                "type": int,
                "default": 1
            },
            ("--y-conv", "-yc"): {
                "help": "Conversion factor for y data",
                "type": float,
                "default": 1
            },
            ("--y-unit", "-yu"): {
                "help": "Unit of y data",
                "type": str,
                "default": "unit"
            },
            ("--y-label", "-yl"): {
                "help": "Label of y data",
                "type": str,
                "default": ""
            },
            ("--y-domain", "-yd"): {
                "help": "Group for sharing y axis range",
                "type": int,
                "default": 2
            },
            ("--z-conv", "-zc"): {
                "help": "Conversion factor for z data",
                "type": float,
                "default": 1
            },
            ("--z-unit", "-zu"): {
                "help": "Unit of z data",
                "type": str,
                "default": "unit"
            },
            ("--z-label", "-zl"): {
                "help": "Label of z data",
                "type": str,
                "default": ""
            },
            ("--z-domain", "-zd"): {
                "help": "Group for sharing y axis range",
                "type": int,
                "default": 3
            },
            ("--w-conv", "-wc"): {
                "help": "Conversion factor for w data",
                "type": float,
                "default": 1
            },
            ("--w-unit", "-wu"): {
                "help": "Unit of w data",
                "type": str,
                "default": "unit"
            },
            ("--w-label", "-wl"): {
                "help": "Label of w data",
                "type": str,
                "default": ""
            },
            ("--w-domain", "-wd"): {
                "help": "Group for sharing y axis range",
                "type": int,
                "default": 4
            },
        }
    },
    "graph-t": {
        "message": "Graph the T values.",
        "main": graph_t.main,
        "output": graph_t.default_output_file,
        "arguments": {
            "csv": {
                "help": "Path to csv result files",
            },
        },
        "options": {
            ("--reverse", "-r"): {
                "help": "Reverse data",
                "type": bool,
                "default": False
            },
            ("--graph-by-shift", "-gs"): {
                "help": "Use shift as x axis instead of time",
                "type": bool,
                "default": False
            },
            ("--t-conv", "-tc"): {
                "help": "Conversion factor for time",
                "type": float,
                "default": 1
            },
            ("--t-unit", "-tu"): {
                "help": "Unit of time",
                "type": str,
                "default": "kyr"
            },
            ("--t-label", "-tl"): {
                "help": "Label time",
                "type": str,
                "default": "Time"
            },
            ("--x-label", "-xl"): {
                "help": "Label for first data series",
                "type": str,
                "default": "1"
            },
            ("--y-label", "-yl"): {
                "help": "Label for second data series",
                "type": str,
                "default": "2"
            },
        }
    }
}


def manual_args():
    args = argparse.Namespace()

    vars(args)["debug"] = False

    print("Select which task to perform:")

    actions_list = list(actions.items())
    for number, (_, details) in enumerate(actions_list):
        print(f"  ({number + 1}) - {details['message']}")

    selected = -1
    while not 0 <= selected < len(actions):
        try:
            selected = int(input("Enter the number corresponding to the task: ")) - 1
        except ValueError:
            pass

    print()

    selected_item = actions_list[selected][1]
    vars(args)["main_function"] = selected_item["main"]
    vars(args)["output_function"] = selected_item["output"]

    for arg, details in selected_item["arguments"].items():
        vars(args)[arg] =  input(f"{details['help']}{'(default: ' + str(details['default']) + ')' if 'default' in details else ''}: ")
        if "default" in details and vars(args)[arg] == "":
            vars(args)[arg] = details["default"]

    for opt, details in selected_item["options"].items():
        opt_name = opt[0].strip("-").replace("-", "_")
        vars(args)[opt_name] = input(f"{details['help']}{f'(default: ' + str(details['default']) + ')' if 'default' in details else ''}: ")
        if "default" in details and vars(args)[opt_name] == "":
            vars(args)[opt_name] = details["default"]
        else:
            vars(args)[opt_name] = details["type"](vars(args)[opt_name])



    vars(args)["output_file"] = input(f"Output file ({args.output_function(args)}): ")
    if vars(args)["output_file"] == "":
        vars(args)["output_file"] = None

    return args



def parse_arguments():
    parser = argparse.ArgumentParser(prog="info-flow")
    parser.add_argument("--output_file", "-o", required=False, type=str, default=None)
    parser.add_argument("--debug", "-d", action=argparse.BooleanOptionalAction, default=False)

    parser.set_defaults(main_function=None)
     
    subparers = parser.add_subparsers()

    for action, details in actions.items():
        subparser = subparers.add_parser(action)
        for argument, arg_details in details["arguments"].items():
            subparser.add_argument(argument, **arg_details)
        for option, opt_details in details["options"].items():
            subparser.add_argument(*option, **opt_details)
        subparser.set_defaults(main_function=details["main"], output_function=details["output"])
    
    args = parser.parse_args()

    if args.main_function is None:
        args = manual_args()

    if args.output_file is None:
        args.output_function(args)

    return args

   
def main():
    args = parse_arguments()

    setup_logging(args)

    args.main_function(args)


if __name__ == "__main__":
    main()

