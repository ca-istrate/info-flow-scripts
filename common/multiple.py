import logging
import os

import pandas as pd
import numpy as np
from scipy.stats import norm

from common.pairing import Pairing
import common.theoretical as theo


logger = logging.getLogger("script")


def default_output_file(args) -> None:
    x, _ = os.path.splitext(os.path.basename(args.x))
    y, _ = os.path.splitext(os.path.basename(args.y))
    window_size = 'full' if args.window_size is None else args.window_size
    window_start = args.window_start
    window_end = 'all' if args.window_end is None else args.window_end
    window_step = args.window_step
    shift_start = args.shift_start
    shift_end = args.shift_end
    shift_step = args.shift_step
    vars(args)["output_file"] = f"{x}-{y}-{window_size}-{window_start}-{window_end}-{window_step}-{shift_start}-{shift_end}-{shift_step}-multiple.out"
    return vars(args)["output_file"]


def main(args) -> None:
    logger.info("Running windowed info-flow calculation.")
    
    logger.debug(f"Input arguments: {args}")
    
    data = Pairing.from_files(args.x, args.y)
    results = []
    
    for shift in range(args.shift_start, args.shift_end+1, args.shift_step):
        logger.debug(f"Testing shift = {shift}")
        shifted_data = data.subset(shift=shift)

        if args.window_size is None and args.window_end is None:
            logger.debug(f"Run single anaysis for entire data")
            single_row, _ = theo.analyze(args, shifted_data)
            single_row["shift"] = [shift]
            results.append(single_row)
            continue

        max_size = len(shifted_data.x)

        if args.window_end is None:
            window_end = max_size - args.window_size
        else:
            window_end = args.window_end
            if args.window_size is not None and window_end + args.window_size > max_size:
                window_end = max_size - args.window_size
                logging.warning(f"Warning: not enough data for the specified window size and window end values. Using window_end of {window_end} instead.")

        if window_end <= 0:
            logging.warning(f"Value of window_end is not positive. No analysies for this data/shift combo are performed. Shift = {shift}")
            continue

        for start in range(args.window_start, window_end, args.window_step):
            logger.debug(f"Run anaysis for window starting at {start}")
            end = None if args.window_size is None else start + args.window_size
            data_subset = shifted_data.subset(start=start, end=end)
            single_row, _ = theo.analyze(args, data_subset)
            single_row["shift"] = [shift]
            results.append(single_row)

    results = pd.concat(results, ignore_index=True)
    results.to_csv(args.output_file+".csv", index=False)



