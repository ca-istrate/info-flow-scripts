import logging
import os

import pandas as pd
import matplotlib.pyplot as plt
from statsmodels.formula.api import ols

from common.pairing import Pairing


logger = logging.getLogger("script")


def default_output_file(args) -> None:
    x, _ = os.path.splitext(os.path.basename(args.x))
    y, _ = os.path.splitext(os.path.basename(args.y))
    vars(args)["output_file"] = f"{x}-{y}-{args.start}-{args.end}-{args.shift}-linearity.out"
    return vars(args)["output_file"]


def fit_bic(data, start, end):
    logger.debug(f"Performing linear check from {start} to {end}")

    df = pd.DataFrame({
        "x": data.x[start:end],
        "y": data.y[start:end]
    })

    model = ols("y ~ x", data=df)
    results = model.fit()

    logger.debug("Fitted linear model with BIC {results.bic}.")
    return (start, end, results.bic)


def fit_multiple(data, args):

    fits = []

    l = max(len(data.x), len(data.y))

    logger.info(f"Fitting linear models of size {args.fit_size} every {args.fit_increment} steps.")
    for start in range(0, l + 1 - args.fit_size, args.fit_increment):
        fits.append(fit_bic(data, start, start + args.fit_size))

    logger.info(f"Results based on Bayesian Information Critertion")
    
    count = args.fit_count if len(fits) > args.fit_count else -1
    fits = sorted(fits, key=lambda x: x[2])[:count]
    for fit in fits:
        logger.info(f"[{fit[0]:>6} - {fit[1]:>6}]  - {fit[2]:>9.2f}" )


def fit_single(data, args):

    df = pd.DataFrame({
        "x": data.x,
        "y": data.y
    })

    logger.info("Fitting linear model.")
    model = ols("y ~ x", data=df)
    results = model.fit()

    logger.info("Results:")
    summary = results.summary()
    logger.info(summary)

    if args.plot:
        min_x = min(df["x"]) - 1
        max_x = max(df["x"]) + 1
        xs = pd.DataFrame({
            "x": [min_x, max_x]
        })

        plt.scatter(df["x"], df["y"])
        plt.plot(xs, results.predict(xs))
        plt.xlim([-5, 5])
        plt.ylim([-5, 5])
        plt.show()


def main(args) -> None:
    logger.info("Running linearity check.")
    
    logger.debug(f"Input arguments: {args}")
    
    data = Pairing.from_files(args.x, args.y, start=args.start, end=args.end, shift=args.shift)

    if args.multiple:
        fit_multiple(data, args)
    else:
        fit_single(data, args)
    

