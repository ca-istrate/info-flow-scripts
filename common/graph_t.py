import os

import logging

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


logger = logging.getLogger("script")


def default_output_file(args) -> None:
    inp = args.csv
    vars(args)["output_file"] = f"{inp}-graph-t.out"
    return vars(args)["output_file"]


def get_data(args, csv: str):
    data = pd.read_csv(csv) 

    data.loc[:, ['x_time_start', 'x_time_end']] *= args.t_conv

    return data


def plot_one(args, data: pd.DataFrame):
    logger.debug(data)

    ax = plt.subplot(211)

    if args.graph_by_shift:
        midpoint = np.array(data.loc[:, 'shift'])
        labels = list(map(lambda x: f"{x:.1f}", midpoint * args.t_conv)) 
        offset = 0.1
    else:
        midpoint = np.array(data.loc[:, 'x_time_start'] + data.loc[:, 'x_time_end']) / 2
        labels = list(map(lambda x: f"{x[0]:.1f}-{x[1]:.1f}", zip(data.loc[:, 'x_time_start'], data.loc[:, 'x_time_end'])))
        offset = args.t_conv

    logger.debug(midpoint)

    ax.axhline(y=0, ls="--", lw=0.5, c="black")

    ax.errorbar(midpoint + offset, data.loc[:, 'T12_theo'], yerr=np.abs(data.loc[:, 'T12_int_dist']), fmt="o", label=r"$T_{" +f"{args.x_label}" + "-" + f"{args.y_label}" + r"}$")
    ax.errorbar(midpoint - offset, data.loc[:, 'T21_theo'], yerr=np.abs(data.loc[:, 'T21_int_dist']), fmt="o", label=r"$T_{" +f"{args.y_label}" + "-" + f"{args.x_label}" + r"}$")

    ax.set_ylabel("Absolute information flow (nats/ut)")
    ax.spines['bottom'].set_visible(False)
    plt.setp(ax.get_xticklabels(), visible=False)
    plt.setp(ax.get_xticklines(), visible=False)

    ax.legend()

    ax2 = plt.subplot(212, sharex=ax)
    ax2.axhline(y=0, ls="--", lw=0.5, c="black")
    ax2.errorbar(midpoint + offset, data.loc[:, 'tau12'], yerr=np.abs(data.loc[:, 'error_tau12']) * 1.96, fmt="o", label=r"$\tau_{" +f"{args.x_label}" + "-" + f"{args.y_label}" + r"}$")
    ax2.errorbar(midpoint - offset, data.loc[:, 'tau21'], yerr=np.abs(data.loc[:, 'error_tau21']) * 1.96, fmt="o", label=r"$\tau_{" +f"{args.y_label}" + "-" + f"{args.x_label}" + r"}$")
    ax2.set_ylabel("Relative information flow (%)")
    ax2.set_xticks(midpoint[::2], labels[::2])
    ax2.set_xlabel(f"{args.t_label} ({args.t_unit})")
    ax2.spines['top'].set_visible(False)
    ax2.legend()

    plt.subplots_adjust(hspace=0)
    if args.reverse:
        plt.gca().invert_xaxis()

    plt.show()



def main(args) -> None:
    logger.info("Graphing T values.")
    
    logger.debug(f"Input arguments: {args}")

    plt.rcParams.update({
        "font.family": "DejaVu Sans"
    })

    data = get_data(args, args.csv)

    plot_one(args, data)


