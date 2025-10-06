import os

import logging

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from common.pairing import Pairing

logger = logging.getLogger("script")

def default_output_file(args) -> None:
    inp = os.path.splitext(os.path.basename(args.x))[0]
    inp += "-" + os.path.splitext(os.path.basename(args.y))[0] if args.y is not None else ""
    inp += "-" + os.path.splitext(os.path.basename(args.z))[0] if args.z is not None else ""
    inp += "-" + os.path.splitext(os.path.basename(args.w))[0] if args.w is not None else ""
    vars(args)["output_file"] = f"{inp}-{args.start}-{args.end}-graph.out"
    return vars(args)["output_file"]


def plot_data(args, data1: Pairing, data2: Pairing, count):

    plot_number = count * 100 + 10

    data_pairs = [
        (data1.x, data1.x_time, args.x_conv, args.x_unit, args.x_label, args.x_domain),
        (data1.y, data1.y_time, args.y_conv, args.y_unit, args.y_label, args.y_domain),
        (data2.x, data2.x_time, args.z_conv, args.z_unit, args.z_label, args.z_domain),
        (data2.y, data2.y_time, args.w_conv, args.w_unit, args.w_label, args.w_domain)
    ]

    axes = {}

    for i, (data, time, conv, unit, label, domain) in enumerate(data_pairs):
        if i + 1 > count:
            break

        axes[i] = (plt.subplot(plot_number+i+1, sharex=axes[0][0] if i != 0 else None), domain)

        axes[i][0].plot(time * args.t_conv, data * conv)

        axes[i][0].set_ylabel(f"{label} ({unit})")
        axes[i][0].grid(axis='x')

        if i != 0:
            axes[i][0].spines['top'].set_visible(False)

        if i + 1 == count:
            axes[i][0].set_xlabel(f"{args.t_label} ({args.t_unit})")
            break
        else:
            axes[i][0].spines['bottom'].set_visible(False)
            plt.setp(axes[i][0].get_xticklabels(), visible=False)
            plt.setp(axes[i][0].get_xticklines(), visible=False)

    dom_lim = {}
    _, doms = zip(*axes.values())
    for dom in doms:
        lows = [x[0].get_ylim()[0] for x in axes.values() if x[1] == dom]
        highs = [x[0].get_ylim()[1] for x in axes.values() if x[1] == dom]
        low = min(lows)
        high = max(highs)
        dom_lim[dom] = (low, high)

    for ax, dom in axes.values():
        ax.set_ylim(*dom_lim[dom])


    plt.subplots_adjust(hspace=0)

    if args.reverse:
        plt.gca().invert_xaxis()

    plt.show()
    



def main(args) -> None:
    logger.info("Graphing data series.")
    
    logger.debug(f"Input arguments: {args}")

    y = args.y
    z = args.z
    w = args.w

    if args.y is None:
        y = args.x
        z = y
        w = z
        count = 1
    elif args.z is None:
        z = args.x
        w = z
        count = 2
    elif args.w is None:
        w = args.x
        count = 3
    else:
        count = 4

    data1 = Pairing.from_files(args.x, y, start=args.start, end=args.end, shift=0)
    data2 = Pairing.from_files(z, w, start=args.start, end=args.end, shift=0)

    plot_data(args, data1, data2, count)

    

