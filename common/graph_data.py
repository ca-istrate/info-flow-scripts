import os

import logging

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from common.pairing import Pairing

logger = logging.getLogger("script")

INCH_TO_CM = 2.54 

def default_output_file(args) -> None:
    inp = os.path.splitext(os.path.basename(args.x))[0]
    inp += "-" + os.path.splitext(os.path.basename(args.y))[0] if args.y is not None else ""
    inp += "-" + os.path.splitext(os.path.basename(args.z))[0] if args.z is not None else ""
    inp += "-" + os.path.splitext(os.path.basename(args.w))[0] if args.w is not None else ""
    vars(args)["output_file"] = f"{inp}-{args.start}-{args.end}-graph.out"
    return vars(args)["output_file"]


def plot_data(args, data1: Pairing, data2: Pairing, count):

    plt.rcParams.update({
        # "text.usetex": True,
        "font.size": 11,
        "font.family": "serif",
        "font.serif": ["Times New Roman"]
    })

    plot_number = count * 100 + 10

    data_pairs = [
        (data1.x, data1.x_time, args.x_conv, args.x_unit, args.x_label, args.x_domain, args.x_color),
        (data1.y, data1.y_time, args.y_conv, args.y_unit, args.y_label, args.y_domain, args.y_color),
        (data2.x, data2.x_time, args.z_conv, args.z_unit, args.z_label, args.z_domain, args.z_color),
        (data2.y, data2.y_time, args.w_conv, args.w_unit, args.w_label, args.w_domain, args.w_color)
    ]


    fig = plt.figure(1)
    axes = {}

    for i, (data, time, conv, unit, label, domain, clr) in enumerate(data_pairs):
        if i + 1 > count:
            break

        axes[i] = (plt.subplot(plot_number+i+1, sharex=axes[0][0] if i != 0 else None), domain)

        axes[i][0].plot(time * args.t_conv, data * conv, color=clr, label=f"{label} ({unit})")

        # axes[i][0].set_ylabel(f"{label} ({unit})")
        axes[i][0].grid(axis='x', which="both")

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

    locs = axes[0][0].get_xticks()[1:-1]
    new_locs = np.linspace(locs[0], locs[-1], len(locs) * 2 - 1)
    axes[0][0].set_xticks(new_locs)

    plt.figlegend(loc="upper center", ncols=4, columnspacing=1, handleheight=1, edgecolor="white")

    # fig = plt.gcf()
    # fig.set_size_inches(16.5 / INCH_TO_CM, 8.5 / INCH_TO_CM)
    fig.set_size_inches(args.fig_width / INCH_TO_CM, args.fig_height / INCH_TO_CM)
    fig.subplots_adjust(bottom=args.fig_margin_bottom, top=args.fig_margin_top, hspace=args.fig_margin_hspace,
                        left=args.fig_margin_left, right=args.fig_margin_right)
    # fig.set_layout_engine("constrained")


    # fig.subplots_adjust(left=0.1, bottom=0.15, right=0.975, top=0.8, hspace=0)

    if args.reverse:
        plt.gca().invert_xaxis()

    fig.savefig(f'{args.output_file}.svg', dpi=200)
    if args.plot_show is "True":
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

    

