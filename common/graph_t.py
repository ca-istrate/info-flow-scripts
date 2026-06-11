import os

import logging

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


logger = logging.getLogger("script")


INCH_TO_CM = 2.54 


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

    plt.rcParams.update({
        # "text.usetex": True,
        "mathtext.fontset": "stix",
        "mathtext.default": "bf",
        "font.size": 11,
        "font.family": "serif",
        "font.serif": ["Times New Roman"],
        "text.usetex": True,
        "text.latex.preamble": r"\usepackage{amsmath}\usepackage{lmodern}\usepackage{mhchem}\usepackage{textcomp}"
    })


    fig = plt.figure(1)
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

    shift_alpha = 0.4

    # if args.graph_by_shift:
    t1_pairing = list(zip(data.loc[:, 'shift'], midpoint + offset, data.loc[:, 'T12_theo'], np.abs(data.loc[:, 'T12_int_dist'])))
    positive_t1_pairing = [point for point in t1_pairing if point[0] <= 0]
    negative_t1_pairing = [point for point in t1_pairing if point[0] > 0]

    if positive_t1_pairing:
        _, t, x, x_err = zip(*positive_t1_pairing)
        ax.errorbar(np.array(t), np.array(x), np.array(x_err), fmt="o", label=r"$\mathrm{T_{" +f"{args.x_label}" + r"\rightarrow " + f"{args.y_label}" + r"}}$", color=args.x_color, markersize=3)

    if negative_t1_pairing:
        _, t, x, x_err = zip(*negative_t1_pairing)
        ax.errorbar(np.array(t), np.array(x), np.array(x_err), fmt="o", color=args.x_color, markersize=3, alpha=shift_alpha)


    t2_pairing = list(zip(data.loc[:, 'shift'], midpoint - offset, data.loc[:, 'T21_theo'], np.abs(data.loc[:, 'T21_int_dist'])))
    positive_t2_pairing = [point for point in t2_pairing if point[0] < 0]
    negative_t2_pairing = [point for point in t2_pairing if point[0] >= 0]

    if negative_t2_pairing:
        _, t, y, y_err = zip(*negative_t2_pairing)
        ax.errorbar(np.array(t), np.array(y), np.array(y_err), fmt="o", label=r"$\mathrm{T_{" +f"{args.y_label}" + r"\rightarrow " + f"{args.x_label}" + r"}}$", color=args.y_color, markersize=3)

    if positive_t2_pairing:
        _, t, y, y_err = zip(*positive_t2_pairing)
        ax.errorbar(np.array(t), np.array(y), np.array(y_err), fmt="o", color=args.y_color, markersize=3, alpha=shift_alpha)

    # else:
    #     ax.errorbar(midpoint + offset, data.loc[:, 'T12_theo'], yerr=np.abs(data.loc[:, 'T12_int_dist']), fmt="o", label=r"$\mathrm{T_{" +f"{args.x_label}" + r"\rightarrow " + f"{args.y_label}" + r"}}$", color=args.x_color, markersize=3)
    #     ax.errorbar(midpoint - offset, data.loc[:, 'T21_theo'], yerr=np.abs(data.loc[:, 'T21_int_dist']), fmt="o", label=r"$\mathrm{T_{" +f"{args.y_label}" + r"\rightarrow " + f"{args.x_label}" + r"}}$", color=args.y_color, markersize=3)

    ax.set_ylabel("T (nats/ut)")
    ax.spines['bottom'].set_visible(False)
    plt.setp(ax.get_xticklabels(), visible=False)
    plt.setp(ax.get_xticklines(), visible=False)

    # ax.legend()

    ax2 = plt.subplot(212, sharex=ax)
    ax2.axhline(y=0, ls="--", lw=0.5, c="black")

    # if args.graph_by_shift:
    tau1_pairing = list(zip(data.loc[:, 'shift'], midpoint + offset, data.loc[:, 'tau12'], np.abs(data.loc[:, 'error_tau12']) * 1.96 ))
    positive_tau1_pairing = [point for point in tau1_pairing if point[0] <= 0]
    negative_tau1_pairing = [point for point in tau1_pairing if point[0] > 0]

    if positive_tau1_pairing:
        _, t, x, x_err = zip(*positive_tau1_pairing)
        ax2.errorbar(np.array(t), np.array(x), np.array(x_err), fmt="s", color=args.x_color, markersize=3)

    if negative_tau1_pairing:
        _, t, x, x_err = zip(*negative_tau1_pairing)
        ax2.errorbar(np.array(t), np.array(x), np.array(x_err), fmt="s", color=args.x_color, markersize=3, alpha=shift_alpha)

    tau2_pairing = list(zip(data.loc[:, 'shift'], midpoint - offset, data.loc[:, 'tau21'], np.abs(data.loc[:, 'error_tau21']) * 1.96 ))
    positive_tau2_pairing = [point for point in tau2_pairing if point[0] < 0]
    negative_tau2_pairing = [point for point in tau2_pairing if point[0] >= 0]

    if negative_tau2_pairing:
        _, t, y, y_err = zip(*negative_tau2_pairing)
        ax2.errorbar(np.array(t), np.array(y), np.array(y_err), fmt="s", color=args.y_color, markersize=3)

    if positive_tau2_pairing:
        _, t, y, y_err = zip(*positive_tau2_pairing)
        ax2.errorbar(np.array(t), np.array(y), np.array(y_err), fmt="s", color=args.y_color, markersize=3, alpha=shift_alpha)

    # else:
    #     ax2.errorbar(midpoint + offset, data.loc[:, 'tau12'], yerr=np.abs(data.loc[:, 'error_tau12']) * 1.96, fmt="s", label=r"$\mathrm{\tau_{" +f"{args.x_label}" + r"\rightarrow " + f"{args.y_label}" + r"}}$", color=args.x_color, markersize=3)
    #     ax2.errorbar(midpoint - offset, data.loc[:, 'tau21'], yerr=np.abs(data.loc[:, 'error_tau21']) * 1.96, fmt="s", label=r"$\mathrm{\tau_{" +f"{args.y_label}" + r"\rightarrow " + f"{args.x_label}" + r"}}$", color=args.y_color, markersize=3)


    ax2.set_ylabel(r"$\mathrm{\tau}$ (\%)")
    ax2.set_xticks(midpoint[::2], labels[::2])
    ax2.set_xlabel(f"{args.t_label} ({args.t_unit})")
    ax2.spines['top'].set_visible(False)

    lines3, labels3 = [], []
    if args.plot_r:
        ax3 = ax2.twinx()
        ax3.scatter(midpoint, data.loc[:, 'R12'], marker="x", label="R", color="black", s=args.marker_size)
        ax3.yaxis.set_major_formatter(plt.FormatStrFormatter('%.2f'))
        ax3.set_ylabel("Correlation Coefficient (R)")
        ax3.spines['top'].set_visible(False)
        lines3, labels3 = ax3.get_legend_handles_labels()


    lines = [
        Line2D([], [], marker='_', color=args.x_color, markersize=args.marker_size, linestyle='None', label=f"${args.x_label}" + r"\rightarrow " + f"{args.y_label}$" ),
        Line2D([], [], marker='_', color=args.y_color, markersize=args.marker_size, linestyle='None', label=f"${args.y_label}" + r"\rightarrow " + f"{args.x_label}$"),
    ]
    labels = [
        "placeholder",
        "placeholder"
    ]
    # lines, labels = ax.get_legend_handles_labels()
    # lines2, labels2 = ax2.get_legend_handles_labels()
    # lines += lines2
    # labels += labels2
    lines += lines3
    labels += labels3

    if args.plot_legend == "True":
        plt.figlegend(handles=lines, loc="upper center", ncols=len(lines), columnspacing=0.5, edgecolor="white")
        # plt.figlegend(lines2, labels2, loc="upper right", ncols=len(lines2), columnspacing=0.2)
    # plt.figlegend(lines1+lines2, labels1+labels2, loc="upper center", ncols=len(lines1+lines2), columnspacing=0.2)
    # plt.figlegend(lines1+lines2, labels1+labels2, loc="center right", ncols=1)

    plt.subplots_adjust(hspace=0)
    if args.reverse:
        plt.gca().invert_xaxis()

    # fig = plt.gcf()
    fig.set_size_inches(args.fig_width / INCH_TO_CM, args.fig_height / INCH_TO_CM)
    fig.subplots_adjust(bottom=args.fig_margin_bottom, top=args.fig_margin_top, hspace=args.fig_margin_hspace)
    # fig.get_layout_engine().set(hspace=0)

    fig.savefig(f'{args.output_file}.svg', dpi=200)
    if args.plot_show == "True":
        plt.show()



def main(args) -> None:
    logger.info("Graphing T values.")
    
    logger.debug(f"Input arguments: {args}")

    plt.rcParams.update({
        "font.family": "DejaVu Sans"
    })

    data = get_data(args, args.csv)

    plot_one(args, data)


