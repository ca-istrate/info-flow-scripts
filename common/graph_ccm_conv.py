import os
import matplotlib.pyplot as plt
import pandas as pd

#"lib_size","cross_map_skill","significance_ebs","significance_bts"
INCH_TO_CM = 2.54 


def default_output_file(args) -> None:
    inp = os.path.splitext(os.path.basename(args.csv))[0]
    vars(args)["output_file"] = f"{inp}-graph-ccm-conv.out"
    return vars(args)["output_file"]


def get_data(args, csv: str):
    data = pd.read_csv(csv) 

    return data


def plot_graph(args, data):
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

    fig, ax = plt.subplots()

    ax.fill_between(data.loc[:, "lib_size"], data.loc[:, "significance_ebs"], color=args.ebisuzaki_color, label="Ebisuzaki")
    ax.fill_between(data.loc[:, "lib_size"], data.loc[:, "significance_bts"], color=args.bootstrap_color, label="Bootstrap")
    ax.plot(data.loc[:, "lib_size"], data.loc[:, "cross_map_skill"], c=args.x_color)

    ax.set_xlabel("Library Size")
    ax.set_ylabel("Cross Map Skill")
    ax.set_ylim(-0.05, 1.05)

    fig.set_size_inches(args.fig_width / INCH_TO_CM, args.fig_height / INCH_TO_CM)
    fig.subplots_adjust(bottom=args.fig_margin_bottom, top=args.fig_margin_top, hspace=args.fig_margin_hspace,
                        left=args.fig_margin_left, right=args.fig_margin_right)

    fig.savefig(f'{args.output_file}.svg', dpi=200)
    if args.plot_show == "True":
        plt.show()


def main(args):
    data = get_data(args, args.csv)

    plot_graph(args, data)
