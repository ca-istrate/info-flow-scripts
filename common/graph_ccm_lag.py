from matplotlib.lines import Line2D
import matplotlib.pyplot as plt
import os
import pandas as pd


# E,tau,tp,nn,lib_column,target_column,lib_size,num_pred,rho,mae,rmse,direction
INCH_TO_CM = 2.54 

def default_output_file(args) -> None:
    inp = os.path.splitext(os.path.basename(args.csv))[0]
    vars(args)["output_file"] = f"{inp}-graph-ccm-conv.out"
    return vars(args)["output_file"]


def get_data(args, csv: str):
    data = pd.read_csv(csv) 

    data.loc[:, ['tp']] *= args.t_conv
    data['rho'] = data['rho'].clip(lower=0)

    return data


def plot_graph(args, data_x, data_y):
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

    ax.axvline(x=0, color="black")

    ax.plot(data_x.loc[:, "tp"], data_x.loc[:, "rho"], c=args.x_color, label=fr"${args.y_label}~\text"'{xmap}'fr"~{args.x_label}$")
    ax.axvline(x=data_x.loc[data_x['rho'].idxmax(), 'tp'], ls=":", color=args.x_color)

    if args.y_plot == "True":
        ax.plot(data_y.loc[:, "tp"], data_y.loc[:, "rho"], c=args.y_color, label=fr"${args.x_label}~\text"'{xmap}'fr"~{args.y_label}$")
        ax.axvline(x=data_y.loc[data_y['rho'].idxmax(), 'tp'], ls=":", color=args.y_color)


    ax.set_xlabel(f"{args.t_label} ({args.t_unit})")
    ax.set_ylabel("Cross Map Skill")
    ax.set_ylim(-0.05, 1.05)


    lines = ax.get_legend_handles_labels()[0]
    lines += [
        Line2D([], [], marker='_', color=args.bootstrap_color, markersize=args.marker_size, linestyle='None', label="Bootstrap"),
        Line2D([], [], marker='_', color=args.ebisuzaki_color, markersize=args.marker_size, linestyle='None', label="Ebisuzaki"),
    ]

    plt.figlegend(handles=lines, loc="upper center", ncols=len(lines), columnspacing=0.5, edgecolor="white")

    fig.set_size_inches(args.fig_width / INCH_TO_CM, args.fig_height / INCH_TO_CM)
    fig.subplots_adjust(bottom=args.fig_margin_bottom, top=args.fig_margin_top, hspace=args.fig_margin_hspace,
                        left=args.fig_margin_left, right=args.fig_margin_right)

    fig.savefig(f'{args.output_file}.svg', dpi=200)
    if args.plot_show == "True":
        plt.show()


def main(args):
    data = get_data(args, args.csv)

    data_x = data[data['lib_column'] == data['lib_column'].iloc[0]]
    data_y = data[data['lib_column'] != data['lib_column'].iloc[0]]

    plot_graph(args, data_x, data_y)
