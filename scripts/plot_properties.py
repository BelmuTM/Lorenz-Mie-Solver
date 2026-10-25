#################################################################################
#                                                                               #
#    Lorenz-Mie Solver                                                          #
#    Copyright (C) 2026  Belmu                                                  #
#                                                                               #
#    This program is free software: you can redistribute it and/or modify       #
#    it under the terms of the GNU General Public License as published by       #
#    the Free Software Foundation, either version 3 of the License, or          #
#    (at your option) any later version.                                        #
#                                                                               #
#    This program is distributed in the hope that it will be useful,            #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of             #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              #
#    GNU General Public License for more details.                               #
#                                                                               #
#    You should have received a copy of the GNU General Public License          #
#    along with this program.  If not, see <https://www.gnu.org/licenses/>.     #
#                                                                               #
#################################################################################

import argparse

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

OUTPUT_DIR = Path("output")
SAVE_PNG   = True


def parse_args():
    parser = argparse.ArgumentParser(description="Plot of the Lorenz-Mie scattering properties (extinction, scattering, asymmetry).")

    parser.add_argument("--no-save", action="store_true", help="Do not write the PNG files")

    return parser.parse_args()


def plot_curve(title: str, x: float, y: float, x_label: str, y_label: str, file_name: str, color: str, no_save_png: bool):

    fig, ax = plt.subplots(figsize=(8, 4.5))

    ax.plot(x, y, color=color, linewidth=1.5)
    ax.set_title(title)
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.set_xlim(x.min(), x.max())
    ax.grid(True, alpha=0.3)

    fig.tight_layout()

    if not no_save_png:
        fig.savefig(OUTPUT_DIR / file_name, dpi=200)


if __name__ == "__main__":

    args = parse_args()

    cross = pd.read_csv(OUTPUT_DIR / "cross_sections.csv")
    asym  = pd.read_csv(OUTPUT_DIR / "asymmetry.csv")

    plot_curve(
        "Extinction cross section",
        cross["wavelength"], cross["extinction"],
        "Wavelength [nm]", r"Extinction [m$^{-1}$]",
        "plot_extinction.png",
        "tab:blue",
        args.no_save
    )

    plot_curve(
        "Scattering cross section",
        cross["wavelength"], cross["scattering"],
        "Wavelength [nm]", r"Scattering [m$^{-1}$]",
        "plot_scattering.png",
        "tab:orange",
        args.no_save
    )

    plot_curve(
        "Asymmetry parameter",
        asym["wavelength"], asym["asymmetry"],
        "Wavelength [nm]", r"Asymmetry $g$",
        "plot_asymmetry.png",
        "tab:green",
        args.no_save
    )

    plt.show()
