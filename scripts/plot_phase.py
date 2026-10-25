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
import numpy as np
import pandas as pd

OUTPUT_DIR = Path("output")
PHASE_TYPES = ("unpolarized", "senkrecht", "parallel")


def parse_args():
    parser = argparse.ArgumentParser(description="Polar plot of the Lorenz-Mie phase function.")

    parser.add_argument("wavelengths", nargs="*", type=float, default=[440.0], help="Wavelengths in nanometers (nm)")
    parser.add_argument("--type", choices=PHASE_TYPES, default="unpolarized")
    parser.add_argument("--no-save", action="store_true", help="Do not write the PNG file")

    return parser.parse_args()


def closest_wavelength(available, target):
    return available[np.argmin(np.abs(available - target))]


def mirror_to_full_circle(angles_deg, values):
    """Mirroring the phase curve as it's symmetrical about the axis."""
    theta = np.deg2rad(np.concatenate([angles_deg, 360.0 - angles_deg[-2::-1]]))
    r     = np.concatenate([values, values[-2::-1]])
    return theta, r


def plot_phase(phase_data, wavelengths, phase_type):

    available = np.sort(phase_data["wavelength"].unique())

    fig, ax = plt.subplots(subplot_kw={"projection": "polar"}, figsize=(8, 8))

    for target in wavelengths:

        wavelength = closest_wavelength(available, target)
        subtable   = phase_data[phase_data["wavelength"] == wavelength].sort_values("angle")

        theta, r = mirror_to_full_circle(subtable["angle"].to_numpy(), subtable[phase_type].to_numpy())

        ax.plot(theta, r, linewidth=1, label=f"{wavelength:.0f} nm")

    ax.set_rscale("log")
    ax.set_theta_zero_location("E")
    ax.set_theta_direction(1)
    ax.set_title(f"Phase function ({phase_type})")
    ax.legend(loc="upper right", bbox_to_anchor=(1.25, 1.1))

    fig.tight_layout()

    return fig


if __name__ == "__main__":

    args = parse_args()

    phase = pd.read_csv(OUTPUT_DIR / "phase.csv")

    figure = plot_phase(phase, args.wavelengths, args.type)

    if not args.no_save:
        figure.savefig(OUTPUT_DIR / "plot_phase_polar.png", dpi=200)

    plt.show()
