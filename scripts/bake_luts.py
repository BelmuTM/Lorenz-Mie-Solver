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

from pathlib import Path

import numpy as np
import pandas as pd

import phase_cdf

OUTPUT_DIR = Path("output")

cross_sections_input = OUTPUT_DIR / "cross_sections.csv"
asymmetry_input      = OUTPUT_DIR / "asymmetry.csv"
phase_input          = OUTPUT_DIR / "phase.csv"

properties_output = OUTPUT_DIR / "lut_properties.dat"
phase_output      = OUTPUT_DIR / "lut_phase.dat"

TOLERANCE = 1e-6


def bake_lorenz_mie_luts():

    cross_sections = pd.read_csv(cross_sections_input)
    asymmetries    = pd.read_csv(asymmetry_input)
    phases         = pd.read_csv(phase_input)

    wavelengths = cross_sections.iloc[:, 0].to_numpy()
    extinction  = cross_sections.iloc[:, 1].to_numpy()
    scattering  = cross_sections.iloc[:, 2].to_numpy()
    asymmetry   = asymmetries.iloc[:, 1].to_numpy()

    wavelength_count = len(wavelengths)
    angle_count      = phases["angle"].nunique()

    wavelength_resolution = wavelength_count - 1
    angle_resolution      = angle_count - 1

    expected_wavelengths = np.linspace(wavelengths[0], wavelengths[-1], wavelength_count)
    expected_angles      = np.linspace(0.0, 180.0, angle_count)

    # Validation

    if len(asymmetries) != wavelength_count:
        raise ValueError("Asymmetry file and cross sections file have different wavelength counts.")

    if len(phases) != wavelength_count * angle_count:
        raise ValueError("Phase file row count is not wavelength count * angle count.")

    if not np.allclose(asymmetries.iloc[:, 0], expected_wavelengths, atol=TOLERANCE):
        raise ValueError("Asymmetry: wavelength mismatch.")

    if not np.allclose(phases["wavelength"], np.repeat(expected_wavelengths, angle_count), atol=TOLERANCE):
        raise ValueError("Phase: wavelength mismatch.")

    if not np.allclose(phases["angle"], np.tile(expected_angles, wavelength_count), atol=TOLERANCE):
        raise ValueError("Phase: angle mismatch.")

    # Properties: extinction, scattering, asymmetry, padding (16 bytes per wavelength)

    np.column_stack([extinction, scattering, asymmetry, np.zeros(wavelength_count)]) \
        .astype("<f4").tofile(properties_output)

    # Phase

    phase = phases["unpolarized"].to_numpy().reshape(wavelength_count, angle_count).tolist()

    # Normalize phase and compute phase CDF
    phase, cdf = phase_cdf.calculate_phase_cdf_normalized(
        wavelengths.tolist(), expected_angles.tolist(), phase
    )

    np.stack([phase, cdf], axis=-1).astype("<f4").tofile(phase_output)

    print(f"Wavelengths: {wavelengths[0]}-{wavelengths[-1]} nm, resolution {wavelength_resolution} ({wavelength_count} entries)")
    print(f"Angles: 0-180 deg, resolution {angle_resolution} ({angle_count} entries)")


if __name__ == "__main__":
    bake_lorenz_mie_luts()
