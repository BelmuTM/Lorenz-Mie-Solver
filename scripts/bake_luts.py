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

import csv
import struct

OUTPUT_DIR = Path("output")

cross_sections_input = OUTPUT_DIR / "cross_sections.csv"
asymmetry_input      = OUTPUT_DIR / "asymmetry.csv"
phase_input          = OUTPUT_DIR / "phase.csv"

properties_output = OUTPUT_DIR / "lut_properties.dat"
phase_output      = OUTPUT_DIR / "lut_phase.dat"

TOLERANCE = 1e-6


def read_rows(path: str) -> list[list[str]]:

    with open(path, "r", newline="") as file:
        reader = csv.reader(file)
        next(reader)  # Skip header

        return [row for row in reader if row]

    return None


def check_close(actual: float, expected: float, message: str) -> None:

    if abs(actual - expected) > TOLERANCE:
        raise ValueError(message)


def write_properties_file(cross_sections: list[list[str]], asymmetries: list[list[str]]) -> None:

    with open(properties_output, "wb") as output_file:

        for cross_section_row, asymmetry_row in zip(cross_sections, asymmetries):

            extinction = float(cross_section_row[1])
            scattering = float(cross_section_row[2])
            asymmetry  = float(asymmetry_row[1])

            # extinction, scattering, asymmetry, padding (16 bytes per wavelength)
            output_file.write(struct.pack("<ffff", extinction, scattering, asymmetry, 0.0))


def write_phase_file(phases: list[list[str]]) -> None:

    with open(phase_output, "wb") as output_file:

        for row in phases:
            unpolarized = float(row[2])

            output_file.write(struct.pack("<f", unpolarized))


def bake_lorenz_mie_luts():

    cross_sections = read_rows(cross_sections_input)
    asymmetries    = read_rows(asymmetry_input)
    phases         = read_rows(phase_input)

    # Wavelength grid
    wavelengths      = [float(row[0]) for row in cross_sections]
    wavelength_count = len(wavelengths)

    wavelength_lower_bound = wavelengths[0]
    wavelength_upper_bound = wavelengths[-1]
    wavelength_resolution  = wavelength_count - 1

    # Angle grid, number of consecutive rows sharing the first wavelength
    first_wavelength = phases[0][0]
    angle_count      = 0

    for row in phases:
        if row[0] != first_wavelength:
            break
        angle_count += 1

    angle_resolution = angle_count - 1

    # Validation
    if len(asymmetries) != wavelength_count:
        raise ValueError("Asymmetry file and cross sections file have different wavelength counts.")

    if len(phases) != wavelength_count * angle_count:
        raise ValueError("Phase file row count is not wavelength count * angle count.")

    for w in range(wavelength_count):

        expected_wavelength = wavelength_lower_bound + (wavelength_upper_bound - wavelength_lower_bound) * w / wavelength_resolution

        check_close(wavelengths[w],           expected_wavelength, f"Cross sections: non-uniform wavelength at index {w}.")
        check_close(float(asymmetries[w][0]), expected_wavelength, f"Asymmetry: wavelength mismatch at index {w}.")

        for a in range(angle_count):

            row = phases[w * angle_count + a]

            expected_angle = 180.0 * a / angle_resolution

            check_close(float(row[0]), expected_wavelength, f"Phase: wavelength mismatch at wavelength {w}, angle {a}.")
            check_close(float(row[1]), expected_angle,      f"Phase: angle mismatch at wavelength {w}, angle {a}.")

    write_properties_file(cross_sections, asymmetries)
    write_phase_file(phases)

    print(f"Wavelengths: {wavelength_lower_bound}-{wavelength_upper_bound} nm, resolution {wavelength_resolution} ({wavelength_count} entries)")
    print(f"Angles: 0-180 deg, resolution {angle_resolution} ({angle_count} entries)")


if __name__ == "__main__":
    bake_lorenz_mie_luts()
