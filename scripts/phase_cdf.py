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

# Script to calculate the forward CDF of the phase and the normalized phase from it

import numpy as np


def calculate_phase_cdf_normalized(phase_data, phase_type="unpolarized"):

    wavelengths = np.sort(phase_data["wavelength"].unique())  # Nanometers
    angles      = np.sort(phase_data["angle"].unique())       # Degrees

    phase = (
        phase_data
        .pivot(index="wavelength", columns="angle", values=phase_type)
        .loc[wavelengths, angles]
        .to_numpy(dtype=np.float64)
    )

    wavelengths_n = len(wavelengths)
    angles_n      = len(angles)

    # Convert the angles from degrees to cos(angle_radians)
    for a in range(0, angles_n):
        angles[a] = np.cos(np.deg2rad(angles[a]))

    cdf = [[0.0 for _ in range(0, angles_n)] for _ in range(0, wavelengths_n)]

    # Integrate the forward CDF using trapezoids
    for w in range(0, wavelengths_n):

        cdf[w][0] = 0.0

        for a in range(1, angles_n):

            cdf[w][a] = (
                 cdf[w][a - 1]
               + 2.0 * np.pi * (phase[w][a - 1] + phase[w][a]) * 0.5
               * (angles[a - 1] - angles[a])
            )

    # Normalize both the phase and the CDF
    for w in range(0, wavelengths_n):

        max_cdf = cdf[w][-1]

        print(max_cdf)

        for a in range(0, angles_n):

            phase[w][a] /= max_cdf
            cdf[w][a]   /= max_cdf

    return phase, cdf
