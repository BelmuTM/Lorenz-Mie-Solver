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

import numpy as np


def calculate_phase_cdf_normalized(wavelengths, angles, phase):
    """
    :param wavelengths: list of wavelengths in nanometers
    :param angles: list of angles in degrees
    :param phase: list of phase values, indexed by wavelength and angle
    :return: normalized phase and its normalized forward CDF
    """

    wavelengths_n = len(wavelengths)
    angles_n      = len(angles)

    # Convert the angles from degrees to cos(angle_radians)
    angles = np.cos(np.deg2rad(angles))

    cdf = np.zeros((wavelengths_n, angles_n))

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

        for a in range(0, angles_n):

            phase[w][a] /= max_cdf
            cdf[w][a]   /= max_cdf

    return phase, cdf
