/********************************************************************************/
/*                                                                              */
/*    Lorenz-Mie Solver                                                         */
/*    Copyright (C) 2026  Belmu                                                 */
/*                                                                              */
/*    This program is free software: you can redistribute it and/or modify      */
/*    it under the terms of the GNU General Public License as published by      */
/*    the Free Software Foundation, either version 3 of the License, or         */
/*    (at your option) any later version.                                       */
/*                                                                              */
/*    This program is distributed in the hope that it will be useful,           */
/*    but WITHOUT ANY WARRANTY; without even the implied warranty of            */
/*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             */
/*    GNU General Public License for more details.                              */
/*                                                                              */
/*    You should have received a copy of the GNU General Public License         */
/*    along with this program.  If not, see <https://www.gnu.org/licenses/>.    */
/*                                                                              */
/********************************************************************************/

#pragma once

#include <algorithm>
#include <cmath>
#include <complex>
#include <valarray>
#include <vector>

using ComplexDouble      = std::complex<double>;
using ArrayComplexDouble = std::valarray<ComplexDouble>;
using ArrayDouble        = std::valarray<double>;

namespace LorenzMie {

    struct Context {
        double        wavelength;     // Nanometers [nm]
        double        particleRadius; // Nanometers [nm]
        ComplexDouble hostIOR;
        ComplexDouble particleIOR;
    };

    struct BulkContext {
        double        wavelength;        // Nanometers [nm]
        double        particleRadiusMin; // Nanometers [nm]
        double        particleRadiusMax; // Nanometers [nm]
        int           particleRadiusBins;
        ArrayDouble   numberDensities;
        ComplexDouble hostIOR;
        ComplexDouble particleIOR;
    };

    struct Coefficients {
        ArrayComplexDouble electric;
        ArrayComplexDouble magnetic;

        [[nodiscard]] int maxOrders() const {
            return static_cast<int>(electric.size()) - 2;
        }
    };

    struct ScatteringProperties {
        double extinction = 0.0;
        double scattering = 0.0;

        double asymmetry = 0.0;
    };

    struct PhaseFunction {
        double senkrecht   = 0.0;
        double parallel    = 0.0;
        double unpolarized = 0.0;
    };

    // Individual particle properties

    Coefficients calculateMieCoefficients(const Context& context);

    ScatteringProperties calculateScatteringProperties(const Context& context, const Coefficients& coefficients);

    // Scattering angle is in degrees, converted to radians internally.
    PhaseFunction calculatePhaseFunction(const Coefficients& coefficients, double scatteringAngle);

    // Bulk medium properties

    void calculateBulkProperties(
        const BulkContext&          context,
        const ArrayDouble&          scatteringAngles,
        ScatteringProperties&       properties,
        std::vector<PhaseFunction>& phases
    );

    double particleSizeDistributionLogNormal(double particleRadius, double meanRadius, double radiusStandardDeviation);

    void calculateBulkNumberDensities(
        BulkContext& context,
        double       targetVolumeFraction,
        double       meanRadius,
        double       radiusStandardDeviation
    );

}
