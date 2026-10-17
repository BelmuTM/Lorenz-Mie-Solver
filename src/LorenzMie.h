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
#include <numbers>
#include <valarray>

using ComplexDouble      = std::complex<double>;
using ArrayComplexDouble = std::valarray<ComplexDouble>;

namespace LorenzMie {

    struct Context {
        double        particleRadius;
        double        wavelength;
        ComplexDouble hostIOR;
        ComplexDouble particleIOR;
    };

    struct Coefficients {
        ArrayComplexDouble electric;
        ArrayComplexDouble magnetic;

        int maxOrders() const {
            return static_cast<int>(electric.size()) - 2;
        }
    };

    struct ScatteringProperties {
        double extinction;
        double scattering;

        double asymmetry;
    };

    struct PhaseFunction {
        double senkrecht;
        double parallel;
        double unpolarized;
    };

    Coefficients calculateMieCoefficients(const Context& context);

    ScatteringProperties calculateScatteringProperties(const Context& context, const Coefficients& coefficients);

    PhaseFunction calculatePhaseFunction(const Coefficients& coefficients, double scatteringAngle);
}
