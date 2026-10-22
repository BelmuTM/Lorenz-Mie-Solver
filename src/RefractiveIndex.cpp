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

/*
 *  [References]:
 *      Polyanskiy, M. (2017). Ciddor 1996 - air.
 *      https://github.com/polyanskiy/refractiveindex.info-scripts/blob/master/scripts/Ciddor%201996%20-%20air.py
 */

#include "RefractiveIndex.h"

#include <cmath>

namespace RefractiveIndex {

    constexpr double ABSOLUTE_ZERO = 273.15; // Kelvins

    double water(double lambda) {

        lambda *= 1e-3; // To micrometers

        const double lambdaSq = lambda * lambda;

        double dispersion  = 1.0;
               dispersion += 5.689093832e-1 * lambdaSq / (lambdaSq - 5.110301794e-3);
               dispersion += 1.719708856e-1 * lambdaSq / (lambdaSq - 1.825180155e-2);
               dispersion += 2.062501582e-2 * lambdaSq / (lambdaSq - 2.624158904e-2);
               dispersion += 1.123965424e-1 * lambdaSq / (lambdaSq - 1.067505178e1 );

        return dispersion;
    }

    static double compressibility(double T, const double p, const double xw) {

        T -= ABSOLUTE_ZERO;

        constexpr double a0 =  1.58123e-6;
        constexpr double a1 = -2.9331e-8;
        constexpr double a2 =  1.1043e-10;
        constexpr double b0 =  5.707e-6;
        constexpr double b1 = -2.051e-8;
        constexpr double c0 =  1.9898e-4;
        constexpr double c1 = -2.376e-6;
        constexpr double d  =  1.83e-11;
        constexpr double e  = -0.765e-8;

        return 1.0 - p / T
            * (a0 + a1 * T + a2 * T * T + (b0 + b1 * T) * xw + (c0 + c1 * T) * xw * xw)
            + std::pow(p / T, 2.0) * (d + e * xw * xw);
    }

    // T: °C
    // p: Pa
    // h: %
    // xc: ppm
    double air(double lambda, double T, const double p, const double h, const double xc) {

        lambda *= 1e-3; // To micrometers

        T += ABSOLUTE_ZERO;

        constexpr double R     =  8.314510;
        constexpr double k0    =  238.0185;
        constexpr double k1    =  5792105.0;
        constexpr double k2    =  57.362;
        constexpr double k3    =  167917.0;
        constexpr double w0    =  295.235;
        constexpr double w1    =  2.6422;
        constexpr double w2    = -0.032380;
        constexpr double w3    =  0.004028;
        constexpr double alpha =  1.00062;
        constexpr double beta  =  3.14e-8;
        constexpr double gamma =  5.6e-7;

        const double sigma  = 1.0 / lambda;
        const double sigma2 = sigma * sigma;

        double svp = 0.0;

        // Saturation vapor pressure of water vapor in air at temperature T
        if (T >= 0.0) {

            constexpr double A =  1.2378847e-5;
            constexpr double B = -1.9121316e-2;
            constexpr double C =  33.93711047;
            constexpr double D =  6.3431645e3;

            svp = std::exp(A * std::pow(T, 2.0) + B * T + C + D / T);

        } else {

            svp = std::pow(10.0, -2663.5 / T + 12.537);

        }

        const double f = alpha + beta * p + gamma * std::pow(T, 2.0); // Enhancement factor of water vapor in air

        const double xw = f * h * svp / p; // Molar fraction of water vapor in moist air

        // Refractive index of standard air at 15 °C, 101325 Pa, 0% humidity, 450 ppm CO2
        const double nas = 1.0 + (k1 / (k0 - sigma2) + k3 / (k2 - sigma2)) * 1e-8;
        // Refractive index of standard air at 15 °C, 101325 Pa, 0% humidity, xc ppm CO2
        const double naxs = 1.0 + (nas - 1.0) * (1.0 + 0.534e-6 * (xc - 450.0));
        // Refractive index of water vapor at standard conditions (20 °C, 1333 Pa)
        const double nws = 1.0 + 1.022 * (w0 + w1 * sigma2 + w2 * std::pow(sigma, 4.0) + w3 * std::pow(sigma, 6.0)) * 1e-8;

        const     double Ma = 1e-3 * (28.9635 + 12.011e-6 * (xc - 400.0)); // Molar mass of dry air
        constexpr double Mw = 0.018015;                                    // Molar mass of water vapor

        const double Za = compressibility(288.15, 101325.0, 0.0); // Compressibility of dry air
        const double Zw = compressibility(293.15, 1333.0  , 1.0); // Compressibility of pure water vapor

        const double paxs = 101325.0 * Ma / (Za * R * 288.15); // Density of standard air
        const double pws  = 1333.0   * Mw / (Zw * R * 293.15); // Density of standard water vapor

        const double pa = p * Ma / (compressibility(T, p, xw) * R * T) * (1.0 - xw); // Density of the dry component of the moist air
        const double pw = p * Mw / (compressibility(T, p, xw) * R * T) * xw;         // Density of the water vapor component

        return 1.0 + pa / paxs * (naxs - 1.0) + pw / pws * (nws - 1.0);
    }
}
