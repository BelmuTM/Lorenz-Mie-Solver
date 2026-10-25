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
 *      Bohren, C. F., & Huffman, D. R. (1983). Absorption and Scattering of Light by Small Particles.
 *      https://onlinelibrary.wiley.com/doi/book/10.1002/9783527618156
 *      Frisvad et al. (2007). Computing the Scattering Properties of Participating Media Using Lorenz-Mie Theory.
 *      https://cseweb.ucsd.edu//~henrik/papers/lorenz_mie_theory/
 *      Wikipedia. (2026). Legendre polynomials. https://en.wikipedia.org/wiki/Legendre_polynomials
 *
 *  [Notes]:
 *      In the Lorenz-Mie solution paper:
 *
 *      r is the particle size
 *
 *      x is the size parameter of the host medium.
 *      y is the size parameter of the particle medium.
 *      Both are calculated by multiplying the medium's wavenumber by its particle radius.
 *
 *      A_n is the logarithmic derivative of the Riccati-Bessel function.
 *      B_n is the logarithmic derivative of the Riccati-Hankel function.
 *      Both originate from the spherical Bessel functions.
 *
 *      S1_n is the "senkrecht" (perpendicular) scattering amplitude.
 *      S2_n is the parallel scattering amplitude.
 *
 *      a_n are the electric multipole coefficients (TM modes).
 *      b_n are the magnetic multipole coefficients (TE modes).
 *
 *      The paper's equations all rely on relations of recurrence, hence the need for a running index n.
 */

#include "LorenzMie.h"

using namespace std::complex_literals;

namespace {

    constexpr double PI = std::numbers::pi;

    constexpr double DEGREES_TO_RADIANS = PI / 180.0;

    double complexAbsSquared(const ComplexDouble x) {
        return std::abs(x) * std::abs(x);
    }

        void calculateRiccatiBesselLogDerivatives(
        ArrayComplexDouble& besselDerivatives,
        const ComplexDouble sizeParameter,
        const int           startOrder
    ) {

        const int maxOrder = static_cast<int>(besselDerivatives.size()) - 1;

        ComplexDouble bessel(0.0, 0.0);

        // Downward recurrence
        for (int n = startOrder; n >= 0; n--) {

            const ComplexDouble np1_over_z = ComplexDouble(n + 1, 0.0) / sizeParameter;

            bessel = np1_over_z - ComplexDouble(1.0, 0.0) / (np1_over_z + bessel);

            if (n <= maxOrder)
                besselDerivatives[n] = bessel;
        }
    }

    void calculateLorenzCoefficients(
        LorenzMie::Coefficients&  coefficients,
        const ComplexDouble       hostIOR,
        const ComplexDouble       hostSizeParameter,
        const ComplexDouble       particleIOR,
        const ArrayComplexDouble& hostBesselValues,
        const ArrayComplexDouble& particleBesselValues
    ) {

        const int maxOrder = static_cast<int>(hostBesselValues.size()) - 1;

        // psi_n(x) * zeta_n(x)
        ComplexDouble psiZetaProduct          = 0.5 * (1.0 - std::exp( 2i * hostSizeParameter));
        // psi_n(x) / zeta_n(x)
        ComplexDouble incidentToOutgoingRatio = 0.5 * (1.0 - std::exp(-2i * hostSizeParameter));

        // B_n(x)
        auto hostHankel = ComplexDouble(1i);

        // Modelling the scattering field described by the host medium

        for (int n = 1; n <= maxOrder; n++) {

            const ComplexDouble n_over_z = ComplexDouble(n, 0.0) / hostSizeParameter;

            psiZetaProduct *= (n_over_z - hostBesselValues[n - 1]) * (n_over_z - hostHankel);

            const ComplexDouble hostBessel     = hostBesselValues[n];
            const ComplexDouble particleBessel = particleBesselValues[n];

            hostHankel = hostBessel + 1i / psiZetaProduct;

            incidentToOutgoingRatio *= (hostHankel + n_over_z) / (hostBessel + n_over_z);

            // Boundary impedance matching (c.f. Maxwell's boundary conditions)

            // Electric (TM modes)
            const ComplexDouble electricNumerator   = hostIOR * particleBessel - particleIOR * hostBessel;
            const ComplexDouble electricDenominator = hostIOR * particleBessel - particleIOR * hostHankel;

            const ComplexDouble electricBoundaryMatching = electricNumerator / electricDenominator;

            // Magnetic (TE modes)
            const ComplexDouble magneticNumerator   = particleIOR * particleBessel - hostIOR * hostBessel;
            const ComplexDouble magneticDenominator = particleIOR * particleBessel - hostIOR * hostHankel;

            const ComplexDouble magneticBoundaryMatching = magneticNumerator / magneticDenominator;

            // Radial scattering amplitudes for electric and magnetic fields
            coefficients.electric[n] = incidentToOutgoingRatio * electricBoundaryMatching;
            coefficients.magnetic[n] = incidentToOutgoingRatio * magneticBoundaryMatching;
        }
    }
}

namespace LorenzMie {

    Coefficients calculateMieCoefficients(const Context& context) {

        const ComplexDouble hostWaveNumber     = 2.0 * PI * context.hostIOR     / context.wavelength;
        const ComplexDouble particleWaveNumber = 2.0 * PI * context.particleIOR / context.wavelength;

        const ComplexDouble hostSizeParameter     = hostWaveNumber     * context.particleRadius;
        const ComplexDouble particleSizeParameter = particleWaveNumber * context.particleRadius;

        constexpr double p = 4.3;
        // ceil(|x| + p|x|^1/3 + 1)
        const int maxOrders = static_cast<int>(
            std::ceil(std::abs(hostSizeParameter) + p * std::cbrt(std::abs(hostSizeParameter)) + 1)
        );

        constexpr int startOffset = 60;

        const int hostStartOrder =
            std::max(maxOrders, static_cast<int>(std::ceil(std::abs(hostSizeParameter)))) + startOffset;

        const int particleStartOrder =
            std::max(maxOrders, static_cast<int>(std::ceil(std::abs(particleSizeParameter)))) + startOffset;

        // Bessel derivatives

        ArrayComplexDouble hostBesselValues    (maxOrders + 2);
        ArrayComplexDouble particleBesselValues(maxOrders + 2);

        calculateRiccatiBesselLogDerivatives(hostBesselValues,     hostSizeParameter,     hostStartOrder);
        calculateRiccatiBesselLogDerivatives(particleBesselValues, particleSizeParameter, particleStartOrder);

        // Lorenz electromagnetic multipole coefficients

        Coefficients coefficients{
            ArrayComplexDouble(maxOrders + 2),
            ArrayComplexDouble(maxOrders + 2)
        };

        calculateLorenzCoefficients(
            coefficients,
            context.hostIOR,
            hostSizeParameter,
            context.particleIOR,
            hostBesselValues,
            particleBesselValues
        );

        return coefficients;
    }

    ScatteringProperties calculateScatteringProperties(const Context& context, const Coefficients& coefficients) {

        const ComplexDouble hostIORSq = context.hostIOR * context.hostIOR;

        double extinctionSum = 0.0;
        double scatteringSum = 0.0;

        double asymmetrySum = 0.0;

        for (int n = 1; n <= coefficients.maxOrders(); n++) {

            const ComplexDouble electricCoeff_n = coefficients.electric[n];
            const ComplexDouble magneticCoeff_n = coefficients.magnetic[n];

            const ComplexDouble electricCoeff_n1 = coefficients.electric[n + 1];
            const ComplexDouble magneticCoeff_n1 = coefficients.magnetic[n + 1];

            // Extinction, scattering cross-sections

            extinctionSum += (2.0 * n + 1.0) * std::real((electricCoeff_n + magneticCoeff_n) / hostIORSq);
            scatteringSum += (2.0 * n + 1.0) * (complexAbsSquared(electricCoeff_n) + complexAbsSquared(magneticCoeff_n));

            // Asymmetry

            // Adjacent order, same type interference
            const double interferenceAdjacentOrder =
                n * (n + 2.0) / (n + 1.0) *
                std::real(electricCoeff_n * std::conj(electricCoeff_n1) + magneticCoeff_n * std::conj(magneticCoeff_n1));

            // Same order, electric-magnetic interference
            const double interferenceElectricMagnetic =
                (2.0 * n + 1.0) / (n * (n + 1.0)) * std::real(electricCoeff_n * std::conj(magneticCoeff_n));

            asymmetrySum += interferenceAdjacentOrder + interferenceElectricMagnetic;
        }

        // Host medium radial optical depth
        const double alpha = 4.0 * PI * context.particleRadius * std::imag(context.hostIOR) / context.wavelength;
        // Corrects for incident-field attenuation across the particle surface
        const double gamma = alpha < 1e-3
            ? 1.0 + 2.0 / 3.0 * alpha + 0.25 * alpha * alpha
            : 2.0 * (1.0 + (alpha - 1.0) * std::exp(alpha)) / (alpha * alpha);

        const double wavelengthSq = context.wavelength * context.wavelength;
        const double hostIORAbsSq = complexAbsSquared(context.hostIOR);

        ScatteringProperties properties{};

        properties.extinction = wavelengthSq / (2.0 * PI) * extinctionSum;
        properties.scattering = wavelengthSq * std::exp(-alpha) / (2.0 * PI * gamma * hostIORAbsSq) * scatteringSum;

        properties.asymmetry = asymmetrySum / (0.5 * scatteringSum);

        return properties;
    }

    PhaseFunction calculatePhaseFunction(const Coefficients& coefficients, const double scatteringAngle) {

        const double cosTheta = std::cos(scatteringAngle * DEGREES_TO_RADIANS);

        double angularPi_m1 = 0.0; // pi_{n-1}
        double angularPi_m2 = 0.0; // pi_{n-2}

        ComplexDouble scatteringAmp_S(0.0, 0.0);
        ComplexDouble scatteringAmp_P(0.0, 0.0);

        double scatteringSum = 0.0;

        for (int n = 1; n <= coefficients.maxOrders(); n++) {

            // Angular functions

            const double angularPi_n = n == 1
                ? 1.0
                : ((2.0 * n - 1.0) * cosTheta * angularPi_m1 - n * angularPi_m2) / (n - 1.0);

            const double angularTau_n = n * cosTheta * angularPi_n - (n + 1.0) * angularPi_m1;

            // Scattering amplitudes (senkrecht, parallel)

            const double amplitudeWeight = (2.0 * n + 1.0) / (n * (n + 1.0));

            const ComplexDouble electricCoeff_n = coefficients.electric[n];
            const ComplexDouble magneticCoeff_n = coefficients.magnetic[n];

            scatteringAmp_S += amplitudeWeight * (electricCoeff_n * angularPi_n + magneticCoeff_n * angularTau_n);
            scatteringAmp_P += amplitudeWeight * (magneticCoeff_n * angularPi_n + electricCoeff_n * angularTau_n);

            scatteringSum += (2.0 * n + 1.0) * (complexAbsSquared(electricCoeff_n) + complexAbsSquared(magneticCoeff_n));

            angularPi_m2 = angularPi_m1;
            angularPi_m1 = angularPi_n;
        }

        const double phaseNormalization = 2.0 * PI * scatteringSum;

        PhaseFunction phase{};

        phase.senkrecht   = complexAbsSquared(scatteringAmp_S) / phaseNormalization;
        phase.parallel    = complexAbsSquared(scatteringAmp_P) / phaseNormalization;
        phase.unpolarized = (phase.senkrecht + phase.parallel) * 0.5;

        return phase;
    }

    void calculateBulkProperties(
        const BulkContext&          context,
        const ArrayDouble&          scatteringAngles,
        ScatteringProperties&       properties,
        std::vector<PhaseFunction>& phases
    ) {

        phases.resize(scatteringAngles.size());

        const double particleRadiusStep =
            (context.particleRadiusMax - context.particleRadiusMin) / context.particleRadiusBins;

        for (int r = 0; r < context.particleRadiusBins; r++) {

            const double numberDensityStep = context.numberDensities[r] * particleRadiusStep;

            const double particleRadius = context.particleRadiusMin + (r + 0.5) * particleRadiusStep;

            const Context mieContext(context.wavelength, particleRadius, context.hostIOR, context.particleIOR);

            const Coefficients mieCoefficients = calculateMieCoefficients(mieContext);

            const auto [extinction, scattering, asymmetry] = calculateScatteringProperties(mieContext, mieCoefficients);

            properties.extinction += extinction * numberDensityStep;
            properties.scattering += scattering * numberDensityStep;

            properties.asymmetry += asymmetry  * scattering * numberDensityStep;

            for (int a = 0; a < scatteringAngles.size(); a++) {

                const auto [senkrecht, parallel, unpolarized] = calculatePhaseFunction(mieCoefficients, scatteringAngles[a]);

                phases[a].senkrecht   += senkrecht   * scattering * numberDensityStep;
                phases[a].parallel    += parallel    * scattering * numberDensityStep;
                phases[a].unpolarized += unpolarized * scattering * numberDensityStep;

            }
        }

        properties.asymmetry /= properties.scattering;

        for (int a = 0; a < scatteringAngles.size(); a++) {

            phases[a].senkrecht   /= properties.scattering;
            phases[a].parallel    /= properties.scattering;
            phases[a].unpolarized /= properties.scattering;
        }
    }

    double particleSizeDistributionLogNormal(
        const double particleRadius,
        const double meanRadius,
        const double radiusStandardDeviation
    ) {
        const double ratioSq = radiusStandardDeviation * radiusStandardDeviation / (meanRadius * meanRadius);

        const double alpha = std::log(meanRadius) - 0.5 * std::log(ratioSq + 1.0);
        const double beta  = std::sqrt(std::log(ratioSq + 1.0));

        // 1 / (r^4) here makes the distribution volume-weighted (proportional to the particle's volume)
        return 1.0 / (std::pow(particleRadius, 4.0) * beta * std::sqrt(2.0 * PI))
            * std::exp(-0.5 * std::pow((std::log(particleRadius) - alpha) / beta, 2.0));
    }

    void calculateBulkNumberDensities(
        BulkContext& context,
        const double targetVolumeFraction,
        const double meanRadius,
        const double radiusStandardDeviation
    ) {

        double volumeMoment = 0.0;

        context.numberDensities.resize(context.particleRadiusBins);

        const double particleRadiusStep =
            (context.particleRadiusMax - context.particleRadiusMin) / context.particleRadiusBins;

        for (int r = 0; r < context.particleRadiusBins; r++) {

            const double particleRadius = context.particleRadiusMin + (r + 0.5) * particleRadiusStep;

            context.numberDensities[r] =
                particleSizeDistributionLogNormal(particleRadius, meanRadius, radiusStandardDeviation);

            volumeMoment += std::pow(particleRadius, 3.0) * context.numberDensities[r] * particleRadiusStep;
        }

        const double numberDensityScale = targetVolumeFraction / (4.0 * PI / 3.0 * volumeMoment);

        for (int r = 0; r < context.particleRadiusBins; r++) {
            context.numberDensities[r] *= numberDensityScale;
        }
    }

}
