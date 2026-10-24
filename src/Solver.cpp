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

#include "LorenzMie.h"

#include "RefractiveIndex.h"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>

#define TO_STR std::to_string

constexpr double PARTICLE_RADIUS = 50'000; // Nanometers [nm]

constexpr int WAVELENGTH_LOWER_BOUND = 390; // Nanometers [nm]
constexpr int WAVELENGTH_UPPER_BOUND = 830; // Nanometers [nm]

constexpr int WAVELENGTH_RESOLUTION = 880; // Amount of wavelengths between the lower and upper bounds
constexpr int ANGLE_RESOLUTION      = 720; // Amount of angles between 0° and 180°

const std::string outputFolder = OUTPUT_DIR;

const std::string crossSectionsOutput = outputFolder + "/cross_sections.csv";
const std::string asymmetryOutput     = outputFolder + "/asymmetry.csv";
const std::string phaseOutput         = outputFolder + "/phase.csv";

int main() {

    const auto startTime = std::chrono::high_resolution_clock::now();

    const std::string line = "------------------------------------------";

    std::cout
        << line
        << "\n"
        << "Lorenz-Mie Solver"
        << "\n\n"
        << "Solving for parameters:"
        << "\n"
        << "    - Particle radius: " << TO_STR(PARTICLE_RADIUS) << "nm"
        << "\n"
        << "    - Wavelength range: " << TO_STR(WAVELENGTH_LOWER_BOUND) << "-" << TO_STR(WAVELENGTH_UPPER_BOUND) << "nm"
        << "\n"
        << "    - Wavelengths: " << TO_STR(WAVELENGTH_RESOLUTION)
        << "\n"
        << "    - Angles: " << TO_STR(ANGLE_RESOLUTION)
        << "\n"
        << line
        << "\n";

    std::filesystem::create_directories(outputFolder);

    std::ofstream crossSectionsFile(crossSectionsOutput);
    std::ofstream asymmetryFile(asymmetryOutput);
    std::ofstream phaseFile(phaseOutput);

    if (!crossSectionsFile || !asymmetryFile || !phaseFile) {
        std::cerr << "Failed to open output files." << "\n";
        return 1;
    }

    crossSectionsFile << "wavelength,extinction,scattering" << "\n";

    asymmetryFile << "wavelength,asymmetry" << "\n";

    phaseFile << "wavelength,angle,unpolarized,senkrecht,parallel" << "\n";

    for (int w = 0; w <= WAVELENGTH_RESOLUTION; w++) {

        const double wavelength = std::lerp(
            WAVELENGTH_LOWER_BOUND,
            WAVELENGTH_UPPER_BOUND,
            static_cast<double>(w) / WAVELENGTH_RESOLUTION
        );

        const ComplexDouble airIOR   = {RefractiveIndex::air(wavelength, 15.0, 101350.0, 0.0, 420.0), 0.0};
        const ComplexDouble waterIOR = {RefractiveIndex::water(wavelength), 0.0};

        const LorenzMie::Context mieContext(wavelength, PARTICLE_RADIUS, airIOR, waterIOR);

        const LorenzMie::Coefficients mieCoefficients = LorenzMie::calculateMieCoefficients(mieContext);

        const auto [extinction, scattering, asymmetry] =
            LorenzMie::calculateScatteringProperties(mieContext, mieCoefficients);

        crossSectionsFile
            << TO_STR(wavelength)
            << ","
            << TO_STR(extinction)
            << ","
            << TO_STR(scattering)
            << "\n";

        asymmetryFile
            << TO_STR(wavelength)
            << ","
            << TO_STR(asymmetry)
            << "\n";

        for (int a = 0; a <= ANGLE_RESOLUTION; a++) {

            const double angle = 180.0 * (static_cast<double>(a) / ANGLE_RESOLUTION);

            const auto [senkrecht, parallel, unpolarized] = LorenzMie::calculatePhaseFunction(mieCoefficients, angle);

            phaseFile
                << TO_STR(wavelength)
                << ","
                << TO_STR(angle)
                << ","
                << TO_STR(unpolarized)
                << ","
                << TO_STR(senkrecht)
                << ","
                << TO_STR(parallel)
                << "\n";
        }
    }

    crossSectionsFile.close();
    asymmetryFile.close();
    phaseFile.close();

    const auto endTime   = std::chrono::high_resolution_clock::now();
    const auto duration = std::chrono::duration_cast<std::chrono::seconds>(endTime - startTime).count();

    std::cout
        << "Done! Took " << TO_STR(duration) << "s."
        << "\n"
        << "Results written to: " << std::filesystem::absolute(outputFolder) << "."
        << "\n";

    return 0;
}
