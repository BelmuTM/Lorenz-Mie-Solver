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
#include <format>
#include <fstream>
#include <iostream>

constexpr double PARTICLE_RADIUS_LOWER_BOUND =  5'000; // Nanometers [nm]
constexpr double PARTICLE_RADIUS_UPPER_BOUND = 15'000; // Nanometers [nm]

constexpr int PARTICLE_RADIUS_RESOLUTION = 100; // Amount of different particle sizes accounted for in the bulk medium

constexpr int WAVELENGTH_LOWER_BOUND = 390; // Nanometers [nm]
constexpr int WAVELENGTH_UPPER_BOUND = 830; // Nanometers [nm]

constexpr int WAVELENGTH_RESOLUTION = 880; // Amount of wavelengths between the lower and upper bounds
constexpr int ANGLE_RESOLUTION      = 720; // Amount of angles between 0° and 180°

constexpr double NM_TO_M = 1e9;

constexpr double waterVolumeFraction          = 0.01 / 10e6; // LWC (liquid water content) / water density (1000 kg/m3)
constexpr double waterMeanVolume              = 5'000; // Nanometers [nm]
constexpr double waterVolumeStandardDeviation = 2'000; // Nanometers [nm]

const std::string outputFolder = OUTPUT_DIR;

const std::string crossSectionsOutput = outputFolder + "/cross_sections.csv";
const std::string asymmetryOutput     = outputFolder + "/asymmetry.csv";
const std::string phaseOutput         = outputFolder + "/phase.csv";

struct WavelengthResult {
    double wavelength = 0.0;

    LorenzMie::ScatteringProperties properties{};

    std::vector<LorenzMie::PhaseFunction> phases{};
};

std::string formatVal(const double val) {
    return std::format("{:.2f}", val);
}

std::string formatValPrecise(const double val) {
    return std::format("{:.15f}", val);
}


int main() {

    const auto startTime = std::chrono::high_resolution_clock::now();

    const std::string line = "------------------------------------------------------------";

    std::cout
        << line
        << "\n"
        << "Lorenz-Mie Solver"
        << "\n\n"
        << "Solving for parameters:"
        << "\n"
        << "    - Particle radius range: " << formatVal(PARTICLE_RADIUS_LOWER_BOUND) << "-" << formatVal(PARTICLE_RADIUS_UPPER_BOUND) << "nm"
        << "\n"
        << "    - Particle radius bins: " << std::to_string(PARTICLE_RADIUS_RESOLUTION)
        << "\n"
        << "    - Wavelength range: " << formatVal(WAVELENGTH_LOWER_BOUND) << "-" << formatVal(WAVELENGTH_UPPER_BOUND) << "nm"
        << "\n"
        << "    - Wavelengths: " << std::to_string(WAVELENGTH_RESOLUTION)
        << "\n"
        << "    - Angles: " << std::to_string(ANGLE_RESOLUTION)
        << "\n"
        << line
        << "\n";

    std::filesystem::create_directories(outputFolder);

    ArrayDouble scatteringAngles{};
    scatteringAngles.resize(ANGLE_RESOLUTION + 1);

    for (int a = 0; a <= ANGLE_RESOLUTION; a++) {
        scatteringAngles[a] = 180.0 * (static_cast<double>(a) / ANGLE_RESOLUTION);
    }

    std::vector<WavelengthResult> results(WAVELENGTH_RESOLUTION + 1);

    #pragma omp parallel for schedule(dynamic)
    for (int w = 0; w <= WAVELENGTH_RESOLUTION; w++) {

        const double wavelength = std::lerp(
            WAVELENGTH_LOWER_BOUND,
            WAVELENGTH_UPPER_BOUND,
            static_cast<double>(w) / WAVELENGTH_RESOLUTION
        );

        const ComplexDouble airIOR   = {RefractiveIndex::air(wavelength, 15.0, 101350.0, 0.0, 420.0), 0.0};
        const ComplexDouble waterIOR = {RefractiveIndex::water(wavelength), 0.0};

        LorenzMie::BulkContext bulkContext(
            wavelength,
            PARTICLE_RADIUS_LOWER_BOUND,
            PARTICLE_RADIUS_UPPER_BOUND,
            PARTICLE_RADIUS_RESOLUTION,
            {},
            airIOR,
            waterIOR
        );

        calculateBulkNumberDensities(bulkContext, waterVolumeFraction, waterMeanVolume, waterVolumeStandardDeviation);

        WavelengthResult& result = results[w];
        result.wavelength = wavelength;

        LorenzMie::calculateBulkProperties(bulkContext, scatteringAngles, result.properties, result.phases);
    }

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

    for (const auto& [wavelength, properties, phases] : results) {

        crossSectionsFile
            << formatValPrecise(wavelength) << ","
            << formatValPrecise(properties.extinction * NM_TO_M) << ","
            << formatValPrecise(properties.scattering * NM_TO_M) << "\n";

        asymmetryFile
            << formatValPrecise(wavelength) << ","
            << formatValPrecise(properties.asymmetry) << "\n";

        for (int a = 0; a <= ANGLE_RESOLUTION; a++) {

            phaseFile
                << formatValPrecise(wavelength) << ","
                << formatValPrecise(scatteringAngles[a]) << ","
                << formatValPrecise(phases[a].unpolarized) << ","
                << formatValPrecise(phases[a].senkrecht) << ","
                << formatValPrecise(phases[a].parallel) << "\n";
        }
    }

    crossSectionsFile.close();
    asymmetryFile.close();
    phaseFile.close();

    const auto endTime  = std::chrono::high_resolution_clock::now();
    const auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime).count();

    std::cout
        << "Done! Took " << formatVal(static_cast<double>(duration) / 1000.0) << "s."
        << "\n"
        << "Results written to: " << std::filesystem::absolute(outputFolder) << "."
        << "\n";

    return 0;
}
