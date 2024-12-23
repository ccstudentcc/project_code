// require5.cpp
#include <iostream>
#include <vector>
#include <fstream>
#include "../packages/BSpline.hpp"

void savePlotPointsToCSV(const std::vector<std::pair<double, double>>& plotPoints, const std::string& filename) {
    std::ofstream outfile(filename);
    for (const auto& point : plotPoints) {
        outfile << point.first << "," << point.second << "\n";
    }
    outfile.close();
    std::cout << "Data saved to " << filename << std::endl;
}

template<size_t Degree>
void generateAndSaveBSpline(const std::vector<double>& knots, const std::vector<double>& coefficients, const std::string& filename) {
    BSpline<Degree> spline;
    spline.setKnots(knots);
    spline.setCoefficients(coefficients);

    // Generate plot points
    std::vector<std::pair<double, double>> plotPoints = spline.generatePlotPoints(100);

    // Save to CSV file
    savePlotPointsToCSV(plotPoints, filename);
}

int main() {
    try {
        // Example 1: Cubic B-spline (degree 3)
        {
            std::vector<double> knots = {0, 1, 2, 3, 4, 5, 6};
            std::vector<double> coefficients = {0, 1, 0, -1, 0, 1, 0, -1, 0};
            generateAndSaveBSpline<3>(knots, coefficients, "../../output/requirement5/cubic_bspline.csv");
        }

        // Example 2: Quintic B-spline (degree 5)
        {
            std::vector<double> knots = {0, 1, 2, 3, 4};
            std::vector<double> coefficients = {0, 2, 3, 2, 0, -1, 2, 3, 2};
            generateAndSaveBSpline<5>(knots, coefficients, "../../output/requirement5/quintic_bspline.csv");
        }
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}

