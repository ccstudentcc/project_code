// C.cpp
#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <string>
#include <Eigen/Dense>
#include <iomanip> 
#include <algorithm>
#include <map>
#include <utility>
#include <sstream>
#include "../packages/BSpline.hpp"

// Define the exact function f(x) = 1 / (1 + x^2)
double f_exact(double x)
{
    return 1.0 / (1.0 + x * x);
}

// First derivative of f(x): f'(x) = -2x / (1 + x^2)^2
double f_prime(double x)
{
    return -2.0 * x / std::pow(1.0 + x * x, 2);
}

// Second derivative of f(x): f''(x) = (6x^2 - 2) / (1 + x^2)^3
double f_double_prime(double x)
{
    return (6.0 * x * x - 2.0) / std::pow(1.0 + x * x, 3);
}

// Convert BoundaryCondition enum to string for file naming
std::string boundaryConditionToString(BoundaryCondition bc)
{
    switch (bc)
    {
    case BoundaryCondition::PERIODIC:
        return "PERIODIC";
    case BoundaryCondition::NATURAL:
        return "NATURAL";
    case BoundaryCondition::COMPLETE:
        return "COMPLETE";
    case BoundaryCondition::NOT_A_KNOT:
        return "NOT_A_KNOT";
    case BoundaryCondition::SECOND:
        return "SECOND";
    default:
        return "UNKNOWN";
    }
}

int main()
{
    // Create 'output/problemC' directory if it doesn't exist
    std::string outputDir = "../../output/problemC";
    int num = 201;

    // Define the boundary conditions to test
    std::vector<BoundaryCondition> bc_list = {
        BoundaryCondition::NATURAL,
        BoundaryCondition::COMPLETE,
        BoundaryCondition::NOT_A_KNOT,
        BoundaryCondition::SECOND};

    int N3 = 11;
    int N2 = 10;

    // Cubic B-spline knots t_i = -5 + i, i = 0,...,N3-1
    std::vector<double> knots_cubic(N3);
    for (int i = 0; i < N3; ++i)
    {
        knots_cubic[i] = -5.0 + i;
    }

    // Quadratic B-spline knots t_i = i - 4.5, i = 0,...,N2-1
    std::vector<double> knots_quadratic(N2+2);
    knots_quadratic[0] = -5.0;
    for (int i = 0; i < N2; ++i)
    {
        knots_quadratic[i+1] = i - 4.5;
    }
    knots_quadratic[N2+1] = 5.0;

    // Compute control points for cubic B-spline
    std::vector<double> f_values_cubic(N3);
    for (int i = 0; i < N3; ++i)
    {
        f_values_cubic[i] = f_exact(knots_cubic[i]);
    }

    // Compute control points for quadratic B-spline
    std::vector<double> f_values_quadratic(N2+2);
    for (int i = 0; i < N2+2; ++i)
    {
        f_values_quadratic[i] = f_exact(knots_quadratic[i]);
    }

    // Compute derivatives if needed
    double f_prime_start = f_prime(knots_cubic[0]);
    double f_prime_end = f_prime(knots_cubic[N3 - 1]);
    double f_double_prime_start = f_double_prime(knots_cubic[0]);
    double f_double_prime_end = f_double_prime(knots_cubic[N3 - 1]);

    // =========================
    // BSpline<2> - quadratic B-spline
    BSpline<2> bspline2(knots_quadratic, f_values_quadratic);
    std::vector<std::pair<double, double>> bspline2_plot_points(num);
    bspline2_plot_points = bspline2.generatePlotPoints(num);
    // Save BSpline<2> data to CSV
    std::ofstream bspline2File(outputDir + "/BSpline2_N" + std::to_string(N2) + ".csv");
    bspline2File << "x,BSpline2,f_exact(x)\n";
    for (int i = 0; i < num ; ++i)
    {
        bspline2File << std::fixed << std::setprecision(8)
                     << bspline2_plot_points[i].first << "," << bspline2_plot_points[i].second << ","
                     << f_exact(bspline2_plot_points[i].first) << "\n";
    }
    bspline2File.close();
    std::cout << "Quatratic B-spline interpolation completed. Results saved in " << outputDir + "/BSpline2_N" + std::to_string(N2) + ".csv" << "\n";


    // =========================
    // BSpline<3> - Cubic B-spline with Various Boundary Conditions
    // =========================
    // Define degree for cubic B-spline

    BSpline<3> bspline3;
    bspline3.setKnots(knots_cubic);
    bspline3.setControlPoints(f_values_cubic);
    std::vector<double> plotPoint3(num);
    double start3 = knots_cubic.front();
    double end3 = knots_cubic.back();
    double step3 = (end3 - start3) / (num - 1);

    for (int i = 0; i < num; ++i)
    {
        plotPoint3[i] = start3 + i * step3;
    }

    // Structure to hold evaluations for all boundary conditions
    std::vector<std::vector<double>> bspline3_evals_list(bc_list.size(), std::vector<double>(num, 0.0));

    for (size_t bc_idx = 0; bc_idx < bc_list.size(); ++bc_idx)
    {
        BoundaryCondition bc = bc_list[bc_idx];

        // Set boundary conditions and derivatives as needed
        if (bc == BoundaryCondition::COMPLETE)
        {
            bspline3.setDerivatives(f_prime_start, f_prime_end);
        }
        if (bc == BoundaryCondition::SECOND)
        {
            bspline3.setSecondDerivatives(f_double_prime_start, f_double_prime_end);
        }

        // Set the boundary condition
        bspline3.setCondition(bc);

        // Compute coefficients
        bspline3.computeCoefficients();

        for (int i = 0; i < num; ++i){
            bspline3_evals_list[bc_idx][i] = bspline3.evaluate(plotPoint3[i]);
        }
    }

    // Save BSpline<3> data to CSV with all boundary conditions
    std::ofstream bspline3File(outputDir + "/BSpline3_N" + std::to_string(N3) + ".csv");
    // Write header
    bspline3File << "x";
    for (auto bc : bc_list)
    {
        bspline3File << ",BSpline3_" << boundaryConditionToString(bc);
    }
    bspline3File << ",f_exact(x)\n";

    // Write data rows
    for (int i = 0; i < num; ++i)
    {
        bspline3File << std::fixed << std::setprecision(8) << plotPoint3[i];
        for (size_t bc_idx = 0; bc_idx < bc_list.size(); ++bc_idx)
        {
            bspline3File << "," << bspline3_evals_list[bc_idx][i];
        }
        bspline3File << "," << f_exact(plotPoint3[i]) << "\n";
    }
    bspline3File.close();

    std::cout << "Cubic B-spline interpolation completed. Results saved in " << outputDir + "/BSpline3_N" + std::to_string(N3) + ".csv" << "\n";

     // =========================
    // BSpline<1> - Linear B-spline
    // -------------------------

    // Define knots for Linear B-spline (degree=1)
    std::vector<double> knots_linear = knots_cubic; // Same as cubic knots
    size_t n_linear_control_points = knots_linear.size() ; 

    // Generate control points for Linear B-spline by sampling exact function at knot positions
    std::vector<double> f_values_linear(n_linear_control_points);
    for (size_t i = 0; i < n_linear_control_points; ++i)
    {
        // For linear spline, control points typically coincide with data points
        f_values_linear[i] = f_exact(knots_linear[i]); // Avoid duplicate points at ends
    }

    // Instantiate BSpline<1>
    BSpline<1> bspline1(knots_linear, f_values_linear);

    // Generate plot points for Linear B-spline
    std::vector<std::pair<double, double>> bspline1_plot_points = bspline1.generatePlotPoints(num);

    // Save BSpline<1> data to CSV
    std::ofstream bspline1File(outputDir + "/BSpline1_N" + std::to_string(knots_linear.size()) + ".csv");
    bspline1File << "x,BSpline1,f_exact(x)\n";
    for (int i = 0; i < num; ++i)
    {
        bspline1File << std::fixed << std::setprecision(8)
                     << bspline1_plot_points[i].first << "," << bspline1_plot_points[i].second << ","
                     << f_exact(bspline1_plot_points[i].first) << "\n";
    }
    bspline1File.close();

    std::cout << "Linear B-spline interpolation completed. Results saved in "
              << outputDir + "/BSpline1_N" + std::to_string(knots_linear.size()) + ".csv" << "\n";

    return 0;
}