// D.cpp
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
    std::string outputDir = "../../output/problemD";

    int num = 1001;
    // Generate uniformly distributed points in [-5, 5]
    std::vector<double> x_plot(num);
    double step_plot = 10.0 / (num - 1);
    for(int i = 0; i < num; ++i) {
        x_plot[i] = -5.0 + i * step_plot;
    }
    
    // Compute exact function values at plot points
    std::vector<double> f_plot_values(num);
    for(int i = 0; i < num; ++i) {
        f_plot_values[i] = f_exact(x_plot[i]);
    }

    // Define the boundary conditions to test
    std::vector<BoundaryCondition> bc_list = {
        BoundaryCondition::NATURAL,
        BoundaryCondition::COMPLETE,
        BoundaryCondition::NOT_A_KNOT,
        BoundaryCondition::SECOND};

    int N3 = 11;
    int N2 = 10;

    std::vector<double> error_x_points = {-3.5, -3.0, -0.5, 0.0, 0.5, 3.0, 3.5};

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
    
    std::cout << "Quadratic B-spline interpolation completed.\n";


    // =========================
    // BSpline<3> - Cubic B-spline with Various Boundary Conditions
    // =========================
    // Define degree for cubic B-spline

    BSpline<3> bspline3;
    bspline3.setKnots(knots_cubic);
    bspline3.setControlPoints(f_values_cubic);
    std::vector<std::vector<double>> bspline3_evals_list(bc_list.size(), std::vector<double>(error_x_points.size(), 0.0));
    std::vector<std::vector<double>> bspline3_maxerror_list(bc_list.size(), std::vector<double>(num, 0.0));

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

        for (int i = 0; i < error_x_points.size(); ++i){
            bspline3_evals_list[bc_idx][i] = bspline3.evaluate(error_x_points[i]);
        }
        for (int i = 0; i < num; ++i){
            bspline3_maxerror_list[bc_idx][i] = bspline3.evaluate(x_plot[i]);
        }
    }
    

    std::cout << "Cubic B-spline interpolation completed. \n";

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

    
    std::cout << "Linear B-spline interpolation completed.\n";

    

    // Create error.csv file and write header
    std::ofstream errorFile(outputDir + "/error.csv");
    errorFile << "x,BSpline1,BSpline2,BSpline3_NATURAL,BSpline3_COMPLETE,BSpline3_NOT_A_KNOT,BSpline3_SECOND\n";
    
    // Calculate E_S(x) for each x point
    for (int i=0; i < error_x_points.size(); i++)
    {
        double x = error_x_points[i];
        double E_BSpline1 = std::abs(bspline1.evaluate(x) - f_exact(x));
        double E_BSpline2 = std::abs(bspline2.evaluate(x) - f_exact(x));

        // E_S(x) for different boundary conditions of cubic B-spline
        std::vector<double> E_BSpline3_list(bc_list.size(),0.0);
        for (size_t bc_idx = 0; bc_idx < bc_list.size(); ++bc_idx)
        {
            E_BSpline3_list[bc_idx]=std::abs(bspline3_evals_list[bc_idx][i] - f_exact(x));
        }

        // Write data row
        errorFile << std::fixed << std::setprecision(4) << x << ",";
        errorFile << std::scientific << std::setprecision(4) << E_BSpline1 << ",";
        errorFile << E_BSpline2 << ",";

        for (int j = 0; j < bc_list.size(); j++)
        {
            if(j == bc_list.size()-1)
                errorFile << E_BSpline3_list[j];
            else{
                errorFile << E_BSpline3_list[j] << ",";
            }
        }
        errorFile << "\n";
    }
    errorFile.close();
    std::cout << "Interpolation errors have been calculated and saved in " << outputDir + "/error.csv" << "\n";



    // Track maximum errors
    double max_error_bspline1 = 0.0;
    double max_error_bspline2 = 0.0;
    std::vector<double> max_error_bspline3(bc_list.size(), 0.0);

    // Save maximum errors to maxerror.csv
    std::ofstream maxErrorFile(outputDir + "/maxerror.csv");

    maxErrorFile << "Method,BSpline1,BSpline2,";

    for(size_t bc_idx = 0; bc_idx < bc_list.size(); ++bc_idx) {
        maxErrorFile << "BSpline3_" << boundaryConditionToString(bc_list[bc_idx]);
        if(bc_idx!=bc_list.size()-1){
            maxErrorFile << ","; 
        }
        else{
            maxErrorFile << "\n";
        }
    }

    // Calculate maximum errors at all plot points
    for(int i = 0; i < num; ++i) {
        // Linear B-spline error
        double error1 = std::abs(bspline1.evaluate(x_plot[i]) - f_plot_values[i]);
        max_error_bspline1 = std::max(max_error_bspline1, error1);

        // Quadratic B-spline error  
        double error2 = std::abs(bspline2.evaluate(x_plot[i]) - f_plot_values[i]);
        max_error_bspline2 = std::max(max_error_bspline2, error2);

        // Cubic B-spline errors for each boundary condition
        for(size_t bc_idx = 0; bc_idx < bc_list.size(); ++bc_idx) {
            double error3 = std::abs(bspline3_maxerror_list[bc_idx][i] - f_plot_values[i]);
            max_error_bspline3[bc_idx] = std::max(max_error_bspline3[bc_idx], error3);
        }
    }

    
    maxErrorFile << "MaxError,";
    maxErrorFile << std::scientific << std::setprecision(4);
    maxErrorFile << max_error_bspline1 << "," << max_error_bspline2;
    for(size_t bc_idx = 0; bc_idx < bc_list.size(); ++bc_idx) {
        maxErrorFile << "," << max_error_bspline3[bc_idx];
    }
    maxErrorFile << "\n";
    maxErrorFile.close();

    std::cout << "Interpolation max errors have been calculated and saved in " << outputDir + "/maxerror.csv" << "\n";

    return 0;
}