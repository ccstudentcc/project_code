// test.cpp
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
#include "../packages/PPForm.hpp"
#include "../packages/BSpline.hpp"

const double Pi = acos(-1.);

// Define the exact function f(x) = sin(x)
double f_exact(double x) {
    return sin(x);
}

// First derivative of f(x): f'(x) = cos(x)
double f_prime(double x) {
    return cos(x);
}

// Second derivative of f(x): f''(x) = -sin(x)
double f_double_prime(double x) {
    return -sin(x);
}

// Convert BoundaryCondition enum to string for file naming
std::string boundaryConditionToString(BoundaryCondition bc) {
    switch (bc) {
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

int main() {
    // Create 'output/test' directory if it doesn't exist
    std::string outputDir = "../../output/test";
    
    // Define the values of N
    std::vector<int> Ns = {11, 41, 91};
    
    // Define the boundary conditions to test
    std::vector<BoundaryCondition> bc_list = {
        BoundaryCondition::PERIODIC,
        BoundaryCondition::NATURAL,
        BoundaryCondition::COMPLETE,
        BoundaryCondition::NOT_A_KNOT,
        BoundaryCondition::SECOND
    };
    
    // Open file to record max errors and convergence rates
    std::ofstream maxErrorConvergenceFile(outputDir + "/A_maxerror_convergence.csv");
    maxErrorConvergenceFile << "N,Method,BoundaryCondition,MaxError,ConvergenceRate\n";
    
    // Structure to hold previous max errors for convergence rate calculation
    std::map<std::pair<std::string, std::string>, double> previousMaxErrors;
    
    // Loop over each N
    for(auto N : Ns) {
        std::cout << "Processing N = " << N << "...\n";
        
        // Generate N non-uniformly spaced nodes in [0, 2\pi]
        std::vector<double> nodes(N);
        for(int i = 0; i < N; ++i){
            // Use a non-linear distribution, e.g., clustering near 0 and 2π
            double t = static_cast<double>(i) / (N - 1);
            nodes[i] = std::pow(t, 2) * 2.0 * Pi; // Example: quadratic distribution
        }
        
        // Compute f(x) at nodes
        std::vector<double> f_values(N);
        for(int i = 0; i < N; ++i){
            f_values[i] = f_exact(nodes[i]);
        }
        
        // =========================
        // Generate uniformly distributed evaluation points
        // =========================
        int evalPoints = 720;
        std::vector<double> uniform_points(evalPoints);
        double step = 2.0 * Pi / evalPoints;
        for(int i = 0; i < evalPoints; ++i){
            uniform_points[i] = i * step;
        }
        
        // Compute f(x) at uniform points
        std::vector<double> f_uniform_values(evalPoints);
        for(int i = 0; i < evalPoints; ++i){
            f_uniform_values[i] = f_exact(uniform_points[i]);
        }
        
        // =========================
        // PPForm<1> - Linear Spline
        // =========================
        PPForm<1> ppform1;
        ppform1.setInterpolationPoints(nodes);
        ppform1.setInterpolationValues(f_values);
        
        // Evaluate at uniform points
        std::vector<double> ppform1_evals(evalPoints);
        double max_error_ppform1 = 0.0;
        for(int i = 0; i < evalPoints; ++i){
            ppform1_evals[i] = ppform1.linearSpline(uniform_points[i]);
            double error = std::abs(ppform1_evals[i] - f_uniform_values[i]);
            if(error > max_error_ppform1){
                max_error_ppform1 = error;
            }
        }
        
        // Save PPForm<1> data to CSV
        std::ofstream ppform1File(outputDir + "/PPForm1_N" + std::to_string(N) + ".csv");
        ppform1File << "x,PPForm1,f_exact(x)\n";
        for(int i = 0; i < evalPoints; ++i){
            ppform1File << std::fixed << std::setprecision(8)
                       << uniform_points[i] << "," << ppform1_evals[i] << "," 
                       << f_uniform_values[i] << "\n";
        }
        ppform1File.close();
        
        // Calculate Convergence Rate for PPForm1
        std::string method_ppform1 = "PPForm1";
        std::string bc_ppform1 = "NA";
        double rate_ppform1 = 0.0;
        std::pair<std::string, std::string> key_ppform1 = {method_ppform1, bc_ppform1};
        
        if(previousMaxErrors.find(key_ppform1) != previousMaxErrors.end()){
            double previousError = previousMaxErrors[key_ppform1];
            if(previousError != 0.0){
                rate_ppform1 = std::log(max_error_ppform1 / previousError) / std::log(
                    (2.0 / (Ns[0] - 1)) / (2.0 / (N - 1))
                );
            }
        }
        previousMaxErrors[key_ppform1] = max_error_ppform1;
        
        // Record max error and convergence rate
        maxErrorConvergenceFile << N << "," << method_ppform1 << "," << bc_ppform1 << "," 
                                << max_error_ppform1 << "," << rate_ppform1 << "\n";
        
        // =========================
        // BSpline<1> - Linear B-spline
        // =========================
        // For linear B-splines, knot vector can be the same as nodes
        std::vector<double> knots_bspline1 = nodes;
        BSpline<1> bspline1;
        bspline1.setKnots(knots_bspline1);
        bspline1.setControlPoints(f_values);
        bspline1.computeCoefficients(); // Ensure coefficients are computed
        
        // Evaluate at uniform points
        std::vector<double> bspline1_evals(evalPoints);
        double max_error_bspline1 = 0.0;
        for(int i = 0; i < evalPoints; ++i){
            bspline1_evals[i] = bspline1.evaluate(uniform_points[i]);
            double error = std::abs(bspline1_evals[i] - f_uniform_values[i]);
            if(error > max_error_bspline1){
                max_error_bspline1 = error;
            }
        }
        
        // Save BSpline<1> data to CSV
        std::ofstream bspline1File(outputDir + "/BSpline1_N" + std::to_string(N) + ".csv");
        bspline1File << "x,BSpline1,f_exact(x)\n";
        for(int i = 0; i < evalPoints; ++i){
            bspline1File << std::fixed << std::setprecision(8)
                        << uniform_points[i] << "," << bspline1_evals[i] << "," 
                        << f_uniform_values[i] << "\n";
        }
        bspline1File.close();
        
        // Calculate Convergence Rate for BSpline1
        std::string method_bspline1 = "BSpline1";
        std::string bc_bspline1 = "NA";
        double rate_bspline1 = 0.0;
        std::pair<std::string, std::string> key_bspline1 = {method_bspline1, bc_bspline1};
        
        if(previousMaxErrors.find(key_bspline1) != previousMaxErrors.end()){
            double previousError = previousMaxErrors[key_bspline1];
            if(previousError != 0.0){
                rate_bspline1 = std::log(max_error_bspline1 / previousError) / std::log(
                    (2.0 / (Ns[0] - 1)) / (2.0 / (N - 1))
                );
            }
        }
        previousMaxErrors[key_bspline1] = max_error_bspline1;
        
        // Record max error and convergence rate
        maxErrorConvergenceFile << N << "," << method_bspline1 << "," << bc_bspline1 << "," 
                                << max_error_bspline1 << "," << rate_bspline1 << "\n";
        
        // =========================
        // PPForm<3> - Cubic Spline with Various Boundary Conditions
        // =========================
        PPForm<3> ppform3;
        ppform3.setInterpolationPoints(nodes);
        ppform3.setInterpolationValues(f_values);
        
        // Compute derivatives if needed
        double f_prime_start = f_prime(nodes[0]);
        double f_prime_end = f_prime(nodes[N-1]);
        double f_double_prime_start = f_double_prime(nodes[0]);
        double f_double_prime_end = f_double_prime(nodes[N-1]);
        
        // Structure to hold evaluations for all boundary conditions
        std::vector<std::vector<double>> ppform3_evals_list(bc_list.size(), std::vector<double>(evalPoints, 0.0));
        std::vector<double> ppform3_max_errors(bc_list.size(), 0.0);
        
        for(size_t bc_idx = 0; bc_idx < bc_list.size(); ++bc_idx){
            BoundaryCondition bc = bc_list[bc_idx];
            
            // Set boundary conditions and derivatives as needed
            if(bc == BoundaryCondition::COMPLETE){
                ppform3.setDx(f_prime_start, f_prime_end);
            }
            if(bc == BoundaryCondition::SECOND){
                ppform3.setDDx(f_double_prime_start, f_double_prime_end);
            }
            
            // Compute coefficients
            ppform3.computeCubicCoefficients(bc);
            
            // Evaluate cubic spline at uniform points
            for(int i = 0; i < evalPoints; ++i){
                try {
                    ppform3_evals_list[bc_idx][i] = ppform3.cubicSpline(uniform_points[i], bc);
                }
                catch (const std::exception &e){
                    std::cerr << "Error evaluating PPForm<3>: " << e.what() << "\n";
                    ppform3_evals_list[bc_idx][i] = std::numeric_limits<double>::quiet_NaN();
                }
                double error = std::abs(ppform3_evals_list[bc_idx][i] - f_uniform_values[i]);
                if(error > ppform3_max_errors[bc_idx]){
                    ppform3_max_errors[bc_idx] = error;
                }
            }
        }
        
        // Save PPForm<3> data to CSV with all boundary conditions
        std::ofstream ppform3File(outputDir + "/PPForm3_N" + std::to_string(N) + ".csv");
        // Write header
        ppform3File << "x";
        for(auto bc : bc_list){
            ppform3File << ",PPForm3_" << boundaryConditionToString(bc);
        }
        ppform3File << ",f_exact(x)\n";
        
        // Write data rows
        for(int i = 0; i < evalPoints; ++i){
            ppform3File << std::fixed << std::setprecision(8) << uniform_points[i];
            for(size_t bc_idx = 0; bc_idx < bc_list.size(); ++bc_idx){
                ppform3File << "," << ppform3_evals_list[bc_idx][i];
            }
            ppform3File << "," << f_uniform_values[i] << "\n";
        }
        ppform3File.close();
        
        // Record max errors and convergence rates for each boundary condition
        for(size_t bc_idx = 0; bc_idx < bc_list.size(); ++bc_idx){
            BoundaryCondition bc = bc_list[bc_idx];
            double max_error = ppform3_max_errors[bc_idx];
            std::string method_ppform3 = "PPForm3";
            std::string bc_str = boundaryConditionToString(bc);
            double rate_ppform3 = 0.0;
            std::pair<std::string, std::string> key_ppform3 = {method_ppform3, bc_str};
            
            if(previousMaxErrors.find(key_ppform3) != previousMaxErrors.end()){
                double previousError = previousMaxErrors[key_ppform3];
                if(previousError != 0.0){
                    rate_ppform3 = std::log(max_error / previousError) / std::log(
                        (2.0 / (Ns[0] - 1)) / (2.0 / (N - 1))
                    );
                }
            }
            previousMaxErrors[key_ppform3] = max_error;
            
            // Record max error and convergence rate
            maxErrorConvergenceFile << N << "," << method_ppform3 << "," << bc_str << "," 
                                    << max_error << "," << rate_ppform3 << "\n";
        }
        
        // =========================
        // BSpline<3> - Cubic B-spline with Various Boundary Conditions
        // =========================
        // Define knot vector for cubic B-spline
        std::vector<double> knots_bspline3 = nodes;
        int degree = 3;
        
        BSpline<3> bspline3;
        bspline3.setKnots(knots_bspline3);
        bspline3.setControlPoints(f_values);
        
        // Structure to hold evaluations for all boundary conditions
        std::vector<std::vector<double>> bspline3_evals_list(bc_list.size(), std::vector<double>(evalPoints, 0.0));
        std::vector<double> bspline3_max_errors(bc_list.size(), 0.0);
        
        for(size_t bc_idx = 0; bc_idx < bc_list.size(); ++bc_idx){
            BoundaryCondition bc = bc_list[bc_idx];
            
            // Set boundary conditions and derivatives as needed
            if(bc == BoundaryCondition::COMPLETE){
                bspline3.setDerivatives(f_prime_start, f_prime_end);
            }
            if(bc == BoundaryCondition::SECOND){
                bspline3.setSecondDerivatives(f_double_prime_start, f_double_prime_end);
            }
            
            // Set the boundary condition
            bspline3.setCondition(bc);
            
            // Compute coefficients
            bspline3.computeCoefficients();
            
            // Evaluate cubic B-spline at uniform points
            for(int i = 0; i < evalPoints; ++i){
                try {
                    bspline3_evals_list[bc_idx][i] = bspline3.evaluate(uniform_points[i]);
                }
                catch(const std::exception &e){
                    std::cerr << "Error evaluating BSpline<3>: " << e.what() << "\n";
                    bspline3_evals_list[bc_idx][i] = std::numeric_limits<double>::quiet_NaN();
                }
                double error = std::abs(bspline3_evals_list[bc_idx][i] - f_uniform_values[i]);
                if(error > bspline3_max_errors[bc_idx]){
                    bspline3_max_errors[bc_idx] = error;
                }
            }
        }
        
        // Save BSpline<3> data to CSV with all boundary conditions
        std::ofstream bspline3File(outputDir + "/BSpline3_N" + std::to_string(N) + ".csv");
        // Write header
        bspline3File << "x";
        for(auto bc : bc_list){
            bspline3File << ",BSpline3_" << boundaryConditionToString(bc);
        }
        bspline3File << ",f_exact(x)\n";
        
        // Write data rows
        for(int i = 0; i < evalPoints; ++i){
            bspline3File << std::fixed << std::setprecision(8) << uniform_points[i];
            for(size_t bc_idx = 0; bc_idx < bc_list.size(); ++bc_idx){
                bspline3File << "," << bspline3_evals_list[bc_idx][i];
            }
            bspline3File << "," << f_uniform_values[i] << "\n";
        }
        bspline3File.close();
        
        // Record max errors and convergence rates for each boundary condition
        for(size_t bc_idx = 0; bc_idx < bc_list.size(); ++bc_idx){
            BoundaryCondition bc = bc_list[bc_idx];
            double max_error = bspline3_max_errors[bc_idx];
            std::string method_bspline3 = "BSpline3";
            std::string bc_str = boundaryConditionToString(bc);
            double rate_bspline3 = 0.0;
            std::pair<std::string, std::string> key_bspline3 = {method_bspline3, bc_str};
            
            if(previousMaxErrors.find(key_bspline3) != previousMaxErrors.end()){
                double previousError = previousMaxErrors[key_bspline3];
                if(previousError != 0.0){
                    rate_bspline3 = std::log(max_error / previousError) / std::log(
                        (2.0 / (Ns[0] - 1)) / (2.0 / (N - 1))
                    );
                }
            }
            previousMaxErrors[key_bspline3] = max_error;
            
            // Record max error and convergence rate
            maxErrorConvergenceFile << N << "," << method_bspline3 << "," << bc_str << "," 
                                    << max_error << "," << rate_bspline3 << "\n";
        }
    }
    
    maxErrorConvergenceFile.close();
    
    std::cout << "Interpolation completed. Results saved in the 'output/test' folder.\n";
    
    return 0;
}
