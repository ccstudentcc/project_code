// A.cpp
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

// Define the exact function f(x) = 1 / (1 + 25x^2)
double f_exact(double x) {
    return 1.0 / (1.0 + 25.0 * x * x);
}

// First derivative of f(x): f'(x) = -50x / (1 + 25x^2)^2
double f_prime(double x) {
    return -50.0 * x / std::pow(1.0 + 25.0 * x * x, 2);
}

// Second derivative of f(x): f''(x) = (-50 + 3750x^2) / (1 + 25x^2)^3
double f_double_prime(double x) {
    return (-50.0 + 3750.0 * x * x) / std::pow(1.0 + 25.0 * x * x, 3);
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

// Helper function to compute max-norm at midpoints
double computeMaxNormAtMidpoints(const std::vector<double>& nodes,
                                 const std::vector<double>& interpolatedValues,
                                 const std::vector<double>& exactValues) {
    double max_norm = 0.0;
    int N = nodes.size();
    for(int i = 0; i < N - 1; ++i){
        double midpoint = (nodes[i] + nodes[i+1]) / 2.0;
        double interpolated = 0.0;
        // Linear interpolation between nodes[i] and nodes[i+1]
        double slope = (interpolatedValues[i+1] - interpolatedValues[i]) / (nodes[i+1] - nodes[i]);
        interpolated = interpolatedValues[i] + slope * (midpoint - nodes[i]);
        double exact = f_exact(midpoint);
        double error = std::abs(interpolated - exact);
        if(error > max_norm){
            max_norm = error;
        }
    }
    return max_norm;
}

int main() {
    // Create 'output/problemA' directory if it doesn't exist
    std::string outputDir = "../../output/problemA";
    
    // Define the values of N
    std::vector<int> Ns = {6, 11, 21, 41, 81};
    
    // Define the boundary conditions to test
    std::vector<BoundaryCondition> bc_list = {
        BoundaryCondition::NATURAL,
        BoundaryCondition::COMPLETE,
        BoundaryCondition::NOT_A_KNOT,
        BoundaryCondition::SECOND
    };
    
    // Open file to record max errors and convergence rates
    std::ofstream maxErrorConvergenceFile(outputDir + "/A_maxerror_convergence.csv");
    maxErrorConvergenceFile << "N,Method,BoundaryCondition,MaxError,ConvergenceRate\n";
    
    // =========================
    // Open file to record max-norm errors at midpoints
    // =========================
    std::ofstream maxNormFile(outputDir + "/A_maxNorm.csv");
    maxNormFile << "N,Method,BoundaryCondition,MaxNorm\n";
    
    // Define the number of plot points
    const int numPlotPoints = 801;
    
    // Generate  uniformly distributed points in [-1, 1]
    std::vector<double> x_plot(numPlotPoints);
    double step_plot = 2.0 / (numPlotPoints - 1);
    for(int i = 0; i < numPlotPoints; ++i){
        x_plot[i] = -1.0 + i * step_plot;
    }
    
    // Compute exact function values at plot points
    std::vector<double> f_plot_values(numPlotPoints);
    for(int i = 0; i < numPlotPoints; ++i){
        f_plot_values[i] = f_exact(x_plot[i]);
    }
    
    // Structure to hold previous max errors for convergence rate calculation
    std::map<std::pair<std::string, std::string>, double> previousMaxErrors;
    
    // Loop over each N
    for(auto N : Ns) {
        std::cout << "Processing N = " << N << "...\n";
        
        // Generate N evenly spaced nodes in [-1, 1]
        std::vector<double> nodes(N);
        double step = 2.0 / (N - 1);
        for(int i = 0; i < N; ++i){
            nodes[i] = -1.0 + i * step;
        }
        
        // Compute f(x) at nodes
        std::vector<double> f_values(N);
        for(int i = 0; i < N; ++i){
            f_values[i] = f_exact(nodes[i]);
        }
        
        // Compute derivatives at endpoints
        double f_prime_start = f_prime(nodes[0]);
        double f_prime_end = f_prime(nodes[N-1]);
        double f_double_prime_start = f_double_prime(nodes[0]);
        double f_double_prime_end = f_double_prime(nodes[N-1]);
        
        // =========================
        // PPForm<1> - Linear Spline
        // =========================
        PPForm<1> ppform1;
        ppform1.setInterpolationPoints(nodes);
        ppform1.setInterpolationValues(f_values);
        
        // Evaluate at plot points
        std::vector<double> ppform1_evals(numPlotPoints, 0.0);
        double max_error_ppform1 = 0.0;
        for(int i = 0; i < numPlotPoints; ++i){
            ppform1_evals[i] = ppform1.linearSpline(x_plot[i]);
            double error = std::abs(ppform1_evals[i] - f_plot_values[i]);
            if(error > max_error_ppform1){
                max_error_ppform1 = error;
            }
        }
        
        // Save PPForm<1> data to CSV
        std::ofstream ppform1File(outputDir + "/PPForm1_N" + std::to_string(N) + ".csv");
        ppform1File << "x,PPForm1,f_exact(x)\n";
        for(int i = 0; i < numPlotPoints; ++i){
            ppform1File << std::fixed << std::setprecision(8)
                       << x_plot[i] << "," << ppform1_evals[i] << "," 
                       << f_plot_values[i] << "\n";
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
        // Compute Max-Norm at Midpoints for PPForm<1>
        // =========================
        double maxNorm_ppform1 = computeMaxNormAtMidpoints(nodes, f_values, f_values); // Note: interpolatedValues = f_values for exact function
        
        // Since PPForm<1> is exact at nodes, but we need interpolated spline values
        // Modify computeMaxNormAtMidpoints to use interpolated spline values
        // Here, ppform1_evals contains interpolated values at plot points, but we need interpolated values at midpoints
        // Therefore, we'll compute interpolated values at midpoints manually
        
        // Generate midpoints
        std::vector<double> midpoints(N-1);
        for(int i = 0; i < N-1; ++i){
            midpoints[i] = (nodes[i] + nodes[i+1]) / 2.0;
        }
        
        // Evaluate PPForm<1> at midpoints and compute errors
        std::vector<double> ppform1_mid_evals(N-1, 0.0);
        double max_norm_ppform1_final = 0.0;
        for(int i = 0; i < N-1; ++i){
            ppform1_mid_evals[i] = ppform1.linearSpline(midpoints[i]);
            double exact = f_exact(midpoints[i]);
            double error = std::abs(ppform1_mid_evals[i] - exact);
            if(error > max_norm_ppform1_final){
                max_norm_ppform1_final = error;
            }
        }
        
        // Record max-norm error to A_maxNorm.csv
        maxNormFile << N << "," << method_ppform1 << "," << bc_ppform1 << "," 
                    << max_norm_ppform1_final << "\n";
        
        // =========================
        // BSpline<1> - Linear B-spline
        // =========================
        // For linear B-splines, knot vector can be the same as nodes
        std::vector<double> knots_bspline1 = nodes;
        BSpline<1> bspline1;
        bspline1.setKnots(knots_bspline1);
        bspline1.setControlPoints(f_values);
        bspline1.computeCoefficients(); // Ensure coefficients are computed
        
        // Evaluate at plot points
        std::vector<double> bspline1_evals(numPlotPoints, 0.0);
        double max_error_bspline1 = 0.0;
        for(int i = 0; i < numPlotPoints; ++i){
            bspline1_evals[i] = bspline1.evaluate(x_plot[i]);
            double error = std::abs(bspline1_evals[i] - f_plot_values[i]);
            if(error > max_error_bspline1){
                max_error_bspline1 = error;
            }
        }
        
        // Save BSpline<1> data to CSV
        std::ofstream bspline1File(outputDir + "/BSpline1_N" + std::to_string(N) + ".csv");
        bspline1File << "x,BSpline1,f_exact(x)\n";
        for(int i = 0; i < numPlotPoints; ++i){
            bspline1File << std::fixed << std::setprecision(8)
                        << x_plot[i] << "," << bspline1_evals[i] << "," 
                        << f_plot_values[i] << "\n";
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
        // Compute Max-Norm at Midpoints for BSpline<1>
        // =========================
        // Generate midpoints
        // Already generated earlier as 'midpoints'
        
        // Evaluate BSpline<1> at midpoints and compute errors
        std::vector<double> bspline1_mid_evals(N-1, 0.0);
        double max_norm_bspline1_final = 0.0;
        for(int i = 0; i < N-1; ++i){
            bspline1_mid_evals[i] = bspline1.evaluate(midpoints[i]);
            double exact = f_exact(midpoints[i]);
            double error = std::abs(bspline1_mid_evals[i] - exact);
            if(error > max_norm_bspline1_final){
                max_norm_bspline1_final = error;
            }
        }
        
        // Record max-norm error to A_maxNorm.csv
        maxNormFile << N << "," << method_bspline1 << "," << bc_bspline1 << "," 
                    << max_norm_bspline1_final << "\n";
        
        // =========================
        // PPForm<3> - Cubic Spline with Various Boundary Conditions
        // =========================
        PPForm<3> ppform3;
        ppform3.setInterpolationPoints(nodes);
        ppform3.setInterpolationValues(f_values);
        
        // Structure to hold evaluations for all boundary conditions
        std::vector<std::vector<double>> ppform3_evals_list(bc_list.size(), std::vector<double>(numPlotPoints, 0.0));
        std::vector<double> ppform3_max_errors(bc_list.size(), 0.0);
        
        // Structure to hold max-norm errors at midpoints for all boundary conditions
        std::vector<double> ppform3_maxNorm_errors(bc_list.size(), 0.0);
        
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
            
            // Evaluate cubic spline at plot points and compute max error
            for(int i = 0; i < numPlotPoints; ++i){
                try {
                    ppform3_evals_list[bc_idx][i] = ppform3.cubicSpline(x_plot[i], bc);
                }
                catch (const std::exception &e){
                    std::cerr << "Error evaluating PPForm<3>: " << e.what() << "\n";
                    ppform3_evals_list[bc_idx][i] = std::numeric_limits<double>::quiet_NaN();
                }
                double error = std::abs(ppform3_evals_list[bc_idx][i] - f_plot_values[i]);
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
        for(int i = 0; i < numPlotPoints; ++i){
            ppform3File << std::fixed << std::setprecision(8) << x_plot[i];
            for(size_t bc_idx = 0; bc_idx < bc_list.size(); ++bc_idx){
                ppform3File << "," << ppform3_evals_list[bc_idx][i];
            }
            ppform3File << "," << f_plot_values[i] << "\n";
        }
        ppform3File.close();
        
        // Record max errors and compute max-norm errors for each boundary condition
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
            
            // =========================
            // Compute Max-Norm at Midpoints for PPForm<3>
            // =========================
            double maxNorm_ppform3 = 0.0;
            for(int i = 0; i < N-1; ++i){
                double midpoint = (nodes[i] + nodes[i+1]) / 2.0;
                double interpolated = ppform3.cubicSpline(midpoint, bc);
                double exact = f_exact(midpoint);
                double error = std::abs(interpolated - exact);
                if(error > maxNorm_ppform3){
                    maxNorm_ppform3 = error;
                }
            }
            ppform3_maxNorm_errors[bc_idx] = maxNorm_ppform3;
            
            // Record max-norm error to A_maxNorm.csv
            maxNormFile << N << "," << method_ppform3 << "," << bc_str << "," 
                        << ppform3_maxNorm_errors[bc_idx] << "\n";
        }
        
        // =========================
        // BSpline<3> - Cubic B-spline with Various Boundary Conditions
        // =========================
        BSpline<3> bspline3;
        bspline3.setKnots(nodes);
        bspline3.setControlPoints(f_values);
        
        // Structure to hold evaluations for all boundary conditions
        std::vector<std::vector<double>> bspline3_evals_list(bc_list.size(), std::vector<double>(numPlotPoints, 0.0));
        std::vector<double> bspline3_max_errors(bc_list.size(), 0.0);
        
        // Structure to hold max-norm errors at midpoints for all boundary conditions
        std::vector<double> bspline3_maxNorm_errors(bc_list.size(), 0.0);
        
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
            
            // Evaluate cubic B-spline at plot points and compute max error
            for(int i = 0; i < numPlotPoints; ++i){
                try {
                    bspline3_evals_list[bc_idx][i] = bspline3.evaluate(x_plot[i]);
                }
                catch(const std::exception &e){
                    std::cerr << "Error evaluating BSpline<3>: " << e.what() << "\n";
                    bspline3_evals_list[bc_idx][i] = std::numeric_limits<double>::quiet_NaN();
                }
                double error = std::abs(bspline3_evals_list[bc_idx][i] - f_plot_values[i]);
                if(error > bspline3_max_errors[bc_idx]){
                    bspline3_max_errors[bc_idx] = error;
                }
            }
            // =========================
            // Compute Max-Norm at Midpoints for BSpline<3>
            // =========================
            double maxNorm_bspline3 = 0.0;
            for(int i = 0; i < N-1; ++i){
                double midpoint = (nodes[i] + nodes[i+1]) / 2.0;
                double interpolated = bspline3.evaluate(midpoint);
                double exact = f_exact(midpoint);
                double error = std::abs(interpolated - exact);
                if(error > maxNorm_bspline3){
                    maxNorm_bspline3 = error;
                }
            }
            bspline3_maxNorm_errors[bc_idx] = maxNorm_bspline3;
            
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
        for(int i = 0; i < numPlotPoints; ++i){
            bspline3File << std::fixed << std::setprecision(8) << x_plot[i];
            for(size_t bc_idx = 0; bc_idx < bc_list.size(); ++bc_idx){
                bspline3File << "," << bspline3_evals_list[bc_idx][i];
            }
            bspline3File << "," << f_plot_values[i] << "\n";
        }
        bspline3File.close();
        
        // Record max errors and compute max-norm errors for each boundary condition
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
            
            // Record max-norm error to A_maxNorm.csv
            maxNormFile << N << "," << method_bspline3 << "," << bc_str << "," 
                        << bspline3_maxNorm_errors[bc_idx] << "\n";
            
        }
    }
    
    maxErrorConvergenceFile.close();
    maxNormFile.close(); 
    
    std::cout << "Interpolation completed. Results saved in the 'output/problemA' folder.\n";
    
    return 0;
}