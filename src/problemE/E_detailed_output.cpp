// E.cpp
#define EIGEN_DONT_ALIGN_STATICALLY
#define EIGEN_DONT_VECTORIZE
#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <string>
#include <algorithm>
#include <Eigen/Dense>
#include "../packages/CurveFitting.hpp"
#include "../packages/NodeType.h"
#include "../packages/SplineType.h"


// Define Pi
const double PI = acos(-1.);

/**
 * @brief Calculates the angle distance between two points on the unit sphere.
 * 
 * @param p1 First point on the unit sphere.
 * @param p2 Second point on the unit sphere.
 * @return double The angle distance in radians.
 */
double sphericalAngleDistance(const Eigen::Vector3d& p1, const Eigen::Vector3d& p2) {
    double dot_product = p1.dot(p2);
    // Clamp the dot product to avoid numerical issues with acos
    dot_product = std::max(-1.0, std::min(1.0, dot_product));
    return std::acos(dot_product);
}

/**
 * @brief Generates a vector of parameter values based on node type and N.
 * 
 * @param N Number of nodes.
 * @param T Parameter range end (e.g., 2*Pi or 6*Pi).
 * @param actualPoints Actual curve points for ChordalLength node generation.
 * @return std::vector<double> Generated parameter values.
 */
std::vector<double> generateParameters(int N, double T, const std::vector<Eigen::Vector2d>& actualPoints) {
    std::vector<double> parameters;
    parameters.reserve(N + 1);
    
    for (double i = 0.0; i <= N; i = i+1.0) {
         parameters.push_back(i * T / N);
    }
    parameters[N]=T;
    
    return parameters;
}

/**
 * @brief Saves fitting results to a CSV file.
 * 
 * @param filename The output CSV filename.
 * @param t_values Parameter values.
 * @param fittingPoints Fitted curve points.
 * @param actualPoints Actual curve points.
 * @param is3D Boolean indicating if the curve is 3D.
 */
void saveFittingResults(const std::string& filename, const std::vector<double>& t_values,
                        const std::vector<Eigen::Vector2d>& fittingPoints2D,
                        const std::vector<Eigen::Vector2d>& actualPoints2D) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Failed to open file for writing: " << filename << std::endl;
        return;
    }
    
    // Write header
    file << "t,fittingPoints,actualPoints\n";
    
    // Write data
    for (size_t i = 0; i < t_values.size(); ++i) {
        file << std::fixed << std::setprecision(6) << t_values[i] << ",";
        file << std::fixed << std::setprecision(6) << "(" << fittingPoints2D[i].x() << "," << fittingPoints2D[i].y() << ")" << ",";
        file << std::fixed << std::setprecision(6) << "(" << actualPoints2D[i].x() << "," << actualPoints2D[i].y() << ")\n";
    }
    
    file.close();
}

/**
 * @brief Saves fitting results to a CSV file for 3D curves.
 * 
 * @param filename The output CSV filename.
 * @param t_values Parameter values.
 * @param fittingPoints3D Fitted curve points.
 * @param actualPoints3D Actual curve points.
 */
void saveFittingResults3D(const std::string& filename, const std::vector<double>& t_values,
                          const std::vector<Eigen::Vector3d>& fittingPoints3D,
                          const std::vector<Eigen::Vector3d>& actualPoints3D) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Failed to open file for writing: " << filename << std::endl;
        return;
    }
    
    // Write header
    file << "t,fittingPoints,actualPoints\n";
    
    // Write data
    for (size_t i = 0; i < t_values.size(); ++i) {
        file << std::fixed << std::setprecision(6) << t_values[i] << ",";
        file << std::fixed << std::setprecision(6) << "(" << fittingPoints3D[i].x() << "," << fittingPoints3D[i].y() << "," << fittingPoints3D[i].z() << ")" << ",";
        file << std::fixed << std::setprecision(6) << "(" << actualPoints3D[i].x() << "," << actualPoints3D[i].y() << "," << actualPoints3D[i].z() << ")\n";
    }
    
    file.close();
}

/**
 * @brief Saves max error results to a CSV file.
 * 
 * @param filename The output CSV filename.
 * @param results The vector of max error entries.
 */
void saveMaxErrors(const std::string& filename, const std::vector<std::tuple<int, std::string, std::string, std::string, std::string, double>>& results) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Failed to open max error file for writing: " << filename << std::endl;
        return;
    }
    
    // Write header
    file << "N,curve,splineType,nodeType,boundary,maxError\n";
    
    // Write data
    for (const auto& entry : results) {
        file << std::get<0>(entry) << "," << std::get<1>(entry) << "," << std::get<2>(entry) << "," << std::get<3>(entry) << "," << std::get<4>(entry) << "," << std::get<5>(entry);
        file << "\n";
    }
    
    file.close();
}

int main() {
    // Define output directory
    std::string outputDir = "../../output/problemE/";
    
    // Define numPoint for evaluation
    int numPoint = 501;
    
    // Define node counts
    std::vector<int> nodeNs = {10, 40, 160};
    
    // Define node types
    std::vector<NodeType> nodeTypes = {NodeType::Equidistant, NodeType::ChordalLength};
    
    // Define spline types
    std::vector<SplineType> splineTypes = {SplineType::PPForm, SplineType::BSpline};
    
    // Define boundary condition
    BoundaryCondition condition = BoundaryCondition::PERIODIC;
    

    
    // Define the three curves
    // Curve r1: Heart curve
    auto r1 = [](double t) -> Eigen::Vector2d {
        double x = std::sqrt(3.0) * std::sin(t);
        double y = (2.0 * std::sqrt(3.0) / 3.0) * std::cos(t) + (2.0 / 3.0) * std::sqrt(std::sqrt(3.0) * std::abs(std::sin(t)));
        return Eigen::Vector2d(x, y);
    };
    
    // Curve r2: Another planar curve
    auto r2 = [](double t) -> Eigen::Vector2d {
        double x = std::sin(t) + t * std::cos(t);
        double y = std::cos(t) - t * std::sin(t);
        return Eigen::Vector2d(x, y);
    };

    // Derivative of Curve r2: Another planar curve
    auto r2_deriv = [](double t) -> Eigen::Vector2d {
        double dx = 2*std::cos(t) + t * (-std::sin(t));
        double dy = -2*std::sin(t) - t * std::cos(t);
        return Eigen::Vector2d(dx, dy);
    };

    // Second derivative of Curve r2: Another planar curve
    auto r2_second_deriv = [](double t) -> Eigen::Vector2d {
        double ddx = -3*std::sin(t) + t * -std::cos(t);
        double ddy = -3*std::cos(t) + t * std::sin(t);
        return Eigen::Vector2d(ddx, ddy);
    };
    
    // Curve r3: Spherical curve
    auto r3 = [&](double t) -> Eigen::Vector3d {
        double u = std::cos(t);
        double v = std::sin(t);
        double x = std::sin(u) * std::cos(v);
        double y = std::sin(u) * std::sin(v);
        double z = std::cos(u);
        return Eigen::Vector3d(x, y, z);
    };
    
    // Define parameter ranges
    struct CurveInfo {
        std::string name;
        std::function<Eigen::Vector2d(double)> func2D;
        std::function<Eigen::Vector3d(double)> func3D;
        bool is3D;
        double t_start;
        double t_end;
        BoundaryCondition condition;
    };
    
    std::vector<CurveInfo> curves = {
        {"r1", r1, nullptr, false, 0, 2.0 * PI, BoundaryCondition::PERIODIC},
        {"r1", r1, nullptr, false, 0, 2.0 * PI, BoundaryCondition::NOT_A_KNOT},
        {"r2", r2, nullptr, false, 0.0, 6.0 * PI, BoundaryCondition::COMPLETE},
        {"r2", r2, nullptr, false, 0.0, 6.0 * PI, BoundaryCondition::NOT_A_KNOT},
        {"r2", r2, nullptr, false, 0.0, 6.0 * PI, BoundaryCondition::SECOND},
        {"r3", nullptr, r3, true, 0.0, 2.0 * PI, BoundaryCondition::NOT_A_KNOT},
        {"r3", nullptr, r3, true, 0.0, 2.0 * PI, BoundaryCondition::PERIODIC}};
    
    std::vector<std::tuple<int, std::string, std::string, std::string, std::string, double>> maxErrorResults;
    
    // Iterate over each curve
    for (const auto& curve : curves) {
        std::cout << "Processing curve: " << curve.name << std::endl;
        
        // Determine if the curve is 3D
        bool is3D = curve.is3D;

        condition = curve.condition;
        
        // Generate actual points for chordal length
        int actualPointsCount = 501; // Enough points for chordal length approximation
        std::vector<double> actual_t_values;
        actual_t_values.reserve(actualPointsCount);
        std::vector<Eigen::Vector2d> actualPoints2D;
        std::vector<Eigen::Vector3d> actualPoints3D;
        
        for (int i = 0; i < actualPointsCount; ++i) {
            double t = curve.t_start + i * (curve.t_end - curve.t_start) / (actualPointsCount - 1);
            actual_t_values.push_back(t);
            if (!is3D) {
                actualPoints2D.push_back(curve.func2D(t));
            }
            else {
                actualPoints3D.push_back(curve.func3D(t));
            }
        }
        
        // Iterate over each node count N
        for (const auto& N : nodeNs) {
            // Iterate over each node type
            for (const auto& nodeType : nodeTypes) {
                // Generate node parameters based on node type
                std::vector<double> nodeParameters = generateParameters(N, curve.t_end, actualPoints2D);

                // Iterate over each spline type
                for (const auto& splineType : splineTypes) {
                    // Initialize CurveFitting object
                    CurveFitting fitter;
                    
                    // Collect node points
                    if (!is3D) {
                        std::vector<Eigen::Vector2d> nodePoints;
                        Eigen::Vector2d start = Eigen::Vector2d::Zero();
                        Eigen::Vector2d end = Eigen::Vector2d::Zero();
                        nodePoints.reserve(N + 1);
                        for (const auto& t : nodeParameters) {
                            nodePoints.push_back(curve.func2D(t));
                        }
                        if (condition == BoundaryCondition::PERIODIC) {
                            nodePoints[N] = nodePoints[0];  
                        }
                        if (condition == BoundaryCondition::COMPLETE) {
                            start = r2_deriv(curve.t_start);
                            end = r2_deriv(curve.t_end);
                        }
                        if (condition == BoundaryCondition::SECOND) {
                            start = r2_second_deriv(curve.t_start);
                            end = r2_second_deriv(curve.t_end);
                        }
                        fitter.planeCurve(nodeParameters, nodePoints, nodeType, splineType, condition, start, end);
                    }
                    else {
                        std::vector<Eigen::Vector3d> nodePoints3D;
                        nodePoints3D.reserve(N + 1);
                        for (const auto& t : nodeParameters) {
                            nodePoints3D.push_back(curve.func3D(t));
                        }
                        if (condition == BoundaryCondition::PERIODIC) {
                            nodePoints3D[N] = nodePoints3D[0];  
                        }
                        fitter.sphereCurve(nodeParameters, nodePoints3D, nodeType, splineType, condition);
                    }
                    
                    // Generate evaluation points
                    std::vector<double> eval_t_values;
                    eval_t_values.reserve(numPoint);
                    std::vector<Eigen::Vector2d> fittedPoints2D;
                    std::vector<Eigen::Vector3d> fittedPoints3D;
                    std::vector<Eigen::Vector2d> actualEvalPoints2D;
                    std::vector<Eigen::Vector3d> actualEvalPoints3D;
                    
                    for (int i = 0; i < numPoint; ++i) {
                        double t_eval = curve.t_start + (double)i * (curve.t_end - curve.t_start) / (numPoint - 1);
                        if (i == numPoint - 1) {
                            t_eval = curve.t_end;
                        }
            
                        eval_t_values.push_back(t_eval);
                        if (!is3D) {
                            Eigen::Vector2d fitted = fitter.plane_evaluate(t_eval);
                            fittedPoints2D.push_back(fitted);
                            Eigen::Vector2d actual = curve.func2D(t_eval);
                            actualEvalPoints2D.push_back(actual);
                        }
                        else {
                            Eigen::Vector3d fitted = fitter.sphere_evaluate(t_eval);
                            fittedPoints3D.push_back(fitted);
                            Eigen::Vector3d actual = curve.func3D(t_eval);
                            actualEvalPoints3D.push_back(actual);
                        }
                    }
                    
                    // Prepare actual points for error calculation
                    std::vector<Eigen::Vector2d> actualPointsFitting2D;
                    std::vector<Eigen::Vector3d> actualPointsFitting3D;
                    
                    if (!is3D) {
                        actualPointsFitting2D.reserve(numPoint);
                        for (const auto& t_eval : eval_t_values) {
                            actualPointsFitting2D.push_back(curve.func2D(t_eval));
                        }
                    }
                    else {
                        actualPointsFitting3D.reserve(numPoint);
                        for (const auto& t_eval : eval_t_values) {
                            actualPointsFitting3D.push_back(curve.func3D(t_eval));
                        }
                    }
                    
                    // Save fitting results to CSV
                    std::string splineTypeStr = (splineType == SplineType::PPForm) ? "PPForm" : "BSpline";
                    std::string nodeTypeStr = (nodeType == NodeType::Equidistant) ? "Equidistant" : "ChordalLength";
                    std::string boundaryStr;
                    switch (condition) {
                        case BoundaryCondition::PERIODIC:
                            boundaryStr = "PERIODIC";
                            break;
                        case BoundaryCondition::COMPLETE:
                            boundaryStr = "COMPLETE";
                            break;
                        case BoundaryCondition::SECOND:
                            boundaryStr = "SECOND";
                            break;
                        case BoundaryCondition::NOT_A_KNOT:
                            boundaryStr = "NOT_A_KNOT";
                            break;
                    }
                    
                    std::string filename;
                    if (!is3D) {
                        filename = outputDir + curve.name + "_N_" + std::to_string(N) + "_" + splineTypeStr + "_" + nodeTypeStr + ".csv";
                        saveFittingResults(filename, eval_t_values, fittedPoints2D, actualEvalPoints2D);
                    }
                    else {
                        filename = outputDir + curve.name + "_N_" + std::to_string(N) + "_" + splineTypeStr + "_" + nodeTypeStr + ".csv";
                        saveFittingResults3D(filename, eval_t_values, fittedPoints3D, actualEvalPoints3D);
                    }
                    
                    // Compute max error
                    double maxError = 0.0;

                    if (nodeType != NodeType::ChordalLength)
                    {
                        if (!is3D)
                        {
                            for (size_t i = 0; i < numPoint; ++i)
                            {
                                double error = (fittedPoints2D[i] - actualEvalPoints2D[i]).norm();
                                if (error > maxError)
                                {
                                    maxError = error;
                                }
                            }
                        }
                        else
                        {
                            for (size_t i = 0; i < numPoint; ++i)
                            {
                                double angleDist = sphericalAngleDistance(fittedPoints3D[i], actualEvalPoints3D[i]);
                                if (angleDist > maxError)
                                {
                                    maxError = angleDist;
                                }
                            }
                        }

                        // Store max error result
                        maxErrorResults.emplace_back(N, curve.name, splineTypeStr, nodeTypeStr, boundaryStr, maxError);
                    }
                }
            }
        }
    }
    
    // Sort maxErrorResults by N ascending, then curve, then splineType, then nodeType, then boundary
    std::sort(maxErrorResults.begin(), maxErrorResults.end(),
        [](const std::tuple<int, std::string, std::string, std::string, std::string, double>& a,
           const std::tuple<int, std::string, std::string, std::string, std::string, double>& b) -> bool {
               if (std::get<0>(a) != std::get<0>(b))
                   return std::get<0>(a) < std::get<0>(b);
               if (std::get<1>(a) != std::get<1>(b))
                   return std::get<1>(a) < std::get<1>(b);
               if (std::get<2>(a) != std::get<2>(b))
                   return std::get<2>(a) < std::get<2>(b);
               if (std::get<3>(a) != std::get<3>(b))
                   return std::get<3>(a) < std::get<3>(b);
               if (std::get<4>(a) != std::get<4>(b))
                   return std::get<4>(a) < std::get<4>(b);
               return std::get<5>(a) < std::get<5>(b);
           });
    
    // Save max error results to CSV
    std::string maxErrorFilename = outputDir + "maxerror.csv";
    saveMaxErrors(maxErrorFilename, maxErrorResults);
    
    std::cout << "Fitting and error computation completed. Results saved in " << outputDir << std::endl;
    
    return 0;
}