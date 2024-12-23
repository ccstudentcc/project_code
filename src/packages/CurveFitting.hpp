// CurveFitting.hpp
#ifndef CURVEFITTING_HPP
#define CURVEFITTING_HPP

#include <vector>
#include <Eigen/Dense>
#include <stdexcept>
#include <cmath>
#include "PPForm.hpp"
#include "BSpline.hpp"
#include "BoundaryCondition.h"
#include "NodeType.h"
#include "SplineType.h"

/**
 * @class CurveFitting
 * @brief A class for curve fitting in different dimensions.
 * 
 * This class specializes in cubic spline fitting.
 * It provides methods for plane and spherical curve fitting.
 */
class CurveFitting {
public:
    /**
     * @brief Default constructor.
     */
    explicit CurveFitting();

    /**
     * @brief Fits a plane spline curve to the given 2D points.
     * 
     * @param parameters A vector of parameter values (not necessarily in [0,1]).
     * @param points A vector of 2D points to fit.
     * @param nodeType The type of node preprocessing to apply (default is Equidistant).
     * @param splineType The spline type to use for fitting (PPForm or BSpline).
     * @param condition The boundary condition to apply.
     * @param start The start information of the curve.
     * @param end The end information of the curve.
     * @throws std::invalid_argument if input sizes are inconsistent or insufficient.
     */
    void planeCurve(const std::vector<double>& parameters,
                const std::vector<Eigen::Vector2d>& points,
                NodeType nodeType = NodeType::Equidistant,
                SplineType splineType = SplineType::BSpline,
                BoundaryCondition condition = BoundaryCondition::PERIODIC,
                Eigen::Vector2d start = Eigen::Vector2d::Zero(),
                Eigen::Vector2d end = Eigen::Vector2d::Zero());

    /**
     * @brief Fits a spherical spline curve to the given 3D points.
     * 
     * @param parameters A vector of parameter values (not necessarily in [0,1]).
     * @param points A vector of 3D points to fit (must be on a sphere).
     * @param nodeType The type of node preprocessing to apply (default is Equidistant).
     * @param splineType The spline type to use for fitting (PPForm or BSpline).
     * @param condition The boundary condition to apply (default is Periodic).
     * @throws std::invalid_argument if input sizes are inconsistent or insufficient.
     */
    void sphereCurve(const std::vector<double>& parameters,
                     const std::vector<Eigen::Vector3d>& points,
                     NodeType nodeType = NodeType::Equidistant,
                     SplineType splineType = SplineType::BSpline,
                     BoundaryCondition condition = BoundaryCondition::PERIODIC);

    /**
     * @brief Returns the type of node preprocessing used.
     * 
     * @return The node preprocessing type.
     */
    bool is3D() const;

    /**
     * @brief Evaluate the plane curve at a given parameter t
     * 
     * @param t Input parameter (should be within the parameter range)
     * @return Corresponding 2D point on the plane
     */
    Eigen::Vector2d plane_evaluate(double t);

    /**
     * @brief Evaluate the spherical curve at a given parameter t
     * 
     * @param t Input parameter (should be within the parameter range)
     * @return Corresponding 3D point on the sphere
     */
    Eigen::Vector3d sphere_evaluate(double t);

    /**
     * @brief Generates a set of points on the fitted plane spline curve.
     * 
     * @param numPoints The number of points to generate along the curve.
     * @return A vector containing the generated 2D points on the plane spline curve.
     */
    std::vector<Eigen::Vector2d> generatePlotPoints_plane(int numPoints);

    /**
     * @brief Generates a set of points on the fitted spherical spline curve.
     * 
     * @param numPoints The number of points to generate along the spherical curve.
     * @return A vector containing the generated 3D points on the spherical spline curve.
     */
    std::vector<Eigen::Vector3d> generatePlotPoints_sphere(int numPoints);

private:
    NodeType nodeType_;           /**< Node preprocessing type */
    SplineType splineType_;       /**< Spline type used for fitting */
    BoundaryCondition condition_; /**< Boundary condition (periodic for cubic splines) */

    std::vector<double> parameters_;                /**< Parameter vector */
    double totalLength_;                            /**< Total length of the curve */
    std::vector<Eigen::Vector2d> planePoints_;      /**< Points for plane fitting */
    std::vector<Eigen::Vector3d> spherePoints_;     /**< Points for sphere fitting */

    Eigen::Vector3d northPole_;  /**< Selected North Pole point for projection */

    // For plane curve fitting using PPForm or BSpline
    PPForm<3> ppformX_;          /**< PPForm spline for X coordinates */
    PPForm<3> ppformY_;          /**< PPForm spline for Y coordinates */
    BSpline<3> bsplineX_;        /**< BSpline for X coordinates */
    BSpline<3> bsplineY_;        /**< BSpline for Y coordinates */

    bool is3D_;                  /**< Flag for 3D curve fitting */

    /**
     * @brief Preprocesses the parameters based on the node type.
     */
    void preprocessParameters();

    /**
     * @brief Calculates a new North Pole point that is relatively far from the data points.
     */
    void calculateNorthPole();

    /**
     * @brief Performs stereographic projection from 3D points on a sphere to a 2D plane using the selected North Pole.
     * 
     * @param point A 3D point on the sphere.
     * @return The corresponding 2D point on the plane.
     */
    Eigen::Vector2d stereographicProjection(const Eigen::Vector3d& point) const;

    /**
     * @brief Maps a 2D point on the plane back to the 3D sphere using inverse stereographic projection with the selected North Pole.
     * 
     * @param point A 2D point on the plane.
     * @return The corresponding 3D point on the sphere.
     */
    Eigen::Vector3d inverseStereographicProjection(const Eigen::Vector2d& point) const;
};

// Implementation of CurveFitting

CurveFitting::CurveFitting()
    : nodeType_(NodeType::Equidistant),
      splineType_(SplineType::BSpline),
      condition_(BoundaryCondition::PERIODIC),
      totalLength_(0.0),
      northPole_(Eigen::Vector3d(0, 0, 1)) {
    parameters_.clear();
    planePoints_.clear();
    spherePoints_.clear();
}
    

void CurveFitting::planeCurve(const std::vector<double>& parameters,
                    const std::vector<Eigen::Vector2d>& points,
                    NodeType nodeType,
                    SplineType splineType ,
                    BoundaryCondition condition,
                    Eigen::Vector2d start, Eigen::Vector2d end){
    parameters_ = parameters;
    planePoints_ = points;
    nodeType_ = nodeType;
    splineType_ = splineType;
    condition_ = condition;
    is3D_ = false;

    if (condition_ != BoundaryCondition::PERIODIC && 
        condition_ != BoundaryCondition::NATURAL && 
        condition_ != BoundaryCondition::NOT_A_KNOT && start == Eigen::Vector2d::Zero() && end == Eigen::Vector2d::Zero()) {
        throw std::invalid_argument("Missing parameter for COMPLETE or SECOND boundary condition.");
    }

    size_t n = points.size();
    if (parameters_.size() != n) {
        throw std::invalid_argument("The size of parameters and points must be equal.");
    }
    if (n < 4) {
        throw std::invalid_argument("At least four points are required for cubic spline fitting.");
        }

    preprocessParameters();

    // Separate X and Y coordinates
    std::vector<double> x_values(n);
    std::vector<double> y_values(n);
    for (size_t i = 0; i < n; ++i) {
        x_values[i] = planePoints_[i].x();
        y_values[i] = planePoints_[i].y();
    }

    if (splineType_ == SplineType::PPForm) {
        // Use PPForm<3> for fitting
        if(condition_ == BoundaryCondition::COMPLETE){
            ppformX_.setDx(start.x(), end.x());
            ppformY_.setDx(start.y(), end.y());
        }
        else if(condition_ == BoundaryCondition::SECOND){
            ppformX_.setDDx(start.x(), end.x());
            ppformY_.setDDx(start.y(), end.y());
        }
        ppformX_.setInterpolationPoints(parameters_);
        ppformX_.setInterpolationValues(x_values);
        ppformX_.computeCubicCoefficients(condition_);

        ppformY_.setInterpolationPoints(parameters_);
        ppformY_.setInterpolationValues(y_values);
        ppformY_.computeCubicCoefficients(condition_);
    } else if(splineType_ == SplineType::BSpline){
        // Use BSpline<3> for fitting
        if(condition_ == BoundaryCondition::COMPLETE){
            bsplineX_.setDerivatives(start.x(), end.x());
            bsplineY_.setDerivatives(start.y(), end.y());
        }
        else if(condition_ == BoundaryCondition::SECOND){
            bsplineX_.setSecondDerivatives(start.x(), end.x());
            bsplineY_.setSecondDerivatives(start.y(), end.y());
        }
        bsplineX_.setKnots(parameters_);
        bsplineX_.setControlPoints(x_values);
        bsplineX_.setCondition(condition_);
        bsplineX_.computeCoefficients();

        bsplineY_.setKnots(parameters_);
        bsplineY_.setControlPoints(y_values);
        bsplineY_.setCondition(condition_);
        bsplineY_.computeCoefficients();
    }
}

void CurveFitting::sphereCurve(const std::vector<double>& parameters,
                               const std::vector<Eigen::Vector3d>& points,
                               NodeType nodeType,
                               SplineType splineType,
                               BoundaryCondition condition) {
    parameters_ = parameters;
    spherePoints_ = points;
    nodeType_ = nodeType;
    splineType_ = splineType;
    condition_ = condition;
    is3D_ = true;

    if (condition_ != BoundaryCondition::PERIODIC && 
        condition_ != BoundaryCondition::NATURAL && 
        condition_ != BoundaryCondition::NOT_A_KNOT) {
        throw std::invalid_argument("Invalid boundary condition. Must be PERIODIC, NATURAL, or NOT_A_KNOT.");
    }

    size_t n = points.size();
    if (parameters_.size() != n) {
        throw std::invalid_argument("The size of parameters and points must be equal.");
    }
    if (n < 4) {
        throw std::invalid_argument("At least four points are required for cubic spline fitting.");
    }
    if (condition_ == BoundaryCondition::PERIODIC && !points.front().isApprox(points.back(), 1e-6)) {
        throw std::invalid_argument("For periodic splines, the first and last points must be approximately the same.");
    }

    // Calculate a new North Pole
    calculateNorthPole();

    // Project points from sphere to plane using the new North Pole
    planePoints_.resize(n);
    for (size_t i = 0; i < n; ++i) {
        planePoints_[i] = stereographicProjection(spherePoints_[i]);
    }

    // Use planeCurve() to perform fitting
    planeCurve(parameters_, planePoints_, nodeType_, splineType_, condition_, Eigen::Vector2d::Zero(), Eigen::Vector2d::Zero());
}

bool CurveFitting::is3D() const {
    return is3D_;
}

Eigen::Vector2d CurveFitting::plane_evaluate(double t) {
    double x = 0.0;
    double y = 0.0;
    if (splineType_ == SplineType::PPForm) {
        x = ppformX_.cubicSpline(t, condition_);
        y = ppformY_.cubicSpline(t, condition_);
    } else {
        x = bsplineX_.evaluate(t);
        y = bsplineY_.evaluate(t);
    }

    return Eigen::Vector2d(x, y);
}

Eigen::Vector3d CurveFitting::sphere_evaluate(double t) {
    Eigen::Vector2d point2D = plane_evaluate(t);
    Eigen::Vector3d point3D = inverseStereographicProjection(point2D);
    return point3D.normalized();  // Ensure the point lies on the unit sphere
}

std::vector<Eigen::Vector2d> CurveFitting::generatePlotPoints_plane(int numPoints){
    std::vector<Eigen::Vector2d> plotPoints;
    plotPoints.reserve(numPoints);

    // Generate parameter values
    double t_start = parameters_.front();
    double t_end = parameters_.back();

    double step = (t_end - t_start) / (numPoints - 1);

    for (int i = 0; i < numPoints; ++i) {
        double t = t_start + i * step;
        double x = 0.0;
        double y = 0.0;

        if (splineType_ == SplineType::PPForm) {
            x = ppformX_.cubicSpline(t, condition_);
            y = ppformY_.cubicSpline(t, condition_);
        } else {
            x = bsplineX_.evaluate(t);
            y = bsplineY_.evaluate(t);
        }

        plotPoints.emplace_back(x, y);
    }

    return plotPoints;
}

std::vector<Eigen::Vector3d> CurveFitting::generatePlotPoints_sphere(int numPoints){
    std::vector<Eigen::Vector3d> points_sphere;
    std::vector<Eigen::Vector2d> points_plane = generatePlotPoints_plane(numPoints);

    for (const auto& point2D : points_plane) {
        points_sphere.push_back(inverseStereographicProjection(point2D));
    }

    return points_sphere;
}

void CurveFitting::preprocessParameters() {
    size_t n = parameters_.size();
    std::vector<double> parameters(n,0.0);

    if (nodeType_ == NodeType::ChordalLength) {
        // Recompute parameters based on cumulative chordal length
        for (size_t i = 1; i < n; ++i) {
            parameters[i] = parameters[i - 1] + (planePoints_[i] - planePoints_[i - 1]).norm();
        }
        // Normalize parameters to start from initial parameter value
        totalLength_ = parameters[n - 1];
        for (size_t i = 1; i < n; ++i) {
            parameters[i] /= totalLength_;
            parameters_[i] = parameters_.front() + parameters[i]*(parameters_.back()-parameters_.front());
        }
    } else {
        // Ensure parameters are sorted and periodic
        if (!std::is_sorted(parameters_.begin(), parameters_.end())) {
            throw std::invalid_argument("Parameters must be sorted in ascending order for equidistant nodes.");
        }
    }
}


void CurveFitting::calculateNorthPole() {
    // Find a point that is relatively far from the data points
    // For simplicity, select the opposite direction of the mean point
    Eigen::Vector3d mean = Eigen::Vector3d::Zero();
    for (const auto& point : spherePoints_) {
        mean += point;
    }
    mean.normalize();

    // Set the North Pole as the opposite of the mean point
    northPole_ = -mean;

    // Check if northPole_ coincides with any data point
    bool coincides = false;
    for (const auto& point : spherePoints_) {
        if (northPole_.isApprox(point)) {
            coincides = true;
            break;
        }
    }

    // If coincides, select another northPole_
    if (coincides) {
        // Choose a default north pole that is not in the data points
        northPole_ = Eigen::Vector3d(0, 0, 1);

        // Check again
        for (const auto& point : spherePoints_) {
            if (northPole_.isApprox(point)) {
                // Slightly perturb the northPole_
                northPole_ += Eigen::Vector3d(0.01, 0.01, 0.01);
                northPole_.normalize();
                break;
            }
        }
    }
}

Eigen::Vector2d CurveFitting::stereographicProjection(const Eigen::Vector3d& point) const {
    // Rotate the sphere so that the North Pole aligns with (0, 0, 1)
    Eigen::Vector3d k = northPole_;
    Eigen::Vector3d z_axis(0, 0, 1);

    Eigen::Quaterniond q = Eigen::Quaterniond::FromTwoVectors(k, z_axis);
    Eigen::Vector3d rotated_point = q * point;

    // Standard stereographic projection onto z = 0 plane
    double denom = 1.0 - rotated_point.z();
    if (std::abs(denom) < 1e-8) {
        throw std::invalid_argument("Point cannot be projected (z = 1).");
    }
    return Eigen::Vector2d(rotated_point.x() / denom, rotated_point.y() / denom);
}

Eigen::Vector3d CurveFitting::inverseStereographicProjection(const Eigen::Vector2d& point) const {
    // Inverse stereographic projection from plane to sphere
    double x = point.x();
    double y = point.y();
    double denom = 1.0 + x * x + y * y;
    Eigen::Vector3d rotated_point(
        2.0 * x / denom,
        2.0 * y / denom,
        1.0 - 2.0 / denom
    );

    // Rotate back to original orientation
    Eigen::Vector3d z_axis(0, 0, 1);
    Eigen::Quaterniond q = Eigen::Quaterniond::FromTwoVectors(z_axis, northPole_);
    Eigen::Vector3d original_point = q * rotated_point;

    return original_point;
}

#endif // CURVEFITTING_HPP