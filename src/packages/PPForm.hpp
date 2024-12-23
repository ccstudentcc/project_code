// PPForm.hpp
#ifndef PPFORM_HPP
#define PPFORM_HPP

#include <vector>
#include <array>
#include <map>
#include <cmath>
#include <stdexcept>
#include <memory>
#include <limits>
#include <Eigen/Dense>
#include "Function.hpp"
#include "BoundaryCondition.h"


/**
 * @brief Template class for Piecewise Polynomial Form (PPForm).
 * 
 * @tparam Dim Dimension of the spline (1D or 3D).
 */
template <size_t Dim>
class PPForm;

/**
 * @brief Specialization of PPForm for one-dimensional linear splines.
 */
template <>
class PPForm<1>
{
public:
    /**
     * @brief Default constructor. Initializes derivatives to maximum double value for detection.
     */
    PPForm() = default;
    /**
     * @brief Constructs a 1D linear spline with given points and values.
     * 
     * @param points Interpolation points.
     * @param values Function values at the interpolation points.
     * 
     * @throws std::invalid_argument if the number of points is less than two or
     *         if the sizes of points and values do not match.
     */
    PPForm(const std::vector<double> &points, const std::vector<double> &values);

    /**
     * @brief Constructs a 1D linear spline with given points and a function.
     * 
     * @param points Interpolation points.
     * @param f Function to evaluate at interpolation points.
     * 
     * @throws std::invalid_argument if the number of points is less than two.
     */
    PPForm(const std::vector<double> &points, const Function &f);

    /**
     * @brief Sets the interpolation points.
     * 
     * @param points New interpolation points.
     * 
     * @throws std::invalid_argument if points are not strictly increasing.
     */
    void setInterpolationPoints(const std::vector<double> &points);

    /**
     * @brief Sets the interpolation values.
     * 
     * @param values New interpolation values.
     * 
     * @throws std::invalid_argument if the size of values does not match the number of points.
     */
    void setInterpolationValues(const std::vector<double> &values);

    /**
     * @brief Evaluates the linear spline at a single point.
     * 
     * @param t The point at which to evaluate the spline.
     * @return The interpolated value.
     * 
     * @throws std::out_of_range if t is outside the range of interpolation points.
     */
    double linearSpline(double t) const;

    /**
     * @brief Evaluates the linear spline at multiple points.
     * 
     * @param t Vector of points at which to evaluate the spline.
     * @return Vector of interpolated values.
     * 
     * @throws std::out_of_range if any t is outside the range of interpolation points.
     */
    std::vector<double> linearSpline(const std::vector<double> &t) const;

private:
    std::vector<double> points_;     /**< @brief Interpolation points */
    std::vector<double> values_;     /**< @brief Function values at interpolation points */
    std::vector<double> h_;          /**< @brief Intervals between points */
    size_t totalsize_;               /**< @brief Number of interpolation points */

    // Implementation details remain unchanged
};



/**
 * @brief Specialization of PPForm for three-dimensional cubic splines.
 */
template <>
class PPForm<3>
{
public:
    /**
     * @brief Default constructor. Initializes derivatives to maximum double value for detection.
     */
    PPForm();

    /**
     * @brief Constructs a 3D cubic spline with given points and values.
     * 
     * @param points Interpolation points.
     * @param values Function values at the interpolation points.
     * 
     * @throws std::invalid_argument if the number of points is less than two or
     *         if the sizes of points and values do not match.
     */
    PPForm(const std::vector<double> &points, const std::vector<double> &values);

    /**
     * @brief Constructs a 3D cubic spline with given points and a function.
     * 
     * @param points Interpolation points.
     * @param f Function to evaluate at interpolation points.
     * 
     * @throws std::invalid_argument if the number of points is less than two.
     */
    PPForm(const std::vector<double> &points, const Function &f);

    /**
     * @brief Sets the interpolation points.
     * 
     * @param points New interpolation points.
     * 
     * @throws std::invalid_argument if points are not strictly increasing.
     */
    void setInterpolationPoints(const std::vector<double> &points);

    /**
     * @brief Sets the interpolation values.
     * 
     * @param values New interpolation values.
     * 
     * @throws std::invalid_argument if the size of values does not match the number of points.
     */
    void setInterpolationValues(const std::vector<double> &values);

    /**
     * @brief Sets the first derivatives at the boundary points.
     * 
     * @param dx_a First derivative at the first point.
     * @param dx_b First derivative at the last point.
     */
    void setDx(const double &dx_a, const double &dx_b);

    /**
     * @brief Sets the second derivatives at the boundary points.
     * 
     * @param ddx_a Second derivative at the first point.
     * @param ddx_b Second derivative at the last point.
     */
    void setDDx(const double &ddx_a, const double &ddx_b);

    /**
     * @brief Evaluates the cubic spline at a single point with a specified boundary condition.
     * 
     * @param t The point at which to evaluate the spline.
     * @param condition The boundary condition to apply.
     * @return The interpolated value.
     * 
     * @throws std::out_of_range if t is outside the range of interpolation points.
     */
    double cubicSpline(double t, BoundaryCondition condition);

    /**
     * @brief Evaluates the cubic spline at multiple points with a specified boundary condition.
     * 
     * @param t Vector of points at which to evaluate the spline.
     * @param condition The boundary condition to apply.
     * @return Vector of interpolated values.
     * 
     * @throws std::out_of_range if any t is outside the range of interpolation points.
     */
    std::vector<double> cubicSpline(const std::vector<double> &t, BoundaryCondition condition);

    /**
     * @brief Clears the coefficients for a specific boundary condition.
     * 
     * @param condition The boundary condition whose coefficients should be cleared.
     */
    void clearCoefficients(BoundaryCondition condition);

    /**
     * @brief Retrieves the coefficients for a specific boundary condition.
     * 
     * @param condition The boundary condition.
     * @return Vector of coefficients. Returns an empty vector if not found.
     */
    std::vector<double> getCoefficients(BoundaryCondition condition) const;

    /**
     * @brief Computes cubic coefficients based on the boundary condition.
     * 
     * @param condition The boundary condition to apply.
     */
    void computeCubicCoefficients(BoundaryCondition condition);

private:
    std::vector<double> points_;     /**< @brief Interpolation points */
    std::vector<double> values_;     /**< @brief Function values at interpolation points */
    std::vector<double> h_;          /**< @brief Intervals between points */
    size_t totalsize_;               /**< @brief Number of interpolation points */
    std::array<double, 2> coeff_dx_; /**< @brief Boundary first derivatives */
    std::array<double, 2> coeff_ddx_; /**< @brief Boundary second derivatives */

    std::map<BoundaryCondition, std::vector<double>> coefficients_map_; /**< @brief Coefficients for different boundary conditions */

    

    /**
     * @brief Builds the coefficient matrix for the linear system based on the boundary condition.
     * 
     * @param condition The boundary condition to apply.
     * @return Eigen::MatrixXd The coefficient matrix.
     */
    Eigen::MatrixXd buildCoefficientMatrix(BoundaryCondition condition);

    /**
     * @brief Builds the right-hand side vector for the linear system based on the boundary condition.
     * 
     * @param condition The boundary condition to apply.
     * @return Eigen::VectorXd The RHS vector.
     */
    Eigen::VectorXd buildBVector(BoundaryCondition condition);
};



// Implementation of PPForm<1>

// PPForm<1> Constructors

PPForm<1>::PPForm(const std::vector<double> &points, const std::vector<double> &values)
    : points_(points), values_(values), totalsize_(points.size())
{
    if (totalsize_ < 2)
    {
        throw std::invalid_argument("At least two points are required for linear interpolation.");
    }
    if (points.size() != values.size())
    {
        throw std::invalid_argument("The number of interpolation points and values must be the same.");
    }

    h_.resize(totalsize_ - 1);
    for (size_t i = 0; i < totalsize_ - 1; ++i)
    {
        h_[i] = points_[i + 1] - points_[i];
        if (h_[i] <= 0)
        {
            throw std::invalid_argument("Interpolation points must be strictly increasing.");
        }
    }
}

PPForm<1>::PPForm(const std::vector<double> &points, const Function &f)
    : points_(points), totalsize_(points.size())
{
    if (totalsize_ < 2)
    {
        throw std::invalid_argument("At least two points are required for linear interpolation.");
    }

    values_.resize(totalsize_);
    for (size_t i = 0; i < totalsize_; ++i)
    {
        values_[i] = f(points_[i]);
    }

    h_.resize(totalsize_ - 1);
    for (size_t i = 0; i < totalsize_ - 1; ++i)
    {
        h_[i] = points_[i + 1] - points_[i];
        if (h_[i] <= 0)
        {
            throw std::invalid_argument("Interpolation points must be strictly increasing.");
        }
    }
}

// PPForm<1> Member Functions

void PPForm<1>::setInterpolationPoints(const std::vector<double> &points)
{
    points_ = points;
    totalsize_ = points.size();
    h_.resize(totalsize_ - 1);
    for (size_t i = 0; i < totalsize_ - 1; ++i)
    {
        h_[i] = points_[i + 1] - points_[i];
        if (h_[i] <= 0)
        {
            throw std::invalid_argument("Interpolation points must be strictly increasing.");
        }
    }
}

void PPForm<1>::setInterpolationValues(const std::vector<double> &values)
{
    values_ = values;
    if (values_.size() != points_.size())
    {
        throw std::invalid_argument("The number of interpolation values must match the number of points.");
    }
}

double PPForm<1>::linearSpline(double t) const
{
    for (size_t i = 0; i < totalsize_ - 1; ++i)
    {
        if (t >= points_[i] && t <= points_[i + 1])
        {
            double slope = (values_[i + 1] - values_[i]) / h_[i];
            return values_[i] + slope * (t - points_[i]);
        }
    }
    throw std::out_of_range("Input t is out of the interpolation range.");
}

std::vector<double> PPForm<1>::linearSpline(const std::vector<double> &t) const
{
    std::vector<double> value(t.size());
    for (size_t j = 0; j < t.size(); ++j)
    {
        bool found = false;
        for (size_t i = 0; i < totalsize_ - 1; ++i)
        {
            if (t[j] >= points_[i] && t[j] <= points_[i + 1])
            {
                double slope = (values_[i + 1] - values_[i]) / h_[i];
                value[j] = values_[i] + slope * (t[j] - points_[i]);
                found = true;
                break;
            }
        }
        if (!found)
        {
            throw std::out_of_range("Input t is out of the interpolation range.");
        }
    }
    return value;
}



// Implementation of PPForm<3>

// PPForm<3> Constructors

PPForm<3>::PPForm(){
    coeff_dx_ = {std::numeric_limits<double>::max(), std::numeric_limits<double>::max()};
    coeff_ddx_ = {std::numeric_limits<double>::max(), std::numeric_limits<double>::max()};
}

PPForm<3>::PPForm(const std::vector<double> &points, const std::vector<double> &values)
    : points_(points), values_(values), totalsize_(points.size())
{
    if (totalsize_ < 2)
    {
        throw std::invalid_argument("At least two points are required for cubic spline interpolation.");
    }
    if (points.size() != values.size())
    {
        throw std::invalid_argument("The number of interpolation points and values must be the same.");
    }

    h_.resize(totalsize_ - 1);
    for (size_t i = 0; i < totalsize_ - 1; ++i)
    {
        h_[i] = points_[i + 1] - points_[i];
        if (h_[i] <= 0)
        {
            throw std::invalid_argument("Interpolation points must be strictly increasing.");
        }
    }
}

PPForm<3>::PPForm(const std::vector<double> &points, const Function &f)
    : points_(points), totalsize_(points.size())
{
    if (totalsize_ < 2)
    {
        throw std::invalid_argument("At least two points are required for cubic spline interpolation.");
    }

    values_.resize(totalsize_);
    for (size_t i = 0; i < totalsize_; ++i)
    {
        values_[i] = f(points_[i]);
    }

    h_.resize(totalsize_ - 1);
    for (size_t i = 0; i < totalsize_ - 1; ++i)
    {
        h_[i] = points_[i + 1] - points_[i];
        if (h_[i] <= 0)
        {
            throw std::invalid_argument("Interpolation points must be strictly increasing.");
        }
    }

    // Initialize boundary derivatives
    coeff_dx_ = {f.nth_derivative(points[0], 1), f.nth_derivative(points[totalsize_ - 1], 1)};
    coeff_ddx_ = {f.nth_derivative(points[0], 2), f.nth_derivative(points[totalsize_ - 1], 2)};

    // Compute coefficients for different boundary conditions
    computeCubicCoefficients(BoundaryCondition::NATURAL);
    computeCubicCoefficients(BoundaryCondition::COMPLETE);
    computeCubicCoefficients(BoundaryCondition::NOT_A_KNOT);
    computeCubicCoefficients(BoundaryCondition::SECOND);
}

// PPForm<3> Member Functions

void PPForm<3>::setInterpolationPoints(const std::vector<double> &points)
{
    points_ = points;
    totalsize_ = points.size();
    h_.resize(totalsize_ - 1);
    for (size_t i = 0; i < totalsize_ - 1; ++i)
    {
        h_[i] = points_[i + 1] - points_[i];
        if (h_[i] <= 0)
        {
            throw std::invalid_argument("Interpolation points must be strictly increasing.");
        }
    }
    coefficients_map_.clear();
}

void PPForm<3>::setInterpolationValues(const std::vector<double> &values)
{
    values_ = values;
    if (values_.size() != points_.size())
    {
        throw std::invalid_argument("The number of interpolation values must match the number of points.");
    }
    coefficients_map_.clear();
}

void PPForm<3>::setDx(const double &dx_a, const double &dx_b)
{
    coeff_dx_[0] = dx_a;
    coeff_dx_[1] = dx_b;
    coefficients_map_.clear();
}

void PPForm<3>::setDDx(const double &ddx_a, const double &ddx_b)
{
    coeff_ddx_[0] = ddx_a;
    coeff_ddx_[1] = ddx_b;
    coefficients_map_.clear();
}

double PPForm<3>::cubicSpline(double t, BoundaryCondition condition)
{
    if (coefficients_map_[condition].empty())
    {
        computeCubicCoefficients(condition);
    }

    if(t < points_[totalsize_-1]+1e-8 && t >= points_[totalsize_-1]){
        return values_[totalsize_-1];
    }
    else if(t > points_[0]-1e-8 && t <= points_[0]){
        return values_[0];
    }

    for (size_t i = 0; i < totalsize_ - 1; ++i)
    {
        if (t > points_[i] && t <= points_[i + 1])
        {
            double x = t - points_[i];
            const std::vector<double> &coeffs = coefficients_map_[condition];
            double a = coeffs[4 * i];
            double b = coeffs[4 * i + 1];
            double c = coeffs[4 * i + 2];
            double d = coeffs[4 * i + 3];
            return a + b * x + c * x * x + d * x * x * x;
        }
    }
    throw std::out_of_range("Input t is out of the interpolation range.");
}

std::vector<double> PPForm<3>::cubicSpline(const std::vector<double> &t, BoundaryCondition condition)
{
    std::vector<double> values(t.size());
    for (size_t j = 0; j < t.size(); ++j)
    {
        values[j]=cubicSpline(t[j], condition);
    }
    return values;
}

void PPForm<3>::clearCoefficients(BoundaryCondition condition)
{
    auto it = coefficients_map_.find(condition);
    if (it != coefficients_map_.end())
    {
        it->second.clear();
        coefficients_map_.erase(it);
    }
}

std::vector<double> PPForm<3>::getCoefficients(BoundaryCondition condition) const
{
    auto it = coefficients_map_.find(condition);
    if (it != coefficients_map_.end())
    {
        return it->second;
    }
    return {}; // Return empty vector if not found
}

void PPForm<3>::computeCubicCoefficients(BoundaryCondition condition)
{
    Eigen::MatrixXd A = buildCoefficientMatrix(condition);
    Eigen::VectorXd b = buildBVector(condition);

    // Solve the linear system
    Eigen::VectorXd c = A.colPivHouseholderQr().solve(b);

    if (condition == BoundaryCondition::PERIODIC)
    {
        c.conservativeResize(totalsize_);
        c(totalsize_ - 1) = c(0);
    }

    // Calculate polynomial coefficients a_i, b_i, c_i, d_i
    std::vector<double> coefficients(4 * (totalsize_ - 1));
    for (size_t i = 0; i < totalsize_ - 1; ++i)
    {
        double a = values_[i];
        double b_coeff = (values_[i + 1] - values_[i]) / h_[i] - h_[i] * (2 * c(i) + c(i + 1)) / 3.0;
        double d = (c(i + 1) - c(i)) / (3 * h_[i]);
        coefficients[4 * i] = a;
        coefficients[4 * i + 1] = b_coeff;
        coefficients[4 * i + 2] = c(i);
        coefficients[4 * i + 3] = d;
    }

    coefficients_map_[condition] = coefficients;
}

Eigen::MatrixXd PPForm<3>::buildCoefficientMatrix(BoundaryCondition condition)
{
    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(totalsize_, totalsize_);

    // Internal nodes
    for (size_t i = 1; i < totalsize_ - 1; ++i)
    {
        A(i, i - 1) = h_[i - 1];
        A(i, i) = 2 * (h_[i - 1] + h_[i]);
        A(i, i + 1) = h_[i];
    }

    // Boundary conditions
    switch (condition)
    {
    case BoundaryCondition::NATURAL:
        A(0, 0) = 1.0;
        A(totalsize_ - 1, totalsize_ - 1) = 1.0;
        break;
    case BoundaryCondition::COMPLETE:
        if(coeff_dx_[0] == std::numeric_limits<double>::max() || coeff_dx_[1] == std::numeric_limits<double>::max()){
            throw std::invalid_argument("Boundary first derivatives must be set for COMPLETE boundary conditions.");
        }
        A(0, 0) = 2 * h_[0];
        A(0, 1) = h_[0];
        A(totalsize_ - 1, totalsize_ - 2) = h_[totalsize_ - 2];
        A(totalsize_ - 1, totalsize_ - 1) = 2 * h_[totalsize_ - 2];
        break;
    case BoundaryCondition::PERIODIC:
        if (std::abs(values_[0] - values_[totalsize_ - 1]) > 1e-8)
        {
            throw std::invalid_argument("For periodic boundary conditions, the first and last values must be the same.");
        }
        A.conservativeResize(totalsize_ - 1, totalsize_ - 1);
        A(0, 0) = 2 * (h_[0] + h_[totalsize_ - 2]);
        A(0, 1) = h_[0];
        A(0, totalsize_ - 2) = h_[totalsize_ - 2];
        for (size_t i = 1; i < totalsize_ - 2; ++i)
        {
            A(i, i - 1) = h_[i];
            A(i, i) = 2 * (h_[i] + h_[i + 1]);
            A(i, i + 1) = h_[i + 1];
        }
        A(totalsize_ - 2, totalsize_ - 3) = h_[totalsize_ - 3];
        A(totalsize_ - 2, totalsize_ - 2) = 2 * (h_[totalsize_ - 3] + h_[totalsize_ - 2]);
        A(totalsize_ - 2, 0) = h_[totalsize_ - 2];
        break;
    case BoundaryCondition::NOT_A_KNOT:
        A(0, 0) = -h_[1];
        A(0, 1) = h_[0] + h_[1];
        A(0, 2) = -h_[0];
        A(totalsize_ - 1, totalsize_ - 3) = -h_[totalsize_ - 2];
        A(totalsize_ - 1, totalsize_ - 2) = h_[totalsize_ - 3] + h_[totalsize_ - 2];
        A(totalsize_ - 1, totalsize_ - 1) = -h_[totalsize_ - 3];
        break;
    case BoundaryCondition::SECOND:
        if(coeff_ddx_[0] == std::numeric_limits<double>::max() || coeff_ddx_[1] == std::numeric_limits<double>::max()){
            throw std::invalid_argument("Boundary second derivatives must be set for SECOND boundary conditions.");
        }
        A(0, 0) = 2.0;
        A(totalsize_ - 1, totalsize_ - 1) = 2.0;
        break;
    default:
        throw std::invalid_argument("Unsupported boundary condition.");
    }

    return A;
}

Eigen::VectorXd PPForm<3>::buildBVector(BoundaryCondition condition)
{
    Eigen::VectorXd b = Eigen::VectorXd::Zero(totalsize_);

    for (size_t i = 1; i < totalsize_ - 1; ++i)
    {
        b(i) = 3 * ((values_[i + 1] - values_[i]) / h_[i] - (values_[i] - values_[i - 1]) / h_[i - 1]);
    }

    switch (condition)
    {
    case BoundaryCondition::NATURAL:
        b(0) = 0.0;
        b(totalsize_ - 1) = 0.0;
        break;
    case BoundaryCondition::COMPLETE:
        if(coeff_dx_[0] == std::numeric_limits<double>::max() || coeff_dx_[1] == std::numeric_limits<double>::max()){
            throw std::invalid_argument("Boundary first derivatives must be set for COMPLETE boundary conditions.");
        }
        b(0) = 3 * ((values_[1] - values_[0]) / h_[0] - coeff_dx_[0]);
        b(totalsize_ - 1) = 3 * (coeff_dx_[1] - (values_[totalsize_ - 1] - values_[totalsize_ - 2]) / h_[totalsize_ - 2]);
        break;
    case BoundaryCondition::PERIODIC:
        b.conservativeResize(totalsize_ - 1);
        b(0) = 3 * ((values_[1] - values_[0]) / h_[0] - (values_[0] - values_[totalsize_ - 2]) / h_[totalsize_ - 2]);
        for (size_t i = 1; i < totalsize_ - 2; ++i)
        {
            b(i) = 3 * ((values_[i + 1] - values_[i]) / h_[i] - (values_[i] - values_[i - 1]) / h_[i - 1]);
        }
        b(totalsize_ - 2) = 3 * ((values_[0] - values_[totalsize_ - 2]) / h_[totalsize_ - 2] - (values_[totalsize_ - 2] - values_[totalsize_ - 3]) / h_[totalsize_ - 3]);
        break;
    case BoundaryCondition::NOT_A_KNOT:
        b(0) = 0.0;
        b(totalsize_ - 1) = 0.0;
        break;
    case BoundaryCondition::SECOND:
        if(coeff_ddx_[0] == std::numeric_limits<double>::max() || coeff_ddx_[1] == std::numeric_limits<double>::max()){
            throw std::invalid_argument("Boundary second derivatives must be set for SECOND boundary conditions.");
        }
        b(0) = coeff_ddx_[0];
        b(totalsize_ - 1) = coeff_ddx_[1];
        break;
    default:
        throw std::invalid_argument("Unsupported boundary condition.");
    }

    return b;
}

#endif // PPFORM_HPP