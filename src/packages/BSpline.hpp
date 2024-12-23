// BSpline.hpp
#ifndef BSPLINE_HPP
#define BSPLINE_HPP

#include <vector>
#include <Eigen/Dense>
#include <stdexcept>
#include <cmath>
#include <limits>
#include <string>
#include "BoundaryCondition.h"

constexpr double EPSILON = 1e-8;

/**
 * @class BSpline
 * @brief A template class for generating B-splines of a specified order.
 *
 * This class provides functionality to create and evaluate B-splines of varying orders.
 * For orders that require specific boundary conditions (e.g., cubic B-splines),
 * specialized implementations handle the computation of coefficients.
 *
 * @tparam Dim The order of the B-spline (Dim=1 corresponds to linear B-spline, Dim=3 to cubic B-spline).
 */
template <size_t Dim>
class BSpline
{
public:
    /**
     * @brief Default constructor. Initializes derivatives to maximum double value for detection.
     */
    BSpline();

    /**
     * @brief Constructs a BSpline with given knots and control points.
     *
     * @param knots A sorted vector of knot values.
     * @param controlPoints A vector of control points.
     * @throws std::invalid_argument if the number of knots is insufficient for the specified Dim.
     */
    BSpline(const std::vector<double> &knots, const std::vector<double> &controlPoints);

    /**
     * @brief Constructs a BSpline with given knots, control points, and boundary condition.
     *
     * @param knots A sorted vector of knot values.
     * @param controlPoints A vector of control points.
     * @param condition The boundary condition to be applied.
     * @throws std::invalid_argument if the knot vector is not sorted or insufficient.
     */
    BSpline(const std::vector<double> &knots, const std::vector<double> &controlPoints, BoundaryCondition condition);

    /**
     * @brief Constructs a BSpline with knots, boundary condition, and specified derivatives.
     *
     * This constructor is specifically for boundary conditions that require derivative specifications.
     *
     * @param knots A sorted vector of knot values.
     * @param controlPoints A vector of control points.
     * @param condition The boundary condition to be applied (COMPLETE or SECOND).
     * @param startParameter The derivative value at the start.
     * @param endParameter The derivative value at the end.
     * @throws std::invalid_argument if derivatives are not set or boundary condition is incompatible.
     */
    BSpline(const std::vector<double> &knots, const std::vector<double> &controlPoints, BoundaryCondition condition, double startParameter, double endParameter);

    /**
     * @brief Sets the knot vector and recomputes the coefficients.
     *
     * @param knots A new sorted vector of knot values.
     * @throws std::invalid_argument if the knot vector is not sorted or insufficient.
     */
    void setKnots(const std::vector<double> &knots);

    /**
     * @brief Sets the control points and recomputes the coefficients.
     *
     * @param controlPoints A new vector of control points.
     * @throws std::invalid_argument if the knot vector and control points sizes are incompatible.
     */
    void setControlPoints(const std::vector<double> &controlPoints);

    /**
     * @brief Sets the boundary condition and recomputes the coefficients.
     *
     * @param condition The boundary condition to be applied.
     * @throws std::invalid_argument if the boundary condition is incompatible with the current setup.
     */
    void setCoefficients(const std::vector<double> &coefficients);

    /**
     * @brief Sets the boundary condition and recomputes the coefficients.
     *
     * @param condition The boundary condition to be applied.
     * @throws std::invalid_argument if the boundary condition is incompatible with the current setup.
     */
    void setCondition(BoundaryCondition condition);

    /**
     * @brief Sets the derivatives at the start and end points for COMPLETE boundary condition and recomputes coefficients.
     *
     * @param startDerivative The derivative at the start.
     * @param endDerivative The derivative at the end.
     * @throws std::invalid_argument if the boundary condition is not COMPLETE.
     */
    void setDerivatives(double startDerivative, double endDerivative);

    /**
     * @brief Sets the second derivatives at the start and end points for SECOND boundary condition and recomputes coefficients.
     *
     * @param startSecondDerivative The second derivative at the start.
     * @param endSecondDerivative The second derivative at the end.
     * @throws std::invalid_argument if the boundary condition is not SECOND.
     */
    void setSecondDerivatives(double startSecondDerivative, double endSecondDerivative);

    /**
     * @brief Evaluates the B-spline at a given parameter value.
     *
     * @param t The parameter value at which to evaluate the B-spline.
     * @return The value of the B-spline at parameter t.
     * @throws std::out_of_range if t is outside the knot vector range.
     */
    double evaluate(double t) const;

    /**
     * @brief Evaluates the B-spline at a given parameter value in every dimension by inputted coefficients and knots.
     *
     * @param t The parameter value at which to evaluate the B-spline.
     * @return The value of the B-spline at parameter t.
     * @throws std::out_of_range if t is outside the knot vector range.
     */
    double value(double t) const;

    /**
     * @brief Generates a set of points on the B-spline curve for plotting.
     *
     * @param numPoints The number of points to generate along the curve.
     * @return A vector of pairs, where each pair contains a parameter value and the corresponding B-spline value.
     * @throws std::invalid_argument if numPoints is non-positive.
     */
    std::vector<std::pair<double, double>> generatePlotPoints(int numPoints) const;

    /**
     * @brief Computes the coefficients of the B-spline based on the current boundary condition.
     *
     * This function is specialized for specific Dims that require particular boundary condition handling.
     *
     * @throws std::runtime_error if coefficient computation fails.
     */
    void computeCoefficients();

protected:
    std::vector<double> originalKnots_;   /**< Original knot vector (user input) */
    std::vector<double> clampedKnots_;    /**< Clamped knot vector with appropriate multiplicities */
    std::vector<double> controlPoints_;   /**< Control points vector */
    std::vector<double> coefficients_;    /**< Coefficients vector */
    BoundaryCondition boundaryCondition_; /**< Current boundary condition */
    double startDerivative_;              /**< Start derivative */
    double endDerivative_;                /**< End derivative */
    double startSecondDerivative_;        /**< Start second derivative */
    double endSecondDerivative_;          /**< End second derivative */

    /**
     * @brief Recursively computes the value of the basis function.
     *
     * @param i The knot interval index.
     * @param k The order of the basis function.
     * @param t The parameter value at which to evaluate the basis function.
     * @return The value of the basis function N_i,k(t).
     */
    double basisFunction(int i, int k, double t) const;

    /**
     * @brief Recursively computes the value of the derivative of the basis function.
     *
     * @param i The knot interval index.
     * @param k The order of the basis function.
     * @param t The parameter value at which to evaluate the basis function.
     * @return The value of the basis function N'_i,k(t).
     */
    double basisFunction_derivative(int i, int k, double x) const;

    /**
     * @brief Recursively computes the value of the second derivative of the basis function.
     *
     * @param i The knot interval index.
     * @param k The order of the basis function.
     * @param t The parameter value at which to evaluate the basis function.
     * @return The value of the basis function N''_i,k(t).
     */
    double basisFunction_second_derivative(int i, int k, double x) const;

    /**
     * @brief Recursively computes the value of the third derivative of the basis function.
     *
     * @param i The knot interval index.
     * @param k The order of the basis function.
     * @param t The parameter value at which to evaluate the basis function.
     * @return The value of the basis function N'''_i,k(t).
     */
    double basisFunction_third_derivative(int i, int k, double x) const;

    /**
     * @brief Computes the value of the cardinal B-spline basis function.
     *
     * @param i The knot interval index.
     * @param degree The degree of the cardinal B-spline.
     * @param t The parameter value at which to evaluate the basis function.
     * @return The value of the cardinal B-spline basis function B_{i,\mathbb{Z}}^n(t).
     */
    double cardinalFunction(int i, int degree, double t) const;

private:
    /**
     * @brief Generates clamped knots by adding knots at the ends.
     *
     * @param knots The original knot vector.
     * @return The clamped knot vector with appropriate knots.
     */
    std::vector<double> generateClampedKnots(const std::vector<double> &knots) const;

    /**
     * @brief Validates and clamps the knot vector.
     *
     * Ensures that the knot vector is sorted and has sufficient multiplicity.
     *
     * @param knots The original knot vector.
     * @return The validated and clamped knot vector.
     * @throws std::invalid_argument if the knot vector is not sorted or has insufficient length.
     */
    std::vector<double> validateAndClampKnots(const std::vector<double> &knots) const;
};

// Implementation

// Default constructor
template <size_t Dim>
BSpline<Dim>::BSpline()
    : originalKnots_(),
      clampedKnots_(),
      controlPoints_(),
      coefficients_(),
      boundaryCondition_(BoundaryCondition::NATURAL),
      startDerivative_(std::numeric_limits<double>::max()),
      endDerivative_(std::numeric_limits<double>::max()),
      startSecondDerivative_(std::numeric_limits<double>::max()),
      endSecondDerivative_(std::numeric_limits<double>::max())
{
}

// Constructor with knots and control points
template <size_t Dim>
BSpline<Dim>::BSpline(const std::vector<double> &knots, const std::vector<double> &controlPoints)
    : originalKnots_(knots),
      clampedKnots_(validateAndClampKnots(knots)),
      controlPoints_(controlPoints),
      coefficients_(),
      boundaryCondition_(BoundaryCondition::NATURAL),
      startDerivative_(std::numeric_limits<double>::max()),
      endDerivative_(std::numeric_limits<double>::max()),
      startSecondDerivative_(std::numeric_limits<double>::max()),
      endSecondDerivative_(std::numeric_limits<double>::max())
{
    // Check if there are enough knots for the given dimension
    if (originalKnots_.size() != controlPoints_.size())
    {
        throw std::invalid_argument("Insufficient number of knots and control points for the specified Dim.");
    }
    computeCoefficients();
}

// Constructor with knots, control points, and boundary condition
template <size_t Dim>
BSpline<Dim>::BSpline(const std::vector<double> &knots, const std::vector<double> &controlPoints, BoundaryCondition condition)
    : originalKnots_(knots),
      clampedKnots_(validateAndClampKnots(knots)),
      controlPoints_(controlPoints),
      coefficients_(),
      boundaryCondition_(condition),
      startDerivative_(std::numeric_limits<double>::max()),
      endDerivative_(std::numeric_limits<double>::max()),
      startSecondDerivative_(std::numeric_limits<double>::max()),
      endSecondDerivative_(std::numeric_limits<double>::max())
{
    // Check if there are enough knots for the given dimension
    if (originalKnots_.size() != controlPoints_.size())
    {
        throw std::invalid_argument("Insufficient number of knots and control points for the specified Dim.");
    }
    computeCoefficients();
}

// Constructor with knots, control points, boundary condition, and derivatives
template <size_t Dim>
BSpline<Dim>::BSpline(const std::vector<double> &knots, const std::vector<double> &controlPoints, BoundaryCondition condition, double startParameter, double endParameter)
    : originalKnots_(knots),
      clampedKnots_(validateAndClampKnots(knots)),
      controlPoints_(controlPoints),
      coefficients_(),
      boundaryCondition_(condition),
      startDerivative_(std::numeric_limits<double>::max()),
      endDerivative_(std::numeric_limits<double>::max()),
      startSecondDerivative_(std::numeric_limits<double>::max()),
      endSecondDerivative_(std::numeric_limits<double>::max())
{
    if (originalKnots_.size() != controlPoints_.size())
    {
        throw std::invalid_argument("Insufficient number of knots and control points for the specified Dim.");
    }
    if (condition == BoundaryCondition::SECOND)
    {
        startSecondDerivative_ = startParameter;
        endSecondDerivative_ = endParameter;
    }
    else if (condition == BoundaryCondition::COMPLETE)
    {
        startDerivative_ = startParameter;
        endDerivative_ = endParameter;
    }
    else
    {
        throw std::invalid_argument("Derivatives can only be set for COMPLETE or SECOND boundary conditions.");
    }

    computeCoefficients();
}

// Setters
template <size_t Dim>
void BSpline<Dim>::setKnots(const std::vector<double> &knots)
{
    originalKnots_ = knots;
    clampedKnots_ = validateAndClampKnots(knots);
}

template <size_t Dim>
void BSpline<Dim>::setControlPoints(const std::vector<double> &controlPoints)
{
    controlPoints_ = controlPoints;
}

template <size_t Dim>
void BSpline<Dim>::setCondition(BoundaryCondition condition)
{
    boundaryCondition_ = condition;
}

template <size_t Dim>
void BSpline<Dim>::setCoefficients(const std::vector<double> &coefficients)
{
    coefficients_ = coefficients;
}

template <size_t Dim>
void BSpline<Dim>::setDerivatives(double startDerivative, double endDerivative)
{
    startDerivative_ = startDerivative;
    endDerivative_ = endDerivative;
}

template <size_t Dim>
void BSpline<Dim>::setSecondDerivatives(double startSecondDerivative, double endSecondDerivative)
{
    startSecondDerivative_ = startSecondDerivative;
    endSecondDerivative_ = endSecondDerivative;
}

// Evaluate method
template <size_t Dim>
double BSpline<Dim>::evaluate(double t) const
{
    if (clampedKnots_.empty())
    {
        throw std::runtime_error("Knots have not been set.");
    }

    // Check if t is within the knot vector range
    if (t <= originalKnots_.front() - EPSILON || t >= originalKnots_.back() + EPSILON)
    {
        throw std::out_of_range("Parameter t is outside the knot vector range.");
    }
    else if (t < originalKnots_.back() + EPSILON && t >= originalKnots_.back())
    {
        return controlPoints_.back();
    }
    else if (t > originalKnots_.front() - EPSILON && t <= originalKnots_.front())
    {
        return controlPoints_.front();
    }

    double result = 0.0;
    for (size_t i = 0; i < coefficients_.size() + Dim - 1; ++i)
    {
        result += coefficients_[i] * basisFunction(i - Dim + 2, Dim, t);
    }

    return result;
}

/**
 * @brief Evaluates a quadratic (degree 2) B-spline at a given parameter value.
 *
 * This specialization for quadratic B-splines uses the cardinal basis function
 * to compute the spline value at the specified parameter.
 *
 * @param t The parameter value at which to evaluate the B-spline.
 * @return The value of the B-spline at parameter t.
 * @throws std::out_of_range if t is outside the knot vector range.
 */
template <>
double BSpline<2>::evaluate(double t) const
{
    int end = originalKnots_.back();
    int start = originalKnots_.front();
    int j = 0;

    // Compute the spline value
    double value = 0.0;

    for (int i = start - 1; i <= end; i++)
    {
        value += coefficients_[j] * cardinalFunction(i, 2, t);
        j++;
    }

    return value;
}

/**
 * @brief Evaluates the cubic B-spline at a given parameter value using B-spline basis functions.
 *
 * This implementation uses the Cox-de Boor recursion formula to calculate the B-spline value
 * as a linear combination of basis functions and control coefficients.
 *
 * @param t The parameter value at which to evaluate the B-spline
 * @return The value of the B-spline at parameter t
 * @throws std::out_of_range If parameter t is outside the knot vector range
 * @see basisFunction
 */
template <>
double BSpline<3>::evaluate(double t) const
{
    // Check if t is within the knot vector range
    if (t <= originalKnots_.front() - EPSILON || t >= originalKnots_.back() + EPSILON)
    {
        throw std::out_of_range("Parameter t is outside the knot vector range.");
    }
    else if (t < originalKnots_.back() + EPSILON && t >= originalKnots_.back())
    {
        return controlPoints_.back();
    }
    else if (t > originalKnots_.front() - EPSILON && t <= originalKnots_.front())
    {
        return controlPoints_.front();
    }

    // Calculate spline value using basis functions
    double result = 0.0;
    for (int i = 0; i < static_cast<int>(originalKnots_.size() + 2); ++i)
    {
        result += coefficients_[i] * basisFunction(i - 1, 3, t);
    }

    return result;
}

template <size_t Dim>
double BSpline<Dim>::value(double t) const
{  
    if(coefficients_.size() != originalKnots_.size() + Dim-1)
    {
        throw std::runtime_error("Coefficients have not been set.");
    }
    if (clampedKnots_.empty())
    {
        throw std::runtime_error("Knots have not been set.");
    }

    // Check if t is within the knot vector range
    if (t <= originalKnots_.front() - EPSILON || t >= originalKnots_.back() + EPSILON)
    {
        throw std::out_of_range("Parameter t is outside the knot vector range.");
    }
    else if (t < originalKnots_.back() + EPSILON && t >= originalKnots_.back())
    {
        t=originalKnots_.back();
    }
    else if (t > originalKnots_.front() - EPSILON && t <= originalKnots_.front())
    {
        t=originalKnots_.front();
    }

    double result = 0.0;
    for (size_t i = 0; i < coefficients_.size() + Dim - 1; ++i)
    {
        result += coefficients_[i] * basisFunction(i - Dim + 2, Dim, t);
    }

    return result;
}

// Generate plot points and return both parameter values and corresponding spline values
template <size_t Dim>
std::vector<std::pair<double, double>> BSpline<Dim>::generatePlotPoints(int numPoints) const
{
    if (clampedKnots_.empty())
    {
        throw std::runtime_error("Knots have not been set.");
    }

    if (numPoints <= 0)
    {
        throw std::invalid_argument("Number of plot points must be positive.");
    }

    std::vector<std::pair<double, double>> plotPoints;
    double start = originalKnots_.front();
    double end = originalKnots_.back();
    double step = (end - start) / (numPoints - 1);

    for (int i = 0; i < numPoints; ++i)
    {
        double t = start + (double)i * step;
        if (controlPoints_.empty())
        {
            plotPoints.emplace_back(t, value(t));
        }
        else
        {
            plotPoints.emplace_back(t, evaluate(t));
        }
    }

    return plotPoints;
}

// Basis function implementation (Cox-de Boor recursion formula)
template <size_t Dim>
double BSpline<Dim>::basisFunction(int i, int k, double t) const
{
    if (k < 0)
    {
        throw std::invalid_argument("Order k cannot be negative.");
    }
    if (t <= clampedKnots_[i - 1 + Dim - 1] || t >= clampedKnots_[i + Dim + Dim - 1])
    {
        return 0.0;
    }
    if (k == 0)
    {
        if (clampedKnots_[i - 1 + Dim - 1] < t && t <= clampedKnots_[i + Dim - 1])
            return 1.0;
        else
            return 0.0;
    }
    else
    {
        double denom1 = clampedKnots_[i + k - 1 + Dim - 1] - clampedKnots_[i - 1 + Dim - 1];
        double term1 = 0.0;
        term1 = (t - clampedKnots_[i - 1 + Dim - 1]) / denom1 * basisFunction(i, k - 1, t);

        double denom2 = clampedKnots_[i + k + Dim - 1] - clampedKnots_[i + Dim - 1];
        double term2 = 0.0;
        term2 = (clampedKnots_[i + k + Dim - 1] - t) / denom2 * basisFunction(i + 1, k - 1, t);

        return term1 + term2;
    }
}

template <size_t Dim>
double BSpline<Dim>::basisFunction_derivative(int i, int k, double x) const
{
    if (k == 0 || x <= clampedKnots_[i - 1 + Dim - 1] || x >= clampedKnots_[i + Dim + Dim - 1])
    {
        return 0.0;
    }
    if (k == 1)
    {
        // Base case: when k = 1, return the definition of the first derivative
        if (x > clampedKnots_[i - 1 + Dim - 1] && x <= clampedKnots_[i + Dim - 1])
        {
            return 1.0 / (clampedKnots_[i + Dim - 1] - clampedKnots_[i - 1 + Dim - 1]);
        }
        else if (x > clampedKnots_[i + Dim - 1] && x <= clampedKnots_[i + 1 + Dim - 1])
        {
            return -1.0 / (clampedKnots_[i + Dim] - clampedKnots_[i + Dim - 1]);
        }
        else
        {
            return 0.0;
        }
    }
    if (k < 0)
    {
        throw std::invalid_argument("Order k cannot be negative.");
    }
    return k * basisFunction(i, k - 1, x) / (clampedKnots_[i + k - 1 + Dim - 1] - clampedKnots_[i - 1 + Dim - 1]) - k * basisFunction(i + 1, k - 1, x) / (clampedKnots_[i + k + Dim - 1] - clampedKnots_[i + Dim - 1]);
}

template <size_t Dim>
double BSpline<Dim>::basisFunction_second_derivative(int i, int k, double x) const
{
    if (k < 0)
    {
        throw std::invalid_argument("Order k cannot be negative.");
    }
    if (k <= 1)
    {
        return 0.0;
    }
    return k * basisFunction_derivative(i, k - 1, x) / (clampedKnots_[i + k - 1 + Dim - 1] - clampedKnots_[i - 1 + Dim - 1]) - k * basisFunction_derivative(i + 1, k - 1, x) / (clampedKnots_[i + k + Dim - 1] - clampedKnots_[i + Dim - 1]);
}

template <size_t Dim>
double BSpline<Dim>::basisFunction_third_derivative(int i, int k, double x) const
{
    if (k < 0)
    {
        throw std::invalid_argument("Order k cannot be negative.");
    }
    if (k <= 2)
    {
        return 0.0;
    }
    return k * basisFunction_second_derivative(i, k - 1, x) / (clampedKnots_[i + k - 1 + Dim - 1] - clampedKnots_[i - 1 + Dim - 1]) - k * basisFunction_second_derivative(i + 1, k - 1, x) / (clampedKnots_[i + k + Dim - 1] - clampedKnots_[i + Dim - 1]);
}

template <size_t Dim>
double BSpline<Dim>::cardinalFunction(int i, int degree, double t) const
{
    double result = 0.0;
    for (int k = -1; k <= degree; ++k)
    {
        double binomialCoeff = std::tgamma(degree + 2) / (std::tgamma(k + 2) * std::tgamma(degree - k + 1));
        double term = std::pow(std::max(0.0, k + i - t), degree);
        result += std::pow(-1, degree - k) * binomialCoeff * term;
    }
    return result / std::tgamma(degree + 1);
}

// Generate clamped knots by adding multiplicities at the ends
template <size_t Dim>
std::vector<double> BSpline<Dim>::generateClampedKnots(const std::vector<double> &knots) const
{
    std::vector<double> clamped;
    // Add Dim multiplicities at the start
    for (size_t i = 0; i < Dim; ++i)
    {
        clamped.push_back(knots.front() - Dim + i);
    }
    // Add interior knots
    for (size_t i = 0; i < knots.size(); ++i)
    {
        clamped.push_back(knots[i]);
    }
    // Add Dim multiplicities at the end
    for (size_t i = 0; i < Dim; ++i)
    {
        clamped.push_back(knots.back() + i + 1);
    }
    return clamped;
}

// Validate and clamp knots
template <size_t Dim>
std::vector<double> BSpline<Dim>::validateAndClampKnots(const std::vector<double> &knots) const
{
    if (knots.size() < 2)
    {
        throw std::invalid_argument("Knot vector must contain at least two knots.");
    }

    // Check if knots are sorted
    for (size_t i = 1; i < knots.size(); ++i)
    {
        if (knots[i] < knots[i - 1])
        {
            throw std::invalid_argument("Knot vector must be sorted in non-decreasing order.");
        }
    }

    // Generate clamped knots
    return generateClampedKnots(knots);
}

/**
 * @brief Specialization of computeCoefficients for Dim = 1 (Linear B-spline).
 *
 * For linear B-splines, the coefficients are directly set from the control points.
 * Boundary conditions are not applicable for linear B-splines.
 */
template <>
void BSpline<1>::computeCoefficients()
{
    // For linear B-spline, coefficients are the same as control points
    coefficients_ = controlPoints_;
}

// Specialization for Degree = 2 (Quadratic B-spline)
template <>
void BSpline<2>::computeCoefficients()
{
    size_t n = controlPoints_.size();
    if (std::floor(originalKnots_.front()) != originalKnots_.front())
    {
        throw std::invalid_argument("The first element of the original knot vector must be an integer.");
    }
    if (n < 3)
        throw std::invalid_argument("At least 3 control points are required for quadratic B-spline.");

    std::vector<double> h(n - 1);
    for (size_t i = 0; i < n - 1; ++i)
    {
        h[i] = originalKnots_[i + 1] - originalKnots_[i];
        if ((i == 0 || i == n - 2) && h[i] != 0.5)
        {
            throw std::invalid_argument("Knot vector must have an interval of 0.5 at the start and end.");
        }
        else if (i != 0 && i != n - 2 && h[i] != 1)
        {
            throw std::invalid_argument("Knot vector must have an interval of 1 between interior knots.");
        }
    }

    coefficients_.resize(n);
    // Construct matrix M
    Eigen::MatrixXd M = Eigen::MatrixXd::Zero(n - 2, n - 2);
    M(0, 0) = 5;
    M(0, 1) = 1;
    for (size_t i = 1; i < n - 3; ++i)
    {
        M(i, i - 1) = 1;
        M(i, i) = 6;
        M(i, i + 1) = 1;
    }
    M(n - 3, n - 4) = 1;
    M(n - 3, n - 3) = 5;

    // Construct vector b
    Eigen::VectorXd b(n - 2);
    b(0) = 8 * controlPoints_[1] - 2 * controlPoints_[0];
    for (size_t i = 1; i < n - 3; ++i)
    {
        b(i) = 8 * controlPoints_[i + 1];
    }
    b(n - 3) = 8 * controlPoints_[n - 2] - 2 * controlPoints_[n - 1];

    // Solve the linear system Ma = b
    Eigen::VectorXd a = M.colPivHouseholderQr().solve(b);

    // Construct coefficients_
    coefficients_.resize(n);
    coefficients_[0] = 2 * controlPoints_[0] - a[0];
    for (size_t i = 0; i < n - 2; ++i)
    {
        coefficients_[i + 1] = a[i];
    }
    coefficients_[n - 1] = 2 * controlPoints_[n - 1] - a[n - 3];
}

/**
 * @brief Specialization of computeCoefficients for cubic B-splines (Dim = 3)
 *
 * This implementation handles different boundary conditions:
 * - PERIODIC: Ensures C2 continuity at the ends by making them periodic
 * - NATURAL: Second derivatives at endpoints are zero
 * - NOT_A_KNOT: Third derivatives are continuous at second and second-to-last knots
 * - COMPLETE: Specified first derivatives at endpoints
 * - SECOND: Specified second derivatives at endpoints
 *
 * @throws std::invalid_argument If insufficient control points or unsupported boundary condition
 */
template <>
void BSpline<3>::computeCoefficients()
{
    size_t n = controlPoints_.size();
    if (n < 4)
    {
        throw std::invalid_argument("At least 4 control points required for cubic B-spline.");
    }

    // Extend coefficient vector to support negative index access
    coefficients_.resize(n + 2); // Add 2 extra points for boundary handling

    // Build linear system
    Eigen::MatrixXd A = Eigen::MatrixXd::Zero(n + 2, n + 2);
    Eigen::VectorXd b = Eigen::VectorXd::Zero(n + 2);

    // Fill matrix based on boundary conditions
    for (size_t i = 1; i <= n; ++i)
    {
        for (size_t j = 0; j < 3; ++j)
        {
            A(i, i + j - 1) = basisFunction(i + j - 2, 3, clampedKnots_[i + 2]);
        }
        b(i) = controlPoints_[i - 1];
    }

    switch (boundaryCondition_)
    {
    case BoundaryCondition::NATURAL:
        // Natural boundary conditions (second derivative at endpoints are zero)
        for (size_t i = 0; i <= 2; ++i)
        {
            A(0, i) = basisFunction_second_derivative(i - 1, 3, clampedKnots_[3]);
            A(n + 1, n - 1 + i) = basisFunction_second_derivative(n - 2 + i, 3, clampedKnots_[n + 2]);
        }
        b(0) = 0.0;
        b(n + 1) = 0.0;
        break;

    case BoundaryCondition::COMPLETE:
        // Complete boundary conditions (specified first derivatives at endpoints)
        for (size_t i = 0; i <= 2; ++i)
        {
            A(0, i) = basisFunction_derivative(i - 1, 3, clampedKnots_[3]);
            A(n + 1, n - 1 + i) = basisFunction_derivative(n - 2 + i, 3, clampedKnots_[n + 2]);
        }
        b(0) = startDerivative_;
        b(n + 1) = endDerivative_;
        break;

    case BoundaryCondition::SECOND:
        // Second derivative boundary conditions (specified second derivatives at endpoints)
        for (size_t i = 0; i <= 2; ++i)
        {
            A(0, i) = basisFunction_second_derivative(i - 1, 3, clampedKnots_[3]);
            A(n + 1, n - 1 + i) = basisFunction_second_derivative(n - 2 + i, 3, clampedKnots_[n + 2]);
        }
        b(0) = startSecondDerivative_;
        b(n + 1) = endSecondDerivative_;
        break;

    case BoundaryCondition::NOT_A_KNOT:
        // Not-a-knot boundary conditions (third derivative exists at second and penultimate knots)
        for (int i = -1; i <= 4; ++i)
        {
            A(0, i + 1) -= basisFunction_third_derivative(i, 3, clampedKnots_[5]);
            A(0, i + 1) += basisFunction_third_derivative(i, 3, clampedKnots_[4]);
        }
        for (int i = -1; i <= 3; ++i)
        {
            A(n + 1, n - 2 + i) -= basisFunction_third_derivative(n - 3 + i, 3, clampedKnots_[n + 2]);
            A(n + 1, n - 2 + i) += basisFunction_third_derivative(n - 3 + i, 3, clampedKnots_[n + 1]);
        }

        b(0) = 0.0;
        b(n + 1) = 0.0;
        break;

    case BoundaryCondition::PERIODIC:
        // Periodic boundary conditions
        if (std::abs(controlPoints_.front() - controlPoints_.back()) > EPSILON)
        {
            throw std::invalid_argument("For periodic boundary conditions, the first and last control points must be equal.");
        }
        controlPoints_.back() = controlPoints_.front();
        for (size_t i = 0; i <= 2; ++i)
        {
            A(n + 1, i) = -basisFunction_second_derivative(i - 1, 3, clampedKnots_[3]);
            A(n + 1, n - 1 + i) = basisFunction_second_derivative(n - 2 + i, 3, clampedKnots_[n + 2]);
            A(0, i) = basisFunction_derivative(i - 1, 3, clampedKnots_[3]);
            A(0, n - 1 + i) = -basisFunction_derivative(n - 2 + i, 3, clampedKnots_[n + 2]);
        }
        b(n + 1) = 0.0;
        b(0) = 0.0;
        b(n) = controlPoints_[0];
        break;

    default:
        throw std::invalid_argument("Unsupported boundary condition.");
    }

    // Solve linear system
    Eigen::VectorXd x = A.colPivHouseholderQr().solve(b);

    // Fill coefficient array including boundary handling
    for (size_t i = 0; i < n + 2; ++i)
    {
        coefficients_[i] = x[i];
    }
}

#endif // BSPLINE_HPP