#ifndef FUNCTION
#define FUNCTION

#include <vector>
#include <cstddef>
#include <stdexcept>
#include <iostream>
#include <iomanip> // For std::setprecision
#include <cmath>

/**
 * @brief Abstract base class for mathematical functions.
 *
 * This class provides a virtual interface for functions that can be evaluated,
 * differentiated, and have their n-th derivatives calculated.
 */
class Function
{
public:
    /**
     * @brief Evaluate the function at a given x.
     *
     * @param x The input value.
     * @return The value of the function at x.
     */
    virtual double operator()(double x) const = 0;

    /**
     * @brief Calculate the first derivative of the function at a given x using central difference.
     *
     * @param x The point at which to evaluate the derivative.
     * @return The value of the first derivative at x.
     */
    virtual double derivative(double x) const
    {
        const double h = 1e-5;                              ///< A small increment
        return ((*this)(x + h) - (*this)(x - h)) / (2 * h); // Central difference formula
    }

    /**
     * @brief Calculate the n-th derivative of the function at a given x.
     *
     * This method uses a recursive central difference approach.
     *
     * @param x The point at which to evaluate the n-th derivative.
     * @param n The order of the derivative.
     * @return The value of the n-th derivative at x.
     * @throws std::invalid_argument if n is negative.
     */
    virtual double nth_derivative(double x, int n) const
    {
        if (n < 0)
        {
            throw std::invalid_argument("Derivative order n must be non-negative.");
        }
        if (n == 0)
        {
            return (*this)(x); // Zeroth derivative is the function itself
        }
        const double h = std::pow(10.0, -5.0 / n); // Adaptive step size for higher derivatives
        return (nth_derivative(x + h, n - 1) - nth_derivative(x - h, n - 1)) / (2 * h);
    }

    /**
     * @brief Virtual destructor.
     */
    virtual ~Function() = default;
};

/**
 * @brief Class for representing polynomial functions.
 *
 * This class implements a polynomial function of the form
 * f(x) = a_n * x^n + a_(n-1) * x^(n-1) + ... + a_1 * x + a_0.
 */
class Polynomial : public Function
{
private:
    std::vector<double> coefficients; ///< Store the coefficients of the polynomial

public:
    /**
     * @brief Default constructor initializes the polynomial to zero.
     */
    Polynomial() : coefficients({0.0}) {}

    /**
     * @brief Constructor with coefficients.
     *
     * @param coeffs A vector of coefficients for the polynomial.
     */
    Polynomial(const std::vector<double> &coeffs) : coefficients(coeffs) {}

    /**
     * @brief Copy constructor.
     *
     * @param other The Polynomial object to copy.
     */
    Polynomial(const Polynomial &other) : coefficients(other.coefficients) {}

    /**
     * @brief Copy assignment operator.
     *
     * @param other The Polynomial object to copy.
     * @return A reference to the current object.
     */
    Polynomial &operator=(const Polynomial &other)
    {
        if (this != &other) // Check for self-assignment
        {
            coefficients = other.coefficients; // Copy coefficients
        }
        return *this; // Return the current object
    }
    
    // Overload the + operator for polynomial addition
    /**
     * @brief Add two polynomials.
     *
     * @param other The polynomial to add.
     * @return A new Polynomial object representing the sum.
     */
    Polynomial operator+(const Polynomial &other) const
    {
        size_t max_size = std::max(coefficients.size(), other.coefficients.size());
        std::vector<double> result_coeffs(max_size, 0.0);

        // Add coefficients from both polynomials
        for (size_t i = 0; i < max_size; ++i)
        {
            if (i < coefficients.size())
            {
                result_coeffs[i] += coefficients[i];
            }
            if (i < other.coefficients.size())
            {
                result_coeffs[i] += other.coefficients[i];
            }
        }

        return Polynomial(result_coeffs); // Return the resulting polynomial
    }

    // Overload the - operator for polynomial subtraction
    /**
     * @brief Subtract one polynomial from another.
     *
     * @param other The polynomial to subtract.
     * @return A new Polynomial object representing the difference.
     */
    Polynomial operator-(const Polynomial &other) const
    {
        size_t max_size = std::max(coefficients.size(), other.coefficients.size());
        std::vector<double> result_coeffs(max_size, 0.0);

        // Subtract coefficients from both polynomials
        for (size_t i = 0; i < max_size; ++i)
        {
            if (i < coefficients.size())
            {
                result_coeffs[i] += coefficients[i];
            }
            if (i < other.coefficients.size())
            {
                result_coeffs[i] -= other.coefficients[i];
            }
        }

        return Polynomial(result_coeffs); // Return the resulting polynomial
    }

    // Overload the * operator for polynomial multiplication
    /**
     * @brief Multiply two polynomials.
     *
     * @param other The polynomial to multiply.
     * @return A new Polynomial object representing the product.
     */
    Polynomial operator*(const Polynomial &other) const
    {
        size_t result_size = coefficients.size() + other.coefficients.size() - 1;
        std::vector<double> result_coeffs(result_size, 0.0);

        // Multiply coefficients of both polynomials
        for (size_t i = 0; i < coefficients.size(); ++i)
        {
            for (size_t j = 0; j < other.coefficients.size(); ++j)
            {
                result_coeffs[i + j] += coefficients[i] * other.coefficients[j];
            }
        }

        return Polynomial(result_coeffs); // Return the resulting polynomial
    }

    /**
     * @brief Evaluate the polynomial at a given x.
     *
     * @param x The input value.
     * @return The value of the polynomial at x.
     */
    virtual double operator()(double x) const override
    {
        double result = 0.0;
        double power = 1.0; // Start from x^0
        for (size_t i = 0; i < coefficients.size(); ++i)
        {
            result += coefficients[i] * power;
            power *= x; // Update the power of x
        }
        return result;
    }

    /**
     * @brief Print the polynomial in human-readable form.
     */
    void print() const
    {
        int n = coefficients.size();

        // Print the first term (constant)
        if (coefficients[0] != 0.0)
        {
            std::cout << std::setprecision(7) << coefficients[0];
        }
        else if(n == 1 && coefficients[0] == 0){
            std::cout << 0 << std::endl;
        }

        // Print the remaining terms
        for (int i = 1; i < n; ++i)
        {
            if (coefficients[i] != 0.0)
            {
                if (coefficients[i] > 0.0)
                {
                    if (i == 1 && coefficients[0] == 0.0)
                    {
                    }
                    else
                    {
                        std::cout << "+";
                    }
                    if (coefficients[i] != 1.0)
                    {
                        std::cout << std::setprecision(7) << coefficients[i];
                    }
                }
                else if (coefficients[i] < 0.0)
                {
                    if (i == 1 && coefficients[0] == 0.0)
                    {
                    }
                    else
                    {
                        std::cout << "-";
                    }
                    if (coefficients[i] != -1.0)
                    {
                        std::cout << std::setprecision(7) << -coefficients[i];
                    }
                }
                if (!(coefficients[i] == 1.0 || coefficients[i] == -1.0))
                {
                    std::cout << "*";
                }
                if (i == 1)
                {
                    std::cout << "x";
                }
                else
                {
                    std::cout << "x^" << i;
                }
            }
        }
        std::cout << std::endl;
    }

    /**
     * @brief Return the first derivative of the polynomial as a new Polynomial.
     *
     * @return A Polynomial object representing the derivative.
     */
    Polynomial derivative_polynomial() const
    {
        std::vector<double> deriv_coeffs;

        // If the polynomial is a constant, the derivative is 0
        if (coefficients.size() <= 1)
        {
            deriv_coeffs.push_back(0);
        }
        else
        {
            // Calculate the coefficients of the derivative
            for (size_t i = 1; i < coefficients.size(); ++i)
            {
                deriv_coeffs.push_back(i * coefficients[i]); // a_i * i corresponds to x^(i-1)
            }
        }

        return Polynomial(deriv_coeffs); // Return the derivative polynomial
    }

    /**
     * @brief Calculate the value of the first derivative at a given x.
     *
     * @param x The point at which to evaluate the derivative.
     * @return The value of the derivative at x.
     */
    virtual double derivative(double x) const override
    {
        // Calculate the derivative value using the derivative polynomial
        Polynomial deriv_poly = this->derivative_polynomial(); // Get the first derivative polynomial
        return deriv_poly(x);                                  // Calculate the value of the derivative polynomial at x
    }

    /**
     * @brief Calculate the n-th derivative of the polynomial and return it as a new Polynomial.
     *
     * This function computes the n-th derivative of the current polynomial
     * and returns it as a new Polynomial object.
     *
     * @param n The order of the derivative to calculate. Must be a non-negative integer.
     * @return A Polynomial object representing the n-th derivative.
     * @throws std::invalid_argument if n is negative.
     */
    Polynomial nth_derivative_polynomial(int n) const
    {
        if (n < 0)
        {
            throw std::invalid_argument("Order must be a non-negative integer");
        }

        Polynomial current(*this); ///< Create a copy of the current polynomial
        for (int i = 0; i < n; ++i)
        {
            current = current.derivative_polynomial(); ///< Get the derivative polynomial of the current polynomial
        }
        return current; ///< Return the n-th derivative polynomial
    }

    /**
     * @brief Calculate the value of the n-th derivative at a given point.
     *
     * This function calculates the n-th derivative of the polynomial
     * and evaluates it at the specified point x.
     *
     * @param x The point at which to evaluate the n-th derivative.
     * @param n The order of the derivative to calculate. Must be a non-negative integer.
     * @return The value of the n-th derivative at point x.
     * @throws std::invalid_argument if n is negative.
     */
    virtual double nth_derivative(double x, int n) const override
    {
        if (n < 0)
        {
            throw std::invalid_argument("Order must be a non-negative integer");
        }

        Polynomial nth_deriv_poly = this->nth_derivative_polynomial(n); ///< Get the n-th derivative polynomial
        return nth_deriv_poly(x);                                       ///< Return the value of the n-th derivative polynomial at x
    }

    /**
     * @brief Find the extrema of the polynomial in the interval [a, b].
     *
     * This function finds the minimum and maximum values of the polynomial
     * in the specified interval, along with their corresponding x values.
     *
     * @param a The left endpoint of the interval.
     * @param b The right endpoint of the interval.
     * @return A pair containing the minimum value and its x coordinate,
     *         and the maximum value and its x coordinate.
     */
    std::pair<std::pair<double, double>, std::pair<double, double>> findExtrema(double a, double b)
    {
        // Get the derivative polynomial
        Polynomial deriv = this->derivative_polynomial();
        
        // Store the critical points
        std::vector<double> critical_points;

        // Simple grid search for critical points
        for (double x = a; x <= b; x += 0.000001) {
            if (std::abs(deriv(x)) < 1e-4) {  // If the derivative is close to zero
                critical_points.push_back(x);
            }
        }

        // Evaluate the polynomial at the endpoints
        double min_value = (*this)(a);
        double min_x = a;
        double max_value = (*this)(a);
        double max_x = a;

        // Check the values at the critical points
        for (double cp : critical_points) {
            double value = (*this)(cp);
            if (value < min_value) {
                min_value = value;
                min_x = cp; // Update min x
            }
            if (value > max_value) {
                max_value = value;
                max_x = cp; // Update max x
            }
        }

        // Check the value at the right endpoint
        double end_value = (*this)(b);
        if (end_value < min_value) {
            min_value = end_value;
            min_x = b; // Update min x
        }
        if (end_value > max_value) {
            max_value = end_value;
            max_x = b; // Update max x
        }

        return {{min_value, min_x}, {max_value, max_x}};
    }
};

#endif // FUNCTION
