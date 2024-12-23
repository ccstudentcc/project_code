// SplineType.h
#ifndef SPLINE_TYPE_H
#define SPLINE_TYPE_H

/**
 * @enum SplineType
 * @brief Specifies the type of spline to use for fitting.
 */
enum class SplineType {
    PPForm,            /**< Piecewise Polynomial Form spline */
    BSpline            /**< B-spline */
};

#endif // SPLINE_TYPE_H