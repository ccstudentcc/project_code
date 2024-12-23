// BoundaryCondition.h
#pragma once

/**
 * @enum BoundaryCondition
 * @brief Specifies the type of boundary condition for the B-spline.
 */
enum class BoundaryCondition { 
    PERIODIC,                     /**< Periodic boundary conditions */
    NATURAL,                      /**< Natural boundary conditions */
    COMPLETE,                     /**< Complete boundary conditions */
    NOT_A_KNOT,                   /**< Not-a-knot boundary conditions */
    SECOND    /**< Specified second derivatives boundary conditions */
};

