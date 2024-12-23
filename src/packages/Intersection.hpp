#pragma once

#include <Eigen/Dense>
#include <vector>
#include <cmath>
#include <limits>
#include <iostream>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "CurveFitting.hpp"

/**
 * @class Intersection
 * @brief A class to detect and store self-intersections of a curve.
 * 
 * This class provides functionalities to detect self-intersections in both
 * 2D and 3D curves. It categorizes intersections into three types:
 * - Closed without self-intersection
 * - Self-intersecting (not just closed)
 * - No self-intersection
 * 
 * It also stores detailed information about each intersection, including the
 * intersecting segments and the exact intersection points.
 */
class Intersection {
public:
    /**
     * @enum IntersectionType
     * @brief Enumeration of possible intersection types.
     */
    enum class IntersectionType {
        CLOSED_NO_INTERSECTION, ///< The curve is closed and does not self-intersect.
        SELF_INTERSECTION,      ///< The curve self-intersects (excluding closure).
        NO_SELF_INTERSECTION    ///< The curve does not self-intersect.
    };

    /**
     * @struct IntersectionInfo
     * @brief Stores information about a single self-intersection.
     */
    struct IntersectionInfo {
        Eigen::VectorXd segmentA_Start;    ///< Start point of the first intersecting segment.
        Eigen::VectorXd segmentA_End;      ///< End point of the first intersecting segment.
        Eigen::VectorXd segmentB_Start;    ///< Start point of the second intersecting segment.
        Eigen::VectorXd segmentB_End;      ///< End point of the second intersecting segment.
        Eigen::VectorXd intersectionPoint; ///< The point of intersection.
    };

    /**
     * @brief Constructor for the Intersection class.
     */
    Intersection();

    /**
     * @brief Detects self-intersections in the given curve by sampling from CurveFitting.
     * 
     * @param fittingCurve The curve object to be checked for self-intersections.
     * @param parameterStart The starting parameter value for sampling the curve.
     * @param parameterEnd The ending parameter value for sampling the curve.
     * @param num_samples The number of sample points to generate along the curve. Default is 2000.
     */
    void detectSelfIntersections(CurveFitting& fittingCurve,
                                 double parameterStart,
                                 double parameterEnd,
                                 int num_samples = 2000);

    /**
     * @brief Detects self-intersections from a direct set of curve points.
     * 
     * This method enables detection using a user-supplied vector of points. It
     * supports parallel execution (via OpenMP if available) to accelerate checks.
     * 
     * @param curvePoints A vector of Eigen::VectorXd representing the points defining the curve.
     */
    void detectSelfIntersectionsFromPoints(const std::vector<Eigen::VectorXd>& curvePoints);

    /**
     * @brief Retrieves the type of intersection detected.
     * 
     * @return IntersectionType The detected intersection type.
     */
    IntersectionType getIntersectionType() const;

    /**
     * @brief Retrieves all intersection information.
     * 
     * @return const std::vector<IntersectionInfo>& A vector containing information about each intersection.
     */
    const std::vector<IntersectionInfo>& getIntersections() const;

    /**
     * @brief Prints a summary of the self-intersection detection results.
     * 
     * This includes the number of intersections and whether the curve is closed.
     */
    void printSummary() const;

private:
    /**
     * @brief Checks if two line segments intersect in N-dimensional space.
     * 
     * If they intersect, stores the intersection info in outInfo and returns true.
     * Otherwise, returns false and outInfo remains unchanged.
     * 
     * @param startA The starting point of the first line segment.
     * @param endA The ending point of the first line segment.
     * @param startB The starting point of the second line segment.
     * @param endB The ending point of the second line segment.
     * @param outInfo IntersectionInfo containing segment data and intersection point if intersection occurs.
     * @param tolerance Threshold for checking parallelism. Default is 1e-9.
     * @return bool True if intersection occurs, false otherwise.
     */
    bool doSegmentsIntersect(const Eigen::VectorXd& startA,
                             const Eigen::VectorXd& endA,
                             const Eigen::VectorXd& startB,
                             const Eigen::VectorXd& endB,
                             IntersectionInfo& outInfo,
                             double tolerance = 1e-9) const;

    /**
     * @brief Clamps a value within a specified range.
     * 
     * @tparam T The data type.
     * @param val The value to be clamped.
     * @param lower The lower bound.
     * @param upper The upper bound.
     * @return T The clamped value.
     */
    template <typename T>
    inline T clamp(const T& val, const T& lower, const T& upper) const {
        return (val < lower) ? lower : (val > upper) ? upper : val;
    }

    /**
     * @brief Determines if two Eigen::VectorXd are approximately equal.
     * 
     * @param a The first vector.
     * @param b The second vector.
     * @param tolerance The tolerance for comparison. Default is 1e-6.
     * @return true If the vectors are approximately equal.
     * @return false Otherwise.
     */
    bool areVectorsEqual(const Eigen::VectorXd& a, const Eigen::VectorXd& b, double tolerance = 1e-6) const;

    IntersectionType intersectionType_;           ///< The detected intersection type.
    bool isClosed_;                               ///< Indicates if the curve is closed.
    std::vector<IntersectionInfo> intersections_; ///< Stores information about each intersection.
};

////////////////////////////////////////////////////////////////////////////////
// Implementation
////////////////////////////////////////////////////////////////////////////////

inline Intersection::Intersection()
    : intersectionType_(IntersectionType::NO_SELF_INTERSECTION),
      isClosed_(false) {}

inline bool Intersection::areVectorsEqual(const Eigen::VectorXd& a,
                                          const Eigen::VectorXd& b,
                                          double tolerance) const {
    if (a.size() != b.size()) return false;
    for (int i = 0; i < a.size(); ++i) {
        if (std::abs(a[i] - b[i]) > tolerance) return false;
    }
    return true;
}

inline bool Intersection::doSegmentsIntersect(const Eigen::VectorXd& startA,
                                             const Eigen::VectorXd& endA,
                                             const Eigen::VectorXd& startB,
                                             const Eigen::VectorXd& endB,
                                             IntersectionInfo& outInfo,
                                             double tolerance) const {
    Eigen::VectorXd directionA = endA - startA;
    Eigen::VectorXd directionB = endB - startB;
    Eigen::VectorXd relativePos = startA - startB;

    double dotAA = directionA.dot(directionA);
    double dotAB = directionA.dot(directionB);
    double dotBB = directionB.dot(directionB);
    double dotAR = directionA.dot(relativePos);
    double dotBR = directionB.dot(relativePos);

    double denominator = dotAA * dotBB - dotAB * dotAB;

    // Non-parallel case
    if (std::fabs(denominator) > tolerance) {
        double paramA = (dotAB * dotBR - dotBB * dotAR) / denominator;
        double paramB = (dotAA * dotBR - dotAB * dotAR) / denominator;

        // Check for intersection within the segment range
        if (paramA >= 0.0 && paramA <= 1.0 && paramB >= 0.0 && paramB <= 1.0) {
            outInfo.segmentA_Start = startA;
            outInfo.segmentA_End   = endA;
            outInfo.segmentB_Start = startB;
            outInfo.segmentB_End   = endB;
            outInfo.intersectionPoint = startA + paramA * directionA;
            return true;
        }
    } else {
        // Parallel or collinear case
        double projR = directionA.dot(relativePos);
        double projQ = directionA.dot(endB - startA);
        // Collinear check
        if (std::fabs(projR) < tolerance && std::fabs(projQ) < tolerance) {
            double lenA = directionA.squaredNorm();
            double lenB = directionB.squaredNorm();

            bool overlap = ((dotAR >= 0 && dotAR <= lenA) ||
                            (dotBR >= 0 && dotBR <= lenB));
            if (overlap) {
                outInfo.segmentA_Start = startA;
                outInfo.segmentA_End   = endA;
                outInfo.segmentB_Start = startB;
                outInfo.segmentB_End   = endB;
                // Use the midpoint of segment A as the "intersection"
                outInfo.intersectionPoint = startA + 0.5 * directionA;
                return true;
            }
        }
    }
    return false;
}

inline void Intersection::detectSelfIntersections(CurveFitting& fittingCurve,
                                                  double parameterStart,
                                                  double parameterEnd,
                                                  int num_samples) {
    intersections_.clear();
    intersectionType_ = IntersectionType::NO_SELF_INTERSECTION;
    isClosed_         = false;

    bool isThreeD = fittingCurve.is3D();
    std::vector<Eigen::VectorXd> sampledPoints;
    sampledPoints.reserve(num_samples + 1);

    // Sample points
    for (int idx = 0; idx <= num_samples; ++idx) {
        double t = clamp<double>(
            parameterStart + idx * (parameterEnd - parameterStart) / num_samples,
            parameterStart, parameterEnd
        );
        Eigen::VectorXd point;
        if (!isThreeD) {
            Eigen::Vector2d val = fittingCurve.plane_evaluate(t);
            point.resize(2);
            point << val.x(), val.y();
        } else {
            Eigen::Vector3d val = fittingCurve.sphere_evaluate(t);
            point.resize(3);
            point << val.x(), val.y(), val.z();
        }
        sampledPoints.push_back(point);
    }

    // Check closure
    if (areVectorsEqual(sampledPoints.front(), sampledPoints.back())) {
        isClosed_ = true;
        intersectionType_ = IntersectionType::CLOSED_NO_INTERSECTION;
    }

    int totalPoints = static_cast<int>(sampledPoints.size());
    bool foundIntersection = false;

#ifdef _OPENMP
#pragma omp parallel
    {
        std::vector<IntersectionInfo> localIntersections;
        localIntersections.reserve(128); // Arbitrary

#pragma omp for nowait
        for (int i = 0; i < totalPoints - 1; ++i) {
            for (int j = i + 2; j < totalPoints - 1; ++j) {
                if (j == i + 1) continue;
                if (isClosed_ && i == 0 && j == totalPoints - 2) continue;

                IntersectionInfo info;
                if (doSegmentsIntersect(sampledPoints[i], sampledPoints[i + 1],
                                        sampledPoints[j], sampledPoints[j + 1],
                                        info)) {
                    localIntersections.push_back(info);
                }
            }
        }

#pragma omp critical
        {
            if (!localIntersections.empty()) {
                foundIntersection = true;
                intersections_.insert(intersections_.end(),
                                      localIntersections.begin(),
                                      localIntersections.end());
            }
        }
    }
#else
    // Serial fallback
    for (int i = 0; i < totalPoints - 1; ++i) {
        for (int j = i + 2; j < totalPoints - 1; ++j) {
            if (j == i + 1) continue;
            if (isClosed_ && i == 0 && j == totalPoints - 2) continue;

            IntersectionInfo info;
            if (doSegmentsIntersect(sampledPoints[i], sampledPoints[i + 1],
                                    sampledPoints[j], sampledPoints[j + 1],
                                    info)) {
                foundIntersection = true;
                intersections_.push_back(info);
            }
        }
    }
#endif

    if (foundIntersection) {
        intersectionType_ = IntersectionType::SELF_INTERSECTION;
    } else if (isClosed_) {
        intersectionType_ = IntersectionType::CLOSED_NO_INTERSECTION;
    } else {
        intersectionType_ = IntersectionType::NO_SELF_INTERSECTION;
    }
}

inline void Intersection::detectSelfIntersectionsFromPoints(const std::vector<Eigen::VectorXd>& curvePoints) {
    intersections_.clear();
    intersectionType_ = IntersectionType::NO_SELF_INTERSECTION;
    isClosed_         = false;

    if (curvePoints.size() < 2) {
        intersectionType_ = IntersectionType::NO_SELF_INTERSECTION;
        return;
    }

    // Check closure
    if (areVectorsEqual(curvePoints.front(), curvePoints.back())) {
        isClosed_ = true;
        intersectionType_ = IntersectionType::CLOSED_NO_INTERSECTION;
    }

    int totalPoints = static_cast<int>(curvePoints.size());
    bool foundIntersection = false;

#ifdef _OPENMP
#pragma omp parallel
    {
        std::vector<IntersectionInfo> localIntersections;
        localIntersections.reserve(128);

#pragma omp for nowait
        for (int i = 0; i < totalPoints - 1; ++i) {
            for (int j = i + 2; j < totalPoints - 1; ++j) {
                if (j == i + 1) continue;
                if (isClosed_ && i == 0 && j == totalPoints - 2) continue;

                IntersectionInfo info;
                if (doSegmentsIntersect(curvePoints[i],   curvePoints[i + 1],
                                        curvePoints[j],   curvePoints[j + 1],
                                        info)) {
                    localIntersections.push_back(info);
                }
            }
        }

#pragma omp critical
        {
            if (!localIntersections.empty()) {
                foundIntersection = true;
                intersections_.insert(intersections_.end(),
                                      localIntersections.begin(),
                                      localIntersections.end());
            }
        }
    }
#else
    // Serial approach
    for (int i = 0; i < totalPoints - 1; ++i) {
        for (int j = i + 2; j < totalPoints - 1; ++j) {
            if (j == i + 1) continue;
            if (isClosed_ && i == 0 && j == totalPoints - 2) continue;

            IntersectionInfo info;
            if (doSegmentsIntersect(curvePoints[i],   curvePoints[i + 1],
                                    curvePoints[j],   curvePoints[j + 1],
                                    info)) {
                foundIntersection = true;
                intersections_.push_back(info);
            }
        }
    }
#endif

    if (foundIntersection) {
        intersectionType_ = IntersectionType::SELF_INTERSECTION;
    } else if (isClosed_) {
        intersectionType_ = IntersectionType::CLOSED_NO_INTERSECTION;
    } else {
        intersectionType_ = IntersectionType::NO_SELF_INTERSECTION;
    }
}

inline Intersection::IntersectionType Intersection::getIntersectionType() const {
    return intersectionType_;
}

inline const std::vector<Intersection::IntersectionInfo>& Intersection::getIntersections() const {
    return intersections_;
}

inline void Intersection::printSummary() const {
    std::cout << "---- Self-Intersection Summary ----\n";
    std::cout << "Curve is " << (isClosed_ ? "closed." : "open.") << "\n";

    switch (intersectionType_) {
    case IntersectionType::CLOSED_NO_INTERSECTION:
        std::cout << "The curve is closed and does not self-intersect.\n";
        break;
    case IntersectionType::SELF_INTERSECTION:
        std::cout << "The curve has self-intersections.\n";
        std::cout << "Number of Intersection Points: " << intersections_.size() << "\n";
        std::cout << "Details of Intersections:\n";
        for (size_t i = 0; i < intersections_.size(); ++i) {
            const auto& info = intersections_[i];
            std::cout << "  Intersection " << (i + 1) << ":\n";
            std::cout << "    Segment A: (" << info.segmentA_Start.transpose()
                      << ") --> (" << info.segmentA_End.transpose() << ")\n";
            std::cout << "    Segment B: (" << info.segmentB_Start.transpose()
                      << ") --> (" << info.segmentB_End.transpose() << ")\n";
            std::cout << "    Intersection Point: ("
                      << info.intersectionPoint.transpose() << ")\n";
        }
        break;
    case IntersectionType::NO_SELF_INTERSECTION:
        std::cout << "The curve does not self-intersect.\n";
        break;
    default:
        std::cout << "Unknown intersection type.\n";
        break;
    }

    std::cout << "-----------------------------------\n";
}