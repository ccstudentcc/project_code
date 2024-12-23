// intersect_test.cpp
#include <fstream>
#include <iostream>
#include <cmath>
#include "../packages/Intersection.hpp"
#include <Eigen/Dense>
#include <vector>
#include <functional>

const double PI = acos(-1.0);
const double step = 0.01;

// Function to generate points
std::vector<Eigen::VectorXd> generatePoints(const std::function<Eigen::VectorXd(double)> &func, double start, double end, double step)
{
    std::vector<Eigen::VectorXd> points;
    double t = start;
    while (t <= end)
    {
        points.push_back(func(t));
        t += step;
    }
    if (t - step < end)
    {
        points.push_back(func(end));
    }
    return points;
}

int main()
{
    // Define r1, r2, r3 curves
    auto r1 = [](double t) -> Eigen::VectorXd
    {
        double x = std::sqrt(3.0) * std::sin(t);
        double y = (2.0 * std::sqrt(3.0) / 3.0) * std::cos(t) + (2.0 / 3.0) * std::sqrt(std::sqrt(3.0) * std::abs(std::sin(t)));
        Eigen::VectorXd pt(2);
        pt << x, y;
        return pt;
    };

    auto r2 = [](double t) -> Eigen::VectorXd
    {
        double x = std::sin(t) + t * std::cos(t);
        double y = std::cos(t) - t * std::sin(t);
        Eigen::VectorXd pt(2);
        pt << x, y;
        return pt;
    };

    auto r3 = [](double t) -> Eigen::VectorXd
    {
        double u = std::cos(t);
        double v = std::sin(t);
        double x = std::sin(u) * std::cos(v);
        double y = std::sin(u) * std::sin(v);
        double z = std::cos(u);
        Eigen::VectorXd pt(3);
        pt << x, y, z;
        return pt;
    };

    // Define r4 curve
    auto r4 = [](double t) -> Eigen::VectorXd
    {
        double x = (1.5) * std::sin(t) + 0.5 * std::sin(1.5 * t);
        double y = (1.5) * std::cos(t) - 0.5 * std::cos(1.5 * t);
        Eigen::VectorXd pt(2);
        pt << x, y;
        return pt;
    };

    // Generate and detect intersections
    auto points1 = generatePoints(r1, 0.0, 2.0 * PI, step);

    Intersection intersection1;
    intersection1.detectSelfIntersectionsFromPoints(points1);

    auto points2 = generatePoints(r2, 0.0, 6.0 * PI, step);
    Intersection intersection2;
    intersection2.detectSelfIntersectionsFromPoints(points2);

    auto points3 = generatePoints(r3, 0.0, 2.0 * PI, step);
    Intersection intersection3;
    intersection3.detectSelfIntersectionsFromPoints(points3);

    // Generate and detect intersections for r4
    auto points4 = generatePoints(r4, 0.0, 4.0 * PI, step);
    Intersection intersection4;
    intersection4.detectSelfIntersectionsFromPoints(points4);

    // Redirect printSummary() output to file
    std::ofstream outFile("../../output/intersect/intersect_curve.txt");
    if (outFile.is_open())
    {
        auto oldBuf = std::cout.rdbuf(outFile.rdbuf());
        std::cout << "Curve r1:" << std::endl;
        intersection1.printSummary();
        std::cout << "Curve r2:" << std::endl;
        intersection2.printSummary();
        std::cout << "Curve r3:" << std::endl;
        intersection3.printSummary();
        std::cout << "Curve r4:" << std::endl;
        intersection4.printSummary();
        std::cout.rdbuf(oldBuf);
        outFile.close();
    }

    return 0;
}
