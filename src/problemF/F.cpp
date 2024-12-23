// F.cpp
#include<iostream>
#include<cmath>
#include<string>
#include<vector>
#include<fstream>
#include <tuple>

// Function to calculate the truncated power function
double truncPowerFunc(double x, int degree) {
    if (x > 0) {
        return pow(x, degree);
    } else {
        return 0.0;
    }
}

// Function to calculate the divided difference of the truncated power function
// using the recursive formula by the definition of divided difference
double divideTruncate(double x, std::vector<double> t, int degree) {
    size_t n = t.size();
    if (n == 0) {
        throw "Knot should not be empty";
    }
    if (n == 1) {
        return truncPowerFunc(t[0] - x, degree);
    }
    if (n == 2) {
        return (truncPowerFunc(t[n - 1] - x, degree) - truncPowerFunc(t[0] - x, degree)) / (t[n - 1] - t[0]);
    } else {
        std::vector<double> t1(t.begin() + 1, t.end());
        std::vector<double> t2(t.begin(), t.end() - 1);
        return (divideTruncate(x, t1, degree) - divideTruncate(x, t2, degree)) / (t[n - 1] - t[0]);
    }
}

// Function to test the truncated power function and return results as a vector of tuples
std::vector<std::tuple<double, double>> testTruncate(double t, int degree, double begin, double end, int numPoints = 100) {
    std::vector<std::tuple<double, double>> results;
    results.reserve(numPoints);
    for (double x = begin; x <= end; x += (end - begin) / numPoints) {
        double value = truncPowerFunc(t - x, degree);
        results.emplace_back(x, value);
    }
    return results;
}

// Function to test the divided difference of the truncated power function and return results as a vector of tuples
std::vector<std::tuple<double, double>> testDividedTruncate(std::vector<double> t, int degree, double begin, double end, int numPoints = 100) {
    std::vector<std::tuple<double, double>> results;
    results.reserve(numPoints);
    for (double x = begin; x <= end; x += (end - begin) / numPoints) {
        double value = divideTruncate(x, t, degree);
        results.emplace_back(x, value);
    }
    return results;
}

// Function to save the truncated power function results to a file
void saveTruncate(double t, int degree, std::string filename, int numpoints = 100) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Failed to open max error file for writing: " << filename << std::endl;
        return;
    }

    std::vector<std::tuple<double, double>> results = testTruncate(t, degree, t - 1, t + 2, numpoints);

    // Write data
    for (const auto &entry : results) {
        file << std::get<0>(entry) << "," << std::get<1>(entry);
        file << "\n";
    }

    file.close();
}

// Function to save the divided difference of the truncated power function results to a file
void saveDivideTruncate(std::vector<double> t, int degree, std::string filename, int numpoints = 100) {
    try {
        if (t.size() != static_cast<size_t>(degree + 2)) {
            throw std::invalid_argument("Knot size should be degree+2");
        }
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
        return;
    }

    std::ofstream outfile(filename, std::ios::trunc);
    if (!outfile.is_open()) {
        std::cerr << "Failed to open max error file for writing: " << filename << std::endl;
        return;
    }

    size_t n = t.size();

    for (int i = 0; i < static_cast<int>(n); i++) {
        for (int j = 0; j <= i; j++) {
            std::vector<double> t1(t.begin() + (i - j), t.begin() + i + 1);
            std::vector<std::tuple<double, double>> results;
            results = testDividedTruncate(t1, degree, t.front() - 1, t.back(), numpoints);
            outfile << "(" << i << "," << j << ")" << std::endl;
            for (const auto &entry : results) {
                outfile << std::get<0>(entry) << "," << std::get<1>(entry) << std::endl;
            }
        }
    }
}

// Main function to execute the save functions and run the plotting scripts
int main() {
    saveTruncate(0, 1, "../../output/problemF/trunc_1.csv");
    std::string command = "python plotF_trunc.py ../../output/problemF/trunc_1.csv";
    int result = system(command.c_str());
	if (result != 0) {
    	std::cerr << "Command failed with exit code: " << result << std::endl;
	}

    saveTruncate(0, 2, "../../output/problemF/trunc_2.csv");
    command = "python plotF_trunc.py ../../output/problemF/trunc_2.csv";
    result = system(command.c_str());
	if (result != 0) {
    	std::cerr << "Command failed with exit code: " << result << std::endl;
	}

    saveDivideTruncate({0, 1, 2}, 1, "../../output/problemF/divided_1.csv");
    command = "python plotF_divide.py ../../output/problemF/divided_1.csv";
    result = system(command.c_str());
    if (result != 0) {
    	std::cerr << "Command failed with exit code: " << result << std::endl;
	}

    saveDivideTruncate({0, 1, 2, 3}, 2, "../../output/problemF/divided_2.csv");
    command = "python plotF_divide.py ../../output/problemF/divided_2.csv";
    result = system(command.c_str());
    if (result != 0) {
    	std::cerr << "Command failed with exit code: " << result << std::endl;
	}
}
