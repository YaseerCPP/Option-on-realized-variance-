#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <algorithm>
#include <numeric>

// Function to simulate asset price paths using the Black-Scholes model
std::vector<std::vector<double>> simulatePricePaths(double S0, double mu, double sigma, double T, int numSteps, int numPaths) {
    std::vector<std::vector<double>> pricePaths(numPaths, std::vector<double>(numSteps + 1, S0));
    double dt = T / numSteps;
    std::random_device rd;
    std::mt19937 generator(rd());
    std::normal_distribution<> distribution(0.0, 1.0);

    for (int i = 0; i < numPaths; ++i) {
        for (int j = 1; j <= numSteps; ++j) {
            double dW = distribution(generator) * sqrt(dt);
            pricePaths[i][j] = pricePaths[i][j-1] * exp((mu - 0.5 * sigma * sigma) * dt + sigma * dW);
        }
    }
    return pricePaths;
}

// Function to calculate realized variance
double calculateRealizedVariance(const std::vector<std::vector<double>>& pricePaths) {
    std::vector<double> returns;
    for (const auto& path : pricePaths) {
        for (size_t i = 1; i < path.size(); ++i) {
            double ret = log(path[i] / path[i-1]);
            returns.push_back(ret);
        }
    }
    double mean = std::accumulate(returns.begin(), returns.end(), 0.0) / returns.size();
    double sq_sum = std::inner_product(returns.begin(), returns.end(), returns.begin(), 0.0);
    return sq_sum / returns.size() - mean * mean;
}

// Main function to price options on realized variance using Monte Carlo simulation
int main() {
    // Parameters
    double S0 = 100.0;    // Initial asset price
    double mu = 0.1;      // Drift
    double sigma = 0.2;   // Volatility
    double T = 1.0;       // Time to maturity in years
    int numSteps = 252;   // Number of time steps (e.g., daily steps)
    int numPaths = 10000; // Number of Monte Carlo paths
    double strikeVariance = 0.04; // Strike variance
    double notional = 1000000;    // Notional amount

    // Simulate price paths
    auto pricePaths = simulatePricePaths(S0, mu, sigma, T, numSteps, numPaths);

    // Calculate the realized variance
    double realizedVariance = calculateRealizedVariance(pricePaths);

    // Calculate payoffs for call and put options on realized variance
    double callPayoff = std::max(realizedVariance - strikeVariance, 0.0) * notional;
    double putPayoff = std::max(strikeVariance - realizedVariance, 0.0) * notional;

    // Output the results
    std::cout << "Realized Variance: " << realizedVariance << std::endl;
    std::cout << "Call Option Payoff: " << callPayoff << std::endl;
    std::cout << "Put Option Payoff: " << putPayoff << std::endl;

    return 0;
}
