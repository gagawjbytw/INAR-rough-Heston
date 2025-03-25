#include <iostream>
#include <vector>
#include <random>
#include <chrono>
#include <thread>
#include <future>
#include <iomanip>
#include <cmath>
#include <map>
#include <string>
#include <algorithm>

using namespace std;

// Parameter validation function
void validate_parameters(double S0, double V0, double mu, double r, double gamma,
                       double theta, double sigma, double rho, double dt) {
    if (S0 <= 0) throw runtime_error("Initial stock price must be positive");
    if (V0 < 0) throw runtime_error("Initial volatility must be non-negative");
    if (gamma <= 0) throw runtime_error("Mean reversion speed must be positive");
    if (theta <= 0) throw runtime_error("Long-term volatility must be positive");
    if (sigma <= 0) throw runtime_error("Volatility of volatility must be positive");
    if (rho <= -1 || rho >= 1) throw runtime_error("Correlation coefficient must be in [-1,1]");
    if (dt <= 0) throw runtime_error("Time step must be positive");
    
    // Check Feller condition
    if (2 * gamma * theta <= sigma * sigma) {
        cout << "Warning: Feller condition not satisfied, volatility process may touch zero" << endl;
    }
}

// Normal distribution random number generator
class NormalGenerator {
private:
    mt19937 gen;
    normal_distribution<double> dist;

public:
    NormalGenerator(int seed = 9999999) : gen(seed), dist(0.0, 1.0) {}  // Using fixed seed 42
    double operator()() { return dist(gen); }
};

// Heston model path simulation
struct PathResult {
    vector<double> S;
    vector<double> V;
    vector<double> S_anti;
    vector<double> V_anti;
};

PathResult simulate_heston(double S0, double V0, double mu, double r, double gamma,
                         double theta, double sigma, double rho, double T,
                         double dt, int seed) {
    int M = static_cast<int>(T / dt);
    PathResult result;
    result.S.resize(M + 1);
    result.V.resize(M + 1);
    result.S_anti.resize(M + 1);
    result.V_anti.resize(M + 1);
    
    NormalGenerator rng(seed);
    
    result.S[0] = result.S_anti[0] = S0;
    result.V[0] = result.V_anti[0] = V0;
    
    for (int i = 0; i < M; ++i) {
        double Z1 = rng();
        double Z2 = rho * Z1 + sqrt(1 - rho * rho) * rng();
        
        // 正向路径
        double dV = gamma * (theta - result.V[i]) * dt + 
                   sigma * sqrt(max(result.V[i], 0.0)) * sqrt(dt) * Z2;
        result.V[i + 1] = max(result.V[i] + dV, 0.0);
        double dS = r * result.S[i] * dt + 
                   sqrt(max(result.V[i], 0.0)) * result.S[i] * sqrt(dt) * Z1;
        result.S[i + 1] = max(result.S[i] + dS, 0.0);
        
        // 对偶路径
        double dV_anti = gamma * (theta - result.V_anti[i]) * dt - 
                        sigma * sqrt(max(result.V_anti[i], 0.0)) * sqrt(dt) * Z2;
        result.V_anti[i + 1] = max(result.V_anti[i] + dV_anti, 0.0);
        double dS_anti = r * result.S_anti[i] * dt - 
                        sqrt(max(result.V_anti[i], 0.0)) * result.S_anti[i] * sqrt(dt) * Z1;
        result.S_anti[i + 1] = max(result.S_anti[i] + dS_anti, 0.0);
    }
    
    return result;
}

// Option payoff calculation
double compute_payoff(const vector<double>& S_path, double K, const string& option_type,
                     double barrier = -1) {
    double S_T = S_path.back();
    
    if (option_type == "european_call") {
        return max(S_T - K, 0.0);
    }
    else if (option_type == "european_put") {
        return max(K - S_T, 0.0);
    }
    else if (option_type == "asian_call") {
        double avg_S = 0;
        for (double s : S_path) avg_S += s;
        avg_S /= S_path.size();
        return max(avg_S - K, 0.0);
    }
    else if (option_type == "asian_put") {
        double avg_S = 0;
        for (double s : S_path) avg_S += s;
        avg_S /= S_path.size();
        return max(K - avg_S, 0.0);
    }
    else if (option_type == "lookback_call") {
        double max_S = *max_element(S_path.begin(), S_path.end());
        return max(max_S - K, 0.0);
    }
    else if (option_type == "lookback_put") {
        double min_S = *min_element(S_path.begin(), S_path.end());
        return max(K - min_S, 0.0);
    }
    else if (option_type == "up_in_call") {
        double max_S = *max_element(S_path.begin(), S_path.end());
        return (max_S >= barrier) ? max(S_T - K, 0.0) : 0.0;
    }
    else if (option_type == "down_in_call") {
        double min_S = *min_element(S_path.begin(), S_path.end());
        return (min_S <= barrier) ? max(S_T - K, 0.0) : 0.0;
    }
    else if (option_type == "up_out_put") {
        double max_S = *max_element(S_path.begin(), S_path.end());
        return (max_S < barrier) ? max(K - S_T, 0.0) : 0.0;
    }
    else if (option_type == "down_out_put") {
        double min_S = *min_element(S_path.begin(), S_path.end());
        return (min_S > barrier) ? max(K - S_T, 0.0) : 0.0;
    }
    
    throw runtime_error("不支持的期权类型: " + option_type);
}

// Monte Carlo simulation function
vector<double> monte_carlo_simulation(double S0, double V0, double mu, double r,
                                    double gamma, double theta, double sigma,
                                    double rho, double T, double dt, double K,
                                    const string& option_type, double barrier,
                                    int n_sims, int n_threads) {
    validate_parameters(S0, V0, mu, r, gamma, theta, sigma, rho, dt);
    
    vector<double> payoffs;
    payoffs.reserve(n_sims);
    
    auto worker = [&](int start, int end) {
        vector<double> local_payoffs;
        local_payoffs.reserve(end - start);
        
        for (int i = start; i < end; ++i) {
            auto paths = simulate_heston(S0, V0, mu, r, gamma, theta, sigma, rho, T, dt, i);
            double payoff = compute_payoff(paths.S, K, option_type, barrier);
            double anti_payoff = compute_payoff(paths.S_anti, K, option_type, barrier);
            local_payoffs.push_back((payoff + anti_payoff) / 2);
        }
        
        return local_payoffs;
    };
    
    vector<future<vector<double>>> futures;
    int batch_size = n_sims / n_threads;
    
    for (int i = 0; i < n_threads; ++i) {
        int start = i * batch_size;
        int end = (i == n_threads - 1) ? n_sims : (i + 1) * batch_size;
        futures.push_back(async(launch::async, worker, start, end));
    }
    
    for (auto& future : futures) {
        auto result = future.get();
        payoffs.insert(payoffs.end(), result.begin(), result.end());
    }
    
    return payoffs;
}

// Calculate price and confidence interval
tuple<double, double, double> compute_price_and_ci(const vector<double>& payoffs,
                                                 double r, double T) {
    double mean_payoff = 0;
    for (double p : payoffs) mean_payoff += p;
    mean_payoff /= payoffs.size();
    
    double variance = 0;
    for (double p : payoffs) {
        double diff = p - mean_payoff;
        variance += diff * diff;
    }
    variance /= (payoffs.size() - 1);
    
    double std_payoff = sqrt(variance);
    double discount = exp(-r * T);
    double price = mean_payoff * discount;
    double ci_width = 1.96 * std_payoff / sqrt(payoffs.size()) * discount;
    
    return {price, price - ci_width, price + ci_width};
}

int main() {
    // Parameter settings
    double S0 = 100;      // Initial price
    double V0 = 0.0392;   // Initial volatility
    double mu = 0;        // Drift rate
    double r = 0;         // Risk-free rate
    double gamma = 0.1;   // Mean reversion speed
    double theta = 0.3156;// Long-term volatility
    double nu = 0.331;
    double sigma = gamma * nu;// Volatility of volatility
    double rho = -0.681;  // Price-volatility correlation
    double T = 1;         // Time to maturity
    double dt = 0.004;    // Time step
    int n_sims = 100000;  // Number of simulation paths
    int n_threads = 10;   // Number of threads
    vector<double> strikes = {80, 90, 100, 110, 120};  // Strike prices
    map<string, double> barriers = {
        {"up_in_call", 110.0},    // Up-and-in barrier
        {"down_out_put", 90.0}    // Down-and-out barrier
    };
    // Note: Using fixed seed 42 for simulation, results will be identical for each run

    // Option types in desired order
    vector<pair<string, vector<string>>> option_types = {
        {"European", {"european_call", "european_put"}},
        {"Asian", {"asian_call", "asian_put"}},
        {"Lookback", {"lookback_call", "lookback_put"}},
        {"Barrier", {"up_in_call", "down_out_put", "down_in_call", "up_out_put"}}
    };

    // Record start time
    auto start_time = chrono::high_resolution_clock::now();
    
    cout << "Starting Heston Model Option Pricing..." << endl;
    cout << fixed << setprecision(4);
    cout << "Parameters: S0=" << S0 << ", V0=" << V0 << ", r=" << r 
         << ", gamma=" << gamma << ", theta=" << theta 
         << ", sigma=" << sigma << ", rho=" << rho << endl;
    cout << "Simulation settings: Paths=" << n_sims << ", Time step=" << dt 
         << ", Threads=" << n_threads << endl;

    // Store results
    vector<pair<string, vector<tuple<double, string, double, tuple<double, double>>>>> results;
    results.reserve(4); // Reserve space for 4 option types

    try {
        // Run simulation
        for (const auto& [opt_group, types] : option_types) {
            cout << "\nCalculating " << opt_group << " options..." << endl;
            vector<tuple<double, string, double, tuple<double, double>>> group_results;
            
            for (double K : strikes) {
                for (const string& opt_type : types) {
                    double current_barrier = (opt_group == "Barrier") ? 
                        barriers[opt_type] : -1;  // Use different barriers for different options
                    vector<double> payoffs = monte_carlo_simulation(
                        S0, V0, mu, r, gamma, theta, sigma, rho, T, dt, K,
                        opt_type, current_barrier, n_sims, n_threads
                    );
                    
                    auto [price, ci_lower, ci_upper] = compute_price_and_ci(payoffs, r, T);
                    group_results.push_back({K, opt_type, price, {ci_lower, ci_upper}});
                }
            }
            results.push_back({opt_group, group_results});
        }

        // Print results
        for (const auto& [opt_group, data] : results) {
            cout << "\n" << opt_group << " Options:" << endl;
            if (opt_group == "Barrier") {
                cout << "Strike    Up-In-Call     95% CI      Down-Out-Put    95% CI" << endl;
                cout << "-------------------------------------------------------------" << endl;
                for (size_t i = 0; i < data.size(); i += 4) {
                    auto [K1, _, price1, ci1] = data[i];
                    auto [K2, __, price2, ci2] = data[i + 1];
                    cout << fixed << setprecision(4)
                         << setw(8) << K1 << "  "
                         << setw(8) << price1 << "  ["
                         << setw(8) << get<0>(ci1) << ", "
                         << setw(8) << get<1>(ci1) << "]  "
                         << setw(8) << price2 << "  ["
                         << setw(8) << get<0>(ci2) << ", "
                         << setw(8) << get<1>(ci2) << "]" << endl;
                }
            } else {
                cout << "Strike    Call        95% CI         Put         95% CI" << endl;
                cout << "--------------------------------------------------------" << endl;
                for (size_t i = 0; i < data.size(); i += 2) {
                    auto [K, _, call_price, call_ci] = data[i];
                    auto [__, ___, put_price, put_ci] = data[i + 1];
                    cout << fixed << setprecision(4)
                         << setw(8) << K << "  "
                         << setw(8) << call_price << "  ["
                         << setw(8) << get<0>(call_ci) << ", "
                         << setw(8) << get<1>(call_ci) << "]  "
                         << setw(8) << put_price << "  ["
                         << setw(8) << get<0>(put_ci) << ", "
                         << setw(8) << get<1>(put_ci) << "]" << endl;
                }
            }
        }
    }
    catch (const exception& e) {
        cerr << "Error occurred: " << e.what() << endl;
        return 1;
    }

    // Calculate runtime
    auto end_time = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::milliseconds>(end_time - start_time);
    cout << "\nTotal time used: " << duration.count() / 1000.0 << " seconds" << endl;

    return 0;
} 