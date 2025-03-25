#include <vector>
#include <cmath>
#include <random>
#include <iostream>
#include <chrono>
#include <thread>
#include <future>
#include <iomanip>
#include <algorithm>

// Function to solve beta from rho
double solve_beta(double rho) {
    // Equation: rho = (1-beta)/sqrt(2(1+beta^2))
    // Rearranging: 2*rho^2*(1+beta^2) = (1-beta)^2
    // 2*rho^2 + 2*rho^2*beta^2 = 1 - 2*beta + beta^2
    // (2*rho^2 - 1)*beta^2 + 2*beta + (2*rho^2 - 1) = 0
    // Quadratic equation: a*beta^2 + b*beta + c = 0
    double a = 2 * rho * rho - 1;
    double b = 2;
    double c = 2 * rho * rho - 1;
    
    // Using quadratic formula to find both roots
    double discriminant = b * b - 4 * a * c;
    double root1 = (-b + std::sqrt(discriminant)) / (2 * a);
    double root2 = (-b - std::sqrt(discriminant)) / (2 * a);
    
    // Return the root that is greater than 1
    if (root1 > 1) return root1;
    if (root2 > 1) return root2;
    
    // If neither root is greater than 1 (shouldn't happen with valid rho)
    throw std::runtime_error("No valid beta found (beta should be > 1)");
}

// Function to solve mu from nu
double solve_mu(double nu, double theta, double beta, double gamma) {
    // Equation: nu = sqrt(theta*(1+beta^2)/(gamma*mu*(1+beta)^2))
    // Rearranging: mu = theta*(1+beta^2)/(gamma*nu^2*(1+beta)^2)
    return theta * (1 + beta * beta) / (gamma * nu * nu * (1 + beta) * (1 + beta));
}

// Function to solve xi from V0
double solve_xi(double V0, double theta) {
    // Equation: V0 = xi*theta
    return V0 / theta;
}

// Helper function to compute fractional kernel
double f_alpha_one(int t, double alpha) {
    if (t <= 0) return 0.0;
    if (t == 1) {
        return 1.0 - 1.0 / (std::tgamma(1.0 - alpha) * std::pow(2.0, alpha));
    }
    return 1.0 / std::tgamma(1.0 - alpha) * (1.0 / std::pow(t, alpha) - 1.0 / std::pow(t + 1, alpha));
}

// Structure to hold path simulation results
struct PathSimResult {
    std::vector<double> path;      // Full price path
    double final_price;            // S_T
    double avg_price;              // Average price for Asian options
    double max_price;              // Maximum price for lookback options
    double min_price;              // Minimum price for lookback options
    bool barrier_hit_110;          // Whether upper barrier (110) was hit
    bool barrier_hit_90;           // Whether lower barrier (90) was hit
};

// Simulate a single path and compute all necessary metrics
PathSimResult simulate_inar_path(
    double alpha, double gamma_val, int T_steps,
    double beta_param, double mu, double theta, double xi, double S0, 
    double barrier_up, double barrier_down,  // Modified to accept two barriers
    std::mt19937& gen
) {
    double a_T = 1.0 - gamma_val * std::pow(T_steps, -alpha);
    double mu_T = mu * std::pow(T_steps, alpha - 1.0);

    std::vector<double> phi(T_steps);
    for (int t = 0; t < T_steps; ++t) {
        phi[t] = a_T * f_alpha_one(t + 1, alpha);
    }
    std::vector<double> phi_cum(T_steps);
    phi_cum[0] = phi[0];
    for (int t = 1; t < T_steps; ++t) {
        phi_cum[t] = phi_cum[t - 1] + phi[t];
    }

    std::vector<int> X_plus(T_steps + 1, 0), X_minus(T_steps + 1, 0);
    std::vector<int> N_plus(T_steps + 1, 0), N_minus(T_steps + 1, 0);
    std::vector<double> price_path(T_steps + 1, S0);
    double sum_prices = S0;  // For Asian options

    for (int t_idx = 1; t_idx <= T_steps; ++t_idx) {
        double sum_phi_prev = (t_idx == 1) ? 0.0 : phi_cum[t_idx - 2];
        double hat_mu = mu_T + xi * mu_T * ((1.0 / (1.0 - a_T)) * (1.0 - sum_phi_prev) - sum_phi_prev);
        double lam_plus = hat_mu;
        double lam_minus = hat_mu;

        for (int s_idx = 1; s_idx < t_idx; ++s_idx) {
            double val = a_T * f_alpha_one(t_idx - s_idx, alpha) / (1.0 + beta_param);
            double increment = val * X_plus[s_idx] + val * beta_param * X_minus[s_idx];
            lam_plus += increment;
            lam_minus += increment;
        }

        std::poisson_distribution<int> dist_plus(lam_plus > 0 ? lam_plus : 0);
        std::poisson_distribution<int> dist_minus(lam_minus > 0 ? lam_minus : 0);
        X_plus[t_idx] = dist_plus(gen);
        X_minus[t_idx] = dist_minus(gen);
        N_plus[t_idx] = N_plus[t_idx - 1] + X_plus[t_idx];
        N_minus[t_idx] = N_minus[t_idx - 1] + X_minus[t_idx];

        double coef = std::sqrt(theta / 2.0 * (1.0 - a_T) / (std::pow(T_steps, alpha) * mu));
        double diff_pm = N_plus[t_idx] - N_minus[t_idx];
        double P_t = coef * diff_pm - (theta / 2.0) * (1.0 - a_T) / (std::pow(T_steps, alpha) * mu) * N_plus[t_idx];
        price_path[t_idx] = S0 * std::exp(P_t);
        sum_prices += price_path[t_idx];
    }

    PathSimResult result;
    result.path = price_path;
    result.final_price = price_path[T_steps];
    result.avg_price = sum_prices / (T_steps + 1);
    result.max_price = *std::max_element(price_path.begin(), price_path.end());
    result.min_price = *std::min_element(price_path.begin(), price_path.end());
    result.barrier_hit_110 = (result.max_price >= barrier_up);
    result.barrier_hit_90 = (result.min_price <= barrier_down);
    return result;
}

// Option payoff calculations
double european_call_payoff(double S_T, double strike) {
    return std::max(S_T - strike, 0.0);
}

double european_put_payoff(double S_T, double strike) {
    return std::max(strike - S_T, 0.0);
}

double asian_call_payoff(double avg_price, double strike) {
    return std::max(avg_price - strike, 0.0);
}

double asian_put_payoff(double avg_price, double strike) {
    return std::max(strike - avg_price, 0.0);
}

double lookback_call_payoff(double max_price, double strike) {
    return std::max(max_price - strike, 0.0);
}

double lookback_put_payoff(double min_price, double strike) {
    return std::max(strike - min_price, 0.0);
}

double barrier_up_in_call_payoff(const PathSimResult& result, double strike, double barrier_up) {
    if (result.barrier_hit_110) {
        return std::max(result.final_price - strike, 0.0);
    }
    return 0.0;
}

double barrier_down_out_put_payoff(const PathSimResult& result, double strike, double barrier_down) {
    if (!result.barrier_hit_90) {
        return std::max(strike - result.final_price, 0.0);
    }
    return 0.0;
}

// Structure to store simulation results
struct SimResults {
    std::vector<std::vector<double>> european_call_prices;
    std::vector<std::vector<double>> european_put_prices;
    std::vector<std::vector<double>> asian_call_prices;
    std::vector<std::vector<double>> asian_put_prices;
    std::vector<std::vector<double>> lookback_call_prices;
    std::vector<std::vector<double>> lookback_put_prices;
    std::vector<std::vector<double>> barrier_up_in_call_prices;
    std::vector<std::vector<double>> barrier_down_out_put_prices;

    SimResults(size_t n_strikes, size_t n_sims) {
        european_call_prices.resize(n_strikes, std::vector<double>(n_sims));
        european_put_prices.resize(n_strikes, std::vector<double>(n_sims));
        asian_call_prices.resize(n_strikes, std::vector<double>(n_sims));
        asian_put_prices.resize(n_strikes, std::vector<double>(n_sims));
        lookback_call_prices.resize(n_strikes, std::vector<double>(n_sims));
        lookback_put_prices.resize(n_strikes, std::vector<double>(n_sims));
        barrier_up_in_call_prices.resize(n_strikes, std::vector<double>(n_sims));
        barrier_down_out_put_prices.resize(n_strikes, std::vector<double>(n_sims));
    }
};

// Function to calculate mean and confidence intervals
struct StatResult {
    double mean;
    double ci_lower;
    double ci_upper;
};

StatResult calculate_stats(const std::vector<double>& prices) {
    double sum = 0.0;
    for (double price : prices) {
        sum += price;
    }
    double mean = sum / prices.size();
    
    double variance = 0.0;
    for (double price : prices) {
        variance += (price - mean) * (price - mean);
    }
    variance /= (prices.size() - 1);  // Using n-1 for sample variance
    
    // Calculate standard error
    double std_error = std::sqrt(variance / prices.size());
    
    // Calculate 95% confidence interval (using 1.96 for 95% CI)
    double ci_lower = mean - 1.96 * std_error;
    double ci_upper = mean + 1.96 * std_error;
    
    return {mean, ci_lower, ci_upper};
}

// Modified parallel Monte Carlo function
std::vector<std::vector<StatResult>> parallel_monte_carlo_all_prices(
    int n_sims, const std::vector<double>& strikes, 
    double barrier_up, double barrier_down,  // Modified to accept two barriers
    double alpha, double gamma_val, int T_steps,
    double beta_param, double mu, double theta, double xi, double S0,
    unsigned int base_seed = 12345
) {
    const unsigned int num_threads = std::thread::hardware_concurrency();
    SimResults results(strikes.size(), n_sims);
    std::vector<std::future<void>> futures;
    
    int sims_per_thread = n_sims / num_threads;
    for (unsigned int t = 0; t < num_threads; ++t) {
        int start_sim = t * sims_per_thread;
        int end_sim = (t == num_threads - 1) ? n_sims : (t + 1) * sims_per_thread;
        
        futures.push_back(std::async(std::launch::async, [&, start_sim, end_sim, t, base_seed]() {
            std::mt19937 gen(base_seed + t);
            
            for (int sim = start_sim; sim < end_sim; ++sim) {
                PathSimResult path_result = simulate_inar_path(
                    alpha, gamma_val, T_steps, beta_param, mu, theta, xi, S0, 
                    barrier_up, barrier_down, gen  // Updated to pass both barriers
                );
                
                for (size_t k = 0; k < strikes.size(); ++k) {
                    double strike = strikes[k];
                    
                    results.european_call_prices[k][sim] = european_call_payoff(path_result.final_price, strike);
                    results.european_put_prices[k][sim] = european_put_payoff(path_result.final_price, strike);
                    results.asian_call_prices[k][sim] = asian_call_payoff(path_result.avg_price, strike);
                    results.asian_put_prices[k][sim] = asian_put_payoff(path_result.avg_price, strike);
                    results.lookback_call_prices[k][sim] = lookback_call_payoff(path_result.max_price, strike);
                    results.lookback_put_prices[k][sim] = lookback_put_payoff(path_result.min_price, strike);
                    results.barrier_up_in_call_prices[k][sim] = barrier_up_in_call_payoff(path_result, strike, barrier_up);
                    results.barrier_down_out_put_prices[k][sim] = barrier_down_out_put_payoff(path_result, strike, barrier_down);
                }
            }
        }));
    }
    
    for (auto& future : futures) {
        future.get();
    }
    
    // Calculate statistics for each option type and strike
    std::vector<std::vector<StatResult>> all_stats(8, std::vector<StatResult>(strikes.size()));
    
    for (size_t k = 0; k < strikes.size(); ++k) {
        all_stats[0][k] = calculate_stats(results.european_call_prices[k]);
        all_stats[1][k] = calculate_stats(results.european_put_prices[k]);
        all_stats[2][k] = calculate_stats(results.asian_call_prices[k]);
        all_stats[3][k] = calculate_stats(results.asian_put_prices[k]);
        all_stats[4][k] = calculate_stats(results.lookback_call_prices[k]);
        all_stats[5][k] = calculate_stats(results.lookback_put_prices[k]);
        all_stats[6][k] = calculate_stats(results.barrier_up_in_call_prices[k]);
        all_stats[7][k] = calculate_stats(results.barrier_down_out_put_prices[k]);
    }
    
    return all_stats;
}

int main() {
    // 设置固定的基础随机种子
    unsigned int base_seed = 100000;
    
    int n_sims = 100000;
    std::vector<double> strikes;
    for (double K = 80; K <= 120; K += 10) {
        strikes.push_back(K);
    }
    
    // Basic parameters
    double barrier_up = 110.0;    // Upper barrier for up-and-in call
    double barrier_down = 90.0;   // Lower barrier for down-and-out put
    double r = 0;
    int T_steps = 250;
    double alpha = 0.62;
    double gamma_val = 0.1;
    double theta = 0.3156;
    double S0 = 100;
    
    // Model parameters that need to be solved
    double rho = -0.681;
    double nu = 0.331;
    double V0 = 0.0392;

    // Solve for derived parameters
    double beta_param = solve_beta(rho);
    double mu = solve_mu(nu, theta, beta_param, gamma_val);
    double xi = solve_xi(V0, theta);
    
    // Print derived parameters
    std::cout << "Derived Parameters:\n";
    std::cout << "beta = " << beta_param << " (from rho = " << rho << ")\n";
    std::cout << "mu = " << mu << " (from nu = " << nu << ")\n";
    std::cout << "xi = " << xi << " (from V0 = " << V0 << ")\n\n";

    auto t0 = std::chrono::high_resolution_clock::now();

    auto all_stats = parallel_monte_carlo_all_prices(
        n_sims, strikes, barrier_up, barrier_down, alpha, gamma_val, T_steps,
        beta_param, mu, theta, xi, S0, base_seed
    );

    auto t1 = std::chrono::high_resolution_clock::now();
    double time_used = std::chrono::duration<double>(t1 - t0).count();

    std::cout << std::fixed << std::setprecision(4);
    
    std::cout << "\nEuropean Options:\n";
    std::cout << "Strike    Call        95% CI         Put         95% CI\n";
    std::cout << "--------------------------------------------------------\n";
    for (size_t i = 0; i < strikes.size(); ++i) {
        std::cout << std::setw(6) << strikes[i] << "  "
                  << std::setw(8) << all_stats[0][i].mean << "  "
                  << "[" << std::setw(6) << all_stats[0][i].ci_lower << ", " 
                  << std::setw(6) << all_stats[0][i].ci_upper << "]  "
                  << std::setw(8) << all_stats[1][i].mean << "  "
                  << "[" << std::setw(6) << all_stats[1][i].ci_lower << ", " 
                  << std::setw(6) << all_stats[1][i].ci_upper << "]\n";
    }

    std::cout << "\nAsian Options:\n";
    std::cout << "Strike    Call        95% CI         Put         95% CI\n";
    std::cout << "--------------------------------------------------------\n";
    for (size_t i = 0; i < strikes.size(); ++i) {
        std::cout << std::setw(6) << strikes[i] << "  "
                  << std::setw(8) << all_stats[2][i].mean << "  "
                  << "[" << std::setw(6) << all_stats[2][i].ci_lower << ", " 
                  << std::setw(6) << all_stats[2][i].ci_upper << "]  "
                  << std::setw(8) << all_stats[3][i].mean << "  "
                  << "[" << std::setw(6) << all_stats[3][i].ci_lower << ", " 
                  << std::setw(6) << all_stats[3][i].ci_upper << "]\n";
    }

    std::cout << "\nLookback Options:\n";
    std::cout << "Strike    Call        95% CI         Put         95% CI\n";
    std::cout << "--------------------------------------------------------\n";
    for (size_t i = 0; i < strikes.size(); ++i) {
        std::cout << std::setw(6) << strikes[i] << "  "
                  << std::setw(8) << all_stats[4][i].mean << "  "
                  << "[" << std::setw(6) << all_stats[4][i].ci_lower << ", " 
                  << std::setw(6) << all_stats[4][i].ci_upper << "]  "
                  << std::setw(8) << all_stats[5][i].mean << "  "
                  << "[" << std::setw(6) << all_stats[5][i].ci_lower << ", " 
                  << std::setw(6) << all_stats[5][i].ci_upper << "]\n";
    }

    std::cout << "\nBarrier Options (Up-barrier=110, Down-barrier=90):\n";
    std::cout << "Strike    Up-In-Call     95% CI      Down-Out-Put    95% CI\n";
    std::cout << "-------------------------------------------------------------\n";
    for (size_t i = 0; i < strikes.size(); ++i) {
        std::cout << std::setw(6) << strikes[i] << "  "
                  << std::setw(11) << all_stats[6][i].mean << "  "
                  << "[" << std::setw(6) << all_stats[6][i].ci_lower << ", " 
                  << std::setw(6) << all_stats[6][i].ci_upper << "]  "
                  << std::setw(11) << all_stats[7][i].mean << "  "
                  << "[" << std::setw(6) << all_stats[7][i].ci_lower << ", " 
                  << std::setw(6) << all_stats[7][i].ci_upper << "]\n";
    }

    std::cout << "\nTotal time used: " << time_used << " seconds\n";

    return 0;
} 
