#include <vector>
#include <cmath>
#include <random>
#include <iostream>
#include <chrono>
#include <thread>
#include <future>
#include <iomanip>

double f_alpha_one(int t, double alpha) {
    if (t <= 0) return 0.0;
    if (t == 1) {
        return 1.0 - 1.0 / (std::tgamma(1.0 - alpha) * std::pow(2.0, alpha));
    }
    return 1.0 / std::tgamma(1.0 - alpha) * (1.0 / std::pow(t, alpha) - 1.0 / std::pow(t + 1, alpha));
}

double simulate_inar_path(
    double alpha, double gamma_val, int T_steps,
    double beta_param, double mu, double theta, double xi, double S0,
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
    }

    double coef = std::sqrt(theta / 2.0 * (1.0 - a_T) / (std::pow(T_steps, alpha) * mu));
    double diff_pm = N_plus[T_steps] - N_minus[T_steps];
    double P_T = coef * diff_pm - (theta / 2.0) * (1.0 - a_T) / (std::pow(T_steps, alpha) * mu) * N_plus[T_steps];
    double S_T = S0 * std::exp(P_T);
    return S_T;
}

double european_call_payoff(double S_T, double strike) {
    return std::max(S_T - strike, 0.0);
}

std::vector<double> simulate_call_payoffs_once(
    const std::vector<double>& strikes, double alpha, double gamma_val, int T_steps,
    double beta_param, double mu, double theta, double xi, double S0,
    std::mt19937& gen
) {
    double S_T = simulate_inar_path(alpha, gamma_val, T_steps, beta_param, mu, theta, xi, S0, gen);
    std::vector<double> payoffs(strikes.size());
    for (size_t i = 0; i < strikes.size(); ++i) {
        payoffs[i] = european_call_payoff(S_T, strikes[i]);
    }
    return payoffs;
}

std::vector<double> parallel_monte_carlo_call_prices(
    int n_sims, const std::vector<double>& strikes, double alpha, double gamma_val, int T_steps,
    double beta_param, double mu, double theta, double xi, double S0
) {
    const unsigned int num_threads = std::thread::hardware_concurrency();
    std::vector<std::future<std::vector<double>>> futures(num_threads);
    
    int sims_per_thread = n_sims / num_threads;
    
    for (unsigned int i = 0; i < num_threads; ++i) {
        futures[i] = std::async(std::launch::async, [=]() {
            std::random_device rd;
            std::mt19937 gen(rd() + i); 
            std::vector<double> sums(strikes.size(), 0.0);
            
            int local_sims = (i == num_threads - 1) ? 
                            sims_per_thread + (n_sims % num_threads) : 
                            sims_per_thread;
                            
            for (int j = 0; j < local_sims; ++j) {
                std::vector<double> payoffs = simulate_call_payoffs_once(
                    strikes, alpha, gamma_val, T_steps, beta_param, mu, theta, xi, S0, gen
                );
                for (size_t k = 0; k < strikes.size(); ++k) {
                    sums[k] += payoffs[k];
                }
            }
            return sums;
        });
    }
    
    std::vector<double> total_sums(strikes.size(), 0.0);
    for (auto& future : futures) {
        std::vector<double> thread_sums = future.get();
        for (size_t i = 0; i < strikes.size(); ++i) {
            total_sums[i] += thread_sums[i];
        }
    }
    
    for (double& sum : total_sums) {
        sum /= n_sims;
    }
    
    return total_sums;
}

int main() {
    int n_sims = 100000;       
    std::vector<double> strikes;
    for (double K = 80.0; K <= 120.0; K += 5.0) {
        strikes.push_back(K);
    }
    double r = 0.0;            
    int T_steps = 1000;        
    double alpha = 0.62;       
    double gamma_val = 0.1;    
    double beta_param = 27.558;
    double mu = 26.859152;     
    double theta = 0.3156;     
    double xi = 0.1242;        
    double S0 = 100.0;         

    auto t0 = std::chrono::high_resolution_clock::now();

    std::vector<double> call_prices = parallel_monte_carlo_call_prices(
        n_sims, strikes, alpha, gamma_val, T_steps, beta_param, mu, theta, xi, S0
    );

    auto t1 = std::chrono::high_resolution_clock::now();
    double time_used = std::chrono::duration<double>(t1 - t0).count();

    std::cout << std::fixed << std::setprecision(4);
    std::cout << "European Call Option Prices:\n";
    std::cout << "Strike    Price\n";
    std::cout << "----------------\n";
    for (size_t i = 0; i < strikes.size(); ++i) {
        std::cout << std::setw(6) << strikes[i] << "    " << std::setw(8) << call_prices[i] << "\n";
    }
    std::cout << "\nTime used: " << time_used << " seconds.\n";

    return 0;
} 