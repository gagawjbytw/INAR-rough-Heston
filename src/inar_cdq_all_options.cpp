#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <chrono>
#include <thread>
#include <future>
#include <iomanip>
#include <algorithm>
#include <complex>
#include <stdexcept>

// =================================================================
// START: Self-contained FFT implementation (keeping your original style)
// =================================================================
namespace FFT {
    using Complex = std::complex<double>;
    const double PI = acos(-1.0);

    void fft_recursive(std::vector<Complex>& a, bool invert) {
        int n = (int)a.size();
        if (n <= 1) return;

        std::vector<Complex> a0(n / 2), a1(n / 2);
        for (int i = 0; 2 * i < n; i++) {
            a0[i] = a[2 * i];
            a1[i] = a[2 * i + 1];
        }

        fft_recursive(a0, invert);
        fft_recursive(a1, invert);

        double ang = 2 * PI / n * (invert ? -1 : 1);
        Complex w(1), wn(cos(ang), sin(ang));

        for (int i = 0; 2 * i < n; i++) {
            a[i] = a0[i] + w * a1[i];
            a[i + n / 2] = a0[i] - w * a1[i];
            w *= wn;
        }

        if (invert) {
            for (auto& x : a) x /= 2.0;
        }
    }

    // Pack the real sequence into a complex array for the FFT
    std::vector<Complex> fft(const std::vector<double>& input, int fft_size) {
        std::vector<Complex> a(fft_size, 0.0);
        for (size_t i = 0; i < input.size(); ++i) a[i].real(input[i]);
        fft_recursive(a, false);
        return a;
    }

    // In-place inverse transform, returning the real part
    std::vector<double> ifft(std::vector<Complex>& a) {
        fft_recursive(a, true);
        std::vector<double> result(a.size());
        for (size_t i = 0; i < a.size(); ++i) result[i] = a[i].real();
        return result;
    }
}
// =================================================================
// END: Self-contained FFT implementation
// =================================================================

// Helper: next power of two
static inline int next_pow2(int x) { int p = 1; while (p < x) p <<= 1; return p; }

// Parameter solver (same as your original version)
double solve_beta(double rho) {
    double a = 2 * rho * rho - 1;
    double b = 2;
    double c = 2 * rho * rho - 1;
    double disc = b * b - 4 * a * c;
    if (disc < 0) throw std::runtime_error("Invalid rho: discriminant < 0");
    double r1 = (-b + std::sqrt(disc)) / (2 * a);
    double r2 = (-b - std::sqrt(disc)) / (2 * a);
    if (r1 > 1) return r1;
    if (r2 > 1) return r2;
    throw std::runtime_error("No valid beta found (beta should be > 1)");
}
double solve_mu(double nu, double theta, double beta, double gamma) {
    return theta * (1 + beta * beta) / (gamma * nu * nu * (1 + beta) * (1 + beta));
}
double solve_xi(double V0, double theta) { return V0 / theta; }

// Fractional kernel (same as your original version)
double f_alpha_one(int t, double alpha) {
    if (t <= 0) return 0.0;
    // Small safeguard to avoid Gamma(0) blowing up when alpha = 1
    double a = std::min(alpha, 1.0 - 1e-9);
    if (t == 1) {
        return 1.0 - 1.0 / std::tgamma(1.0 - a);
    }
    return (1.0 / std::tgamma(1.0 - a)) *
           (1.0 / std::pow(t - 1.0, a) - 1.0 / std::pow((double)t, a));
}

// ===================== Path simulation result struct (interface unchanged) =====================
struct PathSimResult {
    std::vector<double> path;
    double final_price;
    double avg_price;
    double max_price;
    double min_price;
    bool barrier_hit_110;
    bool barrier_hit_90;
};

// ===================== CDQ online convolution simulation (core modification) =====================
//
// Generate X_t^+, X_t^- and Y_t = X_t^+ + beta X_t^- on the fly.
// Recursively process [L,R]: handle the left half [L,M], convolve its contribution into lam_conv for the right half, then descend into [M+1,R].
//
struct CDQSim {
    // Inputs/constants
    const std::vector<double>& kernel;   // K[1..T_steps]
    const std::vector<double>& lam_base; // baseline[1..T_steps]
    double beta_param;
    std::mt19937& gen;

    // Mutable state
    std::vector<double>& lam_conv; // Accumulated convolution contribution
    std::vector<int>& X_plus;
    std::vector<int>& X_minus;
    std::vector<double>& Y;        // Y_t

    // Convolution: A = Y[L..M], B = K[1..(R-L)] -> conv
    // For the right half t ∈ [M+1,R], add lam_conv[t] += conv[t - L - 1]
    void convolve_left_to_right(int L, int M, int R) {
        int nA = M - L + 1;
        int nB = R - L;         // K[1..R-L]
        if (nA <= 0 || nB <= 0) return;

        int N = next_pow2(nA + nB - 1);
        std::vector<double> A(N, 0.0), B(N, 0.0);
        for (int i = 0; i < nA; ++i) A[i] = Y[L + i];
        for (int j = 0; j < nB; ++j) B[j] = kernel[1 + j];

        auto FA = FFT::fft(A, N);
        auto FB = FFT::fft(B, N);
        std::vector<FFT::Complex> FC(N);
        for (int i = 0; i < N; ++i) FC[i] = FA[i] * FB[i];
        auto conv = FFT::ifft(FC);

        for (int t = M + 1; t <= R; ++t) {
            int k = t - L - 1;                  // take the required segment
            if (k >= 0 && k < (int)conv.size()) lam_conv[t] += conv[k];
        }
    }

    // Leaf: sample X^+, X^- online and write Y
    void leaf(int t) {
        double lam = std::max(0.0, lam_base[t] + lam_conv[t]);
        std::poisson_distribution<int> dist;
        X_plus[t]  = dist(gen, decltype(dist)::param_type(lam));
        X_minus[t] = dist(gen, decltype(dist)::param_type(lam));
        Y[t] = (double)X_plus[t] + beta_param * (double)X_minus[t];
    }

    void run(int L, int R) {
        if (L == R) { leaf(L); return; }
        int M = (L + R) >> 1;
        run(L, M);
        convolve_left_to_right(L, M, R);
        run(M + 1, R);
    }
};

PathSimResult simulate_inar_path_fft( // Keep the original name (implementation switched to CDQ)
    double alpha, double gamma_val, int T_steps,
    double beta_param, double mu, double theta, double xi, double S0,
    double barrier_up, double barrier_down,
    std::mt19937& gen
) {
    // Same unit-maturity discretization as the original: dt = 1/T_steps
    // a_T = 1 - gamma * dt^alpha = 1 - gamma * T_steps^{-alpha}
    // mu_T = mu * dt^{1-alpha}  = mu * T_steps^{alpha-1}
    double a_T  = 1.0 - gamma_val * std::pow((double)T_steps, -alpha);
    double mu_T = mu * std::pow((double)T_steps,  alpha - 1.0);

    // Precompute kernel K[1..T_steps]
    std::vector<double> kernel(T_steps + 1, 0.0);
    for (int t = 1; t <= T_steps; ++t) {
        kernel[t] = a_T * f_alpha_one(t, alpha) / (1.0 + beta_param);
    }

    // Baseline matches your original version
    std::vector<double> phi(T_steps);
    for (int t = 0; t < T_steps; ++t) phi[t] = a_T * f_alpha_one(t + 1, alpha);

    std::vector<double> phi_cum(T_steps, 0.0);
    if (T_steps > 0) phi_cum[0] = phi[0];
    for (int t = 1; t < T_steps; ++t) phi_cum[t] = phi_cum[t - 1] + phi[t];

    std::vector<double> lam_base(T_steps + 1, 0.0);
    for (int t_idx = 1; t_idx <= T_steps; ++t_idx) {
        double sum_phi_prev = (t_idx == 1) ? 0.0 : phi_cum[t_idx - 2];
        lam_base[t_idx] = mu_T
            + xi * mu_T * ( (1.0 / (1.0 - a_T)) * (1.0 - sum_phi_prev) - sum_phi_prev );
    }

    // CDQ online convolution & sampling
    std::vector<double> lam_conv(T_steps + 2, 0.0), Y(T_steps + 2, 0.0);
    std::vector<int>    X_plus(T_steps + 2, 0), X_minus(T_steps + 2, 0);

    CDQSim sim{kernel, lam_base, beta_param, gen, lam_conv, X_plus, X_minus, Y};
    sim.run(1, T_steps);

    // Compute the full price path and statistics (single linear pass)
    std::vector<int> N_plus(T_steps + 1, 0), N_minus(T_steps + 1, 0);
    std::vector<double> price_path(T_steps + 1, S0);

    double sum_prices = S0;
    // denom = mu * T_steps^alpha; coefficients match the original version
    double denom      = mu * std::pow((double)T_steps, alpha);
    double coef       = std::sqrt( (theta / 2.0) * (1.0 - a_T) / denom );
    double drift_coef = (theta / 2.0) * (1.0 - a_T) / denom;

    for (int t_idx = 1; t_idx <= T_steps; ++t_idx) {
        N_plus[t_idx]  = N_plus[t_idx - 1]  + X_plus[t_idx];
        N_minus[t_idx] = N_minus[t_idx - 1] + X_minus[t_idx];

        double diff_pm = (double)(N_plus[t_idx] - N_minus[t_idx]);
        double P_t     = coef * diff_pm - drift_coef * (double)N_plus[t_idx];

        price_path[t_idx] = S0 * std::exp(P_t);
        sum_prices += price_path[t_idx];
    }

    PathSimResult result;
    result.path          = price_path;
    result.final_price   = price_path[T_steps];
    result.avg_price     = sum_prices / (T_steps + 1);
    result.max_price     = *std::max_element(price_path.begin(), price_path.end());
    result.min_price     = *std::min_element(price_path.begin(), price_path.end());
    result.barrier_hit_110 = (result.max_price >= barrier_up);
    result.barrier_hit_90  = (result.min_price <= barrier_down);
    return result;
}

// ===================== Payoff definitions (unchanged) =====================
double european_call_payoff(double S_T, double strike) { return std::max(S_T - strike, 0.0); }
double european_put_payoff(double S_T, double strike)  { return std::max(strike - S_T, 0.0); }
double asian_call_payoff(double avg_price, double strike) { return std::max(avg_price - strike, 0.0); }
double asian_put_payoff(double avg_price, double strike)  { return std::max(strike - avg_price, 0.0); }
double lookback_call_payoff(double max_price, double strike) { return std::max(max_price - strike, 0.0); }
double lookback_put_payoff(double min_price, double strike)  { return std::max(strike - min_price, 0.0); }
double barrier_up_in_call_payoff(const PathSimResult& result, double strike, double barrier_up) {
    return result.barrier_hit_110 ? std::max(result.final_price - strike, 0.0) : 0.0;
}
double barrier_down_out_put_payoff(const PathSimResult& result, double strike, double barrier_down) {
    return !result.barrier_hit_90 ? std::max(strike - result.final_price, 0.0) : 0.0;
}

// ===================== Statistics helper (unchanged) =====================
struct StatResult { double mean, ci_lower, ci_upper; };
StatResult calculate_stats(const std::vector<double>& prices) {
    double n = (double)prices.size();
    double sum = 0.0;
    for (double p : prices) sum += p;
    double mean = sum / n;
    double var  = 0.0;
    for (double p : prices) var += (p - mean) * (p - mean);
    var /= (n - 1.0);
    double se = std::sqrt(var / n);
    return { mean, mean - 1.96 * se, mean + 1.96 * se };
}

struct SimResults {
    std::vector<std::vector<double>> prices[8];
    SimResults(size_t n_strikes, size_t n_sims) {
        for (int i = 0; i < 8; ++i) prices[i].resize(n_strikes, std::vector<double>(n_sims));
    }
};

// ===================== Parallel Monte Carlo (interface unchanged) =====================
std::vector<std::vector<StatResult>> parallel_monte_carlo_all_prices(
    int n_sims, const std::vector<double>& strikes,
    double barrier_up, double barrier_down,
    double alpha, double gamma_val, int T_steps,
    double beta_param, double mu, double theta, double xi, double S0,
    unsigned int base_seed = 123456789
) {
    const unsigned int hc = std::thread::hardware_concurrency();
    const unsigned int num_threads = std::max(1u, hc);

    SimResults results(strikes.size(), n_sims);
    std::vector<std::future<void>> futures;

    int sims_per_thread = (n_sims + (int)num_threads - 1) / (int)num_threads;
    for (unsigned int t = 0; t < num_threads; ++t) {
        int start_sim = (int)t * sims_per_thread;
        int end_sim   = std::min(n_sims, (int)((t + 1) * sims_per_thread));
        if (start_sim >= end_sim) break;

        futures.push_back(std::async(std::launch::async, [&, start_sim, end_sim, t, base_seed]() {
            std::mt19937 gen(base_seed + (unsigned int)t * 101);
            for (int sim = start_sim; sim < end_sim; ++sim) {
                PathSimResult pr = simulate_inar_path_fft(
                    alpha, gamma_val, T_steps, beta_param, mu, theta, xi, S0,
                    barrier_up, barrier_down, gen
                );
                for (size_t k = 0; k < strikes.size(); ++k) {
                    double K = strikes[k];
                    results.prices[0][k][sim] = european_call_payoff(pr.final_price, K);
                    results.prices[1][k][sim] = european_put_payoff (pr.final_price, K);
                    results.prices[2][k][sim] = asian_call_payoff   (pr.avg_price, K);
                    results.prices[3][k][sim] = asian_put_payoff    (pr.avg_price, K);
                    results.prices[4][k][sim] = lookback_call_payoff(pr.max_price, K);
                    results.prices[5][k][sim] = lookback_put_payoff (pr.min_price, K);
                    results.prices[6][k][sim] = barrier_up_in_call_payoff(pr, K, barrier_up);
                    results.prices[7][k][sim] = barrier_down_out_put_payoff(pr, K, barrier_down);
                }
            }
        }));
    }
    for (auto& f : futures) f.get();

    std::vector<std::vector<StatResult>> all_stats(8, std::vector<StatResult>(strikes.size()));
    for (size_t k = 0; k < strikes.size(); ++k) {
        for (int i = 0; i < 8; ++i) all_stats[i][k] = calculate_stats(results.prices[i][k]);
    }
    return all_stats;
}

int main() {
    unsigned int base_seed = 123456;

    int n_sims = 500000; // start with a smaller number first, e.g., 5000
    std::vector<double> strikes;
    for (double K = 80; K <= 120; K += 10) strikes.push_back(K);

    double barrier_up = 110.0, barrier_down = 90.0;
    int T_steps = 320;
    double alpha = 0.62, gamma_val = 0.1, theta = 0.3156, S0 = 100.0;
    double rho = -0.681, nu = 0.331, V0 = 0.0392;

    double beta_param = solve_beta(rho);
    double mu  = solve_mu(nu, theta, beta_param, gamma_val);
    double xi  = solve_xi(V0, theta);

    std::cout << "Derived Parameters:\n";
    std::cout << "beta = " << beta_param << ", mu = " << mu << ", xi = " << xi << "\n\n";

    auto t0 = std::chrono::high_resolution_clock::now();
    auto all_stats = parallel_monte_carlo_all_prices(
        n_sims, strikes, barrier_up, barrier_down,
        alpha, gamma_val, T_steps, beta_param, mu, theta, xi, S0, base_seed
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
                  << "[" << std::setw(6) << all_stats[0][i].ci_lower << ", " << std::setw(6) << all_stats[0][i].ci_upper << "]  "
                  << std::setw(8) << all_stats[1][i].mean << "  "
                  << "[" << std::setw(6) << all_stats[1][i].ci_lower << ", " << std::setw(6) << all_stats[1][i].ci_upper << "]\n";
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

    std::cout << "\nTotal time used: " << time_used << " seconds (CDQ-FFT version)\n";
    return 0;
}