#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <chrono>
#include <thread>
#include <future>
#include <iomanip>
#include <algorithm>
#include <stdexcept>
#include <numeric>
#include <complex>

#ifndef M_SQRT1_2
#define M_SQRT1_2 0.70710678118654752440
#endif

// ===================================================================
// 1) Parameter calculations
// ===================================================================
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

double f_alpha_one(int t, double alpha) {
    if (t <= 0) return 0.0;
    double a = std::min(alpha, 1.0 - 1e-9); // avoid Γ(1-α) numeric instability
    if (t == 1) {
        // phi_1 = 1 - 1/Gamma(1-alpha)
        return 1.0 - 1.0 / std::tgamma(1.0 - a);
    }
    // phi_n = [1/Gamma(1-alpha)] * (1/(n-1)^alpha - 1/n^alpha), n>=2
    return (1.0 / std::tgamma(1.0 - a)) *
           (1.0 / std::pow(t - 1.0, a) - 1.0 / std::pow((double)t, a));
}

// ===================================================================
// 2) Iterative FFT (Cooley–Tukey), practical and stable
// ===================================================================
namespace FFT {
    using Complex = std::complex<double>;
    const double PI = std::acos(-1.0);

    // Iterative bit-reversal FFT
    void fft(std::vector<Complex>& a, bool invert) {
        int n = (int)a.size();
        // bit-reversal
        for (int i = 1, j = 0; i < n; ++i) {
            int bit = n >> 1;
            for (; j & bit; bit >>= 1) j ^= bit;
            j ^= bit;
            if (i < j) std::swap(a[i], a[j]);
        }
        for (int len = 2; len <= n; len <<= 1) {
            double ang = 2 * PI / len * (invert ? -1.0 : 1.0);
            Complex wlen(std::cos(ang), std::sin(ang));
            for (int i = 0; i < n; i += len) {
                Complex w(1.0, 0.0);
                for (int j = 0; j < len / 2; ++j) {
                    Complex u = a[i + j];
                    Complex v = a[i + j + len / 2] * w;
                    a[i + j] = u + v;
                    a[i + j + len / 2] = u - v;
                    w *= wlen;
                }
            }
        }
        if (invert) {
            for (int i = 0; i < n; ++i) a[i] /= n;
        }
    }

    inline int next_pow2(int x) {
        int p = 1;
        while (p < x) p <<= 1;
        return p;
    }
}

// ===================================================================
// 3) Single-path simulation (CDQ online convolution)
// ===================================================================
struct EuropeanPathResult {
    double final_price;
};

struct CDQSimulator {
    // References to external data
    const std::vector<double>& kernel;     // 1..T_steps, kernel[0] unused
    const std::vector<double>& lam_base;   // 1..T_steps
    std::vector<double>& lam_conv;         // 1..T_steps, accumulated convolution contribution
    std::vector<double>& Y;                // 1..T_steps, Y_t = X^+ + beta X^-
    double beta_param;

    // Counters (used for terminal scaling)
    long long &N_plus_T, &N_minus_T;

    // Random source
    std::mt19937& gen;

    // Convolution: after convolving A = Y[L..M] with B = kernel[1..(R-L)],
    // add the entries with index k = t-L-1 (t ∈ [M+1,R]) to lam_conv[t]
    void convolve_left_to_right(int L, int M, int R) {
        int nA = M - L + 1;
        int nB = R - L;            // corresponds to kernel[1..R-L]
        if (nA <= 0 || nB <= 0) return;

        int n = FFT::next_pow2(nA + nB - 1);
        std::vector<FFT::Complex> fa(n), fb(n);

        for (int i = 0; i < nA; ++i) fa[i] = FFT::Complex(Y[L + i], 0.0);
        for (int j = 0; j < nB; ++j) fb[j] = FFT::Complex(kernel[1 + j], 0.0);

        FFT::fft(fa, false);
        FFT::fft(fb, false);
        for (int i = 0; i < n; ++i) fa[i] *= fb[i];
        FFT::fft(fa, true);

        // Keep only what the right half needs: t ∈ [M+1, R] ⇒ k = t-L-1 ∈ [M-L, R-L-1]
        for (int t = M + 1; t <= R; ++t) {
            int k = t - L - 1;
            if (k >= 0 && k < nA + nB - 1) {
                lam_conv[t] += fa[k].real();
            }
        }
    }

    // Leaf: draw Y[L] on the fly
    void simulate_leaf(int L) {
        double lam = std::max(0.0, lam_base[L] + lam_conv[L]);
        // Use param_type to avoid rebuilding the distribution
        std::poisson_distribution<int> dist;
        int x_plus  = dist(gen, decltype(dist)::param_type(lam));
        int x_minus = dist(gen, decltype(dist)::param_type(lam));
        N_plus_T  += x_plus;
        N_minus_T += x_minus;
        Y[L] = double(x_plus) + beta_param * double(x_minus);
    }

    void run(int L, int R) {
        if (L == R) {
            simulate_leaf(L);
            return;
        }
        int M = (L + R) >> 1;
        // Solve the left half first to obtain Y[L..M]
        run(L, M);
        // Single FFT to inject the left-half contribution into lam_conv
        convolve_left_to_right(L, M, R);
        // Then recurse into the right half
        run(M + 1, R);
    }
};

EuropeanPathResult simulate_path_final_price_fft(
    double T,
    double alpha, double gamma_val, int T_steps,
    double beta_param, double mu, double theta, double xi, double S0,
    std::mt19937& gen
) {
    // --------- Δt and discretization ----------
    double dt       = T / static_cast<double>(T_steps);
    double tau_disc = 1.0 / dt;
    double a_step   = 1.0 - gamma_val * std::pow(dt, alpha);
    double mu_step  = mu * std::pow(dt, 1.0 - alpha);

    // --------- Precompute kernel K (1..T_steps) ----------
    std::vector<double> kernel(T_steps + 1, 0.0);
    for (int t = 1; t <= T_steps; ++t) {
        kernel[t] = a_step * f_alpha_one(t, alpha) / (1.0 + beta_param);
    }

    // --------- Precompute baseline (same as the original code) ----------
    std::vector<double> phi(T_steps);
    for (int t = 0; t < T_steps; ++t) {
        phi[t] = a_step * f_alpha_one(t + 1, alpha);
    }
    std::vector<double> phi_cum(T_steps, 0.0);
    if (T_steps > 0) phi_cum[0] = phi[0];
    for (int t = 1; t < T_steps; ++t) phi_cum[t] = phi_cum[t - 1] + phi[t];

    std::vector<double> lam_base(T_steps + 1, 0.0);
    for (int t_idx = 1; t_idx <= T_steps; ++t_idx) {
        double sum_phi_prev = (t_idx == 1) ? 0.0 : phi_cum[t_idx - 2];
        lam_base[t_idx] = mu_step
            + xi * mu_step * ( (1.0 / (1.0 - a_step)) * (1.0 - sum_phi_prev) - sum_phi_prev );
    }

    // --------- CDQ online convolution + sampling ----------
    std::vector<double> lam_conv(T_steps + 2, 0.0), Y(T_steps + 2, 0.0);
    long long N_plus_T = 0, N_minus_T = 0;

    CDQSimulator sim{
        kernel, lam_base, lam_conv, Y, beta_param, N_plus_T, N_minus_T, gen
    };
    sim.run(1, T_steps);

    // --------- Terminal scaling ----------
    double denom      = mu * std::pow(tau_disc, alpha);
    double coef       = std::sqrt( (theta / 2.0) * (1.0 - a_step) / denom );
    double drift_coef = (theta / 2.0) * (1.0 - a_step) / denom;

    double diff_pm = static_cast<double>(N_plus_T - N_minus_T);
    double P_T     = coef * diff_pm - drift_coef * static_cast<double>(N_plus_T);

    return { S0 * std::exp(P_T) };
}

// ===================================================================
// 4) Black-Scholes and implied volatility
// ===================================================================
static inline double norm_cdf(double x) { return 0.5 * std::erfc(-x * M_SQRT1_2); }

double black_scholes_call(double S, double K, double T, double r, double sigma) {
    if (T <= 1e-12 || sigma <= 1e-12)
        return std::max(S - K * std::exp(-r * T), 0.0);
    double vsqrt = sigma * std::sqrt(T);
    double d1 = (std::log(S / K) + (r + 0.5 * sigma * sigma) * T) / vsqrt;
    double d2 = d1 - vsqrt;
    return S * norm_cdf(d1) - K * std::exp(-r * T) * norm_cdf(d2);
}

inline double clamp_call_price(double price, double S, double K, double T, double r) {
    const double df = std::exp(-r * T);
    const double fwd = S / df;
    const double lower = df * std::max(fwd - K, 0.0); // = max(S - K*df, 0)
    const double upper = S;                           // Call upper bound
    const double eps = 1e-12 * S + 1e-14;            // Tiny buffer

    if (price <= lower) return lower + eps;
    if (price >= upper) return upper - eps;
    return price;
}

double black_scholes_put(double S, double K, double T, double r, double sigma) {
    if (T <= 1e-12 || sigma <= 1e-12) {
        return std::max(K * std::exp(-r * T) - S, 0.0);
    }
    double vsqrt = sigma * std::sqrt(T);
    double d1 = (std::log(S / K) + (r + 0.5 * sigma * sigma) * T) / vsqrt;
    double d2 = d1 - vsqrt;
    return K * std::exp(-r * T) * norm_cdf(-d2) - S * norm_cdf(-d1);
}

// === Added: put price clamp (avoid bounds issues / instability) ===
inline double clamp_put_price(double price, double S, double K, double T, double r) {
    const double df = std::exp(-r * T);
    const double lower = std::max(K * df - S, 0.0);
    const double upper = K * df;
    const double eps = 1e-12 * K + 1e-14;
    if (price <= lower) return lower + eps;
    if (price >= upper) return upper - eps;
    return price;
}

// === Renamed the original implied_volatility to the call inverter ===
double implied_volatility_call(double price, double S, double K, double T, double r) {
    if (T <= 0.0) return 0.0;
    price = clamp_call_price(price, S, K, T, r);
    auto C = [&](double sig){ return black_scholes_call(S, K, T, r, sig); };
    double lo = 0.0, hi = 0.2;
    while (C(hi) < price && hi < 100.0) hi *= 2.0;
    for (int it = 0; it < 100; ++it) {
        double mid = 0.5 * (lo + hi);
        if (C(mid) < price) lo = mid; else hi = mid;
    }
    return 0.5 * (lo + hi);
}

// === Added: put inverter using the same bisection skeleton as the call inverter ===
double implied_volatility_put(double price, double S, double K, double T, double r) {
    if (T <= 0.0) return 0.0;
    price = clamp_put_price(price, S, K, T, r);
    auto P = [&](double sig){ return black_scholes_put(S, K, T, r, sig); };
    double lo = 0.0, hi = 0.2;
    while (P(hi) < price && hi < 100.0) hi *= 2.0;
    for (int it = 0; it < 100; ++it) {
        double mid = 0.5 * (lo + hi);
        if (P(mid) < price) lo = mid; else hi = mid;
    }
    return 0.5 * (lo + hi);
}

// === Added: smart chooser — use the put inverter when the call is ITM ===
// Logic: use the forward F = S * exp(rT) for moneyness; if K < F the call is ITM so switch to the put inverter.
// Obtain the put price from parity: P = C - S + K e^{-rT}
double implied_volatility_smart(double call_price, double S, double K, double T, double r) {
    const double F = S * std::exp(r * T);
    if (K < F) {
        // Call ITM: convert to a put then invert
        double put_price = call_price - S + K * std::exp(-r * T);
        return implied_volatility_put(put_price, S, K, T, r);
    } else {
        // Call OTM or near ATM: invert the call directly
        return implied_volatility_call(call_price, S, K, T, r);
    }
}

inline void clamp_call_strip(std::vector<double>& prices,
                             const std::vector<double>& strikes,
                             double S, double T, double r) {
    const double df = std::exp(-r * T);
    const double eps = 1e-12 * S + 1e-14;
    for (std::size_t i = 0; i < strikes.size(); ++i) {
        double lower = std::max(S - strikes[i] * df, 0.0);
        double upper = S;
        prices[i] = std::min(upper - eps, std::max(lower + eps, prices[i]));
    }
}

void enforce_monotone_convex(std::vector<double>& prices,
                             const std::vector<double>& strikes,
                             double S, double T, double r) {
    const std::size_t n = strikes.size();
    if (n < 3) {
        clamp_call_strip(prices, strikes, S, T, r);
        return;
    }

    clamp_call_strip(prices, strikes, S, T, r);
    const double df = std::exp(-r * T);

    std::vector<double> h(n - 1), slopes(n - 1);
    for (std::size_t i = 0; i + 1 < n; ++i) {
        h[i] = strikes[i + 1] - strikes[i];
        slopes[i] = (prices[i + 1] - prices[i]) / h[i];
        slopes[i] = std::min(0.0, std::max(-df, slopes[i]));
    }

    struct Block { std::size_t l, r; double w, sw; };
    auto avg = [](const Block& b) { return b.sw / std::max(1e-18, b.w); };

    std::vector<Block> stack;
    stack.reserve(n);
    for (std::size_t i = 0; i < n - 1; ++i) {
        Block blk{i, i, h[i], h[i] * slopes[i]};
        stack.push_back(blk);
        while (stack.size() >= 2 && avg(stack[stack.size() - 2]) > avg(stack.back())) {
            Block tail = stack.back();
            stack.pop_back();
            Block& head = stack.back();
            head.r = tail.r;
            head.w += tail.w;
            head.sw += tail.sw;
            double a = avg(head);
            a = std::min(0.0, std::max(-df, a));
            head.sw = a * head.w;
        }
    }

    std::vector<double> adjusted_slopes(n - 1);
    for (const auto& blk : stack) {
        double a = avg(blk);
        for (std::size_t i = blk.l; i <= blk.r; ++i) adjusted_slopes[i] = a;
    }

    std::vector<double> deltas(n, 0.0);
    for (std::size_t i = 1; i < n; ++i) deltas[i] = deltas[i - 1] + adjusted_slopes[i - 1] * h[i - 1];

    double lower_bound = -1e300;
    double upper_bound =  1e300;
    for (std::size_t i = 0; i < n; ++i) {
        double lower = std::max(S - strikes[i] * df, 0.0) + 1e-12 * S;
        double upper = S - 1e-12 * S;
        lower_bound = std::max(lower_bound, lower - deltas[i]);
        upper_bound = std::min(upper_bound, upper - deltas[i]);
    }

    double base = std::min(upper_bound, std::max(lower_bound, prices.front()));
    for (std::size_t i = 0; i < n; ++i) prices[i] = base + deltas[i];
}

// ===================================================================
// 5) Monte Carlo (path-dimension parallelism)
// ===================================================================
double price_european_call_mc(
    double T, double K, int n_sims,
    double alpha, double gamma_val,
    double beta_param, double mu, double theta, double xi, double S0,
    unsigned int base_seed = 12345
) {
    const unsigned int num_threads = std::max(1u, std::thread::hardware_concurrency());
    std::vector<double> payoffs(n_sims, 0.0);
    std::vector<std::future<void>> futures;
    int T_steps = std::max(10, (int)std::round(320.0 * T));   // adjust the number of timesteps as needed
    int sims_per_thread = (n_sims + (int)num_threads - 1) / (int)num_threads;

    for (unsigned int t = 0; t < num_threads; ++t) {
        int start_sim = (int)t * sims_per_thread;
        int end_sim = std::min(n_sims, (int)((t + 1) * sims_per_thread));
        if (start_sim >= end_sim) break;
        futures.push_back(std::async(std::launch::async, [=, &payoffs]() {
            std::mt19937 gen(base_seed + (unsigned int)t * 101);
            for (int i = start_sim; i < end_sim; ++i) {
                auto result = simulate_path_final_price_fft(
                    T, alpha, gamma_val, T_steps,
                    beta_param, mu, theta, xi, S0, gen
                );
                payoffs[i] = std::max(result.final_price - K, 0.0);
            }
        }));
    }
    for (auto& f : futures) f.get();
    double sum = std::accumulate(payoffs.begin(), payoffs.end(), 0.0);
    return sum / (double)n_sims;
}

std::vector<double> price_call_strip_mc(
    double T,
    const std::vector<double>& strikes,
    int n_sims,
    double alpha, double gamma_val,
    double beta_param, double mu, double theta, double xi, double S0,
    unsigned int base_seed = 12345
) {
    const unsigned int num_threads = std::max(1u, std::thread::hardware_concurrency());
    const int sims_per_thread = (n_sims + static_cast<int>(num_threads) - 1) / static_cast<int>(num_threads);

    const double dt_target = 1e-3;
    const int T_steps = std::max(static_cast<int>(std::ceil(T / dt_target)), 50);

    std::vector<std::vector<double>> locals(num_threads, std::vector<double>(strikes.size(), 0.0));
    std::vector<std::future<void>> futures;
    futures.reserve(num_threads);

    for (unsigned int t = 0; t < num_threads; ++t) {
        int start_sim = static_cast<int>(t) * sims_per_thread;
        int end_sim = std::min(n_sims, static_cast<int>((t + 1) * sims_per_thread));
        if (start_sim >= end_sim) break;

        futures.push_back(std::async(std::launch::async, [=, &locals]() {
            std::mt19937 gen(base_seed + static_cast<unsigned int>(t) * 101);
            auto& acc = locals[t];
            for (int i = start_sim; i < end_sim; ++i) {
                auto result = simulate_path_final_price_fft(
                    T, alpha, gamma_val, T_steps,
                    beta_param, mu, theta, xi, S0, gen
                );
                const double ST = result.final_price;
                for (std::size_t j = 0; j < strikes.size(); ++j) {
                    acc[j] += std::max(ST - strikes[j], 0.0);
                }
            }
        }));
    }
    for (auto& f : futures) f.get();

    std::vector<double> sums(strikes.size(), 0.0);
    for (unsigned int t = 0; t < num_threads; ++t)
        for (std::size_t j = 0; j < strikes.size(); ++j)
            sums[j] += locals[t][j];

    std::vector<double> prices(strikes.size(), 0.0);
    for (std::size_t j = 0; j < strikes.size(); ++j)
        prices[j] = sums[j] / static_cast<double>(n_sims);

    return prices;
}

struct StripPrices {
    std::vector<double> calls;
    std::vector<double> puts;
};

StripPrices price_call_put_strip_mc(
    double T,
    const std::vector<double>& strikes,
    int n_sims,
    double alpha, double gamma_val,
    double beta_param, double mu, double theta, double xi, double S0,
    unsigned int base_seed = 12345,
    double r = 0.0
) {
    const unsigned int num_threads = std::max(1u, std::thread::hardware_concurrency());
    const int sims_per_thread = (n_sims + static_cast<int>(num_threads) - 1) / static_cast<int>(num_threads);

    // Same timestep choice as your original function
    const double dt_target = 1e-3;
    const int T_steps = std::max(static_cast<int>(std::ceil(T / dt_target)), 50);

    std::vector<std::vector<double>> locals_call(num_threads, std::vector<double>(strikes.size(), 0.0));
    std::vector<std::vector<double>> locals_put (num_threads, std::vector<double>(strikes.size(), 0.0));

    std::vector<std::future<void>> futures;
    futures.reserve(num_threads);

    for (unsigned int t = 0; t < num_threads; ++t) {
        int start_sim = static_cast<int>(t) * sims_per_thread;
        int end_sim = std::min(n_sims, static_cast<int>((t + 1) * sims_per_thread));
        if (start_sim >= end_sim) break;

        futures.push_back(std::async(std::launch::async, [=, &locals_call, &locals_put]() {
            std::mt19937 gen(base_seed + static_cast<unsigned int>(t) * 101);
            auto& acc_c = locals_call[t];
            auto& acc_p = locals_put[t];

            for (int i = start_sim; i < end_sim; ++i) {
                auto result = simulate_path_final_price_fft(
                    T, alpha, gamma_val, T_steps,
                    beta_param, mu, theta, xi, S0, gen
                );
                const double ST = result.final_price;
                for (std::size_t j = 0; j < strikes.size(); ++j) {
                    const double K = strikes[j];
                    acc_c[j] += std::max(ST - K, 0.0);
                    acc_p[j] += std::max(K - ST, 0.0);
                }
            }
        }));
    }
    for (auto& f : futures) f.get();

    std::vector<double> sums_c(strikes.size(), 0.0), sums_p(strikes.size(), 0.0);
    for (unsigned int t = 0; t < num_threads; ++t) {
        for (std::size_t j = 0; j < strikes.size(); ++j) {
            sums_c[j] += locals_call[t][j];
            sums_p[j] += locals_put [t][j];
        }
    }

    const double df = std::exp(-r * T);  // Discount factor to stay consistent with BS/bounds
    StripPrices out;
    out.calls.resize(strikes.size());
    out.puts .resize(strikes.size());
    for (std::size_t j = 0; j < strikes.size(); ++j) {
        out.calls[j] = (sums_c[j] / static_cast<double>(n_sims)) * df;
        out.puts [j] = (sums_p[j] / static_cast<double>(n_sims)) * df;
    }
    return out;
}

inline void clamp_put_strip(std::vector<double>& prices,
                            const std::vector<double>& strikes,
                            double S, double T, double r) {
    const double df = std::exp(-r * T);
    const double eps = 1e-12 * std::max(S, 1.0) + 1e-14;
    for (std::size_t i = 0; i < strikes.size(); ++i) {
        double lower = std::max(strikes[i] * df - S, 0.0);
        double upper = strikes[i] * df;
        prices[i] = std::min(upper - eps, std::max(lower + eps, prices[i]));
    }
}

void enforce_monotone_convex_put(std::vector<double>& prices,
                                 const std::vector<double>& strikes,
                                 double S, double T, double r) {
    const std::size_t n = strikes.size();
    if (n < 3) {
        clamp_put_strip(prices, strikes, S, T, r);
        return;
    }

    clamp_put_strip(prices, strikes, S, T, r);
    const double df = std::exp(-r * T);

    std::vector<double> h(n - 1), slopes(n - 1);
    for (std::size_t i = 0; i + 1 < n; ++i) {
        h[i] = strikes[i + 1] - strikes[i];
        slopes[i] = (prices[i + 1] - prices[i]) / h[i];
        // Put slopes should stay within [0, df]
        slopes[i] = std::max(0.0, std::min(df, slopes[i]));
    }

    struct Block { std::size_t l, r; double w, sw; };
    auto avg = [](const Block& b) {
        return b.sw / std::max(1e-18, b.w);
    };

    std::vector<Block> stack;
    stack.reserve(n);
    for (std::size_t i = 0; i < n - 1; ++i) {
        Block blk{i, i, h[i], h[i] * slopes[i]};
        stack.push_back(blk);
        // Same as the call version: use PAV so slopes are non-decreasing ⇒ price is convex in K
        while (stack.size() >= 2 && avg(stack[stack.size() - 2]) > avg(stack.back())) {
            Block tail = stack.back();
            stack.pop_back();
            Block& head = stack.back();
            head.r  = tail.r;
            head.w += tail.w;
            head.sw += tail.sw;
            double a = avg(head);
            a = std::max(0.0, std::min(df, a)); // keep it within [0, df]
            head.sw = a * head.w;
        }
    }

    std::vector<double> adjusted(n - 1);
    for (const auto& blk : stack) {
        double a = avg(blk);
        for (std::size_t i = blk.l; i <= blk.r; ++i) adjusted[i] = a;
    }

    std::vector<double> deltas(n, 0.0);
    for (std::size_t i = 1; i < n; ++i) deltas[i] = deltas[i - 1] + adjusted[i - 1] * h[i - 1];

    double lower_bound = -1e300;
    double upper_bound =  1e300;
    for (std::size_t i = 0; i < n; ++i) {
        double lower = std::max(strikes[i] * df - S, 0.0) + 1e-12 * std::max(S, 1.0);
        double upper = strikes[i] * df - 1e-12 * std::max(S, 1.0);
        lower_bound = std::max(lower_bound, lower - deltas[i]);
        upper_bound = std::min(upper_bound, upper - deltas[i]);
    }

    double base = std::min(upper_bound, std::max(lower_bound, prices.front()));
    for (std::size_t i = 0; i < n; ++i) prices[i] = base + deltas[i];
}

double implied_volatility_smart_from_both(double call_price, double put_price,
                                          double S, double K, double T, double r) {
    const double F = S * std::exp(r * T);
    if (K < F) {
        // Call ITM: invert via the put price (more stable)
        return implied_volatility_put(put_price, S, K, T, r);
    } else {
        // Call OTM or near ATM: still use the call inverter
        return implied_volatility_call(call_price, S, K, T, r);
    }
}
// ===================================================================
// 6) Main routine: output the IV surface (CSV)
// ===================================================================
int main() {
    auto t_start_total = std::chrono::high_resolution_clock::now();

    // ---- Model parameters ----
    double alpha_rough = 0.62;
    double alpha_classical = 1.0;
    double S0 = 100.0;
    double gamma_val = 0.1;
    double theta = 0.3156;
    double rho = -0.681;
    //double rho = 0.8;
    double nu = 0.331;
    double V0 = 0.0392;

    // ---- Derived parameters ----
    double beta_param = solve_beta(rho);
    double mu = solve_mu(nu, theta, beta_param, gamma_val);
    double xi = solve_xi(V0, theta);

    // ---- Grid and Monte Carlo settings ----
    int n_sims = 1000000;  // start smaller first, e.g., 2000
    double r = 0.0;

    // 1) Maturities
    std::vector<double> maturities = {
        1.0/12.0, 2.0/12.0, 3.0/12.0, 4.0/12.0, 5.0/12.0, 6.0/12.0,
        //3.0/12.0, 4.0/12.0, 5.0/12.0, 6.0/12.0,
        7.0/12.0, 8.0/12.0, 9.0/12.0, 10.0/12.0, 11.0/12.0, 12.0/12.0
    };

    //std::vector<double> maturities = {
    //    1.0/12.0
    //};

    // 2) Log moneyness
    std::vector<double> log_moneyness;
    for (int i = -10; i <= 10; ++i) log_moneyness.push_back(i * 0.02); // -20% ~ +20%

    std::cout << "Starting Implied Volatility Surface generation...\n";
    std::cout << "Maturity,LogMoneyness,IV_Rough,IV_Classical\n";
    std::cout << std::fixed << std::setprecision(6);

    for (double T : maturities) {
        std::vector<double> strikes;
        strikes.reserve(log_moneyness.size());
        for (double k : log_moneyness) strikes.push_back(S0 * std::exp(k));

        auto strip_start = std::chrono::high_resolution_clock::now();
        // Previously:
        // auto prices_r = price_call_strip_mc(...);
        // auto prices_c = price_call_strip_mc(...);
        // enforce_monotone_convex(prices_r, strikes, S0, T, r);
        // enforce_monotone_convex(prices_c, strikes, S0, T, r);

        // Changed to:
        auto strip_r = price_call_put_strip_mc(T, strikes, n_sims,
                                            alpha_rough, gamma_val,
                                            beta_param, mu, theta, xi, S0,
                                            12345u, r);

        auto strip_c = price_call_put_strip_mc(T, strikes, n_sims,
                                            alpha_classical, gamma_val,
                                            beta_param, mu, theta, xi, S0,
                                            54321u, r);

        // Apply the monotone-convex projection separately to call/put (put projection optional but recommended)
        enforce_monotone_convex(strip_r.calls, strikes, S0, T, r);
        enforce_monotone_convex_put(strip_r.puts,  strikes, S0, T, r);

        enforce_monotone_convex(strip_c.calls, strikes, S0, T, r);
        enforce_monotone_convex_put(strip_c.puts,  strikes, S0, T, r);
        auto strip_end = std::chrono::high_resolution_clock::now();

        for (std::size_t idx = 0; idx < log_moneyness.size(); ++idx) {
            double k  = log_moneyness[idx];
            double K  = strikes[idx];

            double iv_r = implied_volatility_smart_from_both(
                strip_r.calls[idx], strip_r.puts[idx], S0, K, T, r);
            double iv_c = implied_volatility_smart_from_both(
                strip_c.calls[idx], strip_c.puts[idx], S0, K, T, r);

            std::cout << T << "," << k << "," << iv_r << "," << iv_c << std::endl;
        }
        std::chrono::duration<double> strip_dt = strip_end - strip_start;
        std::cerr << "Finished all strikes for maturity T=" << T
                  << " (pricing strip in " << std::fixed << std::setprecision(2)
                  << strip_dt.count() << " s).\n";
    }

    auto t_end_total = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> dt_total = t_end_total - t_start_total;
    std::cerr << "\nTotal computation finished in " << dt_total.count() << " s.\n";
    return 0;
}
