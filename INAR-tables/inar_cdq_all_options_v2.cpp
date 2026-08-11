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

    struct Plan {
        int n = 0;
        std::vector<int> rev;
        std::vector<Complex> wlen_fwd;
        std::vector<Complex> wlen_inv;
    };

    Plan make_plan(int n) {
        Plan plan;
        plan.n = n;
        plan.rev.resize(n, 0);

        int levels = 0;
        while ((1 << levels) < n) {
            ++levels;
        }
        for (int i = 1; i < n; ++i) {
            plan.rev[i] = (plan.rev[i >> 1] >> 1) | ((i & 1) << (levels - 1));
        }

        plan.wlen_fwd.reserve(levels + 1);
        plan.wlen_inv.reserve(levels + 1);
        plan.wlen_fwd.push_back(Complex(1.0, 0.0));
        plan.wlen_inv.push_back(Complex(1.0, 0.0));
        for (int len = 2; len <= n; len <<= 1) {
            const double ang = 2.0 * PI / len;
            plan.wlen_fwd.emplace_back(std::cos(ang), std::sin(ang));
            plan.wlen_inv.emplace_back(std::cos(-ang), std::sin(-ang));
        }
        return plan;
    }

    void fft_inplace(Complex* a, const Plan& plan, bool invert) {
        const int n = plan.n;
        for (int i = 0; i < n; ++i) {
            const int j = plan.rev[i];
            if (i < j) std::swap(a[i], a[j]);
        }

        int stage = 1;
        for (int len = 2; len <= n; len <<= 1, ++stage) {
            const Complex wlen = invert ? plan.wlen_inv[stage] : plan.wlen_fwd[stage];
            for (int i = 0; i < n; i += len) {
                Complex w(1.0, 0.0);
                const int half = len >> 1;
                for (int j = 0; j < half; ++j) {
                    const Complex u = a[i + j];
                    const Complex v = a[i + j + half] * w;
                    a[i + j] = u + v;
                    a[i + j + half] = u - v;
                    w *= wlen;
                }
            }
        }

        if (invert) {
            const double inv_n = 1.0 / static_cast<double>(n);
            for (int i = 0; i < n; ++i) a[i] *= inv_n;
        }
    }

    void fft_inplace(Complex* a, int n, bool invert) {
        const Plan plan = make_plan(n);
        fft_inplace(a, plan, invert);
    }

    // Pack the real sequence into a complex array for the FFT
    std::vector<Complex> fft(const std::vector<double>& input, int fft_size) {
        std::vector<Complex> a(fft_size, 0.0);
        for (size_t i = 0; i < input.size(); ++i) a[i].real(input[i]);
        const Plan plan = make_plan(fft_size);
        fft_inplace(a.data(), plan, false);
        return a;
    }

    // In-place inverse transform, returning the real part
    std::vector<double> ifft(std::vector<Complex>& a) {
        const Plan plan = make_plan(static_cast<int>(a.size()));
        fft_inplace(a.data(), plan, true);
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

// Exact one-step martingale compensation for equal independent Poisson intensities.
// The log1p/expm1 form avoids cancellation when c_tau is small.
static double risk_neutral_drift(double c_tau) {
    if (!(std::isfinite(c_tau) && c_tau >= 0.0)) {
        throw std::runtime_error("c_tau must be finite and nonnegative");
    }
    const double one_minus_exp_neg_c = -std::expm1(-c_tau);
    if (c_tau < 0.5) {
        return -std::log1p(-(one_minus_exp_neg_c * one_minus_exp_neg_c));
    }
    return c_tau - std::log1p(one_minus_exp_neg_c);
}

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
    if (std::abs(alpha - 1.0) <= 1e-14) {
        return t == 1 ? 1.0 : 0.0;
    }
    double a = alpha;
    if (t == 1) {
        return 1.0 - 1.0 / std::tgamma(1.0 - a);
    }
    return (1.0 / std::tgamma(1.0 - a)) *
           (1.0 / std::pow(t - 1.0, a) - 1.0 / std::pow((double)t, a));
}

// ===================== Path simulation result struct (interface unchanged) =====================
struct PathSimResult {
    double final_price;
    double avg_price;
    double max_price;
    double min_price;
    bool barrier_hit_110;
    bool barrier_hit_90;
};

constexpr int INAR_DIRECT_CONV_THRESHOLD = 1024;

struct InarPrecomputed {
    std::vector<double> kernel;
    std::vector<double> lam_base;
    std::vector<std::vector<FFT::Complex>> kernel_fft_cache;
    std::vector<FFT::Plan> fft_plan_cache;
    struct Op {
        bool is_leaf = false;
        bool direct = false;
        int t = 0;
        int L = 0;
        int M = 0;
        int R = 0;
        int nA = 0;
        int nB = 0;
        int N = 0;
    };
    std::vector<Op> ops;
    double beta_param = 0.0;
    double coef = 0.0;
    double drift_coef = 0.0;
};

struct InarWorkspace {
    std::vector<double> lam_conv;
    std::vector<double> Y;
    std::vector<FFT::Complex> fft_a;
    int total_plus = 0;
    int total_minus = 0;
    double sum_prices = 0.0;
    double max_price = 0.0;
    double min_price = 0.0;
    double final_price = 0.0;
    bool barrier_hit_110 = false;
    bool barrier_hit_90 = false;

    explicit InarWorkspace(int T_steps)
        : lam_conv(T_steps + 2, 0.0),
          Y(T_steps + 2, 0.0),
          fft_a(next_pow2(2 * T_steps), FFT::Complex(0.0, 0.0)) {}

    void reset(double S0, double barrier_up, double barrier_down) {
        std::fill(lam_conv.begin(), lam_conv.end(), 0.0);
        total_plus = 0;
        total_minus = 0;
        sum_prices = S0;
        max_price = S0;
        min_price = S0;
        final_price = S0;
        barrier_hit_110 = S0 >= barrier_up;
        barrier_hit_90 = S0 <= barrier_down;
    }
};

void mark_cdq_kernel_lengths(int L, int R, std::vector<char>& needed) {
    if (L == R) {
        return;
    }
    needed[R - L] = 1;
    const int M = (L + R) >> 1;
    mark_cdq_kernel_lengths(L, M, needed);
    mark_cdq_kernel_lengths(M + 1, R, needed);
}

void build_cdq_ops(int L, int R, std::vector<InarPrecomputed::Op>& ops) {
    if (L == R) {
        InarPrecomputed::Op op;
        op.is_leaf = true;
        op.t = L;
        ops.push_back(op);
        return;
    }

    const int M = (L + R) >> 1;
    build_cdq_ops(L, M, ops);

    InarPrecomputed::Op op;
    op.is_leaf = false;
    op.L = L;
    op.M = M;
    op.R = R;
    op.nA = M - L + 1;
    op.nB = R - L;
    op.direct = op.nA * op.nB <= INAR_DIRECT_CONV_THRESHOLD;
    op.N = next_pow2(op.nA + op.nB - 1);
    ops.push_back(op);

    build_cdq_ops(M + 1, R, ops);
}

InarPrecomputed build_inar_precomputed(
    double alpha,
    double gamma_val,
    int T_steps,
    double beta_param,
    double mu,
    double theta,
    double xi
) {
    InarPrecomputed pre;
    pre.kernel.resize(T_steps + 1, 0.0);
    pre.lam_base.resize(T_steps + 1, 0.0);
    pre.kernel_fft_cache.resize(T_steps + 1);
    pre.fft_plan_cache.resize(next_pow2(2 * T_steps) + 1);
    pre.beta_param = beta_param;

    const double t_steps_d = static_cast<double>(T_steps);
    const double a_T = 1.0 - gamma_val * std::pow(t_steps_d, -alpha);
    const double mu_T = mu * std::pow(t_steps_d, alpha - 1.0);
    const double inv_one_minus_a = 1.0 / (1.0 - a_T);
    const double kernel_scale = a_T / (1.0 + beta_param);

    double phi_prefix = 0.0;
    for (int t_idx = 1; t_idx <= T_steps; ++t_idx) {
        const double f_val = f_alpha_one(t_idx, alpha);
        pre.kernel[t_idx] = kernel_scale * f_val;
        pre.lam_base[t_idx] = mu_T
            + xi * mu_T * (inv_one_minus_a * (1.0 - phi_prefix) - phi_prefix);
        phi_prefix += a_T * f_val;
    }

    const double denom = mu * std::pow(t_steps_d, alpha);
    pre.coef = std::sqrt((theta / 2.0) * (1.0 - a_T) / denom);
    pre.drift_coef = risk_neutral_drift(pre.coef);

    std::vector<char> needed(T_steps + 1, 0);
    mark_cdq_kernel_lengths(1, T_steps, needed);
    auto ensure_plan = [&](int N) -> const FFT::Plan& {
        if (pre.fft_plan_cache[N].n == 0) {
            pre.fft_plan_cache[N] = FFT::make_plan(N);
        }
        return pre.fft_plan_cache[N];
    };
    for (int nB = 1; nB <= T_steps; ++nB) {
        if (!needed[nB]) {
            continue;
        }
        const int nA = (nB + 2) / 2;
        const int N = next_pow2(nA + nB - 1);
        std::vector<FFT::Complex> cache(N, FFT::Complex(0.0, 0.0));
        for (int j = 0; j < nB; ++j) {
            cache[j] = pre.kernel[1 + j];
        }
        FFT::fft_inplace(cache.data(), ensure_plan(N), false);
        pre.kernel_fft_cache[nB] = std::move(cache);
    }

    pre.ops.reserve(2 * T_steps - 1);
    build_cdq_ops(1, T_steps, pre.ops);

    return pre;
}

PathSimResult simulate_inar_path_fft(
    const InarPrecomputed& pre,
    double S0,
    double barrier_up, double barrier_down,
    std::mt19937& gen,
    InarWorkspace& ws,
    std::poisson_distribution<int>& dist
) {
    const int T_steps = static_cast<int>(pre.lam_base.size()) - 1;
    ws.reset(S0, barrier_up, barrier_down);

    for (const InarPrecomputed::Op& op : pre.ops) {
        if (op.is_leaf) {
            const double lam = std::max(0.0, pre.lam_base[op.t] + ws.lam_conv[op.t]);
            const int x_plus = dist(gen, std::poisson_distribution<int>::param_type(lam));
            const int x_minus = dist(gen, std::poisson_distribution<int>::param_type(lam));
            ws.total_plus += x_plus;
            ws.total_minus += x_minus;
            ws.Y[op.t] = static_cast<double>(x_plus) + pre.beta_param * static_cast<double>(x_minus);

            const double diff_pm = static_cast<double>(ws.total_plus - ws.total_minus);
            const double P_t = pre.coef * diff_pm - pre.drift_coef * static_cast<double>(ws.total_plus);
            const double price_t = S0 * std::exp(P_t);
            ws.final_price = price_t;
            ws.sum_prices += price_t;
            ws.max_price = std::max(ws.max_price, price_t);
            ws.min_price = std::min(ws.min_price, price_t);
            ws.barrier_hit_110 = ws.barrier_hit_110 || (price_t >= barrier_up);
            ws.barrier_hit_90  = ws.barrier_hit_90 || (price_t <= barrier_down);
            continue;
        }

        if (op.direct) {
            for (int i = 0; i < op.nA; ++i) {
                const double yi = ws.Y[op.L + i];
                const int base = op.L + i;
                for (int t = op.M + 1; t <= op.R; ++t) {
                    ws.lam_conv[t] += yi * pre.kernel[t - base];
                }
            }
            continue;
        }

        for (int i = 0; i < op.nA; ++i) {
            ws.fft_a[i] = ws.Y[op.L + i];
        }
        std::fill(ws.fft_a.begin() + op.nA, ws.fft_a.begin() + op.N, FFT::Complex(0.0, 0.0));

        const FFT::Plan& plan = pre.fft_plan_cache[op.N];
        FFT::fft_inplace(ws.fft_a.data(), plan, false);
        const std::vector<FFT::Complex>& fft_b = pre.kernel_fft_cache[op.nB];
        for (int i = 0; i < op.N; ++i) {
            ws.fft_a[i] *= fft_b[i];
        }
        FFT::fft_inplace(ws.fft_a.data(), plan, true);

        for (int t = op.M + 1; t <= op.R; ++t) {
            const int k = t - op.L - 1;
            if (k >= 0 && k < op.N) {
                ws.lam_conv[t] += ws.fft_a[k].real();
            }
        }
    }

    PathSimResult result;
    result.final_price = ws.final_price;
    result.avg_price = ws.sum_prices / static_cast<double>(T_steps + 1);
    result.max_price = ws.max_price;
    result.min_price = ws.min_price;
    result.barrier_hit_110 = ws.barrier_hit_110;
    result.barrier_hit_90 = ws.barrier_hit_90;
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

// ===================== Statistics helper =====================
struct StatResult { double mean, ci_lower, ci_upper; };
struct PartialStats {
    std::vector<double> sum;
    std::vector<double> sum_sq;
    int count = 0;

    explicit PartialStats(size_t n_payoffs = 0)
        : sum(n_payoffs, 0.0), sum_sq(n_payoffs, 0.0) {}
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
    const InarPrecomputed pre = build_inar_precomputed(
        alpha, gamma_val, T_steps, beta_param, mu, theta, xi
    );

    const size_t n_strikes = strikes.size();
    const size_t n_payoffs = 8 * n_strikes;
    std::vector<std::future<void>> futures;
    std::vector<PartialStats> partials(num_threads, PartialStats(n_payoffs));

    int sims_per_thread = (n_sims + (int)num_threads - 1) / (int)num_threads;
    for (unsigned int t = 0; t < num_threads; ++t) {
        int start_sim = (int)t * sims_per_thread;
        int end_sim   = std::min(n_sims, (int)((t + 1) * sims_per_thread));
        if (start_sim >= end_sim) break;

        futures.push_back(std::async(std::launch::async, [&, start_sim, end_sim, t, base_seed]() {
            std::mt19937 gen(base_seed + (unsigned int)t * 101);
            std::poisson_distribution<int> dist;
            InarWorkspace ws(T_steps);
            PartialStats local(n_payoffs);
            for (int sim = start_sim; sim < end_sim; ++sim) {
                PathSimResult pr = simulate_inar_path_fft(
                    pre, S0, barrier_up, barrier_down, gen, ws, dist
                );
                for (size_t k = 0; k < n_strikes; ++k) {
                    const double K = strikes[k];
                    const double payoffs[8] = {
                        european_call_payoff(pr.final_price, K),
                        european_put_payoff(pr.final_price, K),
                        asian_call_payoff(pr.avg_price, K),
                        asian_put_payoff(pr.avg_price, K),
                        lookback_call_payoff(pr.max_price, K),
                        lookback_put_payoff(pr.min_price, K),
                        barrier_up_in_call_payoff(pr, K, barrier_up),
                        barrier_down_out_put_payoff(pr, K, barrier_down)
                    };
                    for (int opt = 0; opt < 8; ++opt) {
                        const size_t idx = static_cast<size_t>(opt) * n_strikes + k;
                        local.sum[idx] += payoffs[opt];
                        local.sum_sq[idx] += payoffs[opt] * payoffs[opt];
                    }
                }
                local.count += 1;
            }
            partials[t] = std::move(local);
        }));
    }
    for (auto& f : futures) f.get();

    std::vector<std::vector<StatResult>> all_stats(8, std::vector<StatResult>(n_strikes));
    const double n = static_cast<double>(n_sims);
    for (size_t k = 0; k < n_strikes; ++k) {
        for (int opt = 0; opt < 8; ++opt) {
            const size_t idx = static_cast<size_t>(opt) * n_strikes + k;
            double sum = 0.0;
            double sum_sq = 0.0;
            for (const PartialStats& partial : partials) {
                sum += partial.sum[idx];
                sum_sq += partial.sum_sq[idx];
            }

            const double mean = sum / n;
            double var = 0.0;
            if (n_sims > 1) {
                var = (sum_sq - (sum * sum) / n) / (n - 1.0);
                var = std::max(var, 0.0);
            }
            const double se = std::sqrt(var / n);
            all_stats[opt][k] = {mean, mean - 1.96 * se, mean + 1.96 * se};
        }
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
