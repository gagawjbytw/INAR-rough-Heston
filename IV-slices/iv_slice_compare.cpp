#include <algorithm>
#include <chrono>
#include <cmath>
#include <complex>
#include <future>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

struct RHParams {
    double H;
    double lambda;
    double theta;
    double nu;
    double rho;
    double V0;
    double S0;
    double T;
    double r;
};

struct SliceStats {
    std::string method;
    int steps = -1;
    int factors = -1;
    unsigned int threads = 0;
    std::vector<double> prices;
    std::vector<double> ivs;
    double max_iv_error = 0.0;
    double rms_iv_error = 0.0;
    double seconds = 0.0;
    double discounted_spot_mean = std::numeric_limits<double>::quiet_NaN();
    double spot_se = std::numeric_limits<double>::quiet_NaN();
    double martingale_z_score = std::numeric_limits<double>::quiet_NaN();
    double terminal_integrated_variance_mean = std::numeric_limits<double>::quiet_NaN();
    double terminal_noise_bracket_mean = std::numeric_limits<double>::quiet_NaN();
    double integrated_bracket_gap_max = std::numeric_limits<double>::quiet_NaN();
    double integrated_driver_mean = std::numeric_limits<double>::quiet_NaN();
    double integrated_driver_second_moment = std::numeric_limits<double>::quiet_NaN();
    double integrated_driver_bracket_ratio = std::numeric_limits<double>::quiet_NaN();
    double integrated_variance_driver_second_moment = std::numeric_limits<double>::quiet_NaN();
    double integrated_variance_driver_bracket_ratio = std::numeric_limits<double>::quiet_NaN();
    double projection_hit_rate = std::numeric_limits<double>::quiet_NaN();
    double projection_gap_mean = std::numeric_limits<double>::quiet_NaN();
    double projection_gap_max = std::numeric_limits<double>::quiet_NaN();
};

struct TerminalSpotMoments {
    long double sum = 0.0L;
    long double sum_sq = 0.0L;

    void add(double discounted_spot) {
        const long double x = static_cast<long double>(discounted_spot);
        sum += x;
        sum_sq += x * x;
    }

    void merge(const TerminalSpotMoments& other) {
        sum += other.sum;
        sum_sq += other.sum_sq;
    }
};

static unsigned int checked_worker_threads(int n_sims, unsigned int requested_threads) {
    if (n_sims <= 0) {
        throw std::runtime_error("n_sims must be positive.");
    }
    if (requested_threads == 0) {
        throw std::runtime_error("threads must be positive.");
    }
    return std::min(requested_threads, static_cast<unsigned int>(n_sims));
}

static void set_spot_diagnostics(
    SliceStats& out,
    int n_sims,
    unsigned int threads,
    double S0,
    const TerminalSpotMoments& moments
) {
    const long double n = static_cast<long double>(n_sims);
    const long double mean = moments.sum / n;
    long double sample_variance = 0.0L;
    if (n_sims > 1) {
        sample_variance = (moments.sum_sq - moments.sum * moments.sum / n)
            / static_cast<long double>(n_sims - 1);
        sample_variance = std::max(0.0L, sample_variance);
    }
    const long double se = std::sqrt(sample_variance / n);

    out.threads = threads;
    out.discounted_spot_mean = static_cast<double>(mean);
    out.spot_se = static_cast<double>(se);
    if (se > 0.0L) {
        out.martingale_z_score = static_cast<double>(
            (mean - static_cast<long double>(S0)) / se
        );
    }

    std::ostringstream message;
    message << std::setprecision(10)
            << "[spot-diagnostic] method=" << out.method
            << " paths=" << n_sims
            << " threads=" << threads
            << " discounted_spot_mean=" << out.discounted_spot_mean
            << " spot_se=" << out.spot_se
            << " martingale_z_score=" << out.martingale_z_score;
    std::cerr << message.str() << '\n';
}

double european_call_payoff(double S_T, double K) {
    return std::max(S_T - K, 0.0);
}

double norm_cdf(double x) {
    return 0.5 * std::erfc(-x / std::sqrt(2.0));
}

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

double bs_call_price(double S0, double K, double T, double r, double sigma) {
    if (T <= 0.0) {
        return std::max(S0 - K, 0.0);
    }
    if (sigma <= 0.0) {
        return std::max(S0 - K * std::exp(-r * T), 0.0);
    }
    const double sqrtT = std::sqrt(T);
    const double d1 = (std::log(S0 / K) + (r + 0.5 * sigma * sigma) * T) / (sigma * sqrtT);
    const double d2 = d1 - sigma * sqrtT;
    return S0 * norm_cdf(d1) - K * std::exp(-r * T) * norm_cdf(d2);
}

double implied_vol_call(double price, double S0, double K, double T, double r) {
    const double intrinsic = std::max(S0 - K * std::exp(-r * T), 0.0);
    const double upper_bound = S0;
    if (price <= intrinsic + 1e-12) {
        return 0.0;
    }
    if (price >= upper_bound - 1e-12) {
        return 5.0;
    }

    double lo = 1e-8;
    double hi = 5.0;
    for (int iter = 0; iter < 120; ++iter) {
        const double mid = 0.5 * (lo + hi);
        const double mid_price = bs_call_price(S0, K, T, r, mid);
        if (mid_price < price) {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    return 0.5 * (lo + hi);
}

SliceStats finalize_slice_stats(
    const std::string& method,
    int steps,
    int factors,
    const std::vector<double>& mean_prices,
    const std::vector<double>& strikes,
    double S0,
    double T,
    double r,
    const std::vector<double>& reference_ivs,
    double elapsed
) {
    SliceStats out;
    out.method = method;
    out.steps = steps;
    out.factors = factors;
    out.prices = mean_prices;
    out.ivs.resize(strikes.size(), 0.0);
    out.seconds = elapsed;

    double max_err = 0.0;
    double sum_sq = 0.0;
    for (size_t i = 0; i < strikes.size(); ++i) {
        out.ivs[i] = implied_vol_call(mean_prices[i], S0, strikes[i], T, r);
        const double err = std::abs(out.ivs[i] - reference_ivs[i]);
        max_err = std::max(max_err, err);
        sum_sq += err * err;
    }
    out.max_iv_error = max_err;
    out.rms_iv_error = std::sqrt(sum_sq / static_cast<double>(strikes.size()));
    return out;
}

double sample_inverse_gaussian(double mu, double lambda, std::mt19937& gen) {
    if (!(mu > 0.0) || !(lambda > 0.0)) {
        return 0.0;
    }

    static thread_local std::normal_distribution<double> nd(0.0, 1.0);
    static thread_local std::uniform_real_distribution<double> ud(0.0, 1.0);

    const double xi = nd(gen);
    const double y = xi * xi;
    const double mu2 = mu * mu;
    const double sqrt_term = std::sqrt(4.0 * mu * lambda * y + mu2 * y * y);
    const double x = mu + (mu2 * y) / (2.0 * lambda) - (mu / (2.0 * lambda)) * sqrt_term;
    const double u = ud(gen);
    return (u <= mu / (mu + x)) ? x : mu2 / x;
}

// ============================================================
// Integrated multifactor Euler
// ============================================================
struct KernelApprox {
    std::vector<double> nodes;
    std::vector<double> weights;
};

KernelApprox build_loggrid_kernel_approx(double H, double lambda, int n,
                                         double x_min, double x_max) {
    if (!(H > 0.0 && H < 0.5)) throw std::runtime_error("H must be in (0, 0.5).");
    if (n < 1) throw std::runtime_error("n must be >= 1.");
    if (!(x_min > 0.0 && x_max > x_min)) throw std::runtime_error("Need 0 < x_min < x_max.");

    KernelApprox ka;
    ka.nodes.resize(n);
    ka.weights.resize(n);

    const double cH = 1.0 / (std::tgamma(H + 0.5) * std::tgamma(0.5 - H));
    const double expo = 0.5 - H;

    const double log_min = std::log(x_min);
    const double log_max = std::log(x_max);
    const double dlog = (log_max - log_min) / static_cast<double>(n);

    for (int i = 0; i < n; ++i) {
        const double left = log_min + i * dlog;
        const double right = log_min + (i + 1) * dlog;
        const double mid = 0.5 * (left + right);

        const double xL = std::exp(left);
        const double xR = std::exp(right);
        const double xM = std::exp(mid);

        const double w = lambda * cH / expo * (std::pow(xR, expo) - std::pow(xL, expo));

        ka.nodes[i] = xM;
        ka.weights[i] = w;
    }

    return ka;
}

struct IntegratedTerminalResult {
    double spot = 0.0;
    double terminal_integrated_variance = 0.0;
    double terminal_noise_bracket = 0.0;
    double martingale_driver = 0.0;
    double variance_driver = 0.0;
    long double projection_gap_sum = 0.0L;
    double projection_gap_max = 0.0;
    unsigned long long projection_hits = 0;
};

IntegratedTerminalResult simulate_integrated_multifactor_euler_terminal(
    const RHParams& p,
    const KernelApprox& ka,
    int M_steps,
    std::mt19937& gen
) {
    std::normal_distribution<double> nd(0.0, 1.0);

    const int n = static_cast<int>(ka.nodes.size());
    const double dt = p.T / static_cast<double>(M_steps);
    const double rho_perp = std::sqrt(std::max(0.0, 1.0 - p.rho * p.rho));

    std::vector<double> Xi(n, 0.0), Xi_next(n, 0.0), decay(n);
    for (int i = 0; i < n; ++i) decay[i] = std::exp(-ka.nodes[i] * dt);

    double Xbar = 0.0;
    double M = 0.0;
    double Mperp = 0.0;
    double noise_bracket = 0.0;
    long double projection_gap_sum = 0.0L;
    double projection_gap_max = 0.0;
    unsigned long long projection_hits = 0;

    for (int k = 0; k < M_steps; ++k) {
        const double t_k = k * dt;

        // ka.weights already approximate lambda * K_H, so lambda must not be
        // applied a second time.  The integrated scheme feeds back the
        // nondecreasing projection Xbar used to define the Gaussian bracket.
        const double common = (p.theta * t_k - Xbar + p.nu * M) * dt;
        for (int i = 0; i < n; ++i) {
            Xi_next[i] = decay[i] * (Xi[i] + common);
        }

        const double t_next = (k + 1) * dt;
        double X_next = p.V0 * t_next;
        for (int i = 0; i < n; ++i) X_next += ka.weights[i] * Xi_next[i];

        const double Xbar_next = std::max(Xbar, X_next);
        const double dXbar = std::max(0.0, Xbar_next - Xbar);
        const double projection_gap = std::max(0.0, Xbar_next - X_next);
        projection_gap_sum += static_cast<long double>(projection_gap);
        projection_gap_max = std::max(projection_gap_max, projection_gap);
        projection_hits += static_cast<unsigned long long>(X_next < Xbar);

        const double Z1 = nd(gen);
        const double Z2 = nd(gen);

        M += std::sqrt(dXbar) * Z1;
        Mperp += std::sqrt(dXbar) * Z2;
        noise_bracket += dXbar;

        Xi.swap(Xi_next);
        Xbar = Xbar_next;
    }

    const double martingale_driver = p.rho * M + rho_perp * Mperp;
    const double logS_T =
        std::log(p.S0) + p.r * p.T
        - 0.5 * Xbar
        + martingale_driver;

    IntegratedTerminalResult out;
    out.spot = std::exp(logS_T);
    out.terminal_integrated_variance = Xbar;
    out.terminal_noise_bracket = noise_bracket;
    out.martingale_driver = martingale_driver;
    out.variance_driver = M;
    out.projection_gap_sum = projection_gap_sum;
    out.projection_gap_max = projection_gap_max;
    out.projection_hits = projection_hits;
    return out;
}

SliceStats run_integrated_slice_mc(
    int n_sims,
    int M_steps,
    const RHParams& p,
    const KernelApprox& ka,
    const std::vector<double>& strikes,
    const std::vector<double>& reference_ivs,
    int factors,
    unsigned int requested_threads,
    unsigned int base_seed = 123456u
) {
    auto t0 = std::chrono::high_resolution_clock::now();

    struct PartialStats {
        std::vector<double> sum;
        TerminalSpotMoments spot;
        long double terminal_integrated_variance_sum = 0.0L;
        long double terminal_noise_bracket_sum = 0.0L;
        double integrated_bracket_gap_max = 0.0;
        long double integrated_driver_sum = 0.0L;
        long double integrated_driver_sq_sum = 0.0L;
        long double integrated_variance_driver_sq_sum = 0.0L;
        long double projection_gap_sum = 0.0L;
        double projection_gap_max = 0.0;
        unsigned long long projection_hits = 0;
        explicit PartialStats(size_t n = 0) : sum(n, 0.0) {}
    };

    const unsigned int num_threads = checked_worker_threads(n_sims, requested_threads);
    const int sims_per_thread = (n_sims + static_cast<int>(num_threads) - 1) / static_cast<int>(num_threads);
    std::vector<std::future<void>> futures;
    std::vector<PartialStats> partials(num_threads, PartialStats(strikes.size()));

    for (unsigned int tid = 0; tid < num_threads; ++tid) {
        const int L = static_cast<int>(tid) * sims_per_thread;
        const int R = std::min(n_sims, static_cast<int>((tid + 1) * sims_per_thread));
        if (L >= R) break;

        futures.push_back(std::async(std::launch::async, [&, L, R, tid]() {
            std::mt19937 gen(base_seed + 10007u * tid);
            PartialStats local(strikes.size());
            for (int i = L; i < R; ++i) {
                const IntegratedTerminalResult terminal =
                    simulate_integrated_multifactor_euler_terminal(p, ka, M_steps, gen);
                const double discounted_ST = std::exp(-p.r * p.T) * terminal.spot;
                local.spot.add(discounted_ST);
                local.terminal_integrated_variance_sum += terminal.terminal_integrated_variance;
                local.terminal_noise_bracket_sum += terminal.terminal_noise_bracket;
                local.integrated_bracket_gap_max = std::max(
                    local.integrated_bracket_gap_max,
                    std::abs(terminal.terminal_integrated_variance - terminal.terminal_noise_bracket)
                );
                local.integrated_driver_sum += terminal.martingale_driver;
                local.integrated_driver_sq_sum +=
                    static_cast<long double>(terminal.martingale_driver)
                    * static_cast<long double>(terminal.martingale_driver);
                local.integrated_variance_driver_sq_sum +=
                    static_cast<long double>(terminal.variance_driver)
                    * static_cast<long double>(terminal.variance_driver);
                local.projection_gap_sum += terminal.projection_gap_sum;
                local.projection_gap_max = std::max(
                    local.projection_gap_max, terminal.projection_gap_max
                );
                local.projection_hits += terminal.projection_hits;
                for (size_t k = 0; k < strikes.size(); ++k) {
                    local.sum[k] += std::exp(-p.r * p.T)
                        * european_call_payoff(terminal.spot, strikes[k]);
                }
            }
            partials[tid] = std::move(local);
        }));
    }
    for (auto& f : futures) f.get();

    auto t1 = std::chrono::high_resolution_clock::now();
    const double elapsed = std::chrono::duration<double>(t1 - t0).count();

    std::vector<double> prices(strikes.size(), 0.0);
    TerminalSpotMoments spot_moments;
    long double terminal_integrated_variance_sum = 0.0L;
    long double terminal_noise_bracket_sum = 0.0L;
    double integrated_bracket_gap_max = 0.0;
    long double integrated_driver_sum = 0.0L;
    long double integrated_driver_sq_sum = 0.0L;
    long double integrated_variance_driver_sq_sum = 0.0L;
    long double projection_gap_sum = 0.0L;
    double projection_gap_max = 0.0;
    unsigned long long projection_hits = 0;
    for (const PartialStats& partial : partials) {
        for (size_t k = 0; k < strikes.size(); ++k) {
            prices[k] += partial.sum[k];
        }
        spot_moments.merge(partial.spot);
        terminal_integrated_variance_sum += partial.terminal_integrated_variance_sum;
        terminal_noise_bracket_sum += partial.terminal_noise_bracket_sum;
        integrated_bracket_gap_max = std::max(
            integrated_bracket_gap_max, partial.integrated_bracket_gap_max
        );
        integrated_driver_sum += partial.integrated_driver_sum;
        integrated_driver_sq_sum += partial.integrated_driver_sq_sum;
        integrated_variance_driver_sq_sum += partial.integrated_variance_driver_sq_sum;
        projection_gap_sum += partial.projection_gap_sum;
        projection_gap_max = std::max(projection_gap_max, partial.projection_gap_max);
        projection_hits += partial.projection_hits;
    }
    for (double& x : prices) x /= static_cast<double>(n_sims);

    SliceStats out = finalize_slice_stats(
        "Integrated-MultifactorEuler", M_steps, factors, prices,
        strikes, p.S0, p.T, p.r, reference_ivs, elapsed
    );
    set_spot_diagnostics(out, n_sims, num_threads, p.S0, spot_moments);

    const long double n = static_cast<long double>(n_sims);
    out.terminal_integrated_variance_mean = static_cast<double>(
        terminal_integrated_variance_sum / n
    );
    out.terminal_noise_bracket_mean = static_cast<double>(terminal_noise_bracket_sum / n);
    out.integrated_bracket_gap_max = integrated_bracket_gap_max;
    out.integrated_driver_mean = static_cast<double>(integrated_driver_sum / n);
    out.integrated_driver_second_moment = static_cast<double>(integrated_driver_sq_sum / n);
    out.integrated_variance_driver_second_moment = static_cast<double>(
        integrated_variance_driver_sq_sum / n
    );
    if (out.terminal_noise_bracket_mean > 0.0) {
        out.integrated_driver_bracket_ratio =
            out.integrated_driver_second_moment / out.terminal_noise_bracket_mean;
        out.integrated_variance_driver_bracket_ratio =
            out.integrated_variance_driver_second_moment
            / out.terminal_noise_bracket_mean;
    }
    const long double projection_observations = n * static_cast<long double>(M_steps);
    out.projection_hit_rate = static_cast<double>(
        static_cast<long double>(projection_hits) / projection_observations
    );
    out.projection_gap_mean = static_cast<double>(projection_gap_sum / projection_observations);
    out.projection_gap_max = projection_gap_max;

    std::ostringstream message;
    message << std::setprecision(10)
            << "[integrated-diagnostic] terminal_integrated_variance_mean="
            << out.terminal_integrated_variance_mean
            << " terminal_noise_bracket_mean=" << out.terminal_noise_bracket_mean
            << " integrated_bracket_gap_max=" << out.integrated_bracket_gap_max
            << " driver_mean=" << out.integrated_driver_mean
            << " driver_second_moment=" << out.integrated_driver_second_moment
            << " driver_bracket_ratio=" << out.integrated_driver_bracket_ratio
            << " variance_driver_second_moment="
            << out.integrated_variance_driver_second_moment
            << " variance_driver_bracket_ratio="
            << out.integrated_variance_driver_bracket_ratio
            << " projection_hit_rate=" << out.projection_hit_rate
            << " projection_gap_mean=" << out.projection_gap_mean
            << " projection_gap_max=" << out.projection_gap_max;
    std::cerr << message.str() << '\n';
    return out;
}

// ============================================================
// INAR-CDQ-FFT
// ============================================================
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
        while ((1 << levels) < n) ++levels;
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
}

static inline int next_pow2(int x) {
    int p = 1;
    while (p < x) p <<= 1;
    return p;
}

double solve_beta(double rho) {
    const double a = 2 * rho * rho - 1;
    const double b = 2;
    const double c = 2 * rho * rho - 1;
    const double disc = b * b - 4 * a * c;
    if (disc < 0) throw std::runtime_error("Invalid rho: discriminant < 0");

    const double r1 = (-b + std::sqrt(disc)) / (2 * a);
    const double r2 = (-b - std::sqrt(disc)) / (2 * a);

    if (r1 > 1) return r1;
    if (r2 > 1) return r2;
    throw std::runtime_error("No valid beta found.");
}

double solve_mu(double nu, double theta, double beta, double gamma) {
    return theta * (1 + beta * beta) /
           (gamma * nu * nu * (1 + beta) * (1 + beta));
}

double solve_xi(double V0, double theta) {
    return V0 / theta;
}

double f_alpha_one(int t, double alpha) {
    if (t <= 0) return 0.0;
    if (std::abs(alpha - 1.0) <= 1e-14) {
        return t == 1 ? 1.0 : 0.0;
    }
    const double a = alpha;
    if (t == 1) {
        return 1.0 - 1.0 / std::tgamma(1.0 - a);
    }
    return (1.0 / std::tgamma(1.0 - a)) *
           (1.0 / std::pow(t - 1.0, a) - 1.0 / std::pow(static_cast<double>(t), a));
}

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

    explicit InarWorkspace(int T_steps)
        : lam_conv(T_steps + 2, 0.0),
          Y(T_steps + 2, 0.0),
          fft_a(next_pow2(2 * T_steps), FFT::Complex(0.0, 0.0)) {}

    void reset() {
        std::fill(lam_conv.begin(), lam_conv.end(), 0.0);
        total_plus = 0;
        total_minus = 0;
    }
};

void mark_cdq_kernel_lengths(int L, int R, std::vector<char>& needed) {
    if (L == R) return;
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
    double maturity,
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

    if (!(maturity > 0.0)) {
        throw std::runtime_error("INAR maturity must be positive.");
    }
    const double tau_scale = static_cast<double>(T_steps) / maturity;
    const double a_T = 1.0 - gamma_val * std::pow(tau_scale, -alpha);
    const double mu_T = mu * std::pow(tau_scale, alpha - 1.0);
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

    const double denom = mu * std::pow(tau_scale, alpha);
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
        if (!needed[nB]) continue;
        const int nA = (nB + 2) / 2;
        const int N = next_pow2(nA + nB - 1);
        std::vector<FFT::Complex> cache(N, FFT::Complex(0.0, 0.0));
        for (int j = 0; j < nB; ++j) cache[j] = pre.kernel[1 + j];
        FFT::fft_inplace(cache.data(), ensure_plan(N), false);
        pre.kernel_fft_cache[nB] = std::move(cache);
    }

    pre.ops.reserve(2 * T_steps - 1);
    build_cdq_ops(1, T_steps, pre.ops);
    return pre;
}

double simulate_inar_terminal(
    const InarPrecomputed& pre,
    double S0,
    std::mt19937& gen,
    InarWorkspace& ws,
    std::poisson_distribution<int>& dist
) {
    ws.reset();

    for (const InarPrecomputed::Op& op : pre.ops) {
        if (op.is_leaf) {
            const double lam = std::max(0.0, pre.lam_base[op.t] + ws.lam_conv[op.t]);
            const auto poisson_param = std::poisson_distribution<int>::param_type(lam);
            const int x_plus = dist(gen, poisson_param);
            const int x_minus = dist(gen, poisson_param);
            ws.total_plus += x_plus;
            ws.total_minus += x_minus;
            ws.Y[op.t] = static_cast<double>(x_plus) + pre.beta_param * static_cast<double>(x_minus);
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

        for (int i = 0; i < op.nA; ++i) ws.fft_a[i] = ws.Y[op.L + i];
        std::fill(ws.fft_a.begin() + op.nA, ws.fft_a.begin() + op.N, FFT::Complex(0.0, 0.0));

        const FFT::Plan& plan = pre.fft_plan_cache[op.N];
        FFT::fft_inplace(ws.fft_a.data(), plan, false);
        const std::vector<FFT::Complex>& fft_b = pre.kernel_fft_cache[op.nB];
        for (int i = 0; i < op.N; ++i) ws.fft_a[i] *= fft_b[i];
        FFT::fft_inplace(ws.fft_a.data(), plan, true);

        for (int t = op.M + 1; t <= op.R; ++t) {
            const int k = t - op.L - 1;
            if (k >= 0 && k < op.N) ws.lam_conv[t] += ws.fft_a[k].real();
        }
    }

    const double diff_pm = static_cast<double>(ws.total_plus - ws.total_minus);
    const double log_return = pre.coef * diff_pm - pre.drift_coef * static_cast<double>(ws.total_plus);
    return S0 * std::exp(log_return);
}

SliceStats run_inar_slice_mc(
    int n_sims,
    int T_steps,
    const std::vector<double>& strikes,
    double alpha,
    double gamma_val,
    double beta_param,
    double mu,
    double theta,
    double xi,
    double S0,
    double T,
    double r,
    const std::vector<double>& reference_ivs,
    unsigned int requested_threads,
    unsigned int base_seed = 123456789u
) {
    const InarPrecomputed pre = build_inar_precomputed(alpha, gamma_val, T_steps, T, beta_param, mu, theta, xi);
    auto t0 = std::chrono::high_resolution_clock::now();

    struct PartialStats {
        std::vector<double> sum;
        TerminalSpotMoments spot;
        explicit PartialStats(size_t n = 0) : sum(n, 0.0) {}
    };

    const unsigned int num_threads = checked_worker_threads(n_sims, requested_threads);
    const int sims_per_thread = (n_sims + static_cast<int>(num_threads) - 1) / static_cast<int>(num_threads);
    std::vector<std::future<void>> futures;
    std::vector<PartialStats> partials(num_threads, PartialStats(strikes.size()));

    for (unsigned int tid = 0; tid < num_threads; ++tid) {
        const int L = static_cast<int>(tid) * sims_per_thread;
        const int R = std::min(n_sims, static_cast<int>((tid + 1) * sims_per_thread));
        if (L >= R) break;

        futures.push_back(std::async(std::launch::async, [&, L, R, tid]() {
            std::mt19937 gen(base_seed + tid * 101u);
            std::poisson_distribution<int> dist;
            InarWorkspace ws(T_steps);
            PartialStats local(strikes.size());
            for (int i = L; i < R; ++i) {
                // simulate_inar_terminal returns the discounted microscopic spot.
                const double ST = std::exp(r * T)
                    * simulate_inar_terminal(pre, S0, gen, ws, dist);
                const double discounted_ST = std::exp(-r * T) * ST;
                local.spot.add(discounted_ST);
                for (size_t k = 0; k < strikes.size(); ++k) {
                    local.sum[k] += std::exp(-r * T) * european_call_payoff(ST, strikes[k]);
                }
            }
            partials[tid] = std::move(local);
        }));
    }
    for (auto& f : futures) f.get();

    auto t1 = std::chrono::high_resolution_clock::now();
    const double elapsed = std::chrono::duration<double>(t1 - t0).count();

    std::vector<double> prices(strikes.size(), 0.0);
    TerminalSpotMoments spot_moments;
    for (const PartialStats& partial : partials) {
        for (size_t k = 0; k < strikes.size(); ++k) prices[k] += partial.sum[k];
        spot_moments.merge(partial.spot);
    }
    for (double& x : prices) x /= static_cast<double>(n_sims);

    SliceStats out = finalize_slice_stats(
        "INAR-CDQ-FFT", T_steps, -1, prices,
        strikes, S0, T, r, reference_ivs, elapsed
    );
    set_spot_diagnostics(out, n_sims, num_threads, S0, spot_moments);
    return out;
}

// ============================================================
// iVi-InverseGaussian
// ============================================================
struct IVIParams {
    double alpha;
    double gamma_val;
    double theta;
    double nu;
    double rho;
    double V0;
    double S0;
    double T;
    double r;
};

struct IVIPrecomputed {
    std::vector<double> k;
    std::vector<double> g0_step_integrals;
    double b;
    double c;
    double rho;
    double rho_perp;
    double logS0;
    double discount;
};

IVIPrecomputed build_ivi_precomputed(const IVIParams& p, int n_steps) {
    const double H = p.alpha - 0.5;
    if (!(H > 0.0 && H < 0.5)) {
        throw std::runtime_error("iVi benchmark is set up for H in (0, 0.5).");
    }

    const double dt = p.T / static_cast<double>(n_steps);
    const double beta_kernel = 1.0 / std::tgamma(H + 0.5);
    const double a = p.gamma_val * p.theta;
    const double b = -p.gamma_val;
    const double c = p.gamma_val * p.nu;

    IVIPrecomputed pre;
    pre.k.resize(n_steps, 0.0);
    pre.g0_step_integrals.resize(n_steps, 0.0);
    pre.b = b;
    pre.c = c;
    pre.rho = p.rho;
    pre.rho_perp = std::sqrt(std::max(0.0, 1.0 - p.rho * p.rho));
    pre.logS0 = std::log(p.S0) + p.r * p.T;
    pre.discount = std::exp(-p.r * p.T);

    for (int ell = 0; ell < n_steps; ++ell) {
        const double left = static_cast<double>(ell) * dt;
        const double right = static_cast<double>(ell + 1) * dt;
        pre.k[ell] = beta_kernel
            * (std::pow(right, H + 0.5) - std::pow(left, H + 0.5))
            / (H + 0.5);
    }

    for (int i = 0; i < n_steps; ++i) {
        const double ti = static_cast<double>(i) * dt;
        const double ti1 = static_cast<double>(i + 1) * dt;
        pre.g0_step_integrals[i] =
            p.V0 * dt
            + a * beta_kernel
                * (std::pow(ti1, H + 1.5) - std::pow(ti, H + 1.5))
                / ((H + 1.5) * (H + 0.5));
    }

    return pre;
}

double simulate_ivi_terminal_spot(
    const IVIPrecomputed& pre,
    int n_steps,
    std::mt19937& gen,
    std::vector<double>& U_inc,
    std::vector<double>& Z_inc
) {
    static thread_local std::normal_distribution<double> nd(0.0, 1.0);

    double logS = pre.logS0;
    const double one_minus_bk0 = 1.0 - pre.b * pre.k[0];
    const double ck0 = pre.c * pre.k[0];

    for (int i = 0; i < n_steps; ++i) {
        double alpha_i = pre.g0_step_integrals[i];
        for (int j = 0; j < i; ++j) {
            alpha_i += pre.k[i - j] * (pre.b * U_inc[j] + pre.c * Z_inc[j]);
        }
        alpha_i = std::max(alpha_i, 0.0);

        double u_i = 0.0;
        double z_i = 0.0;
        if (alpha_i > 0.0 && ck0 > 0.0) {
            const double mu = alpha_i / one_minus_bk0;
            const double lambda = (alpha_i / ck0) * (alpha_i / ck0);
            u_i = sample_inverse_gaussian(mu, lambda, gen);
            z_i = (one_minus_bk0 * u_i - alpha_i) / ck0;
        }

        U_inc[i] = u_i;
        Z_inc[i] = z_i;

        const double N_i = nd(gen);
        logS += -0.5 * u_i + pre.rho * z_i + pre.rho_perp * std::sqrt(std::max(0.0, u_i)) * N_i;
    }

    return std::exp(logS);
}

SliceStats run_ivi_slice_mc(
    int n_sims,
    int n_steps,
    const IVIParams& p,
    const std::vector<double>& strikes,
    const std::vector<double>& reference_ivs,
    unsigned int requested_threads,
    unsigned int base_seed = 123456u
) {
    const IVIPrecomputed pre = build_ivi_precomputed(p, n_steps);
    auto t0 = std::chrono::high_resolution_clock::now();

    struct PartialStats {
        std::vector<double> sum;
        TerminalSpotMoments spot;
        explicit PartialStats(size_t n = 0) : sum(n, 0.0) {}
    };

    const unsigned int num_threads = checked_worker_threads(n_sims, requested_threads);
    const int sims_per_thread =
        (n_sims + static_cast<int>(num_threads) - 1) / static_cast<int>(num_threads);

    std::vector<PartialStats> partials(num_threads, PartialStats(strikes.size()));
    std::vector<std::future<void>> futures;
    for (unsigned int tid = 0; tid < num_threads; ++tid) {
        const int L = static_cast<int>(tid) * sims_per_thread;
        const int R = std::min(n_sims, static_cast<int>((tid + 1) * sims_per_thread));
        if (L >= R) break;

        futures.push_back(std::async(std::launch::async, [&, L, R, tid]() {
            std::mt19937 gen(base_seed + 7919u * tid);
            std::vector<double> U_inc(n_steps, 0.0), Z_inc(n_steps, 0.0);
            PartialStats local(strikes.size());
            for (int i = L; i < R; ++i) {
                const double ST = simulate_ivi_terminal_spot(pre, n_steps, gen, U_inc, Z_inc);
                local.spot.add(pre.discount * ST);
                for (size_t k = 0; k < strikes.size(); ++k) {
                    local.sum[k] += pre.discount * european_call_payoff(ST, strikes[k]);
                }
            }
            partials[tid] = std::move(local);
        }));
    }
    for (auto& f : futures) f.get();

    auto t1 = std::chrono::high_resolution_clock::now();
    const double elapsed = std::chrono::duration<double>(t1 - t0).count();

    std::vector<double> prices(strikes.size(), 0.0);
    TerminalSpotMoments spot_moments;
    for (const PartialStats& partial : partials) {
        for (size_t k = 0; k < strikes.size(); ++k) prices[k] += partial.sum[k];
        spot_moments.merge(partial.spot);
    }
    for (double& x : prices) x /= static_cast<double>(n_sims);

    SliceStats out = finalize_slice_stats(
        "iVi-InverseGaussian", n_steps, -1, prices,
        strikes, p.S0, p.T, p.r, reference_ivs, elapsed
    );
    set_spot_diagnostics(out, n_sims, num_threads, p.S0, spot_moments);
    return out;
}

// ============================================================
// HQE-RoughHeston
// ============================================================
struct HQEParams {
    double alpha;
    double gamma_val;
    double theta_long;
    double nu;
    double rho;
    double V0;
    double S0;
    double T;
    double r;
};

struct HQEPrecomputed {
    int n_steps;
    double dt;
    double H;
    double alpha;
    double lambda;
    double V0;
    double theta_long;
    double nu_coeff;
    double rho;
    double rho_perp;
    double beta;
    double gamma;
    double logS0;
    double discount;
    std::vector<double> xi_grid;
    std::vector<double> b_star;
};

long double mittag_leffler(long double z, long double alpha, long double beta) {
    long double sum = 0.0L;
    long double z_pow = 1.0L;
    for (int n = 0; n < 400; ++n) {
        const long double term = z_pow / std::tgamma(static_cast<double>(alpha * n + beta));
        sum += term;
        if (std::fabsl(term) < 1e-18L * std::max(1.0L, std::fabsl(sum))) {
            break;
        }
        z_pow *= z;
    }
    return sum;
}

double kernel_value_hqe(double t, const HQEPrecomputed& pre) {
    if (t <= 0.0) return 0.0;
    const long double z = -static_cast<long double>(pre.lambda) * std::powl(t, pre.alpha);
    const long double ml = mittag_leffler(z, pre.alpha, pre.alpha);
    return pre.nu_coeff * std::pow(t, pre.alpha - 1.0) * static_cast<double>(ml);
}

double integrate_kernel_square_hqe(double a, double b, const HQEPrecomputed& pre) {
    const int panels = (a == 0.0 ? 512 : 256);
    const double h = (b - a) / static_cast<double>(panels);
    double sum = 0.0;
    for (int i = 0; i < panels; ++i) {
        const double x = a + (i + 0.5) * h;
        const double k = kernel_value_hqe(x, pre);
        sum += k * k;
    }
    return sum * h;
}

HQEPrecomputed build_hqe_precomputed(const HQEParams& p, int n_steps) {
    HQEPrecomputed pre{};
    pre.n_steps = n_steps;
    pre.dt = p.T / static_cast<double>(n_steps);
    pre.H = p.alpha - 0.5;
    pre.alpha = p.alpha;
    pre.lambda = p.gamma_val;
    pre.V0 = p.V0;
    pre.theta_long = p.theta_long;
    pre.nu_coeff = p.gamma_val * p.nu;
    pre.rho = p.rho;
    pre.rho_perp = std::sqrt(std::max(0.0, 1.0 - p.rho * p.rho));
    pre.logS0 = std::log(p.S0) + p.r * p.T;
    pre.discount = std::exp(-p.r * p.T);
    pre.xi_grid.assign(n_steps + 1, 0.0);
    pre.b_star.assign(n_steps + 1, 0.0);

    if (!(pre.H > 0.0 && pre.H < 0.5)) {
        throw std::runtime_error("HQE benchmark expects H in (0, 0.5).");
    }
    if (!(pre.lambda > 0.0 && pre.nu_coeff > 0.0)) {
        throw std::runtime_error("HQE needs positive lambda and nu_coeff.");
    }

    pre.xi_grid[0] = p.V0;
    for (int i = 1; i <= n_steps; ++i) {
        const double t = i * pre.dt;
        const long double z = -static_cast<long double>(pre.lambda) * std::powl(t, pre.alpha);
        const long double ml = mittag_leffler(z, pre.alpha, 1.0L);
        pre.xi_grid[i] = p.V0 + (p.theta_long - p.V0) * (1.0 - static_cast<double>(ml));
    }

    const double K0 = pre.nu_coeff / pre.lambda
        * (1.0 - static_cast<double>(mittag_leffler(
            -static_cast<long double>(pre.lambda) * std::powl(pre.dt, pre.alpha),
            pre.alpha,
            1.0L
        )));
    pre.beta = K0 / pre.dt;

    for (int j = 1; j <= n_steps; ++j) {
        const double a = (j - 1) * pre.dt;
        const double b = j * pre.dt;
        const double k_diag = integrate_kernel_square_hqe(a, b, pre);
        pre.b_star[j] = std::sqrt(std::max(0.0, k_diag / pre.dt));
        if (j == 1) {
            pre.gamma = k_diag - K0 * K0 / pre.dt;
        }
    }

    if (!(pre.gamma > 0.0)) {
        throw std::runtime_error("HQE gamma parameter must be positive.");
    }
    return pre;
}

double sample_qe(double m, double s2, std::mt19937& gen) {
    if (!(m > 0.0) || !(s2 > 0.0)) {
        return std::max(0.0, m);
    }

    static thread_local std::normal_distribution<double> nd(0.0, 1.0);
    static thread_local std::uniform_real_distribution<double> ud(0.0, 1.0);

    const double psi = s2 / (m * m);
    if (!(psi > 0.0)) {
        return m;
    }

    if (psi >= 1.5) {
        const double u = std::max(1e-15, ud(gen));
        const double p = 2.0 / (1.0 + psi);
        if (u < p) {
            return m * (1.0 + psi) * 0.5 * std::log(p / u);
        }
        return 0.0;
    }

    const double two_over_psi = 2.0 / psi;
    const double beta = std::sqrt(
        two_over_psi - 1.0 + std::sqrt(two_over_psi) * std::sqrt(two_over_psi - 1.0)
    );
    const double z = nd(gen);
    const double x = beta + z;
    return m / (1.0 + beta * beta) * x * x;
}

double simulate_hqe_terminal_spot(
    const HQEPrecomputed& pre,
    std::mt19937& gen,
    std::vector<double>& chi
) {
    std::fill(chi.begin(), chi.end(), 0.0);

    double v_old = pre.xi_grid[0];
    double x_old = pre.logS0;
    const double eps = 1e-6;
    static thread_local std::normal_distribution<double> nd(0.0, 1.0);

    for (int n = 1; n <= pre.n_steps; ++n) {
        double xi_hat = pre.xi_grid[n];
        for (int k = 1; k < n; ++k) {
            xi_hat += pre.b_star[n - k + 1] * chi[k];
        }
        xi_hat = std::max(eps, xi_hat);

        const double v_bar = (xi_hat + 2.0 * pre.H * v_old) / (2.0 * pre.H + 1.0);
        const double m = 0.5 * xi_hat;
        const double bchi_tilde = sample_qe(m, pre.beta * pre.beta * v_bar * pre.dt, gen);
        const double eps_tilde = sample_qe(m, v_bar * pre.gamma, gen);

        const double chi_new = (bchi_tilde - 0.5 * xi_hat) / pre.beta;
        const double v_new = bchi_tilde + eps_tilde;
        const double v_mean = 0.5 * (v_old + v_new);

        x_old +=
            -0.5 * v_mean * pre.dt
            + pre.rho * chi_new
            + pre.rho_perp * std::sqrt(std::max(0.0, v_mean * pre.dt)) * nd(gen);

        chi[n] = chi_new;
        v_old = v_new;
    }

    return std::exp(x_old);
}

SliceStats run_hqe_slice_mc(
    int n_sims,
    int n_steps,
    const HQEParams& p,
    const std::vector<double>& strikes,
    const std::vector<double>& reference_ivs,
    unsigned int requested_threads,
    unsigned int base_seed = 123456u
) {
    const HQEPrecomputed pre = build_hqe_precomputed(p, n_steps);
    auto t0 = std::chrono::high_resolution_clock::now();

    struct PartialStats {
        std::vector<double> sum;
        TerminalSpotMoments spot;
        explicit PartialStats(size_t n = 0) : sum(n, 0.0) {}
    };

    const unsigned int num_threads = checked_worker_threads(n_sims, requested_threads);
    const int sims_per_thread =
        (n_sims + static_cast<int>(num_threads) - 1) / static_cast<int>(num_threads);

    std::vector<PartialStats> partials(num_threads, PartialStats(strikes.size()));
    std::vector<std::future<void>> futures;

    for (unsigned int tid = 0; tid < num_threads; ++tid) {
        const int L = static_cast<int>(tid) * sims_per_thread;
        const int R = std::min(n_sims, static_cast<int>((tid + 1) * sims_per_thread));
        if (L >= R) break;

        futures.push_back(std::async(std::launch::async, [&, L, R, tid]() {
            std::mt19937 gen(base_seed + 7919u * tid);
            std::vector<double> chi(pre.n_steps + 1, 0.0);
            PartialStats local(strikes.size());
            for (int i = L; i < R; ++i) {
                const double ST = simulate_hqe_terminal_spot(pre, gen, chi);
                local.spot.add(pre.discount * ST);
                for (size_t k = 0; k < strikes.size(); ++k) {
                    local.sum[k] += pre.discount * european_call_payoff(ST, strikes[k]);
                }
            }
            partials[tid] = std::move(local);
        }));
    }
    for (auto& f : futures) f.get();

    auto t1 = std::chrono::high_resolution_clock::now();
    const double elapsed = std::chrono::duration<double>(t1 - t0).count();

    std::vector<double> prices(strikes.size(), 0.0);
    TerminalSpotMoments spot_moments;
    for (const PartialStats& partial : partials) {
        for (size_t k = 0; k < strikes.size(); ++k) prices[k] += partial.sum[k];
        spot_moments.merge(partial.spot);
    }
    for (double& x : prices) x /= static_cast<double>(n_sims);

    SliceStats out = finalize_slice_stats(
        "HQE-RoughHeston", n_steps, -1, prices,
        strikes, p.S0, p.T, p.r, reference_ivs, elapsed
    );
    set_spot_diagnostics(out, n_sims, num_threads, p.S0, spot_moments);
    return out;
}

struct ReferenceSlice {
    double alpha;
    std::vector<double> prices;
};

int main(int argc, char** argv) {
    const int n_sims = 1000000;
    bool short_maturity = false;
    unsigned int num_threads = 10;
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "--t-1-12" || arg == "--short") {
            short_maturity = true;
        } else if (arg == "--threads") {
            if (i + 1 >= argc) {
                throw std::runtime_error("Missing value for --threads.");
            }
            const long parsed = std::stol(argv[++i]);
            if (parsed <= 0 || static_cast<unsigned long>(parsed) > std::numeric_limits<unsigned int>::max()) {
                throw std::runtime_error("--threads must be a positive unsigned integer.");
            }
            num_threads = static_cast<unsigned int>(parsed);
        }
    }

    const double T = short_maturity ? (1.0 / 12.0) : 1.0;
    const double r = 0.0;
    const double gamma_val = 0.1;
    const double theta = 0.3156;
    const double nu = 0.331;
    const double rho = -0.681;
    const double V0 = 0.0392;
    const double S0 = 100.0;

    const std::vector<double> strikes = short_maturity
        ? std::vector<double>{90.0, 95.0, 100.0, 105.0, 110.0}
        : std::vector<double>{80.0, 90.0, 100.0, 110.0, 120.0};
    const std::vector<int> steps_grid = {40, 80, 160, 320};
    const std::vector<double> alphas = {0.55, 0.62, 0.80, 0.95};
    const int integrated_factors = 1000;
    const double x_min = 1e-3;
    const double x_max = 1e3;

    const std::vector<ReferenceSlice> references = short_maturity
        ? std::vector<ReferenceSlice>{
            {0.55, {10.1267556887, 5.6989890033, 2.4153669912, 0.6906980177, 0.1210835182}},
            {0.62, {10.1143601668, 5.6723586879, 2.3896783413, 0.6809088652, 0.1205286640}},
            {0.80, {10.0934043360, 5.6241329401, 2.3423581424, 0.6629959068, 0.1195930284}},
            {0.95, {10.0836030770, 5.5998574438, 2.3181367808, 0.6539010658, 0.1191571922}}
        }
        : std::vector<ReferenceSlice>{
            {0.55, {22.1843478272, 15.0311228460, 9.5419413959, 5.6843618720, 3.1895226582}},
            {0.62, {22.1366349471, 14.9672185800, 9.4737132680, 5.6234557195, 3.1424569602}},
            {0.80, {22.0140617545, 14.8010509346, 9.2954508931, 5.4643254351, 3.0198686621}},
            {0.95, {21.9144690238, 14.6637176595, 9.1471029750, 5.3318627192, 2.9182287550}}
        };

    std::cout << "Implied-volatility slice convergence comparison\n\n";
    std::cout << "paths   = " << n_sims << "\n";
    std::cout << "threads = " << num_threads << "\n";
    std::cout << "T       = " << T << "\n";
    std::cout << "S0      = " << S0 << "\n";
    std::cout << "r       = " << r << "\n";
    std::cout << "gamma   = " << gamma_val << "\n";
    std::cout << "theta   = " << theta << "\n";
    std::cout << "nu      = " << nu << "\n";
    std::cout << "rho     = " << rho << "\n";
    std::cout << "V0      = " << V0 << "\n";
    if (short_maturity) {
        std::cout << "strikes = {90, 95, 100, 105, 110}\n";
    } else {
        std::cout << "strikes = {80, 90, 100, 110, 120}\n";
    }
    std::cout << "steps   = {40, 80, 160, 320}\n";
    std::cout << "integrated_factors = " << integrated_factors << "\n";
    std::cout << "reference = full rough-Heston Fourier (precomputed)\n";
    std::cout << std::endl;

    std::cout << std::fixed << std::setprecision(6);

    for (size_t a_idx = 0; a_idx < alphas.size(); ++a_idx) {
        const double alpha = alphas[a_idx];
        const std::vector<double>& ref_prices = references[a_idx].prices;
        std::vector<double> ref_ivs(strikes.size(), 0.0);
        for (size_t k = 0; k < strikes.size(); ++k) {
            ref_ivs[k] = implied_vol_call(ref_prices[k], S0, strikes[k], T, r);
        }

        RHParams rh{
            alpha - 0.5,
            gamma_val,
            theta,
            nu,
            rho,
            V0,
            S0,
            T,
            r
        };
        const double beta_param = solve_beta(rho);
        const double mu = solve_mu(nu, theta, beta_param, gamma_val);
        const double xi = solve_xi(V0, theta);
        IVIParams ivi_params{alpha, gamma_val, theta, nu, rho, V0, S0, T, r};
        HQEParams hqe_params{alpha, gamma_val, theta, nu, rho, V0, S0, T, r};
        KernelApprox ka = build_loggrid_kernel_approx(rh.H, rh.lambda, integrated_factors, x_min, x_max);

        std::cout << "alpha = " << alpha << "  (H = " << (alpha - 0.5) << ")\n";
        std::cout << "reference prices:";
        for (double x : ref_prices) std::cout << " " << x;
        std::cout << "\n";
        std::cout << "reference ivs   :";
        for (double x : ref_ivs) std::cout << " " << x;
        std::cout << "\n\n";

        std::cout << std::left
                  << std::setw(28) << "Method"
                  << std::setw(10) << "Steps"
                  << std::setw(10) << "Factors"
                  << std::setw(16) << "MaxIVErr"
                  << std::setw(16) << "RMSIVErr"
                  << std::setw(14) << "ATM_IV"
                  << std::setw(12) << "Time(s)"
                  << "\n";
        std::cout << std::string(106, '-') << "\n";

        for (int steps : steps_grid) {
            const SliceStats st = run_integrated_slice_mc(
                n_sims, steps, rh, ka, strikes, ref_ivs, integrated_factors,
                num_threads,
                100000u + static_cast<unsigned int>(steps) + 1000u * static_cast<unsigned int>(a_idx)
            );
            std::cout << std::left
                      << std::setw(28) << st.method
                      << std::setw(10) << st.steps
                      << std::setw(10) << st.factors
                      << std::setw(16) << st.max_iv_error
                      << std::setw(16) << st.rms_iv_error
                      << std::setw(14) << st.ivs[2]
                      << std::setw(12) << st.seconds
                      << "\n";
        }
        for (int steps : steps_grid) {
            const SliceStats st = run_inar_slice_mc(
                n_sims, steps, strikes, alpha, gamma_val, beta_param, mu, theta, xi, S0, T, r, ref_ivs,
                num_threads,
                200000u + static_cast<unsigned int>(steps) + 1000u * static_cast<unsigned int>(a_idx)
            );
            std::cout << std::left
                      << std::setw(28) << st.method
                      << std::setw(10) << st.steps
                      << std::setw(10) << "-"
                      << std::setw(16) << st.max_iv_error
                      << std::setw(16) << st.rms_iv_error
                      << std::setw(14) << st.ivs[2]
                      << std::setw(12) << st.seconds
                      << "\n";
        }
        for (int steps : steps_grid) {
            const SliceStats st = run_ivi_slice_mc(
                n_sims, steps, ivi_params, strikes, ref_ivs,
                num_threads,
                300000u + static_cast<unsigned int>(steps) + 1000u * static_cast<unsigned int>(a_idx)
            );
            std::cout << std::left
                      << std::setw(28) << st.method
                      << std::setw(10) << st.steps
                      << std::setw(10) << "-"
                      << std::setw(16) << st.max_iv_error
                      << std::setw(16) << st.rms_iv_error
                      << std::setw(14) << st.ivs[2]
                      << std::setw(12) << st.seconds
                      << "\n";
        }
        for (int steps : steps_grid) {
            const SliceStats st = run_hqe_slice_mc(
                n_sims, steps, hqe_params, strikes, ref_ivs,
                num_threads,
                400000u + static_cast<unsigned int>(steps) + 1000u * static_cast<unsigned int>(a_idx)
            );
            std::cout << std::left
                      << std::setw(28) << st.method
                      << std::setw(10) << st.steps
                      << std::setw(10) << "-"
                      << std::setw(16) << st.max_iv_error
                      << std::setw(16) << st.rms_iv_error
                      << std::setw(14) << st.ivs[2]
                      << std::setw(12) << st.seconds
                      << "\n";
        }
        std::cout << "\n";
    }

    return 0;
}
