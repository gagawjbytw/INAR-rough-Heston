// Reproducible implied-volatility surface generator for the paper figures.
//
// This driver deliberately reuses the production CDQ--FFT INAR core from
// inar_cdq_all_options_v2.cpp.  In particular, paths use independent exact
// Poisson draws and the exact one-step risk-neutral compensation
//
//     d_tau = c_tau - log(2 - exp(-c_tau)).
//
// Build:
//   clang++ -O3 -std=c++17 -pthread generate_iv_surface.cpp \
//       -o generate_iv_surface
//
// Paper run:
//   ./generate_iv_surface --paths 1000000 --tau 320 \
//       --output output/surface_exact.csv

#define main inar_all_options_unused_main
#include "../INAR-tables/inar_cdq_all_options_v2.cpp"
#undef main

#include <atomic>
#include <filesystem>
#include <fstream>
#include <limits>
#include <sstream>

namespace {

struct Config {
    std::uint64_t paths = 1000000;
    int tau = 320;
    int threads = static_cast<int>(std::max(1u, std::thread::hardware_concurrency()));
    std::uint64_t seed = 20260811ULL;
    std::string output = "code/iv_surface_output/surface_exact.csv";
    std::string models = "both";

    double rough_alpha = 0.62;
    double gamma = 0.1;
    double theta = 0.3156;
    double rho = -0.681;
    double nu = 0.331;
    double V0 = 0.0392;
    double S0 = 100.0;
};

struct ThreadStats {
    std::vector<long double> payoff_sum;
    std::vector<long double> payoff_sq_sum;
    std::vector<long double> spot_sum;
    std::vector<long double> spot_sq_sum;
    std::uint64_t paths = 0;
    std::uint64_t negative_lambda_clamps = 0;
    double min_raw_lambda = std::numeric_limits<double>::infinity();
    double max_raw_lambda = 0.0;

    ThreadStats(std::size_t cells, std::size_t maturities)
        : payoff_sum(cells, 0.0L),
          payoff_sq_sum(cells, 0.0L),
          spot_sum(maturities, 0.0L),
          spot_sq_sum(maturities, 0.0L) {}
};

struct ModelResult {
    std::string model;
    double alpha = 0.0;
    double beta = 0.0;
    double mu = 0.0;
    double xi = 0.0;
    double c_tau = 0.0;
    double d_tau = 0.0;
    std::uint64_t effective_seed = 0;
    double seconds = 0.0;
    std::uint64_t negative_lambda_clamps = 0;
    double min_raw_lambda = 0.0;
    double max_raw_lambda = 0.0;
    std::vector<double> call_price;
    std::vector<double> call_se;
    std::vector<double> implied_vol;
    std::vector<double> spot_mean;
    std::vector<double> spot_se;
};

double normal_cdf_surface(double x) {
    return 0.5 * std::erfc(-x / std::sqrt(2.0));
}

double bs_call_surface(double S0, double K, double T, double sigma) {
    if (T <= 0.0 || sigma <= 0.0) return std::max(S0 - K, 0.0);
    const double sqrt_T = std::sqrt(T);
    const double d1 = (std::log(S0 / K) + 0.5 * sigma * sigma * T) / (sigma * sqrt_T);
    const double d2 = d1 - sigma * sqrt_T;
    return S0 * normal_cdf_surface(d1) - K * normal_cdf_surface(d2);
}

double implied_vol_surface(double call_price, double S0, double K, double T) {
    const double intrinsic = std::max(S0 - K, 0.0);
    if (!(call_price >= intrinsic && call_price < S0)) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    if (call_price <= intrinsic + 1e-14) return 0.0;
    double lo = 1e-8;
    double hi = 5.0;
    if (bs_call_surface(S0, K, T, hi) < call_price) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    for (int iteration = 0; iteration < 100; ++iteration) {
        const double mid = 0.5 * (lo + hi);
        if (bs_call_surface(S0, K, T, mid) < call_price) lo = mid;
        else hi = mid;
    }
    return 0.5 * (lo + hi);
}

std::vector<double> make_log_moneyness_grid() {
    std::vector<double> grid;
    for (int j = -10; j <= 10; ++j) grid.push_back(0.02 * static_cast<double>(j));
    return grid;
}

std::vector<double> make_maturity_grid() {
    std::vector<double> grid;
    for (int month = 1; month <= 12; ++month) {
        grid.push_back(static_cast<double>(month) / 12.0);
    }
    return grid;
}

std::vector<int> make_checkpoints(const std::vector<double>& maturities, int tau) {
    std::vector<int> checkpoints;
    int previous = 0;
    for (double maturity : maturities) {
        const int step = static_cast<int>(std::floor(maturity * static_cast<double>(tau) + 1e-12));
        if (step <= previous) throw std::runtime_error("tau is too small for distinct monthly checkpoints");
        checkpoints.push_back(step);
        previous = step;
    }
    return checkpoints;
}

void simulate_surface_path(
    const InarPrecomputed& pre,
    double S0,
    const std::vector<int>& checkpoints,
    std::mt19937& generator,
    InarWorkspace& workspace,
    std::poisson_distribution<int>& poisson,
    std::vector<double>& checkpoint_prices,
    ThreadStats& diagnostics
) {
    workspace.reset(S0, std::numeric_limits<double>::infinity(),
                    -std::numeric_limits<double>::infinity());
    std::size_t checkpoint_index = 0;

    for (const InarPrecomputed::Op& op : pre.ops) {
        if (op.is_leaf) {
            const double raw_lambda = pre.lam_base[op.t] + workspace.lam_conv[op.t];
            diagnostics.min_raw_lambda = std::min(diagnostics.min_raw_lambda, raw_lambda);
            diagnostics.max_raw_lambda = std::max(diagnostics.max_raw_lambda, raw_lambda);
            if (raw_lambda < 0.0) {
                ++diagnostics.negative_lambda_clamps;
                throw std::runtime_error("negative raw Poisson intensity encountered");
            }
            const double lambda = raw_lambda;
            const auto poisson_parameter = std::poisson_distribution<int>::param_type(lambda);
            const int x_plus = poisson(generator, poisson_parameter);
            const int x_minus = poisson(generator, poisson_parameter);
            workspace.total_plus += x_plus;
            workspace.total_minus += x_minus;
            workspace.Y[op.t] = static_cast<double>(x_plus)
                + pre.beta_param * static_cast<double>(x_minus);

            if (checkpoint_index < checkpoints.size()
                && op.t == checkpoints[checkpoint_index]) {
                const double count_difference =
                    static_cast<double>(workspace.total_plus - workspace.total_minus);
                const double log_price = pre.coef * count_difference
                    - pre.drift_coef * static_cast<double>(workspace.total_plus);
                checkpoint_prices[checkpoint_index] = S0 * std::exp(log_price);
                if (!std::isfinite(checkpoint_prices[checkpoint_index])) {
                    throw std::runtime_error("non-finite checkpoint spot encountered");
                }
                ++checkpoint_index;
            }
            continue;
        }

        if (op.direct) {
            for (int i = 0; i < op.nA; ++i) {
                const double source = workspace.Y[op.L + i];
                const int source_time = op.L + i;
                for (int t = op.M + 1; t <= op.R; ++t) {
                    workspace.lam_conv[t] += source * pre.kernel[t - source_time];
                }
            }
            continue;
        }

        for (int i = 0; i < op.nA; ++i) {
            workspace.fft_a[i] = workspace.Y[op.L + i];
        }
        std::fill(workspace.fft_a.begin() + op.nA,
                  workspace.fft_a.begin() + op.N,
                  FFT::Complex(0.0, 0.0));
        const FFT::Plan& plan = pre.fft_plan_cache[op.N];
        FFT::fft_inplace(workspace.fft_a.data(), plan, false);
        const std::vector<FFT::Complex>& transformed_kernel = pre.kernel_fft_cache[op.nB];
        for (int i = 0; i < op.N; ++i) workspace.fft_a[i] *= transformed_kernel[i];
        FFT::fft_inplace(workspace.fft_a.data(), plan, true);
        for (int t = op.M + 1; t <= op.R; ++t) {
            const int index = t - op.L - 1;
            if (index >= 0 && index < op.N) {
                workspace.lam_conv[t] += workspace.fft_a[index].real();
            }
        }
    }

    if (checkpoint_index != checkpoints.size()) {
        throw std::runtime_error("not all maturity checkpoints were recorded");
    }
}

ModelResult run_model(
    const Config& config,
    const std::string& model_name,
    double alpha,
    const std::vector<double>& maturities,
    const std::vector<int>& checkpoints,
    const std::vector<double>& log_moneyness,
    const std::vector<double>& strikes,
    std::uint64_t seed_offset
) {
    const double beta = solve_beta(config.rho);
    const double mu = solve_mu(config.nu, config.theta, beta, config.gamma);
    const double xi = solve_xi(config.V0, config.theta);
    const InarPrecomputed pre = build_inar_precomputed(
        alpha, config.gamma, config.tau, beta, mu, config.theta, xi);

    const std::size_t maturity_count = maturities.size();
    const std::size_t strike_count = strikes.size();
    const std::size_t cell_count = maturity_count * strike_count;
    const int thread_count = std::max(1, std::min<int>(
        config.threads, static_cast<int>(config.paths)));
    std::vector<ThreadStats> partials;
    partials.reserve(static_cast<std::size_t>(thread_count));
    for (int thread = 0; thread < thread_count; ++thread) {
        partials.emplace_back(cell_count, maturity_count);
    }

    const auto start_time = std::chrono::steady_clock::now();
    std::vector<std::future<void>> futures;
    for (int thread = 0; thread < thread_count; ++thread) {
        const std::uint64_t begin = config.paths * static_cast<std::uint64_t>(thread)
            / static_cast<std::uint64_t>(thread_count);
        const std::uint64_t end = config.paths * static_cast<std::uint64_t>(thread + 1)
            / static_cast<std::uint64_t>(thread_count);
        futures.push_back(std::async(std::launch::async, [&, thread, begin, end, seed_offset]() {
            ThreadStats& stats = partials[static_cast<std::size_t>(thread)];
            const std::uint64_t thread_seed = config.seed + seed_offset
                + 1000003ULL * static_cast<std::uint64_t>(thread);
            std::seed_seq seed_sequence{
                static_cast<std::uint32_t>(thread_seed),
                static_cast<std::uint32_t>(thread_seed >> 32U),
                static_cast<std::uint32_t>(alpha * 1000000.0)
            };
            std::mt19937 generator(seed_sequence);
            std::poisson_distribution<int> poisson;
            InarWorkspace workspace(config.tau);
            std::vector<double> checkpoint_prices(maturity_count, config.S0);

            for (std::uint64_t path = begin; path < end; ++path) {
                simulate_surface_path(pre, config.S0, checkpoints, generator,
                                      workspace, poisson, checkpoint_prices, stats);
                for (std::size_t maturity_index = 0;
                     maturity_index < maturity_count; ++maturity_index) {
                    const double spot = checkpoint_prices[maturity_index];
                    stats.spot_sum[maturity_index] += spot;
                    stats.spot_sq_sum[maturity_index] += spot * spot;
                    for (std::size_t strike_index = 0;
                         strike_index < strike_count; ++strike_index) {
                        const double strike = strikes[strike_index];
                        // Use the OTM payoff.  For K<S0, exact put--call parity
                        // C=P+S0-K greatly reduces deep-ITM Monte Carlo noise.
                        const double payoff = strike < config.S0
                            ? std::max(strike - spot, 0.0)
                            : std::max(spot - strike, 0.0);
                        const std::size_t index = maturity_index * strike_count + strike_index;
                        stats.payoff_sum[index] += payoff;
                        stats.payoff_sq_sum[index] += payoff * payoff;
                    }
                }
                ++stats.paths;
            }
        }));
    }
    for (auto& future : futures) future.get();
    const double elapsed = std::chrono::duration<double>(
        std::chrono::steady_clock::now() - start_time).count();

    ModelResult result;
    result.model = model_name;
    result.alpha = alpha;
    result.beta = beta;
    result.mu = mu;
    result.xi = xi;
    result.c_tau = pre.coef;
    result.d_tau = pre.drift_coef;
    result.effective_seed = config.seed + seed_offset;
    result.seconds = elapsed;
    result.min_raw_lambda = std::numeric_limits<double>::infinity();
    result.call_price.resize(cell_count);
    result.call_se.resize(cell_count);
    result.implied_vol.resize(cell_count);
    result.spot_mean.resize(maturity_count);
    result.spot_se.resize(maturity_count);

    std::uint64_t total_paths = 0;
    for (const ThreadStats& partial : partials) {
        total_paths += partial.paths;
        result.negative_lambda_clamps += partial.negative_lambda_clamps;
        result.min_raw_lambda = std::min(result.min_raw_lambda, partial.min_raw_lambda);
        result.max_raw_lambda = std::max(result.max_raw_lambda, partial.max_raw_lambda);
    }
    if (total_paths != config.paths) throw std::runtime_error("path accounting mismatch");
    if (result.negative_lambda_clamps != 0) {
        throw std::runtime_error(model_name + ": negative raw Poisson intensity encountered");
    }
    const long double sample_count = static_cast<long double>(total_paths);

    for (std::size_t maturity_index = 0;
         maturity_index < maturity_count; ++maturity_index) {
        long double sum = 0.0L;
        long double sum_sq = 0.0L;
        for (const ThreadStats& partial : partials) {
            sum += partial.spot_sum[maturity_index];
            sum_sq += partial.spot_sq_sum[maturity_index];
        }
        const long double mean = sum / sample_count;
        const long double variance = total_paths > 1
            ? std::max(0.0L, (sum_sq - sum * sum / sample_count)
                / static_cast<long double>(total_paths - 1))
            : 0.0L;
        result.spot_mean[maturity_index] = static_cast<double>(mean);
        result.spot_se[maturity_index] = static_cast<double>(
            std::sqrt(variance / sample_count));

        for (std::size_t strike_index = 0;
             strike_index < strike_count; ++strike_index) {
            const std::size_t index = maturity_index * strike_count + strike_index;
            long double payoff_sum = 0.0L;
            long double payoff_sq_sum = 0.0L;
            for (const ThreadStats& partial : partials) {
                payoff_sum += partial.payoff_sum[index];
                payoff_sq_sum += partial.payoff_sq_sum[index];
            }
            const long double payoff_mean = payoff_sum / sample_count;
            const long double payoff_variance = total_paths > 1
                ? std::max(0.0L, (payoff_sq_sum - payoff_sum * payoff_sum / sample_count)
                    / static_cast<long double>(total_paths - 1))
                : 0.0L;
            const double strike = strikes[strike_index];
            const double call = static_cast<double>(payoff_mean)
                + (strike < config.S0 ? config.S0 - strike : 0.0);
            result.call_price[index] = call;
            result.call_se[index] = static_cast<double>(
                std::sqrt(payoff_variance / sample_count));
            result.implied_vol[index] = implied_vol_surface(
                call, config.S0, strike, maturities[maturity_index]);
            if (!std::isfinite(result.implied_vol[index])) {
                std::ostringstream message;
                message << model_name << ": implied-vol inversion failed at T="
                        << maturities[maturity_index] << ", K=" << strike
                        << ", call=" << call;
                throw std::runtime_error(message.str());
            }
        }
    }
    return result;
}

void write_results(
    const Config& config,
    const std::vector<double>& maturities,
    const std::vector<int>& checkpoints,
    const std::vector<double>& log_moneyness,
    const std::vector<double>& strikes,
    const std::vector<ModelResult>& results
) {
    const std::filesystem::path output_path(config.output);
    if (output_path.has_parent_path()) {
        std::filesystem::create_directories(output_path.parent_path());
    }
    std::ofstream output(config.output);
    if (!output) throw std::runtime_error("cannot open output CSV: " + config.output);
    output << std::setprecision(17);
    output << "model,alpha,gamma,theta,rho,nu,V0,S0,beta,mu,xi,tau,threads,paths,seed,c_tau,d_tau,"
              "maturity,step,effective_grid_time,log_moneyness,strike,call_price,call_se,implied_vol,spot_mean,spot_se,"
              "seconds,negative_lambda_clamps,min_raw_lambda,max_raw_lambda,estimator\n";
    for (const ModelResult& result : results) {
        for (std::size_t maturity_index = 0;
             maturity_index < maturities.size(); ++maturity_index) {
            for (std::size_t strike_index = 0;
                 strike_index < strikes.size(); ++strike_index) {
                const std::size_t index = maturity_index * strikes.size() + strike_index;
                output << result.model << ',' << result.alpha << ',' << config.gamma << ','
                       << config.theta << ',' << config.rho << ',' << config.nu << ','
                       << config.V0 << ',' << config.S0 << ',' << result.beta << ','
                       << result.mu << ',' << result.xi << ',' << config.tau << ','
                       << config.threads << ',' << config.paths << ',' << result.effective_seed << ','
                       << result.c_tau << ','
                       << result.d_tau << ',' << maturities[maturity_index] << ','
                       << checkpoints[maturity_index] << ','
                       << static_cast<double>(checkpoints[maturity_index]) / config.tau << ','
                       << log_moneyness[strike_index] << ','
                       << strikes[strike_index] << ',' << result.call_price[index] << ','
                       << result.call_se[index] << ',' << result.implied_vol[index] << ','
                       << result.spot_mean[maturity_index] << ',' << result.spot_se[maturity_index]
                       << ',' << result.seconds << ',' << result.negative_lambda_clamps << ','
                       << result.min_raw_lambda << ',' << result.max_raw_lambda
                       << ",otm_payoff_plus_exact_parity\n";
            }
        }
    }
}

Config parse_arguments(int argc, char** argv) {
    Config config;
    for (int i = 1; i < argc; ++i) {
        const std::string argument(argv[i]);
        auto require_value = [&](const std::string& option) -> std::string {
            if (i + 1 >= argc) throw std::runtime_error("missing value after " + option);
            return std::string(argv[++i]);
        };
        if (argument == "--paths") config.paths = std::stoull(require_value(argument));
        else if (argument == "--tau") config.tau = std::stoi(require_value(argument));
        else if (argument == "--threads") config.threads = std::stoi(require_value(argument));
        else if (argument == "--seed") config.seed = std::stoull(require_value(argument));
        else if (argument == "--output") config.output = require_value(argument);
        else if (argument == "--models") config.models = require_value(argument);
        else if (argument == "--help") {
            std::cout
                << "Usage: generate_iv_surface [--paths N] [--tau N] [--threads N] "
                   "[--seed N] [--models both|rough|classical] [--output FILE]\n";
            std::exit(0);
        } else {
            throw std::runtime_error("unknown argument: " + argument);
        }
    }
    if (config.paths == 0) throw std::runtime_error("paths must be positive");
    if (config.tau < 12) throw std::runtime_error("tau must be at least 12");
    if (config.threads <= 0) throw std::runtime_error("threads must be positive");
    if (config.models != "both" && config.models != "rough"
        && config.models != "classical") {
        throw std::runtime_error("models must be both, rough, or classical");
    }
    return config;
}

} // namespace

int main(int argc, char** argv) {
    try {
        const Config config = parse_arguments(argc, argv);
        const std::vector<double> maturities = make_maturity_grid();
        const std::vector<int> checkpoints = make_checkpoints(maturities, config.tau);
        const std::vector<double> log_moneyness = make_log_moneyness_grid();
        std::vector<double> strikes;
        for (double k : log_moneyness) strikes.push_back(config.S0 * std::exp(k));

        std::vector<ModelResult> results;
        if (config.models == "both" || config.models == "rough") {
            results.push_back(run_model(config, "rough", config.rough_alpha,
                                        maturities, checkpoints, log_moneyness, strikes, 0ULL));
        }
        if (config.models == "both" || config.models == "classical") {
            results.push_back(run_model(config, "classical", 1.0,
                                        maturities, checkpoints, log_moneyness, strikes,
                                        1000000007ULL));
        }
        write_results(config, maturities, checkpoints, log_moneyness, strikes, results);

        std::cout << std::fixed << std::setprecision(8);
        std::cout << "output=" << config.output << " paths=" << config.paths
                  << " tau=" << config.tau << " threads=" << config.threads << '\n';
        for (const ModelResult& result : results) {
            const std::size_t atm_index = 10;
            const std::size_t last_maturity_offset = 11 * strikes.size();
            std::cout << result.model << ": alpha=" << result.alpha
                      << " beta=" << result.beta << " mu=" << result.mu
                      << " xi=" << result.xi << " c_tau=" << result.c_tau
                      << " d_tau=" << result.d_tau << " seconds=" << result.seconds
                      << " E[S_1]=" << result.spot_mean.back()
                      << " (SE=" << result.spot_se.back() << ')'
                      << " IV_ATM_1=" << result.implied_vol[last_maturity_offset + atm_index]
                      << " clamps=" << result.negative_lambda_clamps
                      << " lambda_range=[" << result.min_raw_lambda << ','
                      << result.max_raw_lambda << "]\n";
        }
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "error: " << error.what() << '\n';
        return 1;
    }
}
