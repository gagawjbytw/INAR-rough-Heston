#include <algorithm>
#include <chrono>
#include <cmath>
#include <future>
#include <iomanip>
#include <iostream>
#include <random>
#include <stdexcept>
#include <thread>
#include <tuple>
#include <vector>

using namespace std;

void validate_parameters(double S0, double V0, double r, double gamma,
                         double theta, double sigma, double rho, double dt) {
    if (S0 <= 0.0) throw runtime_error("Initial stock price must be positive");
    if (V0 < 0.0) throw runtime_error("Initial variance must be non-negative");
    if (gamma <= 0.0) throw runtime_error("Mean reversion speed must be positive");
    if (theta <= 0.0) throw runtime_error("Long-term variance must be positive");
    if (sigma <= 0.0) throw runtime_error("Volatility of volatility must be positive");
    if (rho <= -1.0 || rho >= 1.0) throw runtime_error("Correlation must lie in (-1,1)");
    if (dt <= 0.0) throw runtime_error("Time step must be positive");

    if (2.0 * gamma * theta <= sigma * sigma) {
        cout << "Warning: Feller condition not satisfied, variance may hit zero" << endl;
    }
}

class NormalGenerator {
private:
    mt19937 gen;
    normal_distribution<double> dist;

public:
    explicit NormalGenerator(int seed) : gen(seed), dist(0.0, 1.0) {}
    double operator()() { return dist(gen); }
};

struct PathSummary {
    double final_price = 0.0;
    double avg_price = 0.0;
    double max_price = 0.0;
    double min_price = 0.0;
    bool barrier_hit_110 = false;
    bool barrier_hit_90 = false;
};

struct PathPairSummary {
    PathSummary forward;
    PathSummary antithetic;
};

PathPairSummary simulate_heston_summaries(double S0, double V0, double r, double gamma,
                                          double theta, double sigma, double rho,
                                          double T, double dt, int seed,
                                          double barrier_up, double barrier_down) {
    const int M = static_cast<int>(std::llround(T / dt));
    const double sqrt_dt = std::sqrt(dt);

    NormalGenerator rng(seed);

    double S = S0, V = V0;
    double S_anti = S0, V_anti = V0;

    double sum_S = S0, sum_S_anti = S0;
    double max_S = S0, max_S_anti = S0;
    double min_S = S0, min_S_anti = S0;
    bool hit_up = S0 >= barrier_up;
    bool hit_up_anti = S0 >= barrier_up;
    bool hit_down = S0 <= barrier_down;
    bool hit_down_anti = S0 <= barrier_down;

    for (int i = 0; i < M; ++i) {
        const double Z1 = rng();
        const double Z2 = rho * Z1 + std::sqrt(1.0 - rho * rho) * rng();

        const double sqrtV = std::sqrt(std::max(V, 0.0));
        const double dV = gamma * (theta - V) * dt + sigma * sqrtV * sqrt_dt * Z2;
        V = std::max(V + dV, 0.0);
        const double dS = r * S * dt + sqrtV * S * sqrt_dt * Z1;
        S = std::max(S + dS, 0.0);

        const double sqrtV_anti = std::sqrt(std::max(V_anti, 0.0));
        const double dV_anti = gamma * (theta - V_anti) * dt - sigma * sqrtV_anti * sqrt_dt * Z2;
        V_anti = std::max(V_anti + dV_anti, 0.0);
        const double dS_anti = r * S_anti * dt - sqrtV_anti * S_anti * sqrt_dt * Z1;
        S_anti = std::max(S_anti + dS_anti, 0.0);

        sum_S += S;
        sum_S_anti += S_anti;
        max_S = std::max(max_S, S);
        max_S_anti = std::max(max_S_anti, S_anti);
        min_S = std::min(min_S, S);
        min_S_anti = std::min(min_S_anti, S_anti);
        hit_up = hit_up || (S >= barrier_up);
        hit_up_anti = hit_up_anti || (S_anti >= barrier_up);
        hit_down = hit_down || (S <= barrier_down);
        hit_down_anti = hit_down_anti || (S_anti <= barrier_down);
    }

    PathPairSummary result;
    result.forward = {
        S,
        sum_S / static_cast<double>(M + 1),
        max_S,
        min_S,
        hit_up,
        hit_down
    };
    result.antithetic = {
        S_anti,
        sum_S_anti / static_cast<double>(M + 1),
        max_S_anti,
        min_S_anti,
        hit_up_anti,
        hit_down_anti
    };
    return result;
}

double european_call_payoff(double S_T, double strike) { return max(S_T - strike, 0.0); }
double european_put_payoff(double S_T, double strike) { return max(strike - S_T, 0.0); }
double asian_call_payoff(double avg_price, double strike) { return max(avg_price - strike, 0.0); }
double asian_put_payoff(double avg_price, double strike) { return max(strike - avg_price, 0.0); }
double lookback_call_payoff(double max_price, double strike) { return max(max_price - strike, 0.0); }
double lookback_put_payoff(double min_price, double strike) { return max(strike - min_price, 0.0); }
double barrier_up_in_call_payoff(const PathSummary& result, double strike) {
    return result.barrier_hit_110 ? max(result.final_price - strike, 0.0) : 0.0;
}
double barrier_down_out_put_payoff(const PathSummary& result, double strike) {
    return !result.barrier_hit_90 ? max(strike - result.final_price, 0.0) : 0.0;
}

struct StatResult {
    double mean = 0.0;
    double ci_lower = 0.0;
    double ci_upper = 0.0;
};

struct PartialStats {
    vector<double> sum;
    vector<double> sum_sq;
    int count = 0;

    explicit PartialStats(size_t n_payoffs = 0)
        : sum(n_payoffs, 0.0), sum_sq(n_payoffs, 0.0) {}
};

vector<vector<StatResult>> parallel_monte_carlo_all_prices(
    double S0, double V0, double r, double gamma, double theta, double sigma,
    double rho, double T, double dt, const vector<double>& strikes,
    double barrier_up, double barrier_down, int n_sims, unsigned int base_seed
) {
    validate_parameters(S0, V0, r, gamma, theta, sigma, rho, dt);

    const unsigned int hc = thread::hardware_concurrency();
    const unsigned int num_threads = max(1u, hc);
    const size_t n_strikes = strikes.size();
    const size_t n_payoffs = 8 * n_strikes;

    vector<future<void>> futures;
    vector<PartialStats> partials(num_threads, PartialStats(n_payoffs));

    const int sims_per_thread = (n_sims + static_cast<int>(num_threads) - 1) / static_cast<int>(num_threads);
    for (unsigned int t = 0; t < num_threads; ++t) {
        const int start_sim = static_cast<int>(t) * sims_per_thread;
        const int end_sim = min(n_sims, static_cast<int>((t + 1) * sims_per_thread));
        if (start_sim >= end_sim) break;

        futures.push_back(async(launch::async, [&, start_sim, end_sim, t]() {
            PartialStats local(n_payoffs);
            for (int sim = start_sim; sim < end_sim; ++sim) {
                const PathPairSummary pr = simulate_heston_summaries(
                    S0, V0, r, gamma, theta, sigma, rho, T, dt,
                    static_cast<int>(base_seed + sim), barrier_up, barrier_down
                );

                for (size_t k = 0; k < n_strikes; ++k) {
                    const double K = strikes[k];
                    const double payoffs_fwd[8] = {
                        european_call_payoff(pr.forward.final_price, K),
                        european_put_payoff(pr.forward.final_price, K),
                        asian_call_payoff(pr.forward.avg_price, K),
                        asian_put_payoff(pr.forward.avg_price, K),
                        lookback_call_payoff(pr.forward.max_price, K),
                        lookback_put_payoff(pr.forward.min_price, K),
                        barrier_up_in_call_payoff(pr.forward, K),
                        barrier_down_out_put_payoff(pr.forward, K)
                    };
                    const double payoffs_anti[8] = {
                        european_call_payoff(pr.antithetic.final_price, K),
                        european_put_payoff(pr.antithetic.final_price, K),
                        asian_call_payoff(pr.antithetic.avg_price, K),
                        asian_put_payoff(pr.antithetic.avg_price, K),
                        lookback_call_payoff(pr.antithetic.max_price, K),
                        lookback_put_payoff(pr.antithetic.min_price, K),
                        barrier_up_in_call_payoff(pr.antithetic, K),
                        barrier_down_out_put_payoff(pr.antithetic, K)
                    };

                    for (int opt = 0; opt < 8; ++opt) {
                        const double payoff = 0.5 * (payoffs_fwd[opt] + payoffs_anti[opt]);
                        const size_t idx = static_cast<size_t>(opt) * n_strikes + k;
                        local.sum[idx] += payoff;
                        local.sum_sq[idx] += payoff * payoff;
                    }
                }
                local.count += 1;
            }
            partials[t] = std::move(local);
        }));
    }
    for (auto& f : futures) f.get();

    vector<vector<StatResult>> all_stats(8, vector<StatResult>(n_strikes));
    const double n = static_cast<double>(n_sims);
    const double discount = exp(-r * T);
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
                var = max(var, 0.0);
            }
            const double se = sqrt(var / n);
            all_stats[opt][k] = {
                mean * discount,
                (mean - 1.96 * se) * discount,
                (mean + 1.96 * se) * discount
            };
        }
    }
    return all_stats;
}

int main() {
    const double S0 = 100.0;
    const double V0 = 0.0392;
    const double r = 0.0;
    const double gamma = 0.1;
    const double theta = 0.3156;
    const double nu = 0.331;
    const double sigma = gamma * nu;
    const double rho = -0.681;
    const double T = 1.0;
    const double dt = 1.0 / 320.0;
    const int n_sims = 500000;
    const vector<double> strikes = {80, 90, 100, 110, 120};
    const double barrier_up = 110.0;
    const double barrier_down = 90.0;
    const unsigned int base_seed = 123456u;

    auto start_time = chrono::high_resolution_clock::now();

    cout << "Starting Heston Model Option Pricing..." << endl;
    cout << fixed << setprecision(4);
    cout << "Parameters: S0=" << S0 << ", V0=" << V0 << ", r=" << r
         << ", gamma=" << gamma << ", theta=" << theta
         << ", sigma=" << sigma << ", rho=" << rho << endl;
    cout << "Simulation settings: Paths=" << n_sims << ", Time step=" << dt
         << ", Threads=" << max(1u, thread::hardware_concurrency()) << endl;

    try {
        const auto all_stats = parallel_monte_carlo_all_prices(
            S0, V0, r, gamma, theta, sigma, rho, T, dt, strikes,
            barrier_up, barrier_down, n_sims, base_seed
        );

        cout << "\nEuropean Options:" << endl;
        cout << "Strike    Call        95% CI         Put         95% CI" << endl;
        cout << "--------------------------------------------------------" << endl;
        for (size_t i = 0; i < strikes.size(); ++i) {
            cout << setw(8) << strikes[i] << "  "
                 << setw(8) << all_stats[0][i].mean << "  ["
                 << setw(8) << all_stats[0][i].ci_lower << ", "
                 << setw(8) << all_stats[0][i].ci_upper << "]  "
                 << setw(8) << all_stats[1][i].mean << "  ["
                 << setw(8) << all_stats[1][i].ci_lower << ", "
                 << setw(8) << all_stats[1][i].ci_upper << "]" << endl;
        }

        cout << "\nAsian Options:" << endl;
        cout << "Strike    Call        95% CI         Put         95% CI" << endl;
        cout << "--------------------------------------------------------" << endl;
        for (size_t i = 0; i < strikes.size(); ++i) {
            cout << setw(8) << strikes[i] << "  "
                 << setw(8) << all_stats[2][i].mean << "  ["
                 << setw(8) << all_stats[2][i].ci_lower << ", "
                 << setw(8) << all_stats[2][i].ci_upper << "]  "
                 << setw(8) << all_stats[3][i].mean << "  ["
                 << setw(8) << all_stats[3][i].ci_lower << ", "
                 << setw(8) << all_stats[3][i].ci_upper << "]" << endl;
        }

        cout << "\nLookback Options:" << endl;
        cout << "Strike    Call        95% CI         Put         95% CI" << endl;
        cout << "--------------------------------------------------------" << endl;
        for (size_t i = 0; i < strikes.size(); ++i) {
            cout << setw(8) << strikes[i] << "  "
                 << setw(8) << all_stats[4][i].mean << "  ["
                 << setw(8) << all_stats[4][i].ci_lower << ", "
                 << setw(8) << all_stats[4][i].ci_upper << "]  "
                 << setw(8) << all_stats[5][i].mean << "  ["
                 << setw(8) << all_stats[5][i].ci_lower << ", "
                 << setw(8) << all_stats[5][i].ci_upper << "]" << endl;
        }

        cout << "\nBarrier Options (Up-barrier=110, Down-barrier=90):" << endl;
        cout << "Strike    Up-In-Call     95% CI      Down-Out-Put    95% CI" << endl;
        cout << "-------------------------------------------------------------" << endl;
        for (size_t i = 0; i < strikes.size(); ++i) {
            cout << setw(8) << strikes[i] << "  "
                 << setw(11) << all_stats[6][i].mean << "  ["
                 << setw(8) << all_stats[6][i].ci_lower << ", "
                 << setw(8) << all_stats[6][i].ci_upper << "]  "
                 << setw(11) << all_stats[7][i].mean << "  ["
                 << setw(8) << all_stats[7][i].ci_lower << ", "
                 << setw(8) << all_stats[7][i].ci_upper << "]" << endl;
        }
    } catch (const exception& e) {
        cerr << "Error occurred: " << e.what() << endl;
        return 1;
    }

    const auto end_time = chrono::high_resolution_clock::now();
    const double time_used = chrono::duration<double>(end_time - start_time).count();
    cout << "\nTotal time used: " << time_used << " seconds" << endl;
    return 0;
}
