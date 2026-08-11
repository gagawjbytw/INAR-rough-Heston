#define main inar_all_options_embedded_main
#include "inar_cdq_all_options_v2.cpp"
#undef main

#include <filesystem>
#include <fstream>

namespace {

void require_output(const std::ofstream& out, const std::filesystem::path& path) {
    if (!out) {
        throw std::runtime_error("Failed to write " + path.string());
    }
}

double elapsed_seconds(
    const std::chrono::high_resolution_clock::time_point& start,
    const std::chrono::high_resolution_clock::time_point& end
) {
    return std::chrono::duration<double>(end - start).count();
}

} // namespace

int main(int argc, char** argv) {
    int paths = 500000;
    std::filesystem::path output_dir = "results/new_compensation";
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "--paths" && i + 1 < argc) {
            paths = std::stoi(argv[++i]);
        } else if (arg == "--output-dir" && i + 1 < argc) {
            output_dir = argv[++i];
        } else {
            throw std::runtime_error("Usage: rerun_inar_tables [--paths N] [--output-dir DIR]");
        }
    }
    if (paths <= 1) throw std::runtime_error("paths must exceed one");
    std::filesystem::create_directories(output_dir);

    constexpr double gamma_val = 0.1;
    constexpr double theta = 0.3156;
    constexpr double S0 = 100.0;
    constexpr double rho = -0.681;
    constexpr double nu = 0.331;
    constexpr double V0 = 0.0392;
    constexpr double barrier_up = 110.0;
    constexpr double barrier_down = 90.0;
    constexpr unsigned int seed = 123456u;
    constexpr double benchmark = 9.4737132680;
    const unsigned int threads = std::max(1u, std::thread::hardware_concurrency());

    const double beta_param = solve_beta(rho);
    const double mu = solve_mu(nu, theta, beta_param, gamma_val);
    const double xi = solve_xi(V0, theta);

    const auto tau_path = output_dir / "tau_convergence.csv";
    std::ofstream tau_out(tau_path);
    require_output(tau_out, tau_path);
    tau_out << "alpha,T,tau,paths,threads,seed,price,ci_lower,ci_upper,benchmark,deviation,time_seconds\n";
    tau_out << std::setprecision(12);
    for (const int tau : {40, 80, 160, 320}) {
        const auto start = std::chrono::high_resolution_clock::now();
        const auto stats = parallel_monte_carlo_all_prices(
            paths, {100.0}, barrier_up, barrier_down,
            0.62, gamma_val, tau, beta_param, mu, theta, xi, S0, seed
        );
        const auto end = std::chrono::high_resolution_clock::now();
        const StatResult& call = stats[0][0];
        tau_out << "0.62,1," << tau << ',' << paths << ',' << threads << ',' << seed << ','
                << call.mean << ',' << call.ci_lower << ',' << call.ci_upper << ','
                << benchmark << ',' << (call.mean - benchmark) << ','
                << elapsed_seconds(start, end) << '\n';
        tau_out.flush();
        std::cout << "tau=" << tau << " price=" << std::setprecision(10) << call.mean
                  << " time=" << elapsed_seconds(start, end) << "s\n";
    }

    const std::vector<double> strikes = {80.0, 90.0, 100.0, 110.0, 120.0};
    const auto start = std::chrono::high_resolution_clock::now();
    const auto stats = parallel_monte_carlo_all_prices(
        paths, strikes, barrier_up, barrier_down,
        1.0, gamma_val, 320, beta_param, mu, theta, xi, S0, seed
    );
    const auto end = std::chrono::high_resolution_clock::now();
    const double heston_time = elapsed_seconds(start, end);

    const auto european_path = output_dir / "alpha_one_european.csv";
    std::ofstream european_out(european_path);
    require_output(european_out, european_path);
    european_out << "alpha,T,tau,strike,paths,threads,seed,call,call_ci_lower,call_ci_upper,put,put_ci_lower,put_ci_upper\n";
    european_out << std::setprecision(12);
    for (std::size_t i = 0; i < strikes.size(); ++i) {
        european_out << "1,1,320," << strikes[i] << ',' << paths << ',' << threads << ',' << seed << ','
                     << stats[0][i].mean << ',' << stats[0][i].ci_lower << ',' << stats[0][i].ci_upper << ','
                     << stats[1][i].mean << ',' << stats[1][i].ci_lower << ',' << stats[1][i].ci_upper << '\n';
    }

    const auto pathdep_path = output_dir / "alpha_one_path_dependent.csv";
    std::ofstream pathdep_out(pathdep_path);
    require_output(pathdep_out, pathdep_path);
    pathdep_out << "alpha,T,tau,payoff,strike,paths,threads,seed,estimate,ci_lower,ci_upper,time_seconds\n";
    pathdep_out << std::setprecision(12);
    const std::size_t k100 = 2;
    const std::vector<std::pair<std::string, int>> payoff_rows = {
        {"Asian call", 2}, {"Asian put", 3}, {"Lookback call", 4},
        {"Lookback put", 5}, {"Up-and-in call", 6}, {"Down-and-out put", 7}
    };
    for (const auto& row : payoff_rows) {
        const StatResult& value = stats[row.second][k100];
        pathdep_out << "1,1,320," << row.first << ",100," << paths << ',' << threads << ',' << seed << ','
                    << value.mean << ',' << value.ci_lower << ',' << value.ci_upper << ','
                    << heston_time << '\n';
    }

    std::cout << "alpha=1 time=" << std::setprecision(10) << heston_time << "s\n";
    std::cout << "Saved " << tau_path << '\n'
              << "Saved " << european_path << '\n'
              << "Saved " << pathdep_path << '\n';
    return 0;
}
