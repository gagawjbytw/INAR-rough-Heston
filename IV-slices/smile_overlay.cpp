#define main iv_slice_compare_embedded_main
#include "iv_slice_compare.cpp"
#undef main

#include <fstream>
#include <filesystem>
#include <sstream>

struct OverlayConfig {
    double alpha = 0.62;
    bool short_maturity = true;
    bool paper_smile = false;
    bool inar_richardson = false;
    int paths = 1 << 24;
    int integrated_steps = -1;
    int integrated_factors = 1000;
    int inar_steps = -1;
    int ivi_steps = -1;
    int hqe_steps = -1;
    unsigned int threads = 10;
    std::string output_csv;
};

static bool approx_equal(double a, double b, double tol = 1e-12) {
    return std::abs(a - b) <= tol;
}

static std::vector<double> strike_grid(bool short_maturity) {
    if (short_maturity) {
        return {90.0, 95.0, 100.0, 105.0, 110.0};
    }
    return {
        60.6530659713, 67.0320046036, 74.0818220682, 81.8730753078, 90.4837418036,
        100.0000000000, 110.5170918076, 122.1402758160, 134.9858807576, 149.1824697641,
        164.8721270700
    };
}

static std::vector<double> reference_prices(double alpha, bool short_maturity, bool paper_smile) {
    if (paper_smile) {
        if (short_maturity) throw std::runtime_error("Short-maturity paper-smile reference prices are unavailable.");
        if (approx_equal(alpha, 0.55)) return {39.6118131374, 33.4326468122, 26.7468514170, 19.6393012672, 12.3720695911, 5.6786327138, 1.3348427066, 0.1712953054, 0.0217526183, 0.0030376330, 0.0004540409};
        if (approx_equal(alpha, 0.60)) return {39.6021398314, 33.4207792019, 26.7337105281, 19.6274100675, 12.3661092228, 5.6832200884, 1.3439021790, 0.1744595624, 0.0223039828, 0.0031104505, 0.0004610929};
        if (approx_equal(alpha, 0.62)) return {39.5981667096, 33.4158884404, 26.7282928685, 19.6225261818, 12.3636842230, 5.6850558647, 1.3475771613, 0.1758068976, 0.0225414826, 0.0031412192, 0.0004638614};
        if (approx_equal(alpha, 0.80)) return {39.5597910542, 33.3678581860, 26.6747347368, 19.5746914738, 12.3408308040, 5.7019266801, 1.3817403896, 0.1899284108, 0.0250560253, 0.0034379908, 0.0004821865};
        if (approx_equal(alpha, 0.95)) return {39.5248765681, 33.3220839272, 26.6221568193, 19.5278349956, 12.3203146453, 5.7175545811, 1.4113668327, 0.2041816734, 0.0274990630, 0.0036601771, 0.0004794432};
        throw std::runtime_error("Unsupported alpha/maturity for paper-smile reference prices.");
    }
    if (short_maturity) {
        if (approx_equal(alpha, 0.55)) return {10.1267556887, 5.6989890033, 2.4153669912, 0.6906980177, 0.1210835182};
        if (approx_equal(alpha, 0.62)) return {10.1143601668, 5.6723586879, 2.3896783413, 0.6809088652, 0.1205286640};
        if (approx_equal(alpha, 0.80)) return {10.0934043360, 5.6241329401, 2.3423581424, 0.6629959068, 0.1195930284};
        if (approx_equal(alpha, 0.95)) return {10.0836030770, 5.5998574438, 2.3181367808, 0.6539010658, 0.1191571922};
    } else {
        if (approx_equal(alpha, 0.55)) return {39.5327363218, 33.4585987467, 27.1008888388, 20.7264799125, 14.7263466552, 9.5419413959, 5.5249206452, 2.7978822758, 1.2120586020, 0.4391878311, 0.1301411087};
        if (approx_equal(alpha, 0.62)) return {39.5219924356, 33.4378252870, 27.0657610434, 20.6750634817, 14.6619331666, 9.4737132680, 5.4646233842, 2.7540552515, 1.1862614421, 0.4270938587, 0.1257047612};
        if (approx_equal(alpha, 0.80)) return {39.4960053297, 33.3861946261, 26.9766025591, 20.5425822333, 14.4943784717, 9.2954508931, 5.3070991752, 2.6400102593, 1.1196020243, 0.3961290068, 0.1144634657};
        if (approx_equal(alpha, 0.95)) return {39.4766060391, 33.3461467839, 26.9053773904, 20.4344884107, 14.3558250751, 9.1471029750, 5.1759888405, 2.5455704885, 1.0649100837, 0.3710309023, 0.1054763814};
    }
    throw std::runtime_error("Unsupported alpha for precomputed reference prices.");
}

static void fill_default_steps(OverlayConfig& cfg) {
    if (cfg.paper_smile) {
        if (cfg.integrated_steps < 0) cfg.integrated_steps = 160;
        if (cfg.inar_steps < 0) cfg.inar_steps = 160;
        if (cfg.ivi_steps < 0) cfg.ivi_steps = 160;
        if (cfg.hqe_steps < 0) cfg.hqe_steps = 160;
        return;
    }
    if (cfg.short_maturity) {
        if (cfg.integrated_steps < 0) {
            if (approx_equal(cfg.alpha, 0.55)) cfg.integrated_steps = 40;
            else if (approx_equal(cfg.alpha, 0.62)) cfg.integrated_steps = 40;
            else if (approx_equal(cfg.alpha, 0.80)) cfg.integrated_steps = 160;
            else if (approx_equal(cfg.alpha, 0.95)) cfg.integrated_steps = 40;
        }
        if (cfg.inar_steps < 0) {
            if (approx_equal(cfg.alpha, 0.55)) cfg.inar_steps = 320;
            else if (approx_equal(cfg.alpha, 0.62)) cfg.inar_steps = 160;
            else if (approx_equal(cfg.alpha, 0.80)) cfg.inar_steps = 320;
            else if (approx_equal(cfg.alpha, 0.95)) cfg.inar_steps = 320;
        }
        if (cfg.ivi_steps < 0) {
            if (approx_equal(cfg.alpha, 0.55)) cfg.ivi_steps = 320;
            else if (approx_equal(cfg.alpha, 0.62)) cfg.ivi_steps = 160;
            else if (approx_equal(cfg.alpha, 0.80)) cfg.ivi_steps = 80;
            else if (approx_equal(cfg.alpha, 0.95)) cfg.ivi_steps = 160;
        }
        if (cfg.hqe_steps < 0) {
            if (approx_equal(cfg.alpha, 0.55)) cfg.hqe_steps = 320;
            else if (approx_equal(cfg.alpha, 0.62)) cfg.hqe_steps = 80;
            else if (approx_equal(cfg.alpha, 0.80)) cfg.hqe_steps = 80;
            else if (approx_equal(cfg.alpha, 0.95)) cfg.hqe_steps = 160;
        }
    } else {
        if (cfg.integrated_steps < 0) {
            if (approx_equal(cfg.alpha, 0.55)) cfg.integrated_steps = 40;
            else if (approx_equal(cfg.alpha, 0.62)) cfg.integrated_steps = 40;
            else if (approx_equal(cfg.alpha, 0.80)) cfg.integrated_steps = 320;
            else if (approx_equal(cfg.alpha, 0.95)) cfg.integrated_steps = 40;
        }
        if (cfg.inar_steps < 0) {
            if (approx_equal(cfg.alpha, 0.55)) cfg.inar_steps = 320;
            else if (approx_equal(cfg.alpha, 0.62)) cfg.inar_steps = 160;
            else if (approx_equal(cfg.alpha, 0.80)) cfg.inar_steps = 160;
            else if (approx_equal(cfg.alpha, 0.95)) cfg.inar_steps = 320;
        }
        if (cfg.ivi_steps < 0) {
            if (approx_equal(cfg.alpha, 0.55)) cfg.ivi_steps = 40;
            else if (approx_equal(cfg.alpha, 0.62)) cfg.ivi_steps = 160;
            else if (approx_equal(cfg.alpha, 0.80)) cfg.ivi_steps = 320;
            else if (approx_equal(cfg.alpha, 0.95)) cfg.ivi_steps = 160;
        }
        if (cfg.hqe_steps < 0) {
            if (approx_equal(cfg.alpha, 0.55)) cfg.hqe_steps = 160;
            else if (approx_equal(cfg.alpha, 0.62)) cfg.hqe_steps = 160;
            else if (approx_equal(cfg.alpha, 0.80)) cfg.hqe_steps = 80;
            else if (approx_equal(cfg.alpha, 0.95)) cfg.hqe_steps = 80;
        }
    }

    if (cfg.integrated_steps < 0 || cfg.inar_steps < 0 || cfg.ivi_steps < 0 || cfg.hqe_steps < 0) {
        throw std::runtime_error("Unsupported alpha for default step selection.");
    }
}

static OverlayConfig parse_args(int argc, char** argv) {
    OverlayConfig cfg;
    bool alpha_set = false;
    bool integrated_factors_set = false;
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        auto require_value = [&](const std::string& name) -> std::string {
            if (i + 1 >= argc) throw std::runtime_error("Missing value for " + name);
            return argv[++i];
        };

        if (arg == "--alpha") {
            cfg.alpha = std::stod(require_value(arg));
            alpha_set = true;
        }
        else if (arg == "--paths") cfg.paths = std::stoi(require_value(arg));
        else if (arg == "--t-1") cfg.short_maturity = false;
        else if (arg == "--t-1-12" || arg == "--short") cfg.short_maturity = true;
        else if (arg == "--paper-smile") {
            cfg.paper_smile = true;
            cfg.short_maturity = false;
            if (!alpha_set) cfg.alpha = 0.60;
        }
        else if (arg == "--inar-richardson") cfg.inar_richardson = true;
        else if (arg == "--integrated-steps") cfg.integrated_steps = std::stoi(require_value(arg));
        else if (arg == "--integrated-factors") {
            cfg.integrated_factors = std::stoi(require_value(arg));
            integrated_factors_set = true;
        }
        else if (arg == "--inar-steps") cfg.inar_steps = std::stoi(require_value(arg));
        else if (arg == "--ivi-steps") cfg.ivi_steps = std::stoi(require_value(arg));
        else if (arg == "--hqe-steps") cfg.hqe_steps = std::stoi(require_value(arg));
        else if (arg == "--threads") {
            const long parsed = std::stol(require_value(arg));
            if (parsed <= 0 || static_cast<unsigned long>(parsed) > std::numeric_limits<unsigned int>::max()) {
                throw std::runtime_error("--threads must be a positive unsigned integer.");
            }
            cfg.threads = static_cast<unsigned int>(parsed);
        }
        else if (arg == "--output") cfg.output_csv = require_value(arg);
        else throw std::runtime_error("Unknown argument: " + arg);
    }

    if (cfg.paper_smile && !integrated_factors_set) {
        cfg.integrated_factors = 3;
    }

    if (cfg.output_csv.empty()) {
        std::ostringstream oss;
        if (cfg.paper_smile) {
            oss << "smile_overlay_paper_smile_alpha_" << std::fixed << std::setprecision(2) << cfg.alpha
                << "_t_1_paths_" << cfg.paths << ".csv";
        } else {
            oss << "smile_overlay_alpha_" << std::fixed << std::setprecision(2) << cfg.alpha
                << (cfg.short_maturity ? "_t_1_12" : "_t_1")
                << "_paths_" << cfg.paths << ".csv";
        }
        std::string s = oss.str();
        std::replace(s.begin(), s.end(), '.', '_');
        cfg.output_csv = s;
    }

    fill_default_steps(cfg);
    return cfg;
}

static SliceStats richardson_extrapolate_slice(
    const SliceStats& coarse,
    const SliceStats& fine,
    double p,
    const std::vector<double>& strikes,
    const std::vector<double>& ref_ivs,
    double S0,
    double T,
    double r
) {
    const double a = std::pow(2.0, p);
    std::vector<double> prices(strikes.size(), 0.0);
    for (size_t i = 0; i < strikes.size(); ++i) {
        prices[i] = (a * fine.prices[i] - coarse.prices[i]) / (a - 1.0);
        const double intrinsic = std::max(S0 - strikes[i] * std::exp(-r * T), 0.0);
        prices[i] = std::max(prices[i], intrinsic + 1e-12);
        prices[i] = std::min(prices[i], S0 - 1e-12);
    }

    SliceStats out = finalize_slice_stats(
        "INAR-Richardson(p=1)",
        fine.steps,
        -1,
        prices,
        strikes,
        S0,
        T,
        r,
        ref_ivs,
        coarse.seconds + fine.seconds
    );
    out.threads = fine.threads;
    return out;
}

static void write_method_rows(
    std::ofstream& out,
    const SliceStats& st,
    const std::vector<double>& strikes,
    double S0,
    double alpha,
    double T,
    int paths,
    unsigned int seed
) {
    for (size_t i = 0; i < strikes.size(); ++i) {
        out << st.method << ','
            << alpha << ','
            << T << ','
            << paths << ','
            << st.threads << ','
            << seed << ','
            << st.steps << ','
            << st.factors << ','
            << strikes[i] << ','
            << std::log(strikes[i] / S0) << ','
            << st.prices[i] << ','
            << st.ivs[i] << ','
            << st.seconds << ','
            << st.discounted_spot_mean << ','
            << st.spot_se << ','
            << st.martingale_z_score << ','
            << st.terminal_integrated_variance_mean << ','
            << st.terminal_noise_bracket_mean << ','
            << st.integrated_bracket_gap_max << ','
            << st.integrated_driver_mean << ','
            << st.integrated_driver_second_moment << ','
            << st.integrated_driver_bracket_ratio << ','
            << st.integrated_variance_driver_second_moment << ','
            << st.integrated_variance_driver_bracket_ratio << ','
            << st.projection_hit_rate << ','
            << st.projection_gap_mean << ','
            << st.projection_gap_max << '\n';
    }
}

int main(int argc, char** argv) {
    const OverlayConfig cfg = parse_args(argc, argv);
    const unsigned int threads = cfg.threads;

    const double T = cfg.short_maturity ? (1.0 / 12.0) : 1.0;
    const double r = 0.0;
    const double gamma_val = cfg.paper_smile ? 0.3 : 0.1;
    const double theta = cfg.paper_smile ? (0.02 / 0.3) : 0.3156;
    const double nu = cfg.paper_smile ? 1.0 : 0.331;
    const double rho = cfg.paper_smile ? -0.7 : -0.681;
    const double V0 = cfg.paper_smile ? 0.02 : 0.0392;
    const double S0 = 100.0;
    const double x_min = 1e-3;
    const double x_max = 1e3;

    const std::vector<double> strikes = strike_grid(cfg.short_maturity);
    const std::vector<double> ref_prices = reference_prices(cfg.alpha, cfg.short_maturity, cfg.paper_smile);
    std::vector<double> ref_ivs(strikes.size(), 0.0);
    for (size_t i = 0; i < strikes.size(); ++i) {
        ref_ivs[i] = implied_vol_call(ref_prices[i], S0, strikes[i], T, r);
    }

    RHParams rh{cfg.alpha - 0.5, gamma_val, theta, nu, rho, V0, S0, T, r};
    const double beta_param = solve_beta(rho);
    const double mu = solve_mu(nu, theta, beta_param, gamma_val);
    const double xi = solve_xi(V0, theta);
    IVIParams ivi_params{cfg.alpha, gamma_val, theta, nu, rho, V0, S0, T, r};
    HQEParams hqe_params{cfg.alpha, gamma_val, theta, nu, rho, V0, S0, T, r};
    KernelApprox ka = build_loggrid_kernel_approx(rh.H, rh.lambda, cfg.integrated_factors, x_min, x_max);

    std::cout << "Smile overlay benchmark\n\n";
    std::cout << "alpha   = " << cfg.alpha << "\n";
    std::cout << "T       = " << T << "\n";
    std::cout << "paths   = " << cfg.paths << "\n";
    std::cout << "threads = " << threads << "\n";
    std::cout << "paper_smile = " << (cfg.paper_smile ? "true" : "false") << "\n";
    std::cout << "inar_richardson = " << (cfg.inar_richardson ? "true" : "false") << "\n";
    std::cout << "gamma   = " << gamma_val << "\n";
    std::cout << "theta   = " << theta << "\n";
    std::cout << "nu      = " << nu << "\n";
    std::cout << "rho     = " << rho << "\n";
    std::cout << "V0      = " << V0 << "\n";
    std::cout << "strikes =";
    for (double k : strikes) std::cout << " " << k;
    std::cout << "\n";
    std::cout << "integrated_steps   = " << cfg.integrated_steps << " (factors=" << cfg.integrated_factors << ")\n";
    std::cout << "inar_steps         = " << cfg.inar_steps << "\n";
    std::cout << "ivi_steps          = " << cfg.ivi_steps << "\n";
    std::cout << "hqe_steps          = " << cfg.hqe_steps << "\n";
    std::cout << "output_csv         = " << cfg.output_csv << "\n\n";

    const SliceStats st_integrated = run_integrated_slice_mc(
        cfg.paths, cfg.integrated_steps, rh, ka, strikes, ref_ivs,
        cfg.integrated_factors, threads, 111111u
    );
    const SliceStats st_inar = run_inar_slice_mc(
        cfg.paths, cfg.inar_steps, strikes, cfg.alpha, gamma_val, beta_param,
        mu, theta, xi, S0, T, r, ref_ivs, threads, 222222u
    );
    SliceStats st_inar_richardson;
    bool have_inar_richardson = false;
    if (cfg.inar_richardson) {
        if (cfg.inar_steps <= 1 || (cfg.inar_steps % 2) != 0) {
            throw std::runtime_error("INAR Richardson requires an even inar_steps >= 2.");
        }
        const SliceStats st_inar_coarse = run_inar_slice_mc(
            cfg.paths, cfg.inar_steps / 2, strikes, cfg.alpha, gamma_val,
            beta_param, mu, theta, xi, S0, T, r, ref_ivs, threads, 222222u
        );
        st_inar_richardson = richardson_extrapolate_slice(
            st_inar_coarse, st_inar, 1.0, strikes, ref_ivs, S0, T, r
        );
        have_inar_richardson = true;
    }
    const SliceStats st_ivi = run_ivi_slice_mc(
        cfg.paths, cfg.ivi_steps, ivi_params, strikes, ref_ivs, threads, 333333u
    );
    const SliceStats st_hqe = run_hqe_slice_mc(
        cfg.paths, cfg.hqe_steps, hqe_params, strikes, ref_ivs, threads, 444444u
    );

    const std::filesystem::path out_path(cfg.output_csv);
    if (!out_path.parent_path().empty()) {
        std::filesystem::create_directories(out_path.parent_path());
    }
    std::ofstream out(cfg.output_csv);
    if (!out) {
        throw std::runtime_error("Failed to open output CSV.");
    }
    out << "method,alpha,T,paths,threads,seed,steps,factors,strike,log_moneyness,price,iv,time_seconds,"
        << "discounted_spot_mean,spot_se,martingale_z_score,"
        << "terminal_integrated_variance_mean,terminal_noise_bracket_mean,"
        << "integrated_bracket_gap_max,integrated_driver_mean,"
        << "integrated_driver_second_moment,integrated_driver_bracket_ratio,"
        << "integrated_variance_driver_second_moment,integrated_variance_driver_bracket_ratio,"
        << "projection_hit_rate,projection_gap_mean,projection_gap_max\n";
    for (size_t i = 0; i < strikes.size(); ++i) {
        out << "Reference" << ','
            << cfg.alpha << ','
            << T << ','
            << cfg.paths << ','
            << threads << ','
            << 0 << ','
            << 0 << ','
            << 0 << ','
            << strikes[i] << ','
            << std::log(strikes[i] / S0) << ','
            << ref_prices[i] << ','
            << ref_ivs[i] << ','
            << 0.0;
        for (int diagnostic = 0; diagnostic < 14; ++diagnostic) out << ",nan";
        out << '\n';
    }
    write_method_rows(out, st_integrated, strikes, S0, cfg.alpha, T, cfg.paths, 111111u);
    write_method_rows(out, st_inar, strikes, S0, cfg.alpha, T, cfg.paths, 222222u);
    if (have_inar_richardson) {
        write_method_rows(out, st_inar_richardson, strikes, S0, cfg.alpha, T, cfg.paths, 222222u);
    }
    write_method_rows(out, st_ivi, strikes, S0, cfg.alpha, T, cfg.paths, 333333u);
    write_method_rows(out, st_hqe, strikes, S0, cfg.alpha, T, cfg.paths, 444444u);
    out.close();

    std::cout << std::fixed << std::setprecision(6);
    std::cout << std::left
              << std::setw(28) << "Method"
              << std::setw(10) << "Steps"
              << std::setw(10) << "Factors"
              << std::setw(16) << "MaxIVErr"
              << std::setw(16) << "RMSIVErr"
              << std::setw(12) << "Time(s)"
              << "\n";
    std::cout << std::string(92, '-') << "\n";
    std::vector<SliceStats> methods = {st_integrated, st_inar};
    if (have_inar_richardson) {
        methods.push_back(st_inar_richardson);
    }
    methods.push_back(st_ivi);
    methods.push_back(st_hqe);
    for (const SliceStats& st : methods) {
        std::cout << std::setw(28) << st.method
                  << std::setw(10) << st.steps
                  << std::setw(10) << st.factors
                  << std::setw(16) << st.max_iv_error
                  << std::setw(16) << st.rms_iv_error
                  << std::setw(12) << st.seconds
                  << "\n";
    }
    std::cout << "\nCSV saved to " << cfg.output_csv << "\n";
    return 0;
}
