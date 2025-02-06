#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <complex>
#include <cstdlib>  // Pour system()

constexpr double c = 299792458;  // Vitesse de la lumière (m/s)

// ======================== Classes ========================
class Point {
public:
    double r, theta, x, y;
    Point(double radius, double angle) : r(radius), theta(angle) {
        x = r * std::cos(theta);
        y = r * std::sin(theta);
    }
};

class Emitter {
public:
    double y;
    explicit Emitter(double y_pos) : y(y_pos) {}

    double distance(const Point& p) const {
        return std::sqrt(p.x * p.x + (p.y - y) * (p.y - y));
    }
};

class Antenna {
protected:
    double d, frequency, k;
    int n;
    std::vector<Emitter> emitters;

public:
    Antenna(double spacing, int num, double freq) : d(spacing), n(num), frequency(freq) {
        k = 2 * M_PI * frequency / c;
        for (int i = 0; i < n; i++) {
            emitters.emplace_back((i - (n - 1) / 2.0) * d);
        }
    }
};

// ======================== Fenêtres de Pondération ========================
double blackmanWeight(int i, int N) {
    return 0.42 - 0.5 * std::cos((2 * M_PI * i) / (N - 1)) +
           0.08 * std::cos((4 * M_PI * i) / (N - 1));
}

double chebyshevWeight(int i, int N, double ripple_db) {
    if (ripple_db <= 0) return 1.0; // Eviter log ou arccosh de valeurs invalides
    double ripple_linear = pow(10, ripple_db / 20.0);
    double beta = cosh(acosh(ripple_linear) / (N - 1));

    if (beta < 1.0) return 1.0; // Sécurité supplémentaire
    double x = cos(M_PI * i / (N - 1));
    double Tn = cosh((N - 1) * acosh(std::max(-1.0, std::min(1.0, beta * x))));
    double T0 = cosh((N - 1) * acosh(beta)); // Normalisation
    return Tn / T0;
}

double taylorWeight(int i, int N, double sideLobeLevel) {
    double A = acosh(pow(10, std::abs(sideLobeLevel) / 20.0)) / M_PI;
    double n_bar = A * A + 0.5;
    double sum = 0.0;

    for (int m = 1; m <= floor(n_bar - 1); ++m) {
        sum += pow(-1, m) * cos(2 * M_PI * m * i / (N - 1)) / (1 + 2 * m * m / (A * A));
    }
    return 1 + 2 * sum;
}

double computeEmitterWeight(int i, int N, std::string method, double ripple_or_sll) {
    if (method == "blackman") {
        return blackmanWeight(i, N);
    } else if (method == "chebyshev") {
        return chebyshevWeight(i, N, ripple_or_sll);
    } else if (method == "taylor") {
        return taylorWeight(i, N, ripple_or_sll);
    } else {
        return 1.0; // Pas de pondération par défaut
    }
}

// ======================== Calcul du Diagramme ========================
double farFieldPower(double theta, double delta_phi, double frequency, double d, int N,
                     std::string method, double ripple_or_sll) {
    double k = 2 * M_PI * frequency / c;
    double kd_sin_theta = (k * d) * std::sin(theta);
    double phase_shift = delta_phi;

    std::complex<double> sum = 0.0;
    std::complex<double> i(0, 1);

    for (int n = 0; n < N; ++n) {
        double weight = computeEmitterWeight(n, N, method, ripple_or_sll);
        sum += weight * std::exp(i * static_cast<double>(n) * (kd_sin_theta + phase_shift));
    }

    return std::norm(sum); // Norme au carré pour obtenir la puissance
}

// ======================== Génération Gnuplot ========================
void generate_gnuplot_script(const std::string& filename) {
    std::ofstream gnuplotFile("cartesian_plot.gnu");
    gnuplotFile << "set terminal pngcairo enhanced size 800,600\n";
    gnuplotFile << "set output \"cartesian_plot.png\"\n";
    gnuplotFile << "set xlabel \"Angle (°)\"\n";
    gnuplotFile << "set ylabel \"Gain (dB)\"\n";
    gnuplotFile << "set grid\n";
    gnuplotFile << "set title \"Radiation pattern in dB\"\n";
    gnuplotFile << "set xrange [-100:100]\n";
    gnuplotFile << "set yrange [-85:0]\n";
    gnuplotFile << "plot \"" << filename << "\" using 1:2 with lines lw 2 lc rgb \"blue\" title \" \"\n";
    gnuplotFile.close();
}

int main(int argc, char* argv[]) {
    if (argc < 7) {
        std::cerr << "Usage: " << argv[0] << " <N> <L> <frequency> <delta_phi> <method> <ripple_or_sll>\n";
        std::cerr << "  <method> : 'blackman', 'chebyshev' ou 'taylor'\n";
        return 1;
    }

    int nbs = std::stoi(argv[1]);
    double L = std::stod(argv[2]);
    double frequency = std::stod(argv[3]);
    double delta_phi = std::stod(argv[4]);
    std::string method = argv[5];
    double ripple_or_sll = std::stod(argv[6]);

    double d = L / nbs;
    std::string data_filename = "data_cartesian.dat";
    std::ofstream dataFile(data_filename);

    if (!dataFile) {
        std::cerr << "Erreur : Impossible de créer le fichier '" << data_filename << "'\n";
        return 1;
    }

    std::vector<double> theta_values, power_values_db;
    double max_power = 0;
    std::vector<double> raw_powers;

    for (int i = 0; i < 1000; i++) {
        double theta = -M_PI / 2 + i * (M_PI / 999);
        double power = farFieldPower(theta, delta_phi, frequency, d, nbs, method, ripple_or_sll);
        raw_powers.push_back(power);
        if (power > max_power) {
            max_power = power;
        }
    }

    for (int i = 0; i < 1000; i++) {
        double theta = -M_PI / 2 + i * (M_PI / 999);
        double normalized_power = raw_powers[i] / max_power;
        double power_db = 10 * std::log10(normalized_power);
        power_values_db.push_back(power_db);
        theta_values.push_back(theta * 180.0 / M_PI);
        dataFile << theta_values.back() << " " << power_values_db.back() << "\n";
    }

    dataFile.close();
    std::cout << "Fichier '" << data_filename << "' généré.\n";

    generate_gnuplot_script(data_filename);
    std::cout << "Génération du graphique...\n";
    system("gnuplot cartesian_plot.gnu");

    return 0;
}