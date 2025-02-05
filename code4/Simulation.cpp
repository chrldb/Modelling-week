#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <complex>
#include <cstdlib>  // Pour system()


constexpr double c = 299792458;  // Vitesse de la lumière (m/s)

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

class ElectricField : public Antenna {
    double delta_phi;

public:
    ElectricField(double spacing, int num, double freq, double phase_shift)
        : Antenna(spacing, num, freq), delta_phi(phase_shift) {}

    double shift_angle() const {
        double wavelength = c / frequency;
        double sin_theta = -(wavelength / d) * (delta_phi / (2 * M_PI));
        if (std::abs(sin_theta) > 1) return NAN;
        return std::asin(sin_theta) * 180.0 / M_PI;
    }

    double power(const Point& p, double time) const {
        std::complex<double> sum_complex = 0;
        double aire = 0;
        int count_phase = 0;

        for (const auto& emitter : emitters) {
            double r_n = emitter.distance(p);
            double phase = k * r_n;
            sum_complex += (1.0 / r_n) * std::exp(std::complex<double>(0, -phase))
                           * std::exp(std::complex<double>(0, 2 * M_PI * frequency * time))
                           * std::exp(std::complex<double>(0, count_phase * delta_phi));
            aire += 1.0 / (4 * M_PI * r_n * r_n);
            count_phase++;
        }

        return std::norm(sum_complex) * aire;
    }

    double avgPower(const Point& p) const {
        double sum = 0;
        int samples = 100;
        for (int i = 0; i < samples; i++) {
            double t = i / (double)samples / frequency;
            sum += power(p, t);
        }
        return sum / samples;
    }

};


double farFieldPower(double theta,double delta_phi, double frequency, double d) {
    double k = 2 * M_PI * frequency / c;
    int n = 40;
    double kd_sin_theta = (k * d / 2.0) * std::sin(theta);
    double phase_shift = delta_phi / 2.0;

    double numerator = std::sin(n * (kd_sin_theta + phase_shift));
    double denominator = n * std::sin(kd_sin_theta + phase_shift);

    // Éviter la division par zéro
    if (std::abs(denominator) < 1e-10) {
        denominator = 1e-10;
    }

    double array_factor = std::pow(numerator / denominator, 2);
    return array_factor;
}

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
    double radius = 1000;

      if (argc < 5) {
        std::cerr << "Usage: " << argv[0] << " <N> <d> <frequency> <Delta_Phi>\n";
        return 1;
    }

    // Récupération des paramètres depuis la ligne de commande
    int nbs = std::stoi(argv[1]);
    double d = std::stod(argv[2]);
    double frequency = std::stod(argv[3]);
    double delta_phi = std::stod(argv[4]);
    
    ElectricField eF(4./nbs, nbs, frequency, delta_phi);

    std::string data_filename = "data_cartesian.dat";
    std::ofstream dataFile(data_filename);

    std::vector<double> theta_values, power_values_db;

    // Calcul de la puissance maximale pour normalisation
    double max_power = 0;
    std::vector<double> raw_powers;
    for (int i = 0; i < 1000; i++) {
        double theta = -M_PI/2 + i * ( M_PI / 999);
        // Point p(radius, theta);
         //double power = eF.avgPower(p);
         double power = farFieldPower(theta, delta_phi, frequency, d);
        raw_powers.push_back(power);
        if (power > max_power) {
            max_power = power;
        }
    }

    // Normalisation et conversion en dB
    for (int i = 0; i < 1000; i++) {
        double theta = -M_PI/2 + i * ( M_PI / 999);
        double normalized_power = raw_powers[i] / max_power;  // Normalisation
        double power_db = 10 * std::log10(normalized_power);  // dB avec max = 0
        power_values_db.push_back(power_db);
        theta_values.push_back(theta * 180.0 / M_PI);
        dataFile << theta_values.back() << " " << power_values_db.back() << "\n";
    }

    dataFile.close();
    std::cout << "Fichier '" << data_filename << "' généré.\n";

    // Générer le script Gnuplot
    generate_gnuplot_script(data_filename);
    std::cout << "Script Gnuplot 'cartesian_plot.gnu' généré.\n";

    // Lancer Gnuplot automatiquement
    std::cout << "Génération du graphique...\n";
    system("gnuplot cartesian_plot.gnu");


}