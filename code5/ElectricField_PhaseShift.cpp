#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <complex>
#include <cstdlib>  // Pour system()

constexpr double c = 299792458;  // Vitesse de la lumière (m/s)

// Classe représentant un point dans l'espace (coordonnées polaires et cartésiennes)
class Point {
public:
    double r, theta, x, y;
    Point(double radius, double angle) : r(radius), theta(angle) {
        x = r * std::cos(theta);
        y = r * std::sin(theta);
    }
};

// Classe représentant un émetteur positionné sur l'axe y
class Emitter {
public:
    double y;
    explicit Emitter(double y_pos) : y(y_pos) {}

    // Calcul de la distance entre l'émetteur et un point donné
    double distance(const Point& p) const {
        return std::sqrt(p.x * p.x + (p.y - y) * (p.y - y));
    }
};

// Classe représentant une antenne linéaire avec plusieurs émetteurs
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

// Classe représentant le champ électrique émis par l'antenne
class ElectricField : public Antenna {
    double delta_phi;

public:
    ElectricField(double spacing, int num, double freq, double phase_shift)
        : Antenna(spacing, num, freq), delta_phi(phase_shift) {}

    // Calcul de la puissance moyenne reçue en un point
    double avgPower(const Point& p) const {
        double sum = 0;
        int samples = 100;
        for (int i = 0; i < samples; i++) {
            double t = i / (double)samples / frequency;
            sum += power(p, t);
        }
        return sum / samples;
    }

    // Calcul de la puissance instantanée
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
};

// Fonction générant un script Gnuplot pour tracer le diagramme de rayonnement
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

// Fonction de simulation et de génération des graphiques
int main() {
    double radius = 1000;
    double frequency = 1.2e9; // 1.2 GHz
    double delta_phi = M_PI / 4;
    int nbs = 40;
    ElectricField eF(4. / nbs, nbs, frequency, delta_phi);

    std::string data_filename = "data_cartesian.dat";
    std::ofstream dataFile(data_filename);

    std::vector<double> theta_values, power_values_db;
    double max_power = 0;
    std::vector<double> raw_powers;

    // Calcul de la puissance maximale pour normalisation
    for (int i = 0; i < 1000; i++) {
        double theta = -M_PI / 2 + i * (M_PI / 999);
        Point p(radius, theta);
        double power = eF.avgPower(p);
        raw_powers.push_back(power);
        if (power > max_power) {
            max_power = power;
        }
    }

    // Normalisation et conversion en dB
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

    // Générer le script Gnuplot
    generate_gnuplot_script(data_filename);
    std::cout << "Script Gnuplot 'cartesian_plot.gnu' généré.\n";

    // Lancer Gnuplot automatiquement
    std::cout << "Génération du graphique...\n";
    system("gnuplot cartesian_plot.gnu");

    // Ouvrir l'image générée
    std::cout << "Affichage du graphique...\n";
#ifdef __APPLE__
    system("open cartesian_plot.png");
#elif __linux__
    system("xdg-open cartesian_plot.png");
#elif _WIN32
    system("start cartesian_plot.png");
#endif
}
