#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <complex>
#include <cstdlib>
#include <iomanip>
#include <chrono>

constexpr double c = 299792458;

class Point {
public:
    double r, theta, x, y;
    Point(double radius, double angle) : r(radius), theta(angle) {
        x = r * std::cos(theta);
        y = r * std::sin(theta);
    }
};

class ElectricField {
    double d, frequency, k, delta_phi;
    int n;

public:
    ElectricField(double spacing, int num, double freq, double phase_shift)
        : d(spacing), n(num), frequency(freq), delta_phi(phase_shift) {
        k = 2 * M_PI * frequency / c;
    }

    double farFieldPower(double theta) {
        double kd_sin_theta = (k * d / 2.0) * std::sin(theta);
        double phase_shift = delta_phi / 2.0;

        double numerator = std::sin(n * (kd_sin_theta + phase_shift));
        double denominator = n * std::sin(kd_sin_theta + phase_shift);

        if (std::abs(denominator) < 1e-10) {
            denominator = 1e-10;
        }

        double array_factor = std::pow(numerator / denominator, 2);
        return array_factor;
    }
};

void generate_gnuplot_script(int num_curves) {
    std::ofstream gnuplotFile("plot.gnu");
    gnuplotFile << "set terminal pngcairo enhanced size 800,600\n";
    gnuplotFile << "set output \"radiation_pattern.png\"\n";
    gnuplotFile << "set xlabel \"Theta (radians)\"\n";
    gnuplotFile << "set ylabel \"Power (dB)\"\n";
    gnuplotFile << "set title \"Radiation pattern in dB\"\n";
    gnuplotFile << "set grid\n";

    gnuplotFile << "plot ";
    for (int i = 2; i <= num_curves + 1; i++) {
        gnuplotFile << "\"radiation_pattern.dat\" using 1:" << i << " with lines lw 2 title \"Curve " << (i - 1) << "\"";
        if (i != num_curves + 1) gnuplotFile << ", ";
    }
    gnuplotFile << "\n";
    gnuplotFile.close();
}

double phase_from_angle(double theta_rad, double d, double frequency) {
    double wavelength = c / frequency;
    double sin_theta = std::sin(theta_rad);

    if (std::abs(sin_theta) > 1) {
        std::cout << "Warning: Angle invalide, ajustez theta." << std::endl;
        return NAN;
    }
    return -((2 * M_PI * d) / wavelength) * sin_theta;
}

int main(int argc, char* argv[]) {
    // Vérification des arguments
    if (argc < 4) {
        std::cerr << "Usage: " << argv[0] << " <N> <L> <frequency>\n";
        return 1;
    }

    // Récupération des paramètres depuis la ligne de commande
    int N = std::stoi(argv[1]);
    double L = std::stod(argv[2]);
    double frequency = std::stod(argv[3]);
    double radius = 2000000;
   
    double d = L / N;
    double angle = M_PI/4;
 
    // std::cout << "d = " << d << ", L = " << L << ", frequency = " << frequency << std::endl;

    std::ofstream dataFile("radiation_pattern.dat");

    auto start_time = std::chrono::steady_clock::now();

    std::vector<double> theta_values;
    std::vector<std::vector<double>> power_data;

    int num_theta = 1000;
    for (int i = 0; i < num_theta; i++) {
        theta_values.push_back(-M_PI / 4 + i * (M_PI / (2 * (num_theta - 1))));
    }

    // double Interval_phase = phase_from_angle(angle, d, frequency);

     double numerator = 2*(2 * M_PI * frequency / c) * (d) * (std::sqrt(2) / 2);
    double denominator = 2*(std::sqrt((1 - (1 / std::sqrt(2))) * (24.0 / (N * N))));  // Diviser en `double`
    std::cout << "numerator : " << numerator << std::endl;
    std::cout << "denominator : " << denominator << std::endl;
    double nbs_aims_double = numerator / denominator;
    int nbs_aims = std::round(nbs_aims_double); 
   double aim = -std::abs(numerator);

    // while (aim < std::abs(numerator)){
    //       std::cout << "aim : " << aim << std::endl;
    //     ElectricField eF(d, N, frequency, aim);
    //             aim =aim+ denominator;
    //     std::vector<double> power_curve;
    //     double max_power = 0;
    //     for (double theta : theta_values) {
    //         double power = eF.farFieldPower(theta);
    //         if (power > max_power) {
    //             max_power = power;
    //         }
    //         power_curve.push_back(power);
    //     }

    //     for (double& p : power_curve) {
    //         p = 10 * std::log10(p / max_power);
    //     }

    //     power_data.push_back(power_curve);
    // }

    for (double i = -std::abs(numerator)/2; i <= numerator/2; i +=denominator) {
        std::cout << "aim : " << i << std::endl;
        // double delta_phi = -2 * M_PI * sin(aim) * d / (c / frequency);
        ElectricField eF(d, N, frequency, i);

        std::vector<double> power_curve;
        double max_power = 0;
        for (double theta : theta_values) {
            double power = eF.farFieldPower(theta);
            if (power > max_power) {
                max_power = power;
            }
            power_curve.push_back(power);
        }

        for (double& p : power_curve) {
            p = 10 * std::log10(p / max_power);
        }

        power_data.push_back(power_curve);
    }

    for (size_t i = 0; i < theta_values.size(); i++) {
        dataFile << theta_values[i];
        for (size_t j = 0; j < power_data.size(); j++) {
            dataFile << " " << power_data[j][i];
        }
        dataFile << "\n";
    }

    

    dataFile.close();
    std::cout << "NBS d'aims: " << nbs_aims << "\n";
    std::cout << "Calcul terminé en "
              << std::chrono::duration_cast<std::chrono::seconds>(std::chrono::steady_clock::now() - start_time).count()
              << " secondes !\n";

    generate_gnuplot_script(nbs_aims);

    std::system("gnuplot plot.gnu");

    return nbs_aims;
}