program antenna_simulation
    implicit none
    integer, parameter :: num_theta_points = 1000
    integer, parameter :: num_samples = 100
    integer :: num_elements
    real(8) :: frequency, d, radius
    real(8), parameter :: c = 299792458.0
    real(8) :: wavelength, k
    real(8) :: delta_phi
    character(len=32) :: arg

    ! 🔹 Vérifier que l'utilisateur a passé 5 arguments
    if (command_argument_count() /= 5) then
        print *, "Utilisation : ./antenna_simulation <num_elements> <frequency> <d> <radius> <delta_phi>"
        stop
    end if

    ! 🔹 Lire les arguments en ligne de commande
    call get_command_argument(1, arg)
    read(arg, *) num_elements
    call get_command_argument(2, arg)
    read(arg, *) frequency
    call get_command_argument(3, arg)
    read(arg, *) d
    call get_command_argument(4, arg)
    read(arg, *) radius
    call get_command_argument(5, arg)
    read(arg, *) delta_phi

    ! 🔹 Calculer les valeurs dérivées
    wavelength = c / frequency
    k = 2.0 * 3.141592653589793d0 / wavelength

    print *, "Nombre d'éléments :", num_elements
    print *, "Fréquence (Hz) :", frequency
    print *, "Espacement (m) :", d
    print *, "Rayon :", radius
    print *, "Déphasage delta_phi :", delta_phi

    ! 🔹 Lancer le calcul et la génération du graphique
    call compute_and_plot(num_elements, frequency, d, radius, delta_phi)

contains

    ! 📌 Calcul de la puissance moyenne pour chaque angle theta
    subroutine compute_and_plot(num_elements, frequency, d, radius, delta_phi)
        implicit none
        integer, intent(in) :: num_elements
        real(8), intent(in) :: frequency, d, radius, delta_phi
        real(8) :: theta_values(num_theta_points), power_values_db(num_theta_points)
        real(8) :: theta, max_power
        integer :: i

        ! 🔹 Générer la liste des angles theta
        do i = 1, num_theta_points
            theta_values(i) = -3.141592653589793d0/2 + (i-1) * (2.0 * (3.141592653589793d0/2) / num_theta_points)
        end do

        ! 🔹 Calculer la puissance moyenne
        max_power = 0.0
        do i = 1, num_theta_points
            power_values_db(i) = avgPower(radius, theta_values(i), num_elements, frequency, d, delta_phi)
            if (power_values_db(i) > max_power) then
                max_power = power_values_db(i)
            end if
        end do

        ! 🔹 Normalisation et conversion en dB
        open(10, file="data_cartesian2.dat", status="replace")
        do i = 1, num_theta_points
            power_values_db(i) = 10.0 * log10(power_values_db(i) / max_power)
            write(10, *) theta_values(i) * 180.0 / 3.141592653589793d0, power_values_db(i)
        end do
        close(10)

        print *, "Fichier 'data_cartesian2.dat' généré."

        ! 🔹 Générer et exécuter le script Gnuplot
        call generate_gnuplot_script("data_cartesian2.dat")
        call execute_command_line("gnuplot cartesian_plot2.gnu")

    end subroutine compute_and_plot

    ! 📌 Moyenne de la puissance sur plusieurs échantillons temporels
    function avgPower(r, theta, num_elements, frequency, d, delta_phi) result(avg_pwr)
        implicit none
        real(8), intent(in) :: r, theta, frequency, d, delta_phi
        integer, intent(in) :: num_elements
        real(8) :: avg_pwr, sum_pwr, time_step
        integer :: i

        sum_pwr = 0.0d0
        time_step = 1.0 / (frequency * num_samples)

        do i = 1, num_samples
            sum_pwr = sum_pwr + power(r, theta, i * time_step, num_elements, frequency, d, delta_phi)
        end do

        avg_pwr = sum_pwr / num_samples
    end function avgPower

    ! 📌 Fonction de puissance instantanée
    function power(r, theta, time, num_elements, frequency, d, delta_phi) result(pwr)
        implicit none
        real(8), intent(in) :: r, theta, time, frequency, d, delta_phi
        integer, intent(in) :: num_elements
        real(8) :: pwr, r_n, phase
        complex(8) :: sum_complex
        integer :: n

        sum_complex = (0.0d0, 0.0d0)
        do n = 1, num_elements
            r_n = sqrt(r**2 + ((n - (num_elements - 1) / 2.0) * d - r*sin(theta))**2)
            phase = k * r_n
            sum_complex = sum_complex + (1.0 / r_n) * exp(- (0.0d0, 1.0d0) * phase) &
                * exp((0.0d0, 1.0d0) * (2.0 * 3.141592653589793d0 * frequency * time)) &
                * exp((0.0d0, 1.0d0) * delta_phi * (n-1))
        end do

        pwr = abs(sum_complex)**2
    end function power

    ! 📌 Génère le script Gnuplot
    subroutine generate_gnuplot_script(filename)
        implicit none
        character(len=*), intent(in) :: filename
        integer :: file_id

        open(20, file="cartesian_plot2.gnu", status="replace")
        write(20, '(A)') "set terminal pngcairo enhanced size 800,600"
        write(20, '(A)') "set output 'cartesian_plot2.png'"
        write(20, '(A)') "set xlabel 'Angle (°)'"
        write(20, '(A)') "set ylabel 'Gain (dB)'"
        write(20, '(A)') "set grid"
        write(20, '(A)') "set title 'Radiation pattern in dB'"
        write(20, '(A)') "set xrange [-100:100]"
        write(20, '(A)') "set yrange [-85:0]"
        write(20, '(A, A, A)') "plot '", trim(filename), "' using 1:2 with lines lw 2 lc rgb 'blue' title ' '"
        close(20)

        print *, "Script Gnuplot 'cartesian_plot2.gnu' généré."
    end subroutine generate_gnuplot_script

end program antenna_simulation