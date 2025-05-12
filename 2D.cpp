
#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <string>

const int Nx = 100;
const int Ny = 100;
const int Nt = 300;

// paramètres de grilles
const double Lx = 1.0;
const double Ly = 1.0;
const double dx = Lx / (Nx - 1);
const double dy = Ly / (Ny - 1);
const double g = 9.81;
const double h0 = 1.0;
const double c = std::sqrt(g * h0);
const double dt = 0.4 * std::min(dx, dy) / c; 


int idx(int i, int j) { return i * Ny + j; }

int main() {
    std::vector<double> f_past(Nx * Ny, 0.0);     // t - dt
    std::vector<double> f_now(Nx * Ny, 0.0);    // t
    std::vector<double> f_futur(Nx * Ny, 0.0);   // t + dt

    // Initialisation : bosse gaussienne au centre
    for (int i = 0; i < Nx; ++i) {
        for (int j = 0; j < Ny; ++j) {
            double x = i * dx - Lx / 2;
            double y = j * dy - Ly / 2;
            f_now[idx(i,j)] = std::exp(-100 * (x*x + y*y));
        }
    }

    //
    f_past = f_now;

    // Boucle temporelle
    for (int n = 0; n <= Nt; ++n) {
        for (int i = 1; i < Nx - 1; ++i) {
            for (int j = 1; j < Ny - 1; ++j) {
                double laplacian =
                    (f_now[idx(i+1,j)] + f_now[idx(i-1,j)] - 2 * f_now[idx(i,j)]) / (dx*dx) +
                    (f_now[idx(i,j+1)] + f_now[idx(i,j-1)] - 2 * f_now[idx(i,j)]) / (dy*dy);
                f_futur[idx(i,j)] = 2 * f_now[idx(i,j)] - f_past[idx(i,j)] + c*c * dt*dt * laplacian;
            }
        }

        // Conditions aux bords : 
        for (int i = 0; i < Nx; ++i) {
            f_futur[idx(i,0)] = 0;
            f_futur[idx(i,Ny-1)] = 0;
        }
        for (int j = 0; j < Ny; ++j) {
            f_futur[idx(0,j)] = 0;
            f_futur[idx(Nx-1,j)] = 0;
        }


	// définit l'éccccccriture tous les certains pas de temps
        if (n % 10 == 0) {  // determine la fréquence des images
    std::ofstream out("wave_" + std::to_string(n) + ".dat");
    for (int i = 0; i < Nx; ++i) {
        for (int j = 0; j < Ny; ++j) {
            out << i * dx << " " << j * dy << " " << f_now[idx(i,j)] << "\n";
        }
    }
    out.close();
}

        // Mise à jour temporelle
        f_past.swap(f_now);
        f_now.swap(f_futur);
    }

    return 0;
}

