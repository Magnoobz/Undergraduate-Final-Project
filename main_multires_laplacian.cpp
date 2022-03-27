#define _USE_MATH_DEFINES
#include <cmath>
#include <vector>
#include <fstream>
#include <iostream>
#include <chrono>

#include "src/LSMPS/LSMPS.hpp"
#include "src/Thermal Flux Method/Heat_Flux.hpp"
#include "src/Thermal Flux Method/Surface_Vector_S.hpp"
#include "src/Thermal Flux Method/Calculate_Laplacian.hpp"
#include "src/Time Integration/time_integration.hpp"
# include "src/Neighbor Search/Brute_Force.hpp"

# include "src/Particle/Initialize_Particle.hpp"
# include "src/Particle/Particle_Splitting.hpp"
# include "src/Particle/Multires_Movement.hpp"
# include "src/Particle/Packing_Ratio.hpp"

using namespace std;

int main()
{    
    double x_left   = 0;
    double x_right  = 1;
    double y_bottom = 0;
    double y_top    = 1;

    int nx = 60;
    int ny = 60;

    double dx = (x_right-x_left)/(nx);
    double dy = (y_top-y_bottom)/(ny);

    vector<double> x, y, h, k;
    vector<int> split_index;
    vector<int> not_moving;
    vector<int> is_dummy;

    double n_dummy = 5;

    initialize_particle(x_left, x_right, y_bottom, y_top, nx, ny, n_dummy, not_moving, x, y, h);

    int num_particle = x.size();

    for (int i = 0; i < num_particle; i++)
    {
        if (x[i] < 0.8  && x[i] > 0.2 && y[i] > 0.2 && y[i] < 0.8)
        {
            split_index.push_back(i);
        }
    }

    // for (int i = 0; i < num_particle; i++)
    // {
    //     if (pow(x[i]-0.5,2) + pow (y[i]-0.5,2) < pow(0.3,2))
    //     {
    //         split_index.push_back(i);
    //     }
    // }

    auto splitting_start = chrono::high_resolution_clock::now();

    Particle_Splitting(split_index, not_moving, x, y, h);

    num_particle = x.size();

    vector<double> Laplacian_T_analytic(num_particle), T(num_particle), Laplacian_T_LSMPS(num_particle), T_Lap(num_particle);

    vector<vector<int>> neighbor;
    vector<vector<double>> weight_data;

    double R_e = 1.6;
    vector<vector<double>> Ri_a;
    calc_Ri_a(h, R_e, 7, Ri_a);   

    int loop_count = 0;
    int iter = 100;

    while (loop_count < iter)
    {
        // brute_force(x, y, h, 1, neighbor, weight_data, R_e);
        brute_force_2(x, y, h, 1, neighbor, weight_data, R_e);       

        loop_count++;

        vector<vector<double>> Ni_a;
        calc_Ni(x, y, h, neighbor, weight_data, R_e, Ri_a, Ni_a);

        vector<double> ci;
        calc_ci(x, Ri_a, Ni_a, ci);

        vector<double> delta_x, delta_y;
        calc_DeltaX(x, y, h, R_e, 0.1, neighbor, weight_data, ci, delta_x, delta_y);

        for (int i = 0; i < num_particle; i++)
        {
            if (not_moving[i] == 1)
            {
                continue;
            }

            x[i] = x[i] + delta_x[i];
            y[i] = y[i] + delta_y[i];
        }
    }

    for (int i = 0; i < num_particle; i++)
    {
        k.push_back(200);
        if ((x[i] < x_left) || (x[i] > x_right) || (y[i] < y_bottom) || (y[i] > y_top))
        {
            is_dummy.push_back(0);
        }
        else
        {
            is_dummy.push_back(0);
        }
    }

    for (int i = 0; i < num_particle; i++)
    {       
        if (is_dummy[i] == 1)
        {
            continue;
        }

        // Test Function 1
        T[i] = -1/(8*pow(M_PI,2))*sin(2*M_PI*x[i])*sin(2*M_PI*y[i]);
        Laplacian_T_analytic[i] = sin(2*M_PI*x[i])*sin(2*M_PI*y[i]);

        // Test Function 2
        // T[i] = 100+50*(cos(M_PI*x[i])+cos(M_PI*y[i]));
        // Laplacian_T_analytic[i] = -50*pow(M_PI,2)*(cos(M_PI*x[i])+cos(M_PI*y[i]));

    }

    
    // brute_force(x,y,h,1,neighbor,weight_data,R_e);
    brute_force_2(x,y,h,1,neighbor,weight_data,R_e);

    vector<vector<vector<double>>> LSMPS_eta(num_particle, vector<vector<double>>(5));
    // calc_LSMPS_eta(LSMPS_eta, x, y, h, h, neighbor, weight_data);
    calc_LSMPS_eta_2(LSMPS_eta, x, y, h, h, neighbor, weight_data);

    vector<vector<double>> Sij_x(num_particle), Sij_y(num_particle);
    calculate_Sij(LSMPS_eta, Sij_x, Sij_y, h, h, neighbor);

    vector<double> Bi_x, Bi_y;
    calculate_Bi(Sij_x, Sij_y, x, is_dummy, neighbor, Bi_x, Bi_y);
    // vector<double> Bi_x(num_particle), Bi_y(num_particle);

    vector<vector<double>> Sij_Star_x(num_particle), Sij_Star_y(num_particle);
    // calculate_Sij_Star(Sij_x, Sij_y, Sij_Star_x, Sij_Star_y, x, is_dummy, neighbor);
    calculate_Sij_Star_2(Sij_x, Sij_y, Sij_Star_x, Sij_Star_y, x, is_dummy, neighbor);

    vector<double> kdeltaT_x, kdeltaT_y;
    calc_kdeltaT(x, y, k, T, is_dummy, neighbor, LSMPS_eta, kdeltaT_x, kdeltaT_y);

    vector<vector<double>> kdeltaT_ij;
    calc_MPS_like_value(x, y, k, T, is_dummy, neighbor, kdeltaT_ij);

    vector<vector<double>> alpha_ij, kdeltaT_ij_x, kdeltaT_ij_y;
    calc_hybrid(x, y, neighbor, alpha_ij, kdeltaT_x, kdeltaT_y, is_dummy, kdeltaT_ij, kdeltaT_ij_x, kdeltaT_ij_y);

    calc_Laplacian(k, h, h, neighbor, kdeltaT_ij_x, kdeltaT_ij_y, Sij_Star_x, Sij_Star_y, Bi_x, Bi_y, kdeltaT_x, kdeltaT_y, is_dummy, Laplacian_T_LSMPS);

    vector<double> Lap_Value;
    calc_Laplacian_From_Eta(LSMPS_eta, x, y, T, is_dummy, neighbor, Lap_Value);

    ofstream output;

    output.open("Hasil.csv");

    output << "x" << "," << "y" << "," << "h" << "," << "Analytic" << "," << "LSMPS_Conserved" << "," << "LSMPS_Laplacian\n";

    for (int i = 0; i < num_particle; i++)
    {
        if (x[i] >= x_left && x[i] <= x_right && y[i] >= y_bottom && y[i] <= y_top)
        {
            output << x[i] << "," 
                << y[i] << "," 
                << h[i] << ","
                << Laplacian_T_analytic[i] << ","
                << Laplacian_T_LSMPS[i] << "," 
                << Lap_Value[i] << "\n";
        }
    }
}