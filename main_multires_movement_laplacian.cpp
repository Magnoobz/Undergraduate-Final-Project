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

#include "src/Particle/Initialize_Particle.hpp"
#include "src/Particle/Particle_Splitting.hpp"
#include "src/Particle/Multires_Movement.hpp"
#include "src/Particle/Packing_Ratio.hpp"

#include "src/Neighbor Search/Brute_Force.hpp"
#include "src/Neighbor Search/Spatial_Hash.hpp"

using namespace std;

int main()
{    
    double x_left   = 0;
    double x_right  = 1;
    double y_bottom = 0;
    double y_top    = 1;

    double eax = 2;
    double eay = 1;

    int nx = 80;
    int ny = 40;

    double dx = (x_right-x_left)/(nx)*eax;
    double dy = (y_top-y_bottom)/(ny)*eay;

    vector<double> x, y, h, k;
    vector<int> split_index;
    vector<int> not_moving;
    vector<int> is_dummy;

    double n_dummy = 12;

    initialize_particle(x_left*eax, x_right*eax, y_bottom*eay, y_top*eay, nx, ny, n_dummy, not_moving, x, y, h);

    int num_particle = x.size();

    for (int i = 0; i < num_particle; i++)
    {
        if (pow(x[i]/eax-0.5,2)+pow(y[i]/eay-0.5,2)<pow(0.35,2))
        {
            split_index.push_back(i);
        }
    }

    // for (int i = 0; i < num_particle; i++)
    // {
    //     if (pow(x[i]-0.25,2) + pow (y[i]-0.25,2) < pow(0.15,2))
    //     {
    //         split_index.push_back(i);
    //     }
    //     else if (pow(x[i]-0.25,2) + pow (y[i]-0.75,2) < pow(0.15,2))
    //     {
    //         split_index.push_back(i);
    //     }
    //     else if (pow(x[i]-0.75,2) + pow (y[i]-0.25,2) < pow(0.15,2))
    //     {
    //         split_index.push_back(i);
    //     }
    //     else if (pow(x[i]-0.75,2) + pow (y[i]-0.75,2) < pow(0.15,2))
    //     {
    //         split_index.push_back(i);
    //     }
    // }

    Particle_Splitting(split_index, not_moving, x, y, h);

    num_particle = x.size();

    split_index.clear();
    for (int i = 0; i < num_particle; i++)
    {
        if (x[i] < 0.7*eax  && x[i] > 0.3*eax && y[i] > 0.3*eay && y[i] < 0.7*eay)
        {
            split_index.push_back(i);
        }
    }

    Particle_Splitting(split_index, not_moving, x, y, h);

    num_particle = x.size();

    vector<double> Laplacian_T_analytic(num_particle), T(num_particle), Laplacian_T_LSMPS(num_particle), T_Lap(num_particle);

    vector<vector<int>> neighbor;
    vector<vector<double>> weight_data;

    double R_e = 2.1;
    vector<vector<double>> Ri_a;
    calc_Ri_a(h, R_e, 7, Ri_a);

    int ncell_x, ncell_y, ncell;
    vector<vector<vector<int>>> hash_table;
    vector<int> gridpos_x, gridpos_y;

    int loop_count = 0;
    int iter = 100;

    while (loop_count < iter)
    {
        // brute_force(x, y, h, 1, neighbor, weight_data, R_e);
        // brute_force_2(x, y, h, 1, neighbor, weight_data, R_e);
        hash_grid(x, y, h[0]*R_e, ncell_x, ncell_y, ncell, hash_table, gridpos_x, gridpos_y);
        spatial_hash_neighbor_2(x,y,h,1,R_e,ncell_x,ncell_y,gridpos_x,gridpos_y,hash_table,neighbor,weight_data);

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

    vector<double> hx(num_particle);
    vector<double> hy(num_particle);

    vector<double> x_w(num_particle);
    vector<double> y_w(num_particle);
    
    if ((eax == 1) && (eay == 1)){}
    else
    {
        # pragma omp parallel for
        for (int i = 0; i < num_particle; i++)
        {
            x_w[i] = x[i];
            y_w[i] = y[i];
            
            x[i] = x[i]/eax;
            y[i] = y[i]/eay;

            hx[i] = h[i]/eax;
            hy[i] = h[i]/eay;
        }
    }

    double n_out = 10;
    for (int i = 0; i < num_particle; i++)
    {
        k.push_back(200);
        if ((x[i] < x_left-n_out*dx) || (x[i] > x_right+n_out*dx) || (y[i] < y_bottom-n_out*dy) || (y[i] > y_top+n_out*dy))
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

    R_e = 2.8;
    // brute_force(x,y,h,1,neighbor,weight_data,R_e);
    // brute_force_2(x,y,h,1,neighbor,weight_data,R_e);

    hash_grid(x_w, y_w, h[0]*R_e, ncell_x, ncell_y, ncell, hash_table, gridpos_x, gridpos_y);
    spatial_hash_neighbor(x_w, y_w, h[0]*R_e, ncell_x, ncell_y, gridpos_x, gridpos_y, hash_table, neighbor, weight_data);

    vector<vector<vector<double>>> LSMPS_eta(num_particle, vector<vector<double>>(5));
    // calc_LSMPS_eta(LSMPS_eta, x, y, h, h, neighbor, weight_data);
    calc_LSMPS_eta_2(LSMPS_eta, x, y, hx, hy, neighbor, weight_data);

    vector<vector<double>> Sij_x(num_particle), Sij_y(num_particle);
    calculate_Sij(LSMPS_eta, Sij_x, Sij_y, hx, hy, neighbor);

    vector<double> Bi_x, Bi_y;
    calculate_Bi(Sij_x, Sij_y, x, is_dummy, neighbor, Bi_x, Bi_y);
    // vector<double> Bi_x(num_particle), Bi_y(num_particle);

    vector<vector<double>> Sij_Star_x(num_particle), Sij_Star_y(num_particle);
    calculate_Sij_Star(Sij_x, Sij_y, Sij_Star_x, Sij_Star_y, x, is_dummy, neighbor);
    // calculate_Sij_Star_2(Sij_x, Sij_y, Sij_Star_x, Sij_Star_y, x, is_dummy, neighbor);

    vector<double> kdeltaT_x, kdeltaT_y;
    calc_kdeltaT(x, y, k, T, is_dummy, neighbor, LSMPS_eta, kdeltaT_x, kdeltaT_y);

    vector<vector<double>> kdeltaT_ij;
    calc_MPS_like_value(x, y, k, T, is_dummy, neighbor, kdeltaT_ij);

    vector<vector<double>> alpha_ij, kdeltaT_ij_x, kdeltaT_ij_y;
    calc_hybrid(x, y, neighbor, alpha_ij, kdeltaT_x, kdeltaT_y, is_dummy, kdeltaT_ij, kdeltaT_ij_x, kdeltaT_ij_y);

    calc_Laplacian(k, hx, hy, neighbor, kdeltaT_ij_x, kdeltaT_ij_y, Sij_Star_x, Sij_Star_y, Bi_x, Bi_y, kdeltaT_x, kdeltaT_y, is_dummy, Laplacian_T_LSMPS);

    vector<double> Lap_Value;
    calc_Laplacian_From_Eta(LSMPS_eta, x, y, T, is_dummy, neighbor, Lap_Value);

    ofstream output;

    output.open("output/Test Laplacian/Test Function 1/Hasil_Tripleres_elips.csv");

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