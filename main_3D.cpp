#define _USE_MATH_DEFINES
#include <cmath>
#include <vector>
#include <chrono>
#include <fstream>
#include <iostream>

#include "src/LSMPS/LSMPS_3D.hpp"
#include "src/Thermal Flux Method/Heat_Flux_3D.hpp"
#include "src/Thermal Flux Method/Surface_Vector_S_3D.hpp"
#include "src/Thermal Flux Method/Calculate_Laplacian_3D.hpp"
#include "src/Time Integration/time_integration.hpp"
#include "src/Neighbor Search/Brute_Force.hpp"
#include "src/External Heat Flux/External_Heat_Flux.hpp"

using namespace std;

int main()
{
    double x_left   = 0;
    double x_right  = 1;
    double y_bottom = 0;
    double y_top    = 1;
    double z_front  = 1;
    double z_back   = 0;

    double eax = 1;
    double eay = 1;
    double eaz = 1;

    int n  = 50;
    int nx = n/eax;
    int ny = n/eay;
    int nz = n/eaz;

    // int nx = 100;
    // int ny = 10;
    // int nz = 10;

    // nx += 1;
    // ny += 1;
    // nz += 1;

    double dx = (x_right-x_left)/(nx);
    double dy = (y_top-y_bottom)/(ny);
    double dz = (z_front-z_back)/(nz);

    double dx_out = 0.005;
    double dx_in  = 0.008;

    vector<double> x, y, z, hx, hy, hz, k;
    vector<double> x_w, y_w, z_w;
    vector<double> q_x, q_y, q_z;
    vector<int> is_dummy;

    double n_dummy = 4;
    n_dummy += 0.5;

    double x_par = x_left-n_dummy*dx;
    double y_par = y_bottom-n_dummy*dy;
    double z_par = z_back-n_dummy*dz;

    // x_par = x_left - 6*dx_out;
    // y_par = y_bottom - 6*dx_out;

    // for(int i = 0; i < x_right/dx_out + 11; i++){
    //     x_par += dx_out;
    //     y_par = y_bottom - 6*dx_out;
    //     for(int j = 0 ; j < y_top/dx_out + 11; j++){
    //         y_par += dx_out;
    //         if (x_par <= 0.060001 || x_par >= 0.939999 || y_par <= 0.060001 || y_par >= 0.939999){
    //             x.push_back(x_par);
    //             y.push_back(y_par);
    //             hx.push_back(dx_out);
    //             hy.push_back(dx_out);
    //             k.push_back(200);

    //             x_w.push_back(eay*x_par);
    //             y_w.push_back(eax*y_par);
    //         }
    //     }
    // }

    // x_par = 0.06;
    // y_par = 0.06;

    // for(int i = 0; i < (x_right-0.12)/dx_in - 1; i++){
    //     x_par += dx_in;
    //     y_par = 0.06;
    //     for(int j = 0 ; j < (y_top-0.12)/dx_in - 1; j++){
    //         y_par += dx_in;
    //         x.push_back(x_par);
    //         y.push_back(y_par);
    //         hx.push_back(dx_in);
    //         hy.push_back(dx_in);
    //         k.push_back(200);

    //         x_w.push_back(eay*x_par);
    //         y_w.push_back(eax*y_par);
    //     }
    // }

    for (int i = 0; i < nx + 2*n_dummy - 1; i++)
    {
        x_par += dx;
        y_par = y_bottom-n_dummy*dy;

        for (int j = 0; j < ny + 2*n_dummy - 1; j++)
        {
            y_par += dy;
            z_par = z_back-n_dummy*dz;

            for (int m = 0; m < nz + 2*n_dummy - 1; m++)
            {
                z_par += dz;
                
                x.push_back(x_par);
                y.push_back(y_par);
                z.push_back(z_par);
                hx.push_back(dx);
                hy.push_back(dy);
                hz.push_back(dz);

                q_x.push_back(0);
                q_y.push_back(0);
                q_z.push_back(0);

                k.push_back(200);

                x_w.push_back(eay*eaz*x_par);
                y_w.push_back(eax*eaz*y_par);
                z_w.push_back(eax*eay*z_par);

                if ((x_par < x_left) || (x_par > x_right) || (y_par < y_bottom) || (y_par > y_top) || (z_par < z_back) || (z_par > z_front))
                {
                    is_dummy.push_back(1);
                }
                else
                {
                    is_dummy.push_back(0);
                }
            }
        }
    }

    int num_particle = x.size();

    vector<double> Laplacian_T_analytic(num_particle), T(num_particle), Laplacian_T_LSMPS(num_particle), T_Lap(num_particle);

    vector<vector<int>> neighbor(num_particle);
    vector<vector<double>> weight_data(num_particle);

    double xi,yi,zi;

    for (int i = 0; i < num_particle; i++)
    {
        // if (x[i] < 0)
        // { xi == 0; }
        // else if (x[i] > 1)
        // { xi == 1; }
        // else
        // { xi = x[i]; }

        // if (y[i] < 0)
        // { yi == 0; }
        // else if (y[i] > 1)
        // { yi == 1; }
        // else
        // { yi = y[i]; }

        // if (z[i] < 0)
        // { zi == 0; }
        // else if (z[i] > 1)
        // { zi == 1; }
        // else
        // { zi = z[i]; }

        if (is_dummy[i] == 1)
        {
            continue;
        }
        
        // Test Function 1
        // T[i] = -1/(12*pow(M_PI,2))*sin(2*M_PI*x[i])*sin(2*M_PI*y[i])*sin(2*M_PI*z[i]);
        // Laplacian_T_analytic[i] = sin(2*M_PI*x[i])*sin(2*M_PI*y[i])*sin(2*M_PI*z[i]);

        // Test Function 2
        T[i] = 100+50*(cos(M_PI*x[i])+cos(M_PI*y[i])+cos(M_PI*z[i]));
        Laplacian_T_analytic[i] = -50*M_PI*M_PI*(cos(M_PI*x[i])+cos(M_PI*y[i])+cos(M_PI*z[i]));
    }

    auto start_neighbor_search = chrono::high_resolution_clock::now();

    double R_e = 2.4;
    brute_force_3D(x_w, y_w, z_w, hx, eay, eaz, neighbor, weight_data, R_e);

    auto end_neighbor_search = chrono::high_resolution_clock::now();

    double neighbor_time_ms = std::chrono::duration_cast <std::chrono::milliseconds> (end_neighbor_search-start_neighbor_search).count();
    printf("Neighbor Search Time        : %f second\n\n", neighbor_time_ms/1000);

    vector<vector<vector<double>>> LSMPS_eta(num_particle, vector<vector<double>>(9));
    calc_LSMPS_eta_3D(LSMPS_eta, x, y, z, hx, hy, hz, neighbor, weight_data);

    auto end_eta = chrono::high_resolution_clock::now();

    double eta_time_ms = std::chrono::duration_cast <std::chrono::milliseconds> (end_eta-end_neighbor_search).count();
    printf("Calc Eta Time               : %f second\n\n", eta_time_ms/1000);

    vector<vector<double>> Sij_x(num_particle), Sij_y(num_particle), Sij_z(num_particle);
    calculate_Sij_3D(LSMPS_eta, Sij_x, Sij_y, Sij_z, hx, hy, hz, neighbor);

    auto end_sij = chrono::high_resolution_clock::now();

    double sij_time_ms = std::chrono::duration_cast <std::chrono::milliseconds> (end_sij-end_eta).count();
    printf("Sij Time                : %f second\n", sij_time_ms/1000);

    vector<double> Bi_x, Bi_y, Bi_z;
    calculate_Bi_3D(Sij_x, Sij_y, Sij_z, x, is_dummy, neighbor, Bi_x, Bi_y, Bi_z);

    auto end_bi = chrono::high_resolution_clock::now();

    double bi_time_ms = std::chrono::duration_cast <std::chrono::milliseconds> (end_bi-end_sij).count();
    printf("Bi Time                 : %f second\n", bi_time_ms/1000);

    vector<vector<double>> Sij_Star_x(num_particle), Sij_Star_y(num_particle), Sij_Star_z(num_particle);
    calculate_Sij_Star_3D(Sij_x, Sij_y, Sij_z, Sij_Star_x, Sij_Star_y, Sij_Star_z, x, is_dummy, neighbor);

    auto end_sijstar = chrono::high_resolution_clock::now();

    double sijstar_time_ms = std::chrono::duration_cast <std::chrono::milliseconds> (end_sijstar-end_bi).count();
    printf("Sij Star Time           : %f second\n", sijstar_time_ms/1000);

    vector<double> kdeltaT_x, kdeltaT_y, kdeltaT_z;
    calc_kdeltaT_3D(x, y, z, k, T, is_dummy, neighbor, LSMPS_eta, kdeltaT_x, kdeltaT_y, kdeltaT_z);

    auto end_kdt = chrono::high_resolution_clock::now();

    double kdt_time_ms = std::chrono::duration_cast <std::chrono::milliseconds> (end_kdt-end_sijstar).count();
    printf("kdT Time                : %f second\n", kdt_time_ms/1000);


    vector<vector<double>> kdeltaT_ij;
    calc_MPS_like_value_3D(x, y, z, k, T, is_dummy, neighbor, kdeltaT_ij);

    auto end_mps = chrono::high_resolution_clock::now();

    double mps_time_ms = std::chrono::duration_cast <std::chrono::milliseconds> (end_mps-end_kdt).count();
    printf("MPS like value Time     : %f second\n", mps_time_ms/1000);

    vector<vector<double>> alpha_ij, kdeltaT_ij_x, kdeltaT_ij_y, kdeltaT_ij_z;
    calc_hybrid_3D(x, y, z, neighbor, alpha_ij, kdeltaT_x, kdeltaT_y, kdeltaT_z, is_dummy, kdeltaT_ij, kdeltaT_ij_x, kdeltaT_ij_y, kdeltaT_ij_z);

    auto end_hyb = chrono::high_resolution_clock::now();

    double hyb_time_ms = std::chrono::duration_cast <std::chrono::milliseconds> (end_hyb-end_mps).count();
    printf("Calc Hybrid Time        : %f second\n", hyb_time_ms/1000);

    vector<double> heat_flux;
    External_Heat_Flux_3D(q_x, q_y, q_z, hx, hy, hz, heat_flux);

    auto end_flux = chrono::high_resolution_clock::now();

    double flux_time_ms = std::chrono::duration_cast <std::chrono::milliseconds> (end_flux-end_hyb).count();
    printf("Calc Heat Flux Time     : %f second\n", flux_time_ms/1000);

    calc_Laplacian_3D(k, hx, hy, hz, neighbor, kdeltaT_ij_x, kdeltaT_ij_y, kdeltaT_ij_z, Sij_Star_x, Sij_Star_y, Sij_Star_z, Bi_x, Bi_y, Bi_z, kdeltaT_x, kdeltaT_y, kdeltaT_z, heat_flux, is_dummy, Laplacian_T_LSMPS);

    auto end_lap = chrono::high_resolution_clock::now();

    double calc_lap_time_ms = std::chrono::duration_cast <std::chrono::milliseconds> (end_lap-end_flux).count();
    printf("Calc Laplacian Time     : %f second\n", calc_lap_time_ms/1000);

    vector<double> Lap_Value;
    calc_Laplacian_From_Eta_3D(LSMPS_eta, x, y, z, T, is_dummy, neighbor, Lap_Value);

    ofstream output1, output2;

    output1.open("output/Test Laplacian 3D/Test Function 2/Hasil/Hasil_3D_dummy.csv");

    output1 << "x" << "," << "y" << "," << "z" << "," << "Analytic" << "," << "LSMPS_Conserved" << "," << "LSMPS_Laplacian\n";

    for (int i = 0; i < num_particle; i++)
    {
        if (x[i] >= x_left && x[i] <= x_right && y[i] >= y_bottom && y[i] <= y_top && z[i] >= z_back && z[i] <= z_front)
        {
            output1 << x[i] << "," 
                << y[i] << ","
                << z[i] << ","
                << Laplacian_T_analytic[i] << ","
                << Laplacian_T_LSMPS[i] << "," 
                << Lap_Value[i] << "\n";
        }
    }

    output2.open("output/Test Laplacian 3D/Test Function 2/Summary/Hasil_3D_dummy.csv");
    
    output2  << "Number of Particle," << x.size() <<"\n"
            << "Neighbor Search Time," << neighbor_time_ms/1000 << "\n"
            << "Calc Eta Time," << eta_time_ms/1000 << "\n"
            << "Sij Time," << sij_time_ms/1000 << "\n"
            << "Bi Time," << bi_time_ms/1000 << "\n"
            << "Sij Star Time," << sijstar_time_ms/1000 << "\n"
            << "kdt Time," << kdt_time_ms/1000 << "\n"
            << "MPS Like Time," << mps_time_ms/1000 << "\n"
            << "Hybrid Time," << hyb_time_ms/1000 << "\n"
            << "Heat Flux Time," << flux_time_ms/1000 << "\n"
            << "Calculate Laplacian Time," << calc_lap_time_ms/1000 << "\n";
}