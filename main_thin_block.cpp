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
    auto start_time = chrono::high_resolution_clock::now();

    double power    = 2909;
    double radii    = 0.02;
    double p_den    = power/(M_PI*radii*radii);
    double T_init   = 23.7;
    double k_init   = 0.016*T_init + 10.866;
    double cp_init  = 435;
    double rho_init = 8220;
    double T_liq    = 1336;

    double x_left   = 0;
    double x_right  = 0.1524;
    double y_bottom = 0;
    double y_top    = 0.1524;
    double z_front  = 0.0016;
    double z_back   = 0;

    double eax = 1;
    double eay = 1;
    double eaz = 1/12.7;

    double x_1 = x_left/eax;     // 0
    double x_2 = x_right/eax;    // 0.1524

    double y_1 = y_bottom/eay;   // 0
    double y_2 = y_top/eay;      // 0.1524

    double z_1 = z_back/eaz;     // 0
    double z_2 = z_front/eaz;    // 0.02032

    int nx = 120;
    int ny = 120;
    int nz = 16;

    double dx = (x_right-x_left)/(nx);
    double dy = (y_top-y_bottom)/(ny);
    double dz = (z_front-z_back)/(nz);

    double dx_out = 0.005;
    double dx_in  = 0.008;

    vector<double> x, y, z, hx, hy, hz;
    vector<double> k, cp, rho;
    vector<double> x_w, y_w, z_w;
    vector<double> q_x, q_y, q_z;
    vector<int> is_dummy;

    double n_dummy = 4;
    n_dummy += 0.5;

    double x_par = x_left-n_dummy*dx;
    double y_par = y_bottom-n_dummy*dy;
    double z_par = z_back-n_dummy*dz;

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

                if ((z_par<z_back+dz) && (z_par>z_back) && (pow(x_par-0.0762,2)+pow(y_par-0.0762,2) < pow(0.02,2)))
                {
                    q_x.push_back(0);
                    q_y.push_back(0);
                    q_z.push_back(p_den);
                }
                else
                {
                    q_x.push_back(0);
                    q_y.push_back(0);
                    q_z.push_back(0);
                }

                k.push_back(k_init);
                cp.push_back(cp_init);
                rho.push_back(rho_init);

                x_w.push_back(eay*eaz*x_par);
                y_w.push_back(eax*eaz*y_par);
                z_w.push_back(eax*eay*z_par);

                // x_w.push_back(x_par/eax);
                // y_w.push_back(y_par/eay);
                // z_w.push_back(z_par/eaz);

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

    vector<double> T(num_particle);

    vector<vector<int>> neighbor(num_particle);
    vector<vector<double>> weight_data(num_particle);

    double xi,yi,zi;

    for (int i = 0; i < num_particle; i++)
    {
        T[i] = T_init;
    }

    auto start_neighbor_search = chrono::high_resolution_clock::now();

    double R_e = 2.4;
    // brute_force_3D(x_w, y_w, z_w, hx, eay, eaz, neighbor, weight_data, R_e);
    brute_force_3D_2(x_w, y_w, z_w, hx, eay, eaz, neighbor, weight_data, R_e);

    auto end_neighbor_search = chrono::high_resolution_clock::now();

    double neighbor_time_ms = std::chrono::duration_cast <std::chrono::milliseconds> (end_neighbor_search-start_neighbor_search).count();
    printf("Neighbor Search Time        : %f second\n\n", neighbor_time_ms/1000);

    vector<vector<vector<double>>> LSMPS_eta(num_particle, vector<vector<double>>(9));
    // calc_LSMPS_eta_3D(LSMPS_eta, x, y, z, hx, hy, hz, neighbor, weight_data);
    calc_LSMPS_eta_3D_2(LSMPS_eta, x, y, z, hx, hy, hz, neighbor, weight_data);

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

    int count = 0;
    double t  = 0;
    double dt = 1e-2;

    string name = "output/Thin Block/result/out_" + to_string(count) + ".csv";

    ofstream output1;

    output1.open(name);

    output1 << "x" << "," << "y" << "," << "z" << "," << "q" << "," << "LSMPS_Conserved\n";

    for (int i = 0; i < num_particle; i++)
    {
        if (x[i] >= x_left && x[i] <= x_right && y[i] >= y_bottom && y[i] <= y_top && z[i] >= z_back && z[i] <= z_front)
        {
            output1 << x_w[i] << "," 
                << y_w[i] << ","
                << z_w[i] << ","
                << q_x[i] << ","
                << T[i] << "\n";
        }
    }

    auto start_loop_segment = chrono::high_resolution_clock::now();
    auto end_loop_segment = chrono::high_resolution_clock::now();

    while (t < 15 /*count < 1*/)
    {
        if (count % 10 == 0)
        {
            start_loop_segment = chrono::high_resolution_clock::now();
        }

        vector<double> kdeltaT_x, kdeltaT_y, kdeltaT_z;
        calc_kdeltaT_3D(x, y, z, k, T, is_dummy, neighbor, LSMPS_eta, kdeltaT_x, kdeltaT_y, kdeltaT_z);

        // auto end_kdt = chrono::high_resolution_clock::now();

        // double kdt_time_ms = std::chrono::duration_cast <std::chrono::milliseconds> (end_kdt-end_sijstar).count();
        // printf("kdT Time                : %f second\n", kdt_time_ms/1000);

        vector<vector<double>> kdeltaT_ij;
        calc_MPS_like_value_3D(x, y, z, k, T, is_dummy, neighbor, kdeltaT_ij);

        // auto end_mps = chrono::high_resolution_clock::now();

        // double mps_time_ms = std::chrono::duration_cast <std::chrono::milliseconds> (end_mps-end_kdt).count();
        // printf("MPS like value Time     : %f second\n", mps_time_ms/1000);

        vector<vector<double>> alpha_ij, kdeltaT_ij_x, kdeltaT_ij_y, kdeltaT_ij_z;
        calc_hybrid_3D(x, y, z, neighbor, alpha_ij, kdeltaT_x, kdeltaT_y, kdeltaT_z, is_dummy, kdeltaT_ij, kdeltaT_ij_x, kdeltaT_ij_y, kdeltaT_ij_z);

        // auto end_hyb = chrono::high_resolution_clock::now();

        // double hyb_time_ms = std::chrono::duration_cast <std::chrono::milliseconds> (end_hyb-end_mps).count();
        // printf("Calc Hybrid Time        : %f second\n", hyb_time_ms/1000);

        vector<double> heat_flux;
        External_Heat_Flux_3D(q_x, q_y, q_z, hx, hy, hz, heat_flux);

        // auto end_flux = chrono::high_resolution_clock::now();

        // double flux_time_ms = std::chrono::duration_cast <std::chrono::milliseconds> (end_flux-end_hyb).count();
        // printf("Calc Heat Flux Time     : %f second\n", flux_time_ms/1000);

        vector<double> dTdt;
        calc_dTdt_3D(cp, rho, hx, hy, hz, neighbor, kdeltaT_ij_x, kdeltaT_ij_y, kdeltaT_ij_z, Sij_Star_x, Sij_Star_y, Sij_Star_z, Bi_x, Bi_y, Bi_z, kdeltaT_x, kdeltaT_y, kdeltaT_z, heat_flux, is_dummy, dTdt);

        // auto end_dtdt = chrono::high_resolution_clock::now();

        // double dtdt_time_ms = std::chrono::duration_cast <std::chrono::milliseconds> (end_lap-end_flux).count();
        // printf("dTdt Time               : %f second\n", dtdt_time_ms/1000);

        time_integration(dt, T, dTdt);

        // auto end_int = chrono::high_resolution_clock::now();

        // double int_time_ms = std::chrono::duration_cast <std::chrono::milliseconds> (end_int-end_dtdt).count();
        // printf("Time Integration Time   : %f second\n", int_time_ms/1000);

        k.clear();
        # pragma omp parallel for
        for (int i = 0; i < num_particle; i++)
        {
            if (T[i] > 1336)
            {
                T[i] = 1336;
            }
            k.push_back(0.016*T[i] + 10.866);
        }

        t += dt;
        count += 1;

        if (count % 10 == 0)
        {
            end_loop_segment = chrono::high_resolution_clock::now();
            double loop_time_ms = std::chrono::duration_cast <std::chrono::milliseconds> (end_loop_segment-start_loop_segment).count();
            
            cout << count <<":"<<"\t"<< t << "\t\t Segment time: " << loop_time_ms/1000 << endl;
        }
        
        if (count % 10 == 0)
        {
            string name = "output/Thin Block/result/out_" + to_string(count) + ".csv";

            ofstream output1;

            output1.open(name);

            output1 << "x" << "," << "y" << "," << "z" << "," << "q" << "," << "LSMPS_Conserved\n";

            for (int i = 0; i < num_particle; i++)
            {
                if (x[i] >= x_left && x[i] <= x_right && y[i] >= y_bottom && y[i] <= y_top && z[i] >= z_back && z[i] <= z_front)
                {
                    output1 << x_w[i] << "," 
                        << y_w[i] << ","
                        << z_w[i] << ","
                        << q_z[i] << ","
                        << T[i] << "\n";
                }
            }
        }
    }

    auto end_calculation = chrono::high_resolution_clock::now();
  
    double calc_time_ms = std::chrono::duration_cast <std::chrono::milliseconds> (end_calculation-start_time).count();

    printf("\nNeighbor Search Time        : %f second\n", neighbor_time_ms/1000);
    printf("Calc Eta Time               : %f second\n", eta_time_ms/1000);
    printf("Sij Time                    : %f second\n", sij_time_ms/1000);
    printf("Bi Time                     : %f second\n", bi_time_ms/1000);
    printf("Sij Star Time               : %f second\n", sijstar_time_ms/1000);
    printf("Calculation Time            : %f second\n\n", calc_time_ms/1000);

    ofstream output2;
    output2.open("output/Thin Block/result/summary.csv");
    
    output2  << "Number of Particle," << x.size() <<"\n"
            << "Neighbor Search Time," << neighbor_time_ms/1000 << "\n"
            << "Calc Eta Time," << eta_time_ms/1000 << "\n"
            << "Sij Time," << sij_time_ms/1000 << "\n"
            << "Bi Time," << bi_time_ms/1000 << "\n"
            << "Sij Star Time," << sijstar_time_ms/1000 << "\n"
            << "Calculation Time," << calc_time_ms/1000 << "\n";    
}