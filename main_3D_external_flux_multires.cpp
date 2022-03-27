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

    double x_left   = 0;
    double x_right  = 0.1;
    double y_bottom = 0;
    double y_top    = 0.01;
    double z_front  = 0.01;
    double z_back   = 0;

    double x_change = 0.4;

    double x_multires = x_left+x_change*(x_right-x_left);

    double eax = 1;
    double eay = 1;
    double eaz = 1;

    // int n  = 40;
    // int nx = n/eax + 1;
    // int ny = n/eay + 1;
    // int nz = n/eaz + 1;

    int nx = 100;
    int ny = 10;
    int nz = 10;

    double ratio = 0.8;

    // int nx1 = x_change*nx;
    // int nx2 = (1-x_change)*nx*ratio;
    // int ny1 = ny;
    // int ny2 = ny*ratio;
    // int nz1 = nz;
    // int nz2 = nz*ratio;

    int nx1 = 40;
    int nx2 = 48;
    int ny1 = 10;
    int ny2 = 8;
    int nz1 = 10;
    int nz2 = 8;

    // nx += 1;
    // ny += 1;
    // nz += 1;

    double dx = (x_right-x_left)/(nx);
    double dy = (y_top-y_bottom)/(ny);
    double dz = (z_front-z_back)/(nz);

    // double dx1 = (x_multires-x_left)/(nx1);
    // double dy1 = (y_top-y_bottom)/(ny1);
    // double dz1 = (z_front-z_back)/(nz1);

    // double dx2 = (x_right-x_multires)/(nx2);
    // double dy2 = (y_top-y_bottom)/(ny2);
    // double dz2 = (z_front-z_back)/(nz2);

    double dx1 = 0.001;
    double dy1 = 0.001;
    double dz1 = 0.001;

    double dx2 = 0.00125;
    double dy2 = 0.00125;
    double dz2 = 0.00125;

    vector<double> x, y, z, hx, hy, hz;
    vector<double> k, cp, rho;
    vector<double> x_w, y_w, z_w;
    vector<double> q_x, q_y, q_z;
    vector<int> is_dummy;

    double n_dummy = 4;
    n_dummy += 0.5;

    double x_par = x_left-n_dummy*dx1;
    double y_par = y_bottom-n_dummy*dy1;
    double z_par = z_back-n_dummy*dz1;

    for (int i = 0; i < nx1 + n_dummy - 0.5; i++)
    {
        x_par += dx1;
        y_par = y_bottom-n_dummy*dy1;

        for (int j = 0; j < ny1 + 2*n_dummy - 1; j++)
        {
            y_par += dy1;
            z_par = z_back-n_dummy*dz1;

            for (int m = 0; m < nz1 + 2*n_dummy - 1; m++)
            {
                z_par += dz1;
                
                x.push_back(x_par);
                y.push_back(y_par);
                z.push_back(z_par);
                hx.push_back(dx1);
                hy.push_back(dy1);
                hz.push_back(dz1);

                if ((x_par<x_left+dx1) && (x_par>x_left))
                {
                    q_x.push_back(1e5);
                    q_y.push_back(0);
                    q_z.push_back(0);
                }
                else
                {
                    q_x.push_back(0);
                    q_y.push_back(0);
                    q_z.push_back(0);
                }

                k.push_back(200);
                cp.push_back(900);
                rho.push_back(2700);

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

    x_par = x_par + 0.5*dx1 - 0.5*dx2;

    for (int i = 0; i < nx2 + n_dummy - 0.5; i++)
    {
        x_par += dx2;
        y_par = y_bottom-n_dummy*dy2;

        for (int j = 0; j < ny2 + 2*n_dummy - 1; j++)
        {
            y_par += dy2;
            z_par = z_back-n_dummy*dz2;

            for (int m = 0; m < nz2 + 2*n_dummy - 1; m++)
            {
                z_par += dz2;
                
                x.push_back(x_par);
                y.push_back(y_par);
                z.push_back(z_par);
                hx.push_back(dx2);
                hy.push_back(dy2);
                hz.push_back(dz2);

                q_x.push_back(0);
                q_y.push_back(0);
                q_z.push_back(0);
                

                k.push_back(200);
                cp.push_back(900);
                rho.push_back(2700);

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

    vector<double> T(num_particle);

    vector<vector<int>> neighbor(num_particle);
    vector<vector<double>> weight_data(num_particle);

    double xi,yi,zi;

    for (int i = 0; i < num_particle; i++)
    {
    //     if (x[i] < 0)
    //     { xi == 0; }
    //     else if (x[i] > 1)
    //     { xi == 1; }
    //     else
    //     { xi = x[i]; }

    //     if (y[i] < 0)
    //     { yi == 0; }
    //     else if (y[i] > 1)
    //     { yi == 1; }
    //     else
    //     { yi = y[i]; }

    //     if (z[i] < 0)
    //     { zi == 0; }
    //     else if (z[i] > 1)
    //     { zi == 1; }
    //     else
    //     { zi = z[i]; }

        T[i] = 27;
    }

    auto start_neighbor_search = chrono::high_resolution_clock::now();

    // double R_e = 2.4;
    // brute_force_3D(x_w, y_w, z_w, hx, eay, eaz, neighbor, weight_data, R_e);
    double R_e = 2.9;
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
    double dt = 2e-3;
    
    string name1 = "output/3D External Flux/Multires/result 29/big/out_" + to_string(count) + ".csv";
    string name2 = "output/3D External Flux/Multires/result 29/small/out_" + to_string(count) + ".csv";

    ofstream output1, output2;

    output1.open(name1);
    output2.open(name2);

    output1 << "x" << "," << "y" << "," << "z" << "," << "qx" << "," << "LSMPS_Conserved\n";
    output2 << "x" << "," << "y" << "," << "z" << "," << "qx" << "," << "LSMPS_Conserved\n";

    for (int i = 0; i < num_particle; i++)
    {
        if (x[i] >= x_left && x[i] <= x_right && y[i] >= y_bottom && y[i] <= y_top && z[i] >= z_back && z[i] <= z_front)
        {
            if (hx[i] == dx1)
            {
                output2 << x[i] << "," 
                    << y[i] << ","
                    << z[i] << ","
                    << q_x[i] << ","
                    << T[i] << "\n";
            }
            else
            {
                output1 << x[i] << "," 
                    << y[i] << ","
                    << z[i] << ","
                    << q_x[i] << ","
                    << T[i] << "\n";
            }
            
        }
    }

    auto start_loop_segment = chrono::high_resolution_clock::now();
    auto end_loop_segment = chrono::high_resolution_clock::now();

    while (t < 100 /*count < 1*/)
    {
        if (count % 100 == 0)
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

        t += dt;
        count += 1;

        if (count % 100 == 0)
        {
            end_loop_segment = chrono::high_resolution_clock::now();
            double loop_time_ms = std::chrono::duration_cast <std::chrono::milliseconds> (end_loop_segment-start_loop_segment).count();
            
            cout << count <<":"<<"\t"<< t << "\t\t Segment time: " << loop_time_ms/1000 << endl;
        }
        
        if (count % 500 == 0)
        {
            string name1 = "output/3D External Flux/Multires/result 29/big/out_" + to_string(count) + ".csv";
            string name2 = "output/3D External Flux/Multires/result 29/small/out_" + to_string(count) + ".csv";

            ofstream output1, output2;

            output1.open(name1);
            output2.open(name2);

            output1 << "x" << "," << "y" << "," << "z" << "," << "qx" << "," << "LSMPS_Conserved\n";
            output2 << "x" << "," << "y" << "," << "z" << "," << "qx" << "," << "LSMPS_Conserved\n";

            for (int i = 0; i < num_particle; i++)
            {
                if (x[i] >= x_left && x[i] <= x_right && y[i] >= y_bottom && y[i] <= y_top && z[i] >= z_back && z[i] <= z_front)
                {
                    if (hx[i] == dx1)
                    {
                        output2 << x[i] << "," 
                            << y[i] << ","
                            << z[i] << ","
                            << q_x[i] << ","
                            << T[i] << "\n";
                    }
                    else
                    {
                        output1 << x[i] << "," 
                            << y[i] << ","
                            << z[i] << ","
                            << q_x[i] << ","
                            << T[i] << "\n";
                    }
                    
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

    ofstream output3;
    output3.open("output/3D External Flux/Multires/result 29/summary.csv");
    
    output3  << "Number of Particle," << x.size() <<"\n"
            << "Neighbor Search Time," << neighbor_time_ms/1000 << "\n"
            << "Calc Eta Time," << eta_time_ms/1000 << "\n"
            << "Sij Time," << sij_time_ms/1000 << "\n"
            << "Bi Time," << bi_time_ms/1000 << "\n"
            << "Sij Star Time," << sijstar_time_ms/1000 << "\n"
            << "Calculation Time," << calc_time_ms/1000 << "\n";    
}