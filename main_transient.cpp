#define _USE_MATH_DEFINES
#include <cmath>
#include <vector>
#include <chrono>
#include <iostream>
#include <fstream>

#include "src/LSMPS/LSMPS.hpp"
#include "src/Thermal Flux Method/Heat_Flux.hpp"
#include "src/Thermal Flux Method/Surface_Vector_S.hpp"
#include "src/Thermal Flux Method/Calculate_Laplacian.hpp"
#include "src/Time Integration/time_integration.hpp"
#include "src/Neighbor Search/Brute_Force.hpp"

using namespace std;

int main()
{
    double x_left   = 0;
    double x_right  = 1;
    double y_bottom = 0;
    double y_top    = 1;

    double t  = 0;
    double dt = 1e-2;

    double eax = 1;
    double eay = 1;

    int n  = 80;
    int nx = ceil(n/eax);
    int ny = ceil(n/eay);

    double dx = (x_right-x_left)/(nx);
    double dy = (y_top-y_bottom)/(ny);

    vector<double> x, y, hx, hy, k, cp, rho;
    vector<double> x_w, y_w;
    vector<int> is_dummy;

    double n_dummy = 4;
    n_dummy += 0.5;

    double x_par = x_left-n_dummy*dx;
    double y_par = y_bottom-n_dummy*dy;


    for (int i = 0; i < nx+2*n_dummy - 1; i++)
    {
        x_par += dx;
        y_par = y_bottom-n_dummy*dy;

        for (int j = 0; j < ny+2*n_dummy - 1; j++)
        {
            y_par += dy;

            x.push_back(x_par);
            y.push_back(y_par);
            hx.push_back(dx);
            hy.push_back(dy);
            k.push_back(200);
            cp.push_back(900);
            rho.push_back(2700);

            x_w.push_back(eay*x_par);
            y_w.push_back(eax*y_par);

            if ((x_par < x_left) || (x_par > x_right) || (y_par < y_bottom) || (y_par > y_top))
            {
                is_dummy.push_back(1);
            }
            else
            {
                is_dummy.push_back(0);
            }
        }
    }

    int num_particle = x.size();

    vector<double> Laplacian_T_analytic(num_particle), T(num_particle), Laplacian_T_LSMPS(num_particle), T_Lap(num_particle);

    vector<vector<int>> neighbor(num_particle);
    vector<vector<double>> weight_data(num_particle);

    for (int i = 0; i < num_particle; i++)
    {
        if (is_dummy[i] == 1)
        {
            continue;
        }
        
        T[i] = -1/(8*pow(M_PI,2))*sin(2*M_PI*x[i])*sin(2*M_PI*y[i]);
    }

    auto start_neighbor_search = chrono::high_resolution_clock::now();

    double R_e = 2.1;
    brute_force(x_w,y_w,hx,eay,neighbor,weight_data,R_e);
    
    auto end_neighbor_search = chrono::high_resolution_clock::now();

    double neighbor_time_ms = std::chrono::duration_cast <std::chrono::milliseconds> (end_neighbor_search-start_neighbor_search).count();
    printf("Neighbor Search Time        : %f second\n\n", neighbor_time_ms/1000);

    vector<vector<vector<double>>> LSMPS_eta(num_particle, vector<vector<double>>(5));
    calc_LSMPS_eta(LSMPS_eta, x, y, hx, hy, neighbor, weight_data);

    auto end_eta = chrono::high_resolution_clock::now();

    double eta_time_ms = std::chrono::duration_cast <std::chrono::milliseconds> (end_eta-end_neighbor_search).count();
    printf("Calc Eta Time               : %f second\n\n", eta_time_ms/1000);

    vector<vector<double>> Sij_x(num_particle), Sij_y(num_particle);
    calculate_Sij(LSMPS_eta, Sij_x, Sij_y, hx, hy, neighbor);

    auto end_sij = chrono::high_resolution_clock::now();

    double sij_time_ms = std::chrono::duration_cast <std::chrono::milliseconds> (end_sij-end_eta).count();
    printf("Sij Time                : %f second\n", sij_time_ms/1000);

    vector<double> Bi_x, Bi_y;
    calculate_Bi(Sij_x, Sij_y, x, is_dummy, neighbor, Bi_x, Bi_y);

    auto end_bi = chrono::high_resolution_clock::now();

    double bi_time_ms = std::chrono::duration_cast <std::chrono::milliseconds> (end_bi-end_sij).count();
    printf("Bi Time                 : %f second\n", bi_time_ms/1000);

    vector<vector<double>> Sij_Star_x(num_particle), Sij_Star_y(num_particle);
    calculate_Sij_Star(Sij_x, Sij_y, Sij_Star_x, Sij_Star_y, x, is_dummy, neighbor);

    auto end_sijstar = chrono::high_resolution_clock::now();

    double sijstar_time_ms = std::chrono::duration_cast <std::chrono::milliseconds> (end_sijstar-end_bi).count();
    printf("Sij Star Time           : %f second\n", sijstar_time_ms/1000);

    int count = 0;

    auto start_loop_segment = chrono::high_resolution_clock::now();
    auto end_loop_segment = chrono::high_resolution_clock::now();

    while (t < 200 /*count < 1*/)
    {
        if (count % 10 == 0)
        {
            start_loop_segment = chrono::high_resolution_clock::now();
        }

        vector<double> kdeltaT_x, kdeltaT_y;
        calc_kdeltaT(x, y, k, T, is_dummy, neighbor, LSMPS_eta, kdeltaT_x, kdeltaT_y);

        // auto end_kdt = chrono::high_resolution_clock::now();

        // double kdt_time_ms = std::chrono::duration_cast <std::chrono::milliseconds> (end_kdt-end_sijstar).count();
        // printf("kdT Time                : %f second\n", kdt_time_ms/1000);

        vector<vector<double>> kdeltaT_ij;
        calc_MPS_like_value(x, y, k, T, is_dummy, neighbor, kdeltaT_ij);

        // auto end_mps = chrono::high_resolution_clock::now();

        // double mps_time_ms = std::chrono::duration_cast <std::chrono::milliseconds> (end_mps-end_kdt).count();
        // printf("MPS like value Time     : %f second\n", mps_time_ms/1000);

        vector<vector<double>> alpha_ij, kdeltaT_ij_x, kdeltaT_ij_y;
        calc_hybrid(x, y, neighbor, alpha_ij, kdeltaT_x, kdeltaT_y, is_dummy, kdeltaT_ij, kdeltaT_ij_x, kdeltaT_ij_y);

        // auto end_hyb = chrono::high_resolution_clock::now();

        // double hyb_time_ms = std::chrono::duration_cast <std::chrono::milliseconds> (end_hyb-end_mps).count();
        // printf("Calc Hybrid Time        : %f second\n", hyb_time_ms/1000);

        vector<double> dTdt;
        calc_dTdt(cp, rho, hx, hy, neighbor, kdeltaT_ij_x, kdeltaT_ij_y, Sij_Star_x, Sij_Star_y, Bi_x, Bi_y,kdeltaT_x,kdeltaT_y, is_dummy, dTdt);

        // auto end_dtdt = chrono::high_resolution_clock::now();

        // double dtdt_time_ms = std::chrono::duration_cast <std::chrono::milliseconds> (end_dtdt-end_hyb).count();
        // printf("dTdt Time               : %f second\n", dtdt_time_ms/1000);

        time_integration(dt, T, dTdt);

        // auto end_int = chrono::high_resolution_clock::now();

        // double int_time_ms = std::chrono::duration_cast <std::chrono::milliseconds> (end_int-end_dtdt).count();
        // printf("Time Integration Time   : %f second\n", int_time_ms/1000);

        t += dt;
        count += 1;
        
        if (count % 10 == 0)
        {
            end_loop_segment = chrono::high_resolution_clock::now();
            double loop_time_ms = std::chrono::duration_cast <std::chrono::milliseconds> (end_loop_segment-start_loop_segment).count();
            
            cout << count <<":"<<"\t"<< t << "\t\t Segment time: " << loop_time_ms/1000 << endl;
        }

        if (count == 1 || count % 1000 == 0)
        {
            string name = "output/Laplacian Transient/result/out_" + to_string(count) + ".csv";

            calc_Laplacian(k, hx, hy, neighbor, kdeltaT_ij_x, kdeltaT_ij_y, Sij_Star_x, Sij_Star_y,Bi_x, Bi_y,kdeltaT_x,kdeltaT_y, is_dummy, Laplacian_T_LSMPS);

            for (int i = 0; i < num_particle; i++)
            {
                Laplacian_T_analytic[i] = sin(2*M_PI*x[i])*sin(2*M_PI*y[i])*exp(-8*pow(M_PI,2)*t/rho[i]/cp[i]);
            }

            ofstream output;

            output.open(name);

            output << "x" << "," << "y" << "," << "Analytic" << "," << "LSMPS_Conserved\n";

            for (int i = 0; i < num_particle; i++)
            {
                if (x[i] >= x_left && x[i] <= x_right && y[i] >= y_bottom && y[i] <= y_top)
                {
                    output << x[i] << "," 
                        << y[i] << ","
                        << Laplacian_T_analytic[i] << ","
                        << Laplacian_T_LSMPS[i] << "\n";
                }
            }
        }              
    }

    auto end_calculation = chrono::high_resolution_clock::now();

    
    double calc_time_ms = std::chrono::duration_cast <std::chrono::milliseconds> (end_calculation-end_sijstar).count();

    printf("\nNeighbor Search Time        : %f second\n\n", neighbor_time_ms/1000);
    printf("Calc Eta Time               : %f second\n\n", eta_time_ms/1000);
    printf("Sij Time                    : %f second\n", sij_time_ms/1000);
    printf("Bi Time                     : %f second\n", bi_time_ms/1000);
    printf("Sij Star Time               : %f second\n", sijstar_time_ms/1000);
    printf("Calculation Time            : %f second\n\n", calc_time_ms/1000);

    string name = "output/Laplacian Transient/result/Summary.csv";
    
    ofstream summary;
    summary.open(name);

    summary << "Number of Particle," << x.size() <<"\n\n"
            << "Simulation Time," << t << "\n"
            << "Simulation Step," << count << "\n"
            << "Neighbor Search Time," << neighbor_time_ms/1000 << "\n"
            << "Calc Eta Time," << eta_time_ms/1000 << "\n"
            << "Sij Time," << sij_time_ms/1000 << "\n"
            << "Bi Time," << bi_time_ms/1000 << "\n"
            << "Sij Star Time," << sijstar_time_ms/1000 << "\n"
            << "Calculation Time," << calc_time_ms/1000 << "\n";
}