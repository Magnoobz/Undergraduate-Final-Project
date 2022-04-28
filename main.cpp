#define _USE_MATH_DEFINES
#include <cmath>
#include <vector>
#include <fstream>
#include <iostream>

#include "src/LSMPS/LSMPS.hpp"
#include "src/Thermal Flux Method/Heat_Flux.hpp"
#include "src/Thermal Flux Method/Surface_Vector_S.hpp"
#include "src/Thermal Flux Method/Calculate_Laplacian.hpp"
#include "src/Time Integration/time_integration.hpp"

#include "src/Neighbor Search/Brute_Force.hpp"
#include "src/Neighbor Search/Spatial_Hash.hpp"

using namespace std;

int main()
{    
    double x_left   = 0;
    double x_right  = 1;
    double y_bottom = 0;
    double y_top    = 1;

    double eax = 1;
    double eay = 2;

    int n  = 100;
    int nx = n/eax;
    int ny = n/eay;

    double dx = (x_right-x_left)/(nx);
    double dy = (y_top-y_bottom)/(ny);

    double dx_out = 0.01;
    double dy_out = dx_out*eay/eax;

    double dx_in  = 0.012;
    double dy_in  = dx_in*eay/eax;

    vector<double> x, y, hx, hy, k;
    vector<double> x_w, y_w;
    vector<int> is_dummy;

    double n_dummy = 6;
    n_dummy += 0.5;

    double x_par = x_left-n_dummy*dx_out;
    double y_par = y_bottom-n_dummy*dy_out;

    double n_out = 4;
    for(int i = 0; i < (x_right-x_left)/dx_out + 2*n_dummy-1; i++){
        x_par += dx_out;
        y_par = y_bottom - n_dummy*dy_out;
        for(int j = 0 ; j < (y_top-y_bottom)/dy_out + 2*n_dummy-1; j++){
            y_par += dy_out;
            if (x_par <= 0.2 || x_par >= 0.8 || y_par <= 0.2 || y_par >= 0.8){
                x.push_back(x_par);
                y.push_back(y_par);
                hx.push_back(dx_out);
                hy.push_back(dy_out);
                k.push_back(200);

                x_w.push_back(eay*x_par);
                y_w.push_back(eax*y_par);

                if ((x_par < x_left-n_out*dx_out) || (x_par > x_right+n_out*dx_out) || (y_par < y_bottom-n_out*dy_out) || (y_par > y_top+n_out*dy_out))
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

    x_par = 0.2 - 0.5*dx_in;
    y_par = 0.2 - 0.5*dy_in;

    for(int i = 0; i < (x_right-x_left-0.4)/dx_in; i++){
        x_par += dx_in;
        y_par = 0.2 - 0.5*dy_in;
        for(int j = 0 ; j < (y_top-y_bottom-0.4)/dy_in; j++){
            y_par += dy_in;
            x.push_back(x_par);
            y.push_back(y_par);
            hx.push_back(dx_in);
            hy.push_back(dy_in);
            k.push_back(200);

            x_w.push_back(eay*x_par);
            y_w.push_back(eax*y_par);

            is_dummy.push_back(0);
        }
    }

    // double x_par = x_left-n_dummy*dx;
    // double y_par = y_bottom-n_dummy*dy;

    // for (int i = 0; i < nx+2*n_dummy - 1; i++)
    // {
    //     x_par += dx;
    //     y_par = y_bottom-n_dummy*dy;

    //     for (int j = 0; j < ny+2*n_dummy - 1; j++)
    //     {
    //         y_par += dy;

    //         x.push_back(x_par);
    //         y.push_back(y_par);
    //         hx.push_back(dx);
    //         hy.push_back(dy);
    //         k.push_back(200);

    //         x_w.push_back(eay*x_par);
    //         y_w.push_back(eax*y_par);

    //         if ((x_par < x_left) || (x_par > x_right) || (y_par < y_bottom) || (y_par > y_top))
    //         {
    //             is_dummy.push_back(1);
    //         }
    //         else
    //         {
    //             is_dummy.push_back(0);
    //         }
    //     }
    // }

    int num_particle = x.size();

    vector<double> Laplacian_T_analytic(num_particle), T(num_particle), Laplacian_T_LSMPS(num_particle), T_Lap(num_particle);

    vector<vector<int>> neighbor(num_particle);
    vector<vector<double>> weight_data(num_particle);

    double xi,yi;

    for (int i = 0; i < num_particle; i++)
    {        
        if (is_dummy[i] == 1)
        {
            continue;
        }

        // Test Function 1
        // T[i] = -1/(8*pow(M_PI,2))*sin(2*M_PI*x[i])*sin(2*M_PI*y[i]);
        // Laplacian_T_analytic[i] = sin(2*M_PI*x[i])*sin(2*M_PI*y[i]);

        // Test Function 2
        T[i] = 100+50*(cos(M_PI*x[i])+cos(M_PI*y[i]));
        Laplacian_T_analytic[i] = -50*pow(M_PI,2)*(cos(M_PI*x[i])+cos(M_PI*y[i]));

    }

    double R_e = 4.2;
    // brute_force(x_w,y_w,hx,eay,neighbor,weight_data,R_e);
    // brute_force_2(x_w,y_w,hx,eay,neighbor,weight_data,R_e);

    int ncell_x, ncell_y, ncell;
    vector<vector<vector<int>>> hash_table;
    vector<int> gridpos_x, gridpos_y;
    hash_grid(x_w, y_w, hx[0]*eay*R_e, ncell_x, ncell_y, ncell, hash_table, gridpos_x, gridpos_y);
    spatial_hash_neighbor(x_w, y_w, hx[0]*eay*R_e, ncell_x, ncell_y, gridpos_x, gridpos_y, hash_table, neighbor, weight_data);

    vector<vector<vector<double>>> LSMPS_eta(num_particle, vector<vector<double>>(5));
    // calc_LSMPS_eta(LSMPS_eta, x, y, hx, hy, neighbor, weight_data);
    calc_LSMPS_eta_2(LSMPS_eta, x, y, hx, hy, neighbor, weight_data);

    vector<vector<double>> Sij_x(num_particle), Sij_y(num_particle);
    calculate_Sij(LSMPS_eta, Sij_x, Sij_y, hx, hy, neighbor);

    vector<double> Bi_x, Bi_y;
    calculate_Bi(Sij_x, Sij_y, x, is_dummy, neighbor, Bi_x, Bi_y);
    // vector<double> Bi_x(num_particle), Bi_y(num_particle);

    vector<vector<double>> Sij_Star_x(num_particle), Sij_Star_y(num_particle);
    calculate_Sij_Star(Sij_x, Sij_y, Sij_Star_x, Sij_Star_y, x, is_dummy, neighbor);

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

    string name = "output/Test Laplacian/Test Function 2/Hasil_SpatialHash" + to_string(dx_in/dx_out) + "_" + to_string(eay/eax) + "_" + to_string(R_e) + ".csv";

    output.open(name);
    // output.open("Hasil_tes_ubah_LSMPS_2.csv");

    output << "x" << "," << "y" << "," << "Analytic" << "," << "LSMPS_Conserved" << "," << "hx\n";

    for (int i = 0; i < num_particle; i++)
    {
        if (x[i] >= x_left && x[i] <= x_right && y[i] >= y_bottom && y[i] <= y_top)
        {
            output << x[i] << "," 
                << y[i] << ","
                << Laplacian_T_analytic[i] << ","
                << Laplacian_T_LSMPS[i] << "," 
                << hx[i] << "\n";
        }
    }
}