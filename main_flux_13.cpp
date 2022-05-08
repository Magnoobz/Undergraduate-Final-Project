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
#include "src/External Heat Flux/External_Heat_Flux.hpp"

#include "src/Particle/Initialize_Particle_3D.hpp"
#include "src/Particle/Particle_Splitting_3D.hpp"
#include "src/Particle/Multires_Movement_3D.hpp"
#include "src/Particle/Packing_Ratio_3D.hpp"

#include "src/Neighbor Search/Brute_Force.hpp"
#include "src/Neighbor Search/Spatial_Hash_3D.hpp"

using namespace std;

int main()
{
    auto start_time = chrono::high_resolution_clock::now();

    int option = 14;

    double x_left   = 0;
    double x_right  = 0.1;
    double y_bottom = 0;
    double y_top    = 0.01;
    double z_front  = 0.01;
    double z_back   = 0;

    double eax = 1.5;
    double eay = 1;
    double eaz = 1;

    int nx = 150;
    int ny = 10;
    int nz = 10;

    double dx = (x_right-x_left)/nx;

    vector<double> x, y, z, hx, hy, hz, h_temp;
    vector<double> k, cp, rho;
    vector<double> x_w, y_w, z_w;
    vector<double> q_x, q_y, q_z;
    vector<int> is_dummy, not_moving;

    double n_dummy = 3;
    
    initialize_particle_3D_2(x_left,0.03,y_bottom,y_top,z_back,z_front,45,ny,nz,0,3,3,3,3,3,not_moving,x_w,y_w,z_w,h_temp);
    initialize_particle_3D_2(0.03,0.04,y_bottom,y_top,z_back,z_front,12,8,8,0,0,3,3,3,3,not_moving,x_w,y_w,z_w,h_temp);
    initialize_particle_3D_2(0.04,x_right,y_bottom,y_top,z_back,z_front,54,6,6,3,0,3,3,3,3,not_moving,x_w,y_w,z_w,h_temp);
    int num_particle = x_w.size();

    vector<double> T(num_particle);

    vector<vector<int>> neighbor(num_particle);
    vector<vector<double>> weight_data(num_particle);


    x.resize(num_particle);
    y.resize(num_particle);
    z.resize(num_particle);

    hx.resize(num_particle);
    hy.resize(num_particle);
    hz.resize(num_particle);

    k.resize(num_particle);
    cp.resize(num_particle);
    rho.resize(num_particle);

    q_x.resize(num_particle);
    q_y.resize(num_particle);
    q_z.resize(num_particle);

    is_dummy.resize(num_particle);

    # pragma omp parallel for
    for (int i = 0; i < num_particle; i++)
    {
        x[i] = x_w[i]/eax;
        y[i] = y_w[i]/eay;
        z[i] = z_w[i]/eaz;
        
        hx[i] = h_temp[i]/eax;
        hy[i] = h_temp[i]/eay;
        hz[i] = h_temp[i]/eaz;

        if ((x[i] < x_left)||(x[i] >x_right)||(y[i] < y_bottom)||(y[i] > y_top)||(z[i] < z_back)||(z[i] > z_front))
        {
            is_dummy[i] = 1;
        }        
        else
        {
            is_dummy[i] = 0;
        }

        k[i]   = 200;
        cp[i]  = 900;
        rho[i] = 2700;

        T[i] = 27;
    }
    
    double area_flux = 0;

    # pragma omp parallel for
    for (int i = 0; i < num_particle; i++)
    {
        if (is_dummy[i] == 1){continue;}

        if (x[i] < dx)
        {
            area_flux += hy[i]*hz[i];
        }
    }

    double P_den = 10/area_flux;
    # pragma omp parallel for
    for (int i = 0; i < num_particle; i++)
    {
        if (is_dummy[i] == 1){continue;}

        if (x[i] < dx)
        {
            q_x[i] = P_den;
        }
    }

    double total_energy = 0;
    # pragma omp parallel for
    for (int i = 0; i < num_particle; i++)
    {
        if (is_dummy[i] == 1){continue;}
        
        total_energy += (T[i]+273.15)*cp[i]*rho[i]*hx[i]*hy[i]*hz[i];
    }

    int ncell_x, ncell_y, ncell_z, ncell;
    vector<vector<vector<vector<int>>>> hash_table;
    vector<int> gridpos_x, gridpos_y, gridpos_z;

    double R_e = 2.9;
    hash_grid_3D(x_w, y_w, z_w, h_temp[0]*R_e,ncell_x, ncell_y, ncell_z, ncell, hash_table, gridpos_x, gridpos_y, gridpos_z);
    spatial_hash_neighbor_3D(x_w,y_w, z_w,h_temp[0]*R_e,ncell_x, ncell_y, ncell_z,gridpos_x, gridpos_y, gridpos_z, hash_table, neighbor, weight_data);

    vector<vector<vector<double>>> LSMPS_eta(num_particle, vector<vector<double>>(9));
    calc_LSMPS_eta_3D_2(LSMPS_eta, x, y, z, hx, hy, hz, neighbor, weight_data);

    vector<vector<double>> Sij_x(num_particle), Sij_y(num_particle), Sij_z(num_particle);
    calculate_Sij_3D(LSMPS_eta, Sij_x, Sij_y, Sij_z, hx, hy, hz, neighbor);

    vector<double> Bi_x, Bi_y, Bi_z;
    calculate_Bi_3D(Sij_x, Sij_y, Sij_z, x, is_dummy, neighbor, Bi_x, Bi_y, Bi_z);

    vector<vector<double>> Sij_Star_x(num_particle), Sij_Star_y(num_particle), Sij_Star_z(num_particle);
    calculate_Sij_Star_3D(Sij_x, Sij_y, Sij_z, Sij_Star_x, Sij_Star_y, Sij_Star_z, x, is_dummy, neighbor);

    int count = 0;
    double t  = 0;
    double dt = 2e-3;
    
    ofstream energy;
    energy.open("output/3D/result_" + to_string(option) + "/total_energy.csv");
    energy << "t,energy\n";

    energy << t << "," << total_energy << "\n";
    
    string name0 = "output/3D/result_" + to_string(option) + "/out_" + to_string(count) + ".csv";
    // string name1 = "output/3D/result_" + to_string(option) + "/1_out_" + to_string(count) + ".csv";
    // string name2 = "output/3D/result_" + to_string(option) + "/2_out_" + to_string(count) + ".csv";
    // string name3 = "output/3D/result_" + to_string(option) + "/3_out_" + to_string(count) + ".csv";

    ofstream output0, output1, output2, output3;

    output0.open(name0);
    // output1.open(name1);
    // output2.open(name2);
    // output3.open(name3);

    output0 << "x" << "," << "y" << "," << "z" << "," << "LSMPS_Conserved\n";
    // output1 << "x" << "," << "y" << "," << "z" << "," << "LSMPS_Conserved\n";
    // output2 << "x" << "," << "y" << "," << "z" << "," << "LSMPS_Conserved\n";
    // output3 << "x" << "," << "y" << "," << "z" << "," << "LSMPS_Conserved\n";

    for (int i = 0; i < num_particle; i++)
    {
        if (x[i] >= x_left && x[i] <= x_right && y[i] >= y_bottom && y[i] <= y_top && z[i] >= z_back && z[i] <= z_front)
        {
            output0 << x[i] << "," << y[i] << "," << z[i] << "," << T[i] << "\n";
            
            // if ((h_temp[i] < 1.02*h_temp[0]) && (h_temp[i] > 0.98*h_temp[0]))
            // {
            //     output1 << x[i] << "," << y[i] << "," << z[i] << "," << T[i] << "\n";
            // }
            // else if ((h_temp[i] < 1.14*h_temp[0]) && (h_temp[i] > 1.09*h_temp[0]))
            // {
            //     output2 << x[i] << "," << y[i] << "," << z[i] << "," << T[i] << "\n";
            // }
            // else if ((h_temp[i] < 1.28*h_temp[0]) && (h_temp[i] > 1.22*h_temp[0]))
            // {
            //     output3 << x[i] << "," << y[i] << "," << z[i] << "," << T[i] << "\n";
            // }
        }
    }

    auto start_loop_segment = chrono::high_resolution_clock::now();
    auto end_loop_segment = chrono::high_resolution_clock::now();

    while (t < 100)
    {
        if (count % 100 == 0)
        {
            start_loop_segment = chrono::high_resolution_clock::now();
        }

        vector<double> kdeltaT_x, kdeltaT_y, kdeltaT_z;
        calc_kdeltaT_3D(x, y, z, k, T, is_dummy, neighbor, LSMPS_eta, kdeltaT_x, kdeltaT_y, kdeltaT_z);

        vector<vector<double>> kdeltaT_ij;
        calc_MPS_like_value_3D(x, y, z, k, T, is_dummy, neighbor, kdeltaT_ij);

        vector<vector<double>> alpha_ij, kdeltaT_ij_x, kdeltaT_ij_y, kdeltaT_ij_z;
        calc_hybrid_3D(x, y, z, neighbor, alpha_ij, kdeltaT_x, kdeltaT_y, kdeltaT_z, is_dummy, kdeltaT_ij, kdeltaT_ij_x, kdeltaT_ij_y, kdeltaT_ij_z);

        vector<double> heat_flux;
        External_Heat_Flux_3D(q_x, q_y, q_z, hx, hy, hz, heat_flux);

        vector<double> dTdt;
        calc_dTdt_3D(cp, rho, hx, hy, hz, neighbor, kdeltaT_ij_x, kdeltaT_ij_y, kdeltaT_ij_z, Sij_Star_x, Sij_Star_y, Sij_Star_z, Bi_x, Bi_y, Bi_z, kdeltaT_x, kdeltaT_y, kdeltaT_z, heat_flux, is_dummy, dTdt);

        time_integration(dt, T, dTdt);

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
            total_energy = 0;
            # pragma omp parallel for
            for (int i = 0; i < num_particle; i++)
            {
                if (is_dummy[i] == 1){continue;}
                
                total_energy += (T[i]+273.15)*cp[i]*rho[i]*hx[i]*hy[i]*hz[i];
            }
            energy << t << "," << total_energy << "\n";

            string name10 = "output/3D/result_" + to_string(option) + "/out_" + to_string(count) + ".csv";
            // string name11 = "output/3D/result_" + to_string(option) + "/1_out_" + to_string(count) + ".csv";
            // string name12 = "output/3D/result_" + to_string(option) + "/2_out_" + to_string(count) + ".csv";
            // string name13 = "output/3D/result_" + to_string(option) + "/3_out_" + to_string(count) + ".csv";

            ofstream output10, output11, output12, output13;

            output10.open(name10);
            // output11.open(name11);
            // output12.open(name12);
            // output13.open(name13);

            output10 << "x" << "," << "y" << "," << "z" << "," << "LSMPS_Conserved\n";
            // output11 << "x" << "," << "y" << "," << "z" << "," << "LSMPS_Conserved\n";
            // output12 << "x" << "," << "y" << "," << "z" << "," << "LSMPS_Conserved\n";
            // output13 << "x" << "," << "y" << "," << "z" << "," << "LSMPS_Conserved\n";

            for (int i = 0; i < num_particle; i++)
            {
                if (x[i] >= x_left && x[i] <= x_right && y[i] >= y_bottom && y[i] <= y_top && z[i] >= z_back && z[i] <= z_front)
                {
                    output10 << x[i] << "," << y[i] << "," << z[i] << "," << T[i] << "\n";
                    
                    // if ((h_temp[i] < 1.02*h_temp[0]) && (h_temp[i] > 0.98*h_temp[0]))
                    // {
                    //     output11 << x[i] << "," << y[i] << "," << z[i] << "," << T[i] << "\n";
                    // }
                    // else if ((h_temp[i] < 1.14*h_temp[0]) && (h_temp[i] > 1.09*h_temp[0]))
                    // {
                    //     output12 << x[i] << "," << y[i] << "," << z[i] << "," << T[i] << "\n";
                    // }
                    // else if ((h_temp[i] < 1.28*h_temp[0]) && (h_temp[i] > 1.22*h_temp[0]))
                    // {
                    //     output13 << x[i] << "," << y[i] << "," << z[i] << "," << T[i] << "\n";
                    // }
                }
            }
        }
    }

    auto end_calculation = chrono::high_resolution_clock::now();

    int num_not_dummy = 0;

    for (int i = 0; i < num_particle; i++)
    {
        if (x[i] >= x_left && x[i] <= x_right && y[i] >= y_bottom && y[i] <= y_top && z[i] >= z_back && z[i] <= z_front)
        {
            num_not_dummy += 1;
        }
    }

    double calc_time_ms = std::chrono::duration_cast <std::chrono::milliseconds> (end_calculation-start_time).count();

    printf("Number of Particle            : %d\n", num_particle);
    printf("Not Dummy Particle            : %d\n\n", num_not_dummy);
    printf("Calculation Time              : %f second\n", calc_time_ms/1000);

    ofstream output;
    output.open("output/3D/result_" + to_string(option) + "/Summary.csv");
    
    output  << "Number of Particle," << num_particle <<"\n"
            << "Not Dummy Particle," << num_not_dummy <<"\n"
            << "Calculation Time," << calc_time_ms/1000 << "\n"; 
}