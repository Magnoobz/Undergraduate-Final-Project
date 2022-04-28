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

    vector<double> radius;
    radius.push_back(0.000505);
    radius.push_back(0.00072);
    radius.push_back(0.00102);
    radius.push_back(0.001435);
    vector<double> time;
    time.push_back(0.48);
    time.push_back(0.24);
    time.push_back(0.12);
    time.push_back(0.06);
    vector<double> dts;
    dts.push_back(2e-3);
    dts.push_back(1e-3);
    dts.push_back(5e-4);
    dts.push_back(2.5e-4);

    int option = 2;
    
    double radii    = radius[option-1];
    double t_end    = time[option-1];
    double dt       = dts[option-1];

    double energy   = 295;
    double sigma    = 1*radii/3.;

    double p_den;
    double T_init   = 24.85;
    double k_init   = 42.636;
    double cp_init  = 510;
    double rho_init = 7600;
    double T_liq    = 1500;
    double h_fus    = 2.5e5;
    double T_vap    = 3000;
    double h_vap    = 6.8e6;

    double x_left   = 0;
    double x_right  = 0.01;
    double y_bottom = 0;
    double y_top    = 0.01;
    double z_front  = 0.005;
    double z_back   = 0;

    double eax = 1;
    double eay = 1;
    double eaz = 1;

    int nx = 40;
    int ny = 40;
    int nz = 20;

    double dx = (x_right-x_left)/(nx);
    double dy = (y_top-y_bottom)/(ny);
    double dz = (z_front-z_back)/(nz);

    vector<double> x, y, z, hx, hy, hz;
    vector<double> k, cp, rho;
    vector<double> x_w, y_w, z_w;
    vector<double> q_x, q_y, q_z;
    vector<double> h_fusion, h_vaporization;
    
    vector<int> split_index;
    vector<int> not_moving;    
    vector<int> is_dummy;

    double n_dummy = 5;

    initialize_particle_3D(x_left, x_right, y_bottom, y_top, z_back, z_front, nx, ny, nz, n_dummy, not_moving, x, y, z, hx);
    int num_particle = x.size();
        
    vector<vector<int>> neighbor;
    vector<vector<double>> weight_data;

    int ncell_x, ncell_y, ncell_z, ncell;
    vector<vector<vector<vector<int>>>> hash_table;
    vector<int> gridpos_x, gridpos_y, gridpos_z;
    
    auto start_particle_movement = chrono::high_resolution_clock::now();
    
    double R_e = 2.4;

    for (int i = 0; i < num_particle; i++)
    {
        if (x[i] < 0.008 && x[i] > 0.002 && y[i] < 0.008 && y[i] > 0.002 && z[i] < 0.004)
        {
            split_index.push_back(i);
        }
    }
    Particle_Splitting_3D(split_index, not_moving, x, y, z, hx);
    num_particle = x.size();
    split_index.clear();

    vector<vector<double>> Ri_a;
    calc_Ri_a_3D(hx, R_e, 7, Ri_a);

    int loop_count = 0;
    int iter = 30;

    while (true)
    {
        if(loop_count == iter)
        {
            break;
        }
        
        hash_grid_3D(x, y, z, hx[0]*R_e*eay*eaz,ncell_x, ncell_y, ncell_z, ncell, hash_table, gridpos_x, gridpos_y, gridpos_z);
        spatial_hash_neighbor_3D_2(x,y,z,hx,eay,eaz,R_e,ncell_x, ncell_y, ncell_z,gridpos_x, gridpos_y, gridpos_z, hash_table, neighbor, weight_data);  

        loop_count++;

        vector<vector<double>> Ni_a;
        calc_Ni_3D(x, y, z, hx, neighbor, weight_data, R_e, Ri_a, Ni_a);

        vector<double> ci;
        calc_ci_3D(x, Ri_a, Ni_a, ci);

        vector<double> delta_x, delta_y, delta_z;
        calc_DeltaX_3D(x, y, z, hx, R_e, 0.1, neighbor, weight_data, ci, delta_x, delta_y, delta_z);

        for (int i = 0; i < num_particle; i++)
        {
            if (not_moving[i] == 1)
            {
                continue;
            }

            x[i] = x[i] + delta_x[i];
            y[i] = y[i] + delta_y[i];
            z[i] = z[i] + delta_z[i];
        }
    }

    for (int i = 0; i < num_particle; i++)
    {
        if (x[i] < 0.0065 && x[i] > 0.0035 && y[i] < 0.0065 && y[i] > 0.0035 && z[i] < 0.025)
        {
            split_index.push_back(i);
        }
    }
    Particle_Splitting_3D(split_index, not_moving, x, y, z, hx);
    num_particle = x.size();  
    split_index.clear();

    not_moving.clear();
    not_moving.resize(num_particle);

    #pragma omp parallel for
    for (int i = 0; i < num_particle; i++)
    {
        if (hx[i] == dx)
        {
            not_moving[i] = 1;
        }
        else
        {
            not_moving[i] = 0;
        }
    }

    calc_Ri_a_3D(hx, R_e, 7, Ri_a);   

    loop_count = 0;
    iter = 200;

    while (true)
    {
        if (loop_count = iter)
        {
            break;
        }

        hash_grid_3D(x, y, z, hx[0]*R_e*eay*eaz,ncell_x, ncell_y, ncell_z, ncell, hash_table, gridpos_x, gridpos_y, gridpos_z);
        spatial_hash_neighbor_3D_2(x,y,z,hx,eay,eaz,R_e,ncell_x, ncell_y, ncell_z,gridpos_x, gridpos_y, gridpos_z, hash_table, neighbor, weight_data);  

        loop_count++;

        vector<vector<double>> Ni_a;
        calc_Ni_3D(x, y, z, hx, neighbor, weight_data, R_e, Ri_a, Ni_a);

        vector<double> ci;
        calc_ci_3D(x, Ri_a, Ni_a, ci);

        vector<double> delta_x, delta_y, delta_z;
        calc_DeltaX_3D(x, y, z, hx, R_e, 0.1, neighbor, weight_data, ci, delta_x, delta_y, delta_z);

        for (int i = 0; i < num_particle; i++)
        {
            if (not_moving[i] == 1)
            {
                continue;
            }

            x[i] = x[i] + delta_x[i];
            y[i] = y[i] + delta_y[i];
            z[i] = z[i] + delta_z[i];
        }
    }

    auto end_particle_movement = chrono::high_resolution_clock::now();

    double movement_ms = std::chrono::duration_cast <std::chrono::milliseconds> (end_particle_movement-start_particle_movement).count();
    printf("Particle Movement Time      : %f second\n\n", movement_ms/1000);
    
    hy = hx;
    hz = hx;
    
    x_w.resize(num_particle);
    y_w.resize(num_particle);
    z_w.resize(num_particle);
    
    q_x.resize(num_particle);
    q_y.resize(num_particle);
    q_z.resize(num_particle);
    
    k.resize(num_particle);
    cp.resize(num_particle);
    rho.resize(num_particle);
    h_fusion.resize(num_particle);
    h_vaporization.resize(num_particle);
    
    is_dummy.resize(num_particle);
    
    vector<double> T(num_particle);
    
    double area_temp = 0;
    vector<int> flux;
    
    # pragma omp parallel for
    for (int i = 0; i < num_particle; i++)
    {
        if ((x[i] < x_left) || (x[i] > x_right) || (y[i] < y_bottom) || (y[i] > y_top) || (z[i] < z_back) || (z[i] > z_front))
        {
            is_dummy[i] = 1;
        }
        else
        {
            is_dummy[i] = 0;
        }

        k[i] = (k_init);
        cp[i] = (cp_init);
        rho[i] = (rho_init);
        h_fusion[i] = (h_fus);
        h_vaporization[i] = (h_vap);

        T[i] = T_init;
        
        x_w [i] = x[i]*eay*eaz;
        y_w [i] = y[i]*eax*eaz;
        z_w [i] = z[i]*eax*eay;

        if ((z[i]<z_back+1*dz) && (z[i]>z_back) && (pow(x[i]-0.005,2)+pow(y[i]-0.005,2)+0*pow(z[i],2) < pow(radii,2)))
        {
            double gauss_dist = 1/(2*M_PI*pow(sigma,2))*exp(-(pow(x[i]-0.005,2)+pow(y[i]-0.005,2))/(2*sigma*sigma))*pow(1-(z[i]-z_back)/(1*dz),5);
            area_temp += hx[i]*hy[i]*gauss_dist;

            q_z[i] = gauss_dist;
            flux.push_back(i);
        }
    }
    
    p_den = energy/(t_end*area_temp);
    
    # pragma omp parallel for
    for (int i : flux)
    {
        q_z[i] = q_z[i]*p_den;
    }
    
    auto start_neighbor_search = chrono::high_resolution_clock::now();
    
    R_e = 2;
    hash_grid_3D(x_w,y_w,z_w, hx[0]*R_e*eay*eaz,ncell_x, ncell_y, ncell_z, ncell, hash_table, gridpos_x, gridpos_y, gridpos_z);   
    spatial_hash_neighbor_3D(x_w,y_w,z_w,hx[0]*R_e*eay*eaz,ncell_x, ncell_y, ncell_z,gridpos_x, gridpos_y, gridpos_z, hash_table, neighbor, weight_data);

    
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

    string name11 = "output/Thick Block/Laser "+ to_string(option) +"/result/output_q_1.csv";
    string name12 = "output/Thick Block/Laser "+ to_string(option) +"/result/output_q_2.csv";
    string name13 = "output/Thick Block/Laser "+ to_string(option) +"/result/output_q_3.csv";

    ofstream output11, output12, output13;

    output11.open(name11);
    output12.open(name12);
    output13.open(name13);

    output11 << "x" << "," << "y" << "," << "z" << "," << "q\n";
    output12 << "x" << "," << "y" << "," << "z" << "," << "q\n";
    output13 << "x" << "," << "y" << "," << "z" << "," << "q\n";

    for (int i = 0; i < num_particle; i++)
    {
        if (x[i] >= x_left && x[i] <= x_right && y[i] >= y_bottom && y[i] <= y_top && z[i] >= z_back && z[i] <= z_front)
        {
            if (hx[i] == dx)
            {
                output11 << x_w[i] << "," << y_w[i] << "," << z_w[i] << "," << q_z[i] << "\n";
            }
            else if (hx[i] > 0.7*dx)
            {
                output12 << x_w[i] << "," << y_w[i] << "," << z_w[i] << "," << q_z[i] << "\n";
            }
            else
            {
                output13 << x_w[i] << "," << y_w[i] << "," << z_w[i] << "," << q_z[i] << "\n";
            }
        }
    }
    
    string name1 = "output/Thick Block/Laser "+ to_string(option) +"/result/output_1_" + to_string(count) + ".csv";
    string name2 = "output/Thick Block/Laser "+ to_string(option) +"/result/output_2_" + to_string(count) + ".csv";
    string name3 = "output/Thick Block/Laser "+ to_string(option) +"/result/output_3_" + to_string(count) + ".csv";

    ofstream output1, output2, output3;

    output1.open(name1);
    output2.open(name2);
    output3.open(name3);

    output1 << "x" << "," << "y" << "," << "z" << "," << "LSMPS_Conserved\n";
    output2 << "x" << "," << "y" << "," << "z" << "," << "LSMPS_Conserved\n";
    output3 << "x" << "," << "y" << "," << "z" << "," << "LSMPS_Conserved\n";

    for (int i = 0; i < num_particle; i++)
    {
        if (x[i] >= x_left && x[i] <= x_right && y[i] >= y_bottom && y[i] <= y_top && z[i] >= z_back && z[i] <= z_front)
        {
            if (hx[i] == dx)
            {
                output1 << x_w[i] << "," << y_w[i] << "," << z_w[i] << "," << T[i] << "\n";
            }
            else if (hx[i]>0.7*dx)
            {
                output2 << x_w[i] << "," << y_w[i] << "," << z_w[i] << "," << T[i] << "\n";
            }
            else
            {
                output3 << x_w[i] << "," << y_w[i] << "," << z_w[i] << "," << T[i] << "\n";
            }
        }
    }

    auto start_loop_segment = chrono::high_resolution_clock::now();
    auto end_loop_segment = chrono::high_resolution_clock::now();

    while (t < t_end /*count < 1*/)
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

        # pragma omp parallel for
        for (int i = 0; i < num_particle; i++)
        {
            if (T[i] > T_liq && h_fusion[i] > 0)
            {
                double heat_energy = (T[i]-T_liq)*cp[i];

                if (heat_energy > h_fusion[i])
                {
                    double dT = h_fusion[i]/cp[i];
                    T[i] = T[i] - dT;
                    h_fusion[i] = 0;
                }
                else
                {
                    T[i] = T_liq;
                    h_fusion[i] = h_fusion[i] - heat_energy;
                }
            }
            else if (T[i] > T_vap && h_vaporization[i] > 0)
            {
                double heat_energy = (T[i]-T_vap)*cp[i];

                if (heat_energy > h_vaporization[i])
                {
                    double dT = h_vaporization[i]/cp[i];
                    T[i] = T[i] - dT;
                    // T[i] = T_vap;
                    h_vaporization[i] = 0;
                }
                else
                {
                    T[i] = T_vap;
                    h_vaporization[i] = h_vaporization[i] - heat_energy;
                }
            }
            // else if (h_vaporization[i] == 0){T[i] = T_vap;}
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
            string name1 = "output/Thick Block/Laser "+ to_string(option) +"/result/output_1_" + to_string(count) + ".csv";
            string name2 = "output/Thick Block/Laser "+ to_string(option) +"/result/output_2_" + to_string(count) + ".csv";
            string name3 = "output/Thick Block/Laser "+ to_string(option) +"/result/output_3_" + to_string(count) + ".csv";

            ofstream output4, output5, output6;

            output4.open(name1);
            output5.open(name2);
            output6.open(name3);

            output4 << "x" << "," << "y" << "," << "z" << "," << "LSMPS_Conserved\n";
            output5 << "x" << "," << "y" << "," << "z" << "," << "LSMPS_Conserved\n";
            output6 << "x" << "," << "y" << "," << "z" << "," << "LSMPS_Conserved\n";

            for (int i = 0; i < num_particle; i++)
            {
                if (x[i] >= x_left && x[i] <= x_right && y[i] >= y_bottom && y[i] <= y_top && z[i] >= z_back && z[i] <= z_front)
                {
                    if (hx[i] == dx)
                    {
                        output4 << x_w[i] << "," << y_w[i] << "," << z_w[i] << "," << T[i] << "\n";
                    }
                    else if (hx[i]>0.7*dx)
                    {
                        output5 << x_w[i] << "," << y_w[i] << "," << z_w[i] << "," << T[i] << "\n";
                    }
                    else
                    {
                        output6 << x_w[i] << "," << y_w[i] << "," << z_w[i] << "," << T[i] << "\n";
                    }
                }
            }
        }
    }

    # pragma omp parallel for
    vector<double> melted;
    for (int i = 0; i < num_particle; i++)
    {
        if(T[i] > T_liq)
        {
            melted.push_back(z[i]);
        }
    }

    double depth = *max_element(melted.begin(),melted.end());

    auto end_calculation = chrono::high_resolution_clock::now();
  
    double calc_time_ms = std::chrono::duration_cast <std::chrono::milliseconds> (end_calculation-start_time).count();

    printf("\nNeighbor Search Time        : %f second\n", neighbor_time_ms/1000);
    printf("Calc Eta Time               : %f second\n", eta_time_ms/1000);
    printf("Sij Time                    : %f second\n", sij_time_ms/1000);
    printf("Bi Time                     : %f second\n", bi_time_ms/1000);
    printf("Sij Star Time               : %f second\n", sijstar_time_ms/1000);
    printf("Calculation Time            : %f second\n\n", calc_time_ms/1000);
    printf("Depth                       : %f m\n\n", depth);

    
    ofstream output7;
    output7.open("output/Thick Block/Laser "+ to_string(option) +"/result/Summary.csv");
    
    output7  << "Number of Particle," << x.size() <<"\n"
            << "Movement Time," << movement_ms/1000 << "\n"
            << "Neighbor Search Time," << neighbor_time_ms/1000 << "\n"
            << "Calc Eta Time," << eta_time_ms/1000 << "\n"
            << "Sij Time," << sij_time_ms/1000 << "\n"
            << "Bi Time," << bi_time_ms/1000 << "\n"
            << "Sij Star Time," << sijstar_time_ms/1000 << "\n"
            << "Calculation Time," << calc_time_ms/1000 << "\n\n"
            << "Depth," << depth << "\n";    
}