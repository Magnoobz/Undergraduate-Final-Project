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

    // Solid Domain
    double x0 = 0;
    double x1 = 0.1;
    double y0 = 0;
    double y1 = 0.05;
    double z0 = 0.0;
    double z1 = 0.002;
    
    // Material Properties
    double k = 0.4;         // Thermal Conductivity (W/mC)
    double Q = 1.353e5;     // Volume heating (W/m3)
    double q = 3500;        // Heat Flux (W/m2)
    double T = 25;          // Initial Temperature (C)
    double h = 60;          // Convective Parameter (W/m2C)
    
    double rho = 1200;      // density (kg/m3)
    double cp  = 10;        // heat capacity (J/kgK)

    // Dommain Discretization
    double eax = 1;
    double eay = 1;
    double eaz = 1;
    int    nx  = 50;
    int    ny  = 25;
    int    nz  = 1;

    // Other Parameters
    double Re      = 2.4;
    double n_dummy = 3;
    double t       = 0;
    double dt      = 4e-3;
    
    // Calculation
    double x_length = x1-x0;
    double y_length = y1-y0;
    double z_length = z1-z0;

    double dx = x_length*eax/nx;

    double total_flux = q*y_length*z_length; // W

    // Create Variables Vector
    int num_particle;

    vector<double> x, y, z;
    vector<double> xw, yw, zw;
    vector<double> h_temp, hx, hy, hz;
    vector<double> QQ, kk, Cp, Rho;
    vector<double> hhx, hhy, hhz;
    vector<double> qx, qy, qz;

    vector<int> split_index;
    vector<int> not_moving;
    vector<int> is_dummy;

    vector<vector<int>>    neighbor;
    vector<vector<double>> weight_data;

    int ncell_x, ncell_y, ncell_z, ncell;
    vector<vector<vector<vector<int>>>> hash_table;
    vector<int> gridpos_x, gridpos_y, gridpos_z;

    int loop_count, iter;
    vector<vector<double>> Ri_a, Ni_a;
    vector<double> ci, delta_x, delta_y, delta_z;

    vector<double> Temp;
    
    double area_temp;
    vector<int> flux;

    vector<int> top, bottom, left, right, right1, right2, front, back;

    loop_count = 0;
    iter       = 50000;

    // Initialize Particle
    initialize_particle_3D(x0*eax,x1*eax,y0*eay,y1*eay,z0*eaz,z1*eaz,nx,ny,nz,n_dummy,not_moving,xw,yw,zw,h_temp);
    num_particle = xw.size();

    // Particle Movement
    auto start_movement = chrono::high_resolution_clock::now();

    // for (int i = 0; i < num_particle; i++)
    // {
    //     if (xw[i] < 0.03*eax || xw[i] > 0.07*eax)
    //     {
    //         split_index.push_back(i);
    //     }
    // }
    // Particle_Splitting_3D(split_index,not_moving,xw,yw,zw,h_temp);
    // num_particle = xw.size();
    // split_index.clear();

    calc_Ri_a_3D(h_temp,Re,7,Ri_a);

    loop_count = 0;
    iter       = 0;

    while (true)
    {
        if(loop_count == iter)
        {
            break;
        }

        loop_count++;        
        hash_grid_3D(xw,yw,zw,h_temp[0]*Re,ncell_x,ncell_y,ncell_z,ncell,hash_table,gridpos_x,gridpos_y,gridpos_z);
        spatial_hash_neighbor_3D_2(xw,yw,zw,h_temp,1,1,Re,ncell_x,ncell_y,ncell_z,gridpos_x,gridpos_y,gridpos_z,hash_table,neighbor,weight_data);  
        calc_Ni_3D(xw,yw,zw,h_temp,neighbor,weight_data,Re,Ri_a,Ni_a);
        calc_ci_3D(xw,Ri_a,Ni_a,ci);
        calc_DeltaX_3D(xw,yw,zw,h_temp,Re,0.1,neighbor,weight_data,ci,delta_x,delta_y,delta_z);

        # pragma omp parallel for
        for (int i = 0; i < num_particle; i++)
        {
            if (not_moving[i] == 1)
            {
                continue;
            }

            xw[i] = xw[i] + delta_x[i];
            yw[i] = yw[i] + delta_y[i];
            zw[i] = zw[i] + delta_z[i];
        }
    }

    auto end_movement = chrono::high_resolution_clock::now();
    
    double movement_ms = chrono::duration_cast <chrono::milliseconds> (end_movement-start_movement).count();
    printf("Particle Movement Time      : %f second\n\n", movement_ms/1000);

    // Assign Variables Value
    x.resize(num_particle);
    y.resize(num_particle);
    z.resize(num_particle);

    hx.resize(num_particle);
    hy.resize(num_particle);
    hz.resize(num_particle);
    
    qx.resize(num_particle);
    qy.resize(num_particle);
    qz.resize(num_particle);

    kk.resize(num_particle);
    Cp.resize(num_particle);
    Rho.resize(num_particle);
    QQ.resize(num_particle);
    
    hhx.resize(num_particle);
    hhy.resize(num_particle);
    hhz.resize(num_particle);

    Temp.resize(num_particle);

    is_dummy.resize(num_particle);

    top.resize(num_particle);
    bottom.resize(num_particle);
    front.resize(num_particle);
    back.resize(num_particle);
    left.resize(num_particle);
    right.resize(num_particle);
    right1.resize(num_particle);
    right2.resize(num_particle);

    # pragma omp parallel for
    for (int i = 0; i < num_particle; i++)
    {
        QQ[i] = Q;
        
        x[i] = xw[i]/eax;
        y[i] = yw[i]/eay;
        z[i] = zw[i]/eaz;

        kk[i]   = k;
        Temp[i] = T/*+(0.099-x[i])*5000*/;
        Cp[i]   = cp;
        Rho[i]  = rho;

        hx[i] = h_temp[i]/eax;
        hy[i] = h_temp[i]/eay;
        hz[i] = h_temp[i]/eaz;
        
        if ((x[i] < x0) || (x[i] > x1+0.004) || (y[i] < y0) || (y[i] > y1) || (z[i] < z0) || (z[i] > z1))
        {
            is_dummy[i] = 1;
        }
        else
        {
            is_dummy[i] = 0;
        }
    }

    # pragma omp parallel for
    for (int i = 0; i < num_particle; i++)
    {
        if (is_dummy[i] == 1){continue;}

        if (x[i] < x0 + 0.004){left[i] = 1;}
        else if ((x[i] > x1 - 0.002)&&(x[i] < x1)){right[i] = 1;}
        else if ((x[i] > x1 - 0.004)&&(x[i] < x1-0.0002)){right1[i] = 1;}
        else if ((x[i] > x1 - 0.006)&&(x[i] < x1-0.0004)){right2[i] = 1;}

        if (y[i] < y0 + 0.002){bottom[i] = 1;}
        else if (y[i] > y1 - 0.004){top[i] = 1;}

        if (z[i] < z0 + 0.002){back[i] = 1;}
        else if (z[i] > z1 - 0.002){front[i] = 1;}
    }

    

    # pragma omp parallel for
    for (int i = 0; i < num_particle; i++)
    {
        if (is_dummy[i] == 1){continue;}

        if (left[i] == 1)
        {
            double multiplier = pow(1-(x[i]-0.001)/0.004,1);
            area_temp += hy[i]*hz[i]*multiplier;
            qx[i] = multiplier;
            flux.push_back(i);
        }
        if ((top[i] == 1) || (top[i] == 1))
        {
            if ((right[i] == 1) /*|| (left[i] == 1)*/){continue;}
            
            hhy[i] = -h/2;
        }

        if (right[i] == 1)
        {
            QQ[i] = 0;
        }

        // if (right[i] == 1)
        // {
        //     is_dummy[i] = 1;
        // }
    }

    q = total_flux/area_temp;

    # pragma omp parallel for
    for (int i : flux)
    {
        qx[i] = qx[i]*q;
    }

    vector<int> real, mirror;

    for (int i = 0; i < num_particle; i++)
    {
        if (right1[i] == 1)
        {
            double x_temp = x[i];
            double y_temp = y[i];
            double z_temp = z[i];

            for (int j = 0; j < num_particle; j++)
            {
                if (is_dummy[j] == 1){continue;}

                if (y[j] == y_temp)
                {
                    if(z[j] == z_temp)
                    {
                        if ((x[j] > x1) && (x[j] < x1+0.002))
                        {
                            real.push_back(i);
                            mirror.push_back(j);
                        }
                    }
                }
            }
        }
    }

    for (int i = 0; i < num_particle; i++)
    {
        if (right2[i] == 1)
        {
            double x_temp = x[i];
            double y_temp = y[i];
            double z_temp = z[i];

            for (int j = 0; j < num_particle; j++)
            {
                if (is_dummy[j] == 1){continue;}

                if (y[j] == y_temp)
                {
                    if(z[j] == z_temp)
                    {
                        if ((x[j] > x1+0.002) && (x[j] < x1+0.004))
                        {
                            real.push_back(i);
                            mirror.push_back(j);
                        }
                    }
                }
            }
        }
    }

    int num_mirror=mirror.size();
    # pragma omp parallel for
    for (int i = 0; i < num_mirror; i++)
    {
        Temp[mirror[i]] = 2*T-Temp[real[i]];
    }

    // Neighbor Search
    auto start_neighbor_search = chrono::high_resolution_clock::now();
    
    Re = 2.9;
    hash_grid_3D(xw,yw,zw,h_temp[0]*Re,ncell_x,ncell_y,ncell_z,ncell,hash_table,gridpos_x,gridpos_y,gridpos_z);   
    spatial_hash_neighbor_3D(xw,yw,zw,h_temp[0]*Re,ncell_x,ncell_y,ncell_z,gridpos_x,gridpos_y,gridpos_z, hash_table, neighbor, weight_data);

    auto end_neighbor_search = chrono::high_resolution_clock::now();

    double neighbor_time_ms = std::chrono::duration_cast <std::chrono::milliseconds> (end_neighbor_search-start_neighbor_search).count();
    printf("Neighbor Search Time        : %f second\n\n", neighbor_time_ms/1000);

    // LSMPS, Sij, Sij Star
    vector<vector<vector<double>>> LSMPS_eta(num_particle, vector<vector<double>>(9));
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

    string name1 = "output/Benchmark 1/result/output_1_0.csv";
    string name2 = "output/Benchmark 1/result/output_2_0.csv";
    string name3 = "output/Benchmark 1/result/output_3_0.csv";

    ofstream output1, output2, output3;

    output1.open(name1);
    // output2.open(name2);
    // output3.open(name3);

    output1 << "x" << "," << "y" << "," << "z" << "," << "LSMPS_Conserved\n";
    // output2 << "x" << "," << "y" << "," << "z" << "," << "LSMPS_Conserved\n";
    // output3 << "x" << "," << "y" << "," << "z" << "," << "LSMPS_Conserved\n";

    for (int i = 0; i < num_particle; i++)
    {
        if (x[i] >= x0 && x[i] <= x1+0.004 && y[i] >= y0 && y[i] <= y1 && z[i] >= z0 && z[i] <= z1)
        {
            if (h_temp[i] == dx)
            {
                output1 << xw[i] << "," << yw[i] << "," << zw[i] << "," << Temp[i] << "\n";
            }
            // else if (h_temp[i]>0.7*dx)
            // {
            //     output2 << xw[i] << "," << yw[i] << "," << zw[i] << "," << Temp[i] << "\n";
            // }
            // else
            // {
            //     output3 << xw[i] << "," << yw[i] << "," << zw[i] << "," << Temp[i] << "\n";
            // }
        }
    }
    
    

    auto start_loop_segment = chrono::high_resolution_clock::now();
    auto end_loop_segment = chrono::high_resolution_clock::now();

    while (/*loop_count < iter*/ true)
    {
        if (loop_count % 10 == 0)
        {
            start_loop_segment = chrono::high_resolution_clock::now();
        }

        vector<double> kdeltaT_x, kdeltaT_y, kdeltaT_z;
        calc_kdeltaT_3D(x, y, z, kk, Temp, is_dummy, neighbor, LSMPS_eta, kdeltaT_x, kdeltaT_y, kdeltaT_z);

        vector<vector<double>> kdeltaT_ij;
        calc_MPS_like_value_3D(x, y, z, kk, Temp, is_dummy, neighbor, kdeltaT_ij);

        vector<vector<double>> alpha_ij, kdeltaT_ij_x, kdeltaT_ij_y, kdeltaT_ij_z;
        calc_hybrid_3D(x, y, z, neighbor, alpha_ij, kdeltaT_x, kdeltaT_y, kdeltaT_z, is_dummy, kdeltaT_ij, kdeltaT_ij_x, kdeltaT_ij_y, kdeltaT_ij_z);

        vector<double> heat_flux;
        External_Heat_Flux_3D(qx, qy, qz, hx, hy, hz, heat_flux);
        Internal_Heat_flux_3D(QQ,hx,hy,hz,heat_flux);
        
        vector<double> T_conv(num_particle);
        # pragma omp parallel for
        for (int i = 0; i < num_particle; i++)
        {
            T_conv[i] = Temp[i] - 25;
        }
        Convection_Heating_3D(hhx,hhy,hhz,T_conv,hx,hy,hz,heat_flux);

        vector<double> dTdt;
        calc_dTdt_3D(Cp, Rho, hx, hy, hz, neighbor, kdeltaT_ij_x, kdeltaT_ij_y, kdeltaT_ij_z, Sij_Star_x, Sij_Star_y, Sij_Star_z, Bi_x, Bi_y, Bi_z, kdeltaT_x, kdeltaT_y, kdeltaT_z, heat_flux, is_dummy, dTdt);
        
        // for (int i = 0; i < num_particle; i++)
        // {
        //     if (right[i] == 1)
        //     {
        //         double x_temp = x[i];
        //         double y_temp = y[i];
        //         double z_temp = z[i];
        //         double dTdt_temp = dTdt[i];

        //         # pragma omp parallel for
        //         for (int j = 0; j < num_particle; j++)
        //         {
        //             if ((y[j] == y_temp) && (z[j] == z_temp))
        //             {
        //                 dTdt[j] += (-x_temp+x[j])*dTdt_temp; 
        //             }
        //         }
        //     }
        // }
        vector<double> T_Temp;
        T_Temp = Temp;
        time_integration(dt, Temp, dTdt);

        # pragma omp parallel for
        for (int i = 0; i < num_mirror; i++)
        {
            Temp[mirror[i]] = 2*T-Temp[real[i]];
        }
        
        double T_diff = 0;
        # pragma omp parallel for
        for (int i = 0; i < num_particle; i++)
        {
            if (is_dummy[i] == 1){continue;}

            T_diff += abs(Temp[i] - T_Temp[i]);
        }

        t+=dt;
        loop_count += 1;

        if (loop_count % 10 == 0)
        {
            end_loop_segment = chrono::high_resolution_clock::now();
            double loop_time_ms = std::chrono::duration_cast <std::chrono::milliseconds> (end_loop_segment-start_loop_segment).count();
            
            cout << loop_count <<":"<<"\t"<< t << "\t\t Segment time: " << loop_time_ms/1000 << "\tsecond \t\t" << T_diff << endl;
        }
        
        if (loop_count % 100 == 0)
        {
            string name1 = "output/Benchmark 1/result/output_1_" + to_string(loop_count) + ".csv";
            string name2 = "output/Benchmark 1/result/output_2_" + to_string(loop_count) + ".csv";
            string name3 = "output/Benchmark 1/result/output_3_" + to_string(loop_count) + ".csv";

            ofstream output4, output5, output6;

            output4.open(name1);
            // output5.open(name2);
            // output6.open(name3);

            output4 << "x" << "," << "y" << "," << "z" << "," << "LSMPS_Conserved\n";
            // output5 << "x" << "," << "y" << "," << "z" << "," << "LSMPS_Conserved\n";
            // output6 << "x" << "," << "y" << "," << "z" << "," << "LSMPS_Conserved\n";

            for (int i = 0; i < num_particle; i++)
            {
                if (x[i] >= x0 && x[i] <= x1+0.004 && y[i] >= y0 && y[i] <= y1 && z[i] >= z0 && z[i] <= z1)
                {
                    if (h_temp[i] == dx)
                    {
                        output4 << xw[i] << "," << yw[i] << "," << zw[i] << "," << Temp[i] << "\n";
                    }
                    // else if (h_temp[i]>0.7*dx)
                    // {
                    //     output5 << xw[i] << "," << yw[i] << "," << zw[i] << "," << Temp[i] << "\n";
                    // }
                    // else
                    // {
                    //     output6 << xw[i] << "," << yw[i] << "," << zw[i] << "," << Temp[i] << "\n";
                    // }
                }
            }
        }

        if (T_diff < 1e-2){break;}
    }
    
    auto end_time = chrono::high_resolution_clock::now();

    double calc_time_ms = std::chrono::duration_cast <std::chrono::milliseconds> (end_time-start_time).count();

    printf("\nNeighbor Search Time        : %f second\n", neighbor_time_ms/1000);
    printf("Calc Eta Time               : %f second\n", eta_time_ms/1000);
    printf("Sij Time                    : %f second\n", sij_time_ms/1000);
    printf("Bi Time                     : %f second\n", bi_time_ms/1000);
    printf("Sij Star Time               : %f second\n", sijstar_time_ms/1000);
    printf("Calculation Time            : %f second\n\n", calc_time_ms/1000);

    
    ofstream output7;
    output7.open("output/Benchmark 1/result/Summary.csv");
    
    output7  << "Number of Particle," << x.size() <<"\n"
            << "Movement Time," << movement_ms/1000 << "\n"
            << "Neighbor Search Time," << neighbor_time_ms/1000 << "\n"
            << "Calc Eta Time," << eta_time_ms/1000 << "\n"
            << "Sij Time," << sij_time_ms/1000 << "\n"
            << "Bi Time," << bi_time_ms/1000 << "\n"
            << "Sij Star Time," << sijstar_time_ms/1000 << "\n"
            << "Calculation Time," << calc_time_ms/1000 << "\n";
}