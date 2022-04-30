#define _USE_MATH_DEFINES
#include <cmath>
#include <vector>
#include <chrono>
#include <fstream>
#include <iostream>

#include "src/LSMPS/LSMPS.hpp"
#include "src/Thermal Flux Method/Heat_Flux.hpp"
#include "src/Thermal Flux Method/Surface_Vector_S.hpp"
#include "src/Thermal Flux Method/Calculate_Laplacian.hpp"
#include "src/Time Integration/time_integration.hpp"
#include "src/External Heat Flux/External_Heat_Flux.hpp"

#include "src/Particle/Initialize_Particle.hpp"
#include "src/Particle/Particle_Splitting.hpp"
#include "src/Particle/Multires_Movement.hpp"
#include "src/Particle/Packing_Ratio.hpp"

#include "src/Neighbor Search/Brute_Force.hpp"
#include "src/Neighbor Search/Spatial_Hash.hpp"

using namespace std;

int main()
{
    auto start_time = chrono::high_resolution_clock::now();

    // Solid Domain
    double x0 = 0;
    double x1 = 0.1;
    double y0 = 0;
    double y1 = 0.05;
    
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
    int    nx  = 100;
    int    ny  = 50;

    // Other Parameters
    double Re      = 2.1;
    double n_dummy = 3;
    double t       = 0;
    double dt      = 1e-3;
    
    // Calculation
    double x_length = x1-x0;
    double y_length = y1-y0;

    double dx = x_length*eax/nx;
    double dy = y_length*eay/ny;

    x0 = x0 - dx/eax/2;
    x1 = x1 + dx/eax/2;
    y0 = y0 - dy/eay/2;
    y1 = y1 + dy/eay/2;

    nx = nx+1;
    ny = ny+1;

    double total_flux = q*y_length; // W

    // Create Variables Vector
    int num_particle;

    vector<double> x, y;
    vector<double> xw, yw;
    vector<double> h_temp, hx, hy;
    vector<double> QQ, kk, Cp, Rho;
    vector<double> hhx, hhy;
    vector<double> qx, qy;

    vector<int> split_index;
    vector<int> not_moving;
    vector<int> is_dummy;

    vector<vector<int>>    neighbor;
    vector<vector<double>> weight_data;

    int ncell_x, ncell_y, ncell;
    vector<vector<vector<int>>> hash_table;
    vector<int> gridpos_x, gridpos_y;

    int loop_count, iter;
    vector<vector<double>> Ri_a, Ni_a;
    vector<double> ci, delta_x, delta_y;

    vector<double> Temp;
    
    double area_temp;
    vector<int> flux;

    vector<int> top, bottom, left, right, right1, right2;

    // Initialize Particle
    initialize_particle(x0*eax,x1*eax,y0*eay,y1*eay,nx,ny,n_dummy,not_moving,xw,yw,h_temp);
    num_particle = xw.size();

    // Particle Movement
    auto start_movement = chrono::high_resolution_clock::now();

    // for (int i = 0; i < num_particle; i++)
    // {
    //     if (yw[i] > 0.035*eay || xw[i] > 0.07*eax)
    //     {
    //         split_index.push_back(i);
    //     }
    // }
    // Particle_Splitting(split_index,not_moving,xw,yw,h_temp);
    // num_particle = xw.size();
    // split_index.clear();

    // for (int i = 0; i < num_particle; i++)
    // {
    //     if (yw[i] > 0.025*eay || xw[i] > 0.05*eax)
    //     {
    //         split_index.push_back(i);
    //     }
    // }
    // Particle_Splitting(split_index,not_moving,xw,yw,h_temp);
    // num_particle = xw.size();
    // split_index.clear();

    calc_Ri_a(h_temp,Re,7,Ri_a);

    loop_count = 0;
    iter       = 0;

    while (true)
    {
        if(loop_count == iter)
        {
            break;
        }

        loop_count++;        
        hash_grid(xw,yw,h_temp[0]*Re,ncell_x,ncell_y,ncell,hash_table,gridpos_x,gridpos_y);
        spatial_hash_neighbor_2(xw,yw,h_temp,1,Re,ncell_x,ncell_y,gridpos_x,gridpos_y,hash_table,neighbor,weight_data);  
        calc_Ni(xw,yw,h_temp,neighbor,weight_data,Re,Ri_a,Ni_a);
        calc_ci(xw,Ri_a,Ni_a,ci);
        calc_DeltaX(xw,yw,h_temp,Re,0.1,neighbor,weight_data,ci,delta_x,delta_y);

        # pragma omp parallel for
        for (int i = 0; i < num_particle; i++)
        {
            if (not_moving[i] == 1)
            {
                continue;
            }

            xw[i] = xw[i] + delta_x[i];
            yw[i] = yw[i] + delta_y[i];
        }
    }

    auto end_movement = chrono::high_resolution_clock::now();
    
    double movement_ms = chrono::duration_cast <chrono::milliseconds> (end_movement-start_movement).count();
    printf("Particle Movement Time      : %f second\n\n", movement_ms/1000);

    // Assign Variables Value
    x.resize(num_particle);
    y.resize(num_particle);

    hx.resize(num_particle);
    hy.resize(num_particle);
    
    qx.resize(num_particle);
    qy.resize(num_particle);

    kk.resize(num_particle);
    Cp.resize(num_particle);
    Rho.resize(num_particle);
    QQ.resize(num_particle);
    
    hhx.resize(num_particle);
    hhy.resize(num_particle);

    Temp.resize(num_particle);

    is_dummy.resize(num_particle);

    top.resize(num_particle);
    bottom.resize(num_particle);
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

        kk[i]   = k;
        Temp[i] = T/*+(0.099-x[i])*5000*/;
        Cp[i]   = cp;
        Rho[i]  = rho;

        hx[i] = h_temp[i]/eax;
        hy[i] = h_temp[i]/eay;
        
        if ((x[i] < x0) || (x[i] > x1+2*dx) || (y[i] < y0) || (y[i] > y1))
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

        if (x[i] < x0 + 2*dx){left[i] = 1;}
        else if ((x[i] > x1 - dx)&&(x[i] < x1)){right[i] = 1;}
        else if ((x[i] > x1 - 2*dx)&&(x[i] < x1-dx)){right1[i] = 1;}
        else if ((x[i] > x1 - 3*dx)&&(x[i] < x1-2*dx)){right2[i] = 1;}

        if (y[i] < y0 + dx){bottom[i] = 1;}
        else if (y[i] > y1 - 2*dx){top[i] = 1;}
    }

    

    # pragma omp parallel for
    for (int i = 0; i < num_particle; i++)
    {
        if (is_dummy[i] == 1){continue;}

        if (left[i] == 1)
        {
            double multiplier = pow(1-(x[i]-0.5*dx)/2*dx,1/2.);
            area_temp += hy[i]*multiplier;
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

            for (int j = 0; j < num_particle; j++)
            {
                if (is_dummy[j] == 1){continue;}

                if (y[j] == y_temp)
                {
                    if ((x[j] > x1) && (x[j] < x1+dx))
                    {
                        real.push_back(i);
                        mirror.push_back(j);
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

            for (int j = 0; j < num_particle; j++)
            {
                if (is_dummy[j] == 1){continue;}

                if (y[j] == y_temp)
                {
                    if ((x[j] > x1+dx) && (x[j] < x1+2*dx))
                    {
                        real.push_back(i);
                        mirror.push_back(j);
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
    
    Re = 2.1;
    hash_grid(xw,yw,h_temp[0]*Re,ncell_x,ncell_y,ncell,hash_table,gridpos_x,gridpos_y);   
    spatial_hash_neighbor(xw,yw,h_temp[0]*Re,ncell_x,ncell_y,gridpos_x,gridpos_y,hash_table,neighbor,weight_data);

    auto end_neighbor_search = chrono::high_resolution_clock::now();

    double neighbor_time_ms = std::chrono::duration_cast <std::chrono::milliseconds> (end_neighbor_search-start_neighbor_search).count();
    printf("Neighbor Search Time        : %f second\n\n", neighbor_time_ms/1000);

    // LSMPS, Sij, Sij Star
    vector<vector<vector<double>>> LSMPS_eta(num_particle, vector<vector<double>>(9));
    calc_LSMPS_eta_2(LSMPS_eta, x, y, hx, hy, neighbor, weight_data);
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

    string name1 = "output/Benchmark 1/result/output_1_0.csv";
    string name2 = "output/Benchmark 1/result/output_2_0.csv";
    string name3 = "output/Benchmark 1/result/output_3_0.csv";

    ofstream output1, output2, output3;

    output1.open(name1);
    // output2.open(name2);
    // output3.open(name3);

    output1 << "x" << "," << "y" << "," << "LSMPS_Conserved\n";
    // output2 << "x" << "," << "y" << "," << "LSMPS_Conserved\n";
    // output3 << "x" << "," << "y" << "," << "LSMPS_Conserved\n";

    for (int i = 0; i < num_particle; i++)
    {
        if (x[i] >= x0 && x[i] <= x1 && y[i] >= y0 && y[i] <= y1)
        {
            if (h_temp[i] == dx)
            {
                output1 << xw[i] << "," << yw[i] << "," << Temp[i] << "\n";
            }
            // else if (h_temp[i]>0.7*dx)
            // {
            //     output2 << xw[i] << "," << yw[i] << "," << Temp[i] << "\n";
            // }
            // else
            // {
            //     output3 << xw[i] << "," << yw[i] << "," << Temp[i] << "\n";
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

        vector<double> kdeltaT_x, kdeltaT_y;
        calc_kdeltaT(x, y, kk, Temp, is_dummy, neighbor, LSMPS_eta, kdeltaT_x, kdeltaT_y);

        vector<vector<double>> kdeltaT_ij;
        calc_MPS_like_value(x, y, kk, Temp, is_dummy, neighbor, kdeltaT_ij);

        vector<vector<double>> alpha_ij, kdeltaT_ij_x, kdeltaT_ij_y, kdeltaT_ij_z;
        calc_hybrid(x, y, neighbor, alpha_ij, kdeltaT_x, kdeltaT_y, is_dummy, kdeltaT_ij, kdeltaT_ij_x, kdeltaT_ij_y);

        vector<double> heat_flux;
        External_Heat_Flux(qx, qy, hx, hy, heat_flux);
        Internal_Heat_flux(QQ,hx,hy,heat_flux);
        
        vector<double> T_conv(num_particle);
        # pragma omp parallel for
        for (int i = 0; i < num_particle; i++)
        {
            T_conv[i] = Temp[i] - 25;
        }
        Convection_Heating(hhx,hhy,T_conv,hx,hy,heat_flux);

        vector<double> dTdt;
        calc_dTdt(Cp, Rho, hx, hy, neighbor, kdeltaT_ij_x, kdeltaT_ij_y, Sij_Star_x, Sij_Star_y, Bi_x, Bi_y, kdeltaT_x, kdeltaT_y, heat_flux, is_dummy, dTdt);
        
        vector<double> T_Temp;
        T_Temp = Temp;
        time_integration(dt, Temp, dTdt);

        # pragma omp parallel for
        for (int i = 0; i < num_mirror; i++)
        {
            Temp[mirror[i]] = 2*T-Temp[real[i]];
        }

        // for (int i = 0; i < num_particle; i++)
        // {
        //     if (right[i] == 1)
        //     {
        //         double x_temp = x[i];
        //         double y_temp = y[i];
        //         double Temp_temp = Temp[i]-T;

        //         # pragma omp parallel for
        //         for (int j = 0; j < num_particle; j++)
        //         {
        //             if ((y[j] == y_temp))
        //             {
        //                 Temp[j] += -Temp_temp; 
        //             }
        //         }
        //     }
        // }
        
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

            output4 << "x" << "," << "y" << "," << "LSMPS_Conserved\n";
            // output5 << "x" << "," << "y" << "," << "LSMPS_Conserved\n";
            // output6 << "x" << "," << "y" << "," << "LSMPS_Conserved\n";

            for (int i = 0; i < num_particle; i++)
            {
                if (x[i] >= x0 && x[i] <= x1 && y[i] >= y0 && y[i] <= y1)
                {
                    if (h_temp[i] == dx)
                    {
                        output4 << xw[i] << "," << yw[i] << "," << Temp[i] << "\n";
                    }
                    // else if (h_temp[i]>0.7*dx)
                    // {
                    //     output5 << xw[i] << "," << yw[i] << "," << Temp[i] << "\n";
                    // }
                    // else
                    // {
                    //     output6 << xw[i] << "," << yw[i] << "," << Temp[i] << "\n";
                    // }
                }
            }
        }

        if (T_diff < 5e-3){break;}
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