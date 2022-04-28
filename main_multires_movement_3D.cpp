# include <cmath>
# include <vector>
# include <fstream>
# include <iostream>
# include <chrono>

# include "src/Particle/Initialize_Particle_3D.hpp"
# include "src/Particle/Particle_Splitting_3D.hpp"
# include "src/Particle/Multires_Movement_3D.hpp"
# include "src/Particle/Packing_Ratio_3D.hpp"

# include "src/Neighbor Search/Brute_Force.hpp"
# include "src/Neighbor Search/Spatial_Hash_3D.hpp"

using namespace std;

int main()
{
    double x_left   = 0;
    double x_right  = 1;
    double y_bottom = 0;
    double y_top    = 1;
    double z_back   = 0;
    double z_front  = 1;

    double eax = 1;
    double eay = 1;
    double eaz = 1;

    int nx = 20;
    int ny = 20;
    int nz = 20;

    double dx = (x_right-x_left)/nx;
    double dy = (y_top-y_bottom)/ny;
    double dz = (z_front-z_back)/nz;

    double n_dummy = 5;

    vector<double> x, y, z, h;
    vector<int> split_index;
    vector<int> not_moving;
    vector<int> is_dummy;

    auto initialize_start = chrono::high_resolution_clock::now();

    initialize_particle_3D(x_left, x_right, y_bottom, y_top, z_back, z_front, nx, ny, nz, n_dummy, not_moving, x, y, z, h);

    auto split_index_start = chrono::high_resolution_clock::now();

    int num_particle = x.size();

    // for (int i = 0; i < num_particle; i++)
    // {
    //     if ((x[i] > 0.3) && (x[i] < 0.7) && (y[i] > 0.3) && (y[i] < 0.7))
    //     {
    //         split_index.push_back(i);
    //     }
    // }
    
    for (int i = 0; i < num_particle; i++)
    {
        if (pow(x[i]-0.5,2)+pow(y[i]-0.5,2)<pow(0.35,2))
        {
            split_index.push_back(i);
        }
    }

    auto splitting_start = chrono::high_resolution_clock::now();

    Particle_Splitting_3D(split_index, not_moving, x, y, z, h);

    num_particle = x.size();

    
    split_index.clear();
    for (int i = 0; i < num_particle; i++)
    {
        if (pow(x[i]-0.5,2)+pow(y[i]-0.5,2)<pow(0.2,2))
        {
            split_index.push_back(i);
        }
    }

    Particle_Splitting_3D(split_index, not_moving, x, y, z, h);

    num_particle = x.size();
    
    auto splitting_end = chrono::high_resolution_clock::now();    
    
    vector<vector<double>> weight_data;
    vector<vector<int>> neighbor;

    int ncell_x, ncell_y, ncell_z, ncell;
    vector<vector<vector<vector<int>>>> hash_table;
    vector<int> gridpos_x, gridpos_y, gridpos_z;

    double R_e = 2.4; 

    vector<vector<double>> Ri_a;
    calc_Ri_a_3D(h, R_e, 7, Ri_a);   

    int loop_count = 0;
    int iter = 200;

    while (loop_count < iter)
    {
        // brute_force_3D(x, y, z, h, 1, 1, neighbor, weight_data, R_e);
        hash_grid_3D(x, y, z, h[0]*R_e*eay*eaz,ncell_x, ncell_y, ncell_z, ncell, hash_table, gridpos_x, gridpos_y, gridpos_z);
        spatial_hash_neighbor_3D_2(x,y,z,h,eay,eaz,R_e,ncell_x, ncell_y, ncell_z,gridpos_x, gridpos_y, gridpos_z, hash_table, neighbor, weight_data);  

        loop_count++;

        vector<vector<double>> Ni_a;
        calc_Ni_3D(x, y, z, h, neighbor, weight_data, R_e, Ri_a, Ni_a);

        vector<double> ci;
        calc_ci_3D(x, Ri_a, Ni_a, ci);

        vector<double> delta_x, delta_y, delta_z;
        calc_DeltaX_3D(x, y, z, h, R_e, 0.1, neighbor, weight_data, ci, delta_x, delta_y, delta_z);

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

        ofstream output1, output2, output3;

        string name1 = "output/Tes Distribusi Partikel 3D/Tes Iterasi 4/big/Distribusi_" + to_string(nx) + "Partikel_" + to_string(loop_count) + "_iter.csv";
        string name2 = "output/Tes Distribusi Partikel 3D/Tes Iterasi 4/small/Distribusi_" + to_string(nx) + "Partikel_" + to_string(loop_count) + "_iter.csv";
        string name3 = "output/Tes Distribusi Partikel 3D/Tes Iterasi 4/smaller/Distribusi_" + to_string(nx) + "Partikel_" + to_string(loop_count) + "_iter.csv";

        output1.open(name1);
        output2.open(name2);
        output3.open(name3);

        output1 << "x" << "," << "y" << "," << "z" << "," << "h\n";
        output2 << "x" << "," << "y" << "," << "z" << "," << "h\n";
        output3 << "x" << "," << "y" << "," << "z" << "," << "h\n";

        for (int i = 0; i < x.size(); i++)
        {
            if ((x[i]<0.59) && (y[i]<0.59))
            {
                if (h[i] == dx)
                {
                    output1 << x[i] << "," << y[i] << "," << z[i] << "," << h[i] << "\n";
                }
                else if (h[i]>0.7*dx)
                {
                    output2 << x[i] << "," << y[i] << "," << z[i] << "," << h[i] << "\n";
                }
                else
                {
                    output3 << x[i] << "," << y[i] << "," << z[i] << "," << h[i] << "\n";
                }
            }
        }
    }

    

    
    auto end_calc = chrono::high_resolution_clock::now();
    
    double init_ms = chrono::duration_cast <std::chrono::milliseconds> (split_index_start-initialize_start).count();
    double split_index_ms = chrono::duration_cast <std::chrono::milliseconds> (splitting_start-split_index_start).count();
    double splitting_ms = chrono::duration_cast <std::chrono::milliseconds> (splitting_end-splitting_start).count();
    double calc_ms = chrono::duration_cast <std::chrono::milliseconds> (end_calc-initialize_start).count();
    
       
    
    printf("Init Time            : %f second\n", init_ms/1000);
    printf("Split Index Time     : %f second\n", split_index_ms/1000);
    printf("Splitting Time       : %f second\n", splitting_ms/1000);
    printf("Total Time           : %f second\n", calc_ms/1000);

}