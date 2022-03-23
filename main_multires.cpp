# include <vector>
# include <fstream>
# include <iostream>
# include <chrono>

# include "src\Particle\Initialize_Particle.hpp"
# include "src\Particle\Particle_Splitting.hpp"
# include "src\Particle\Multires_Movement.hpp"
# include "src\Neighbor Search\Brute_Force.hpp"
# include "src\Particle\Packing_Ratio.hpp"

using namespace std;

int main()
{
    double x_left   = 0;
    double x_right  = 1;
    double y_bottom = 0;
    double y_top    = 1;

    int nx = 40;
    int ny = 40;

    double dx = (x_right-x_left)/nx;
    double dy = (y_top-y_bottom)/ny;

    double n_dummy = 5;

    vector<double> x, y, h;
    vector<int> split_index;
    vector<int> not_moving;
    vector<int> is_dummy;

    auto initialize_start = chrono::high_resolution_clock::now();

    initialize_particle(x_left, x_right, y_bottom, y_top, nx, ny, n_dummy, not_moving, x, y, h);

    auto split_index_start = chrono::high_resolution_clock::now();

    int num_particle = x.size();

    for (int i = 0; i < num_particle; i++)
    {
        if ((x[i] > 0.3) && (x[i] < 0.7))
        {
            split_index.push_back(i);
        }
    }

    auto splitting_start = chrono::high_resolution_clock::now();

    Particle_Splitting(split_index, not_moving, x, y, h);


    num_particle = x.size();

    split_index.clear();
    for (int i = 0; i < num_particle; i++)
    {
        if ((x[i] > 0.4) && (x[i] < 0.6))
        {
            split_index.push_back(i);
        }
    }

    Particle_Splitting(split_index, not_moving, x, y, h);

    num_particle = x.size();
    
    auto splitting_end = chrono::high_resolution_clock::now();    
    
    vector<vector<double>> weight_data;
    vector<vector<int>> neighbor;

    double R_e = 2.1; 

    vector<vector<double>> Ri_a;
    calc_Ri_a(h, R_e, 7, Ri_a);   

    int loop_count = 0;
    int iter = 200;

    while (loop_count < iter)
    {
        brute_force(x, y, h, 1, neighbor, weight_data, R_e);       

        loop_count++;

        vector<vector<double>> Ni_a;
        calc_Ni(x, y, h, neighbor, weight_data, R_e, Ri_a, Ni_a);

        vector<double> ci;
        calc_ci(x, Ri_a, Ni_a, ci);

        vector<double> delta_x, delta_y;
        calc_DeltaX(x, y, h, R_e, 0.1, neighbor, weight_data, ci, delta_x, delta_y);

        for (int i = 0; i < num_particle; i++)
        {
            if (not_moving[i] == 1)
            {
                continue;
            }

            x[i] = x[i] + delta_x[i];
            y[i] = y[i] + delta_y[i];
        }

        ofstream output;

        string name = "output/Test Multiresolusi/Tes Iterasi 5/Distribusi_" + to_string(nx) + "Partikel_" + to_string(loop_count) + "_iter.csv";

        output.open(name);

        output << "x" << "," << "y" << "," << "h\n";

        for (int i = 0; i < x.size(); i++)
        {
            output << x[i] << "," << y[i] << "," << h[i] << "\n";
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