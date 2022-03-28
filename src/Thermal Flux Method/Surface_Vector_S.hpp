#include <vector>
#include "../eigen-3.3.7/Eigen/Dense"

using namespace std;
using namespace Eigen;

void calculate_Sij(vector<vector<vector<double>>> LSMPS_Eta,
                   vector<vector<double>> &Sij_x,
                   vector<vector<double>> &Sij_y,
                   vector<double> hx,
                   vector<double> hy,
                   vector<vector<int>> neighbor)
{
    int no_particle = hx.size();

    #pragma omp parallel for num_threads(30)
    for (int i = 0; i < no_particle; i++)
    {
        int no_neighbor = neighbor[i].size();

        double Vi = hx[i]*hy[i];

        for (int j = 0; j < no_neighbor; j++)
        {
            Sij_x[i].push_back(2*Vi*LSMPS_Eta[i][0][j]);
            Sij_y[i].push_back(2*Vi*LSMPS_Eta[i][1][j]);
        }
    }
}

void calculate_Bi(vector<vector<double>> Sij_x,
                  vector<vector<double>> Sij_y,
                  vector<double> x,
                  vector<int> is_dummy,
                  vector<vector<int>> neighbor,
                  vector<double> &Bi_x,
                  vector<double> &Bi_y)
{
    int no_particle = x.size();

    vector<double> temp_x(no_particle);
    vector<double> temp_y(no_particle);

    #pragma omp parallel for num_threads(30)
    for (int i = 0; i < no_particle; i++)
    {
        if (is_dummy[i] == 1)
        {
            continue;
        }
        
        int no_neighbor = neighbor[i].size();

        for (int j = 0; j < no_neighbor; j++)
        {
            temp_x[i] = temp_x[i] - Sij_x[i][j];
            temp_y[i] = temp_y[i] - Sij_y[i][j];
        }
    }

    Bi_x = temp_x;
    Bi_y = temp_y;
}

void calculate_Sij_Star(vector<vector<double>> Sij_x,
                        vector<vector<double>> Sij_y,
                        vector<vector<double>> &Sij_Star_x,
                        vector<vector<double>> &Sij_Star_y,
                        vector<double> x,
                        vector<int> is_dummy,
                        vector<vector<int>> neighbor)
{
    int no_particle = x.size();

    #pragma omp parallel for num_threads(30)
    for (int i = 0; i < no_particle; i++)
    {
        if (is_dummy[i] == 1)
        {
            continue;
        }
        
        int no_neighbor = neighbor[i].size();

        for (int j = 0; j < no_neighbor; j++)
        {
            if (is_dummy[i] == 1)
            {
                Sij_Star_x[i].push_back(0);
                Sij_Star_y[i].push_back(0);
                continue;
            }
            
            for (int k = 0; k < neighbor[neighbor[i][j]].size(); k++)
            {
                if (neighbor[neighbor[i][j]][k] == i)
                {
                    Sij_Star_x[i].push_back(0.5*(Sij_x[i][j] - Sij_x[neighbor[i][j]][k]));
                    Sij_Star_y[i].push_back(0.5*(Sij_y[i][j] - Sij_y[neighbor[i][j]][k]));
                    break;
                }
            }
        }
    }
}

void calculate_Sij_Star_2(vector<vector<double>> Sij_x,
                        vector<vector<double>> Sij_y,
                        vector<vector<double>> &Sij_Star_x,
                        vector<vector<double>> &Sij_Star_y,
                        vector<double> x,
                        vector<int> is_dummy,
                        vector<vector<int>> neighbor)
{
    int no_particle = x.size();

    #pragma omp parallel for num_threads(30)
    for (int i = 0; i < no_particle; i++)
    {
        if (is_dummy[i] == 1)
        {
            continue;
        }
        
        int no_neighbor = neighbor[i].size();

        for (int j = 0; j < no_neighbor; j++)
        {
            if (is_dummy[i] == 1)
            {
                Sij_Star_x[i].push_back(0);
                Sij_Star_y[i].push_back(0);
                continue;
            }
            
            for (int k = 0; k < neighbor[neighbor[i][j]].size(); k++)
            {
                if (neighbor[neighbor[i][j]][k] == i)
                {
                    double par = 0.5;
                    
                    Sij_Star_x[i].push_back(par*Sij_x[i][j] - (1-par)*Sij_x[neighbor[i][j]][k]);
                    Sij_Star_y[i].push_back(par*Sij_y[i][j] - (1-par)*Sij_y[neighbor[i][j]][k]);
                    break;
                }
            }
        }
    }
}