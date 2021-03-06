#include <vector>
#include "../eigen-3.3.7/Eigen/Dense"

using namespace std;
using namespace Eigen;

void calculate_Sij_3D(vector<vector<vector<double>>> LSMPS_Eta,
                      vector<vector<double>> &Sij_x,
                      vector<vector<double>> &Sij_y,
                      vector<vector<double>> &Sij_z,
                      vector<double> hx,
                      vector<double> hy,
                      vector<double> hz,
                      vector<vector<int>> neighbor)
{
    int no_particle = hx.size();

    // #pragma omp parallel for
    for (int i = 0; i < no_particle; i++)
    {
        int no_neighbor = neighbor[i].size();
        Sij_x[i].resize(no_neighbor);
        Sij_y[i].resize(no_neighbor);
        Sij_z[i].resize(no_neighbor);

        double Vi = hx[i]*hy[i]*hz[i];

        #pragma omp parallel for
        for (int j = 0; j < no_neighbor; j++)
        {
            Sij_x[i][j] = (2*Vi*LSMPS_Eta[i][0][j]);
            Sij_y[i][j] = (2*Vi*LSMPS_Eta[i][1][j]);
            Sij_z[i][j] = (2*Vi*LSMPS_Eta[i][2][j]);
        }
    }
}

void calculate_Bi_3D(vector<vector<double>> Sij_x,
                     vector<vector<double>> Sij_y,
                     vector<vector<double>> Sij_z,
                     vector<double> x,
                     vector<int> is_dummy,
                     vector<vector<int>> neighbor,
                     vector<double> &Bi_x,
                     vector<double> &Bi_y,
                     vector<double> &Bi_z)
{
    int no_particle = x.size();

    vector<double> temp_x(no_particle);
    vector<double> temp_y(no_particle);
    vector<double> temp_z(no_particle);

    // #pragma omp parallel for
    for (int i = 0; i < no_particle; i++)
    {
        if (is_dummy[i] == 1)
        {
            continue;
        }
        
        int no_neighbor = neighbor[i].size();

        #pragma omp parallel for
        for (int j = 0; j < no_neighbor; j++)
        {
            temp_x[i] = temp_x[i] - Sij_x[i][j];
            temp_y[i] = temp_y[i] - Sij_y[i][j];
            temp_z[i] = temp_z[i] - Sij_z[i][j];
        }
    }

    Bi_x = temp_x;
    Bi_y = temp_y;
    Bi_z = temp_z;
}

void calculate_Sij_Star_3D(vector<vector<double>> Sij_x,
                           vector<vector<double>> Sij_y,
                           vector<vector<double>> Sij_z,
                           vector<vector<double>> &Sij_Star_x,
                           vector<vector<double>> &Sij_Star_y,
                           vector<vector<double>> &Sij_Star_z,
                           vector<double> x,
                           vector<int> is_dummy,
                           vector<vector<int>> neighbor)
{
    int no_particle = x.size();

    // #pragma omp parallel for
    for (int i = 0; i < no_particle; i++)
    {
        if (is_dummy[i] == 1)
        {
            continue;
        }
        
        int no_neighbor = neighbor[i].size();
        Sij_Star_x[i].resize(no_neighbor);
        Sij_Star_y[i].resize(no_neighbor);
        Sij_Star_z[i].resize(no_neighbor);

        #pragma omp parallel for
        for (int j = 0; j < no_neighbor; j++)
        {
            if (is_dummy[i] == 1)
            {
                Sij_Star_x[i][j] = (0);
                Sij_Star_y[i][j] = (0);
                Sij_Star_z[i][j] = (0);
                continue;
            }
            
            for (int k = 0; k < neighbor[neighbor[i][j]].size(); k++)
            {
                if (neighbor[neighbor[i][j]][k] == i)
                {
                    Sij_Star_x[i][j] = (0.5*(Sij_x[i][j] - Sij_x[neighbor[i][j]][k]));
                    Sij_Star_y[i][j] = (0.5*(Sij_y[i][j] - Sij_y[neighbor[i][j]][k]));
                    Sij_Star_z[i][j] = (0.5*(Sij_z[i][j] - Sij_z[neighbor[i][j]][k]));
                    break;
                }
            }
        }
    }
}