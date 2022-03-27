#include <vector>
#include "../eigen-3.3.7/Eigen/Dense"

using namespace std;
using namespace Eigen;

void calc_Laplacian_3D(vector<double> k,
                       vector<double> hx,
                       vector<double> hy,
                       vector<double> hz,
                       vector<vector<int>> neighbor,
                       vector<vector<double>> kdeltaT_ij_x,
                       vector<vector<double>> kdeltaT_ij_y,
                       vector<vector<double>> kdeltaT_ij_z,
                       vector<vector<double>> Sij_Star_x,
                       vector<vector<double>> Sij_Star_y,
                       vector<vector<double>> Sij_Star_z,
                       vector<double> Bi_x,
                       vector<double> Bi_y,
                       vector<double> Bi_z,
                       vector<double> kdeltaT_x,
                       vector<double> kdeltaT_y,
                       vector<double> kdeltaT_z,
                       vector<double> heat_flux,
                       vector<int> is_dummy,
                       vector<double> &Laplacian)
{
    #pragma omp parallel
    int no_particle = k.size();

    vector<double> temp(no_particle);

    for (int i = 0; i < no_particle; i++)
    {
        if(is_dummy[i] == 1)
        {
            continue;
        }
        
        int no_neighbor = neighbor[i].size();

        for (int j = 0; j < no_neighbor; j++)
        {
            if(is_dummy[neighbor[i][j]]==1)
            {
                continue;
            }
            
            temp[i] += kdeltaT_ij_x[i][j]*Sij_Star_x[i][j] + kdeltaT_ij_y[i][j]*Sij_Star_y[i][j] + + kdeltaT_ij_z[i][j]*Sij_Star_z[i][j];
        }

        temp[i] += Bi_x[i]*kdeltaT_x[i] + Bi_y[i]*kdeltaT_y[i] + Bi_z[i]*kdeltaT_z[i] + heat_flux[i];
        
        temp[i] = temp[i]/(k[i]*hx[i]*hy[i]*hz[i]);
    }

    Laplacian = temp;
}

void calc_dTdt_3D(vector<double> cp,
                  vector<double> rho,
                  vector<double> hx,
                  vector<double> hy,
                  vector<double> hz,
                  vector<vector<int>> neighbor,
                  vector<vector<double>> kdeltaT_ij_x,
                  vector<vector<double>> kdeltaT_ij_y,
                  vector<vector<double>> kdeltaT_ij_z,
                  vector<vector<double>> Sij_Star_x,
                  vector<vector<double>> Sij_Star_y,
                  vector<vector<double>> Sij_Star_z,
                  vector<double> Bi_x,
                  vector<double> Bi_y,
                  vector<double> Bi_z,
                  vector<double> kdeltaT_x,
                  vector<double> kdeltaT_y,
                  vector<double> kdeltaT_z,
                  vector<double> heat_flux,
                  vector<int> is_dummy,
                  vector<double> &dTdt)
{
    #pragma omp parallel
    int no_particle = cp.size();

    vector<double> temp(no_particle);

    for (int i = 0; i < no_particle; i++)
    {
        if(is_dummy[i] == 1)
        {
            continue;
        }
        
        int no_neighbor = neighbor[i].size();

        for (int j = 0; j < no_neighbor; j++)
        {
            if(is_dummy[neighbor[i][j]]==1)
            {
                continue;
            }
            
            temp[i] += kdeltaT_ij_x[i][j]*Sij_Star_x[i][j] + kdeltaT_ij_y[i][j]*Sij_Star_y[i][j] + kdeltaT_ij_z[i][j]*Sij_Star_z[i][j];
        }

        temp[i] += Bi_x[i]*kdeltaT_x[i] + Bi_y[i]*kdeltaT_y[i] + Bi_z[i]*kdeltaT_z[i] + heat_flux[i];
        
        temp[i] = temp[i]/(cp[i]*rho[i]*hx[i]*hy[i]*hz[i]);
    }

    dTdt = temp;
}

void calc_Laplacian_From_Eta_3D(vector<vector<vector<double>>> Eta_LSMPS,
                                vector<double> x,
                                vector<double> y,
                                vector<double> z,
                                vector<double> T,
                                vector<int> is_dummy,
                                vector<vector<int>> neighbor,
                                vector<double> &Lap_Value)
{
    #pragma omp parallel
    int no_particle = x.size();

    vector<double> temp_x;
    vector<double> temp_y;
    vector<double> temp_z;

    vector<double> temp(no_particle);

    for (int i = 0; i < no_particle; i++)
    {
        if(is_dummy[i] == 1)
        {
            continue;
        }
        
        int no_neighbor = neighbor[i].size();

        double xi = x[i];
        double yi = y[i];
        double zi = z[i];
        double Ti = T[i];

        temp_x = Eta_LSMPS[i][6];
        temp_y = Eta_LSMPS[i][7];
        temp_z = Eta_LSMPS[i][8];

        for (int j = 0; j < no_neighbor; j++)
        {
            int idxj = neighbor[i][j];

            if(is_dummy[idxj]==1)
            {
                continue;
            }

            double xj = x[idxj];
            double yj = y[idxj];
            double zj = z[idxj];
            double Tj = T[idxj];

            temp[i] += (Tj-Ti)*(temp_x[j]+temp_y[j]+temp_z[j]);
        }
    }

    Lap_Value = temp;
}