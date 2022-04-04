#include <vector>
#include "../eigen-3.3.7/Eigen/Dense"

using namespace std;
using namespace Eigen;

void calc_kdeltaT_3D(vector<double> x,
                     vector<double> y,
                     vector<double> z,
                     vector<double> k,
                     vector<double> T,
                     vector<int> is_dummy,
                     vector<vector<int>> neighbor,
                     vector<vector<vector<double>>> LSMPS_Eta,
                     vector<double> &kdeltaT_x,
                     vector<double> &kdeltaT_y,
                     vector<double> &kdeltaT_z)
{
    int no_particle = x.size();

    vector<vector<double>> eta_x(no_particle);
    vector<vector<double>> eta_y(no_particle);
    vector<vector<double>> eta_z(no_particle);

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
        
        eta_x[i] = LSMPS_Eta[i][0];
        eta_y[i] = LSMPS_Eta[i][1];
        eta_z[i] = LSMPS_Eta[i][2];

        double xi = x[i];
        double yi = y[i];
        double zi = z[i];
        double ki = k[i];
        double Ti = T[i];

        int no_neighbor = neighbor[i].size();

        #pragma omp parallel for
        for (int j = 0; j < no_neighbor; j++)
        {
            int idxj = neighbor[i][j];

            if (is_dummy[idxj] == 1)
            {
                continue;
            }

            double xj = x[idxj];
            double yj = y[idxj];
            double zj = z[idxj];
            double kj = k[idxj];
            double Tj = T[idxj];

            temp_x[i] += (Tj-Ti)*eta_x[i][j];
            temp_y[i] += (Tj-Ti)*eta_y[i][j];
            temp_z[i] += (Tj-Ti)*eta_z[i][j];
        }

        temp_x[i] = ki*temp_x[i];
        temp_y[i] = ki*temp_y[i];
        temp_z[i] = ki*temp_z[i];
    }

    kdeltaT_x = temp_x;
    kdeltaT_y = temp_y;
    kdeltaT_z = temp_z;
}

void calc_LSMPS_like_3D(vector<double> x,
                        vector<vector<int>> neighbor,
                        vector<double> kdeltaT_x,
                        vector<double> kdeltaT_y,
                        vector<double> kdeltaT_z,
                        vector<int> is_dummy,
                        vector<vector<double>> &kdeltaT_ij_x,
                        vector<vector<double>> &kdeltaT_ij_y,
                        vector<vector<double>> &kdeltaT_ij_z)
{
    int no_particle = x.size();

    vector<vector<double>> temp_x(no_particle);
    vector<vector<double>> temp_y(no_particle);
    vector<vector<double>> temp_z(no_particle);

    // #pragma omp parallel for
    for (int i = 0; i < no_particle; i++)
    {
        if (is_dummy[i] == 0)
        {
            continue;
        }
        
        int no_neighbor = neighbor[i].size();
        temp_x[i].resize(no_neighbor);
        temp_y[i].resize(no_neighbor);
        temp_z[i].resize(no_neighbor);

        #pragma omp parallel for
        for (int j = 0; j < no_neighbor; j++)
        {
            int idxj = neighbor[i][j];

            if(is_dummy[idxj] == 1)
            {
                temp_x[i][j] = (0.5*kdeltaT_x[i]);
                temp_y[i][j] = (0.5*kdeltaT_y[i]);
                temp_z[i][j] = (0.5*kdeltaT_z[i]);
                continue;
            }

            temp_x[i][j] = (0.5*(kdeltaT_x[i]+kdeltaT_x[idxj]));
            temp_y[i][j] = (0.5*(kdeltaT_y[i]+kdeltaT_y[idxj]));
            temp_z[i][j] = (0.5*(kdeltaT_z[i]+kdeltaT_z[idxj]));
        }
    }

    kdeltaT_ij_x = temp_x;
    kdeltaT_ij_y = temp_y;
    kdeltaT_ij_z = temp_z;
}

void calc_MPS_like_value_3D(vector<double> x,
                            vector<double> y,
                            vector<double> z,
                            vector<double> k,
                            vector<double> T,
                            vector<int> is_dummy,
                            vector<vector<int>> neighbor,
                            vector<vector<double>> &kdeltaT_ij)
{
    int no_particle = x.size();

    vector<vector<double>> temp(no_particle);

    // #pragma omp parallel for
    for (int i = 0; i < no_particle; i++)
    {
        if (is_dummy[i] == 1)
        {
            continue;
        }

        double xi = x[i];
        double yi = y[i];
        double zi = z[i];
        double ki = k[i];
        double Ti = T[i];

        int no_neighbor = neighbor[i].size();
        temp[i].resize(no_neighbor);

        #pragma omp parallel for
        for (int j = 0; j < no_neighbor; j++)
        {
            int idxj = neighbor[i][j];

            if (is_dummy[idxj] == 1)
            {
                temp[i][j] = (0);
                continue;
            }

            double xj = x[idxj];
            double yj = y[idxj];
            double zj = z[idxj];
            double kj = k[idxj];
            double Tj = T[idxj];

            double dist = sqrt(pow(xj-xi,2)+pow(yj-yi,2)+pow(zj-zi,2));

            temp[i][j] = (0.5*(ki+kj)*(Tj-Ti)/dist);
        }
    }

    kdeltaT_ij = temp;
}

void calc_MPS_like_3D(vector<double> x,
                      vector<double> y,
                      vector<double> z,
                      vector<int> is_dummy,
                      vector<vector<int>> neighbor,
                      vector<vector<double>> kdeltaT_ij,
                      vector<vector<double>> &kdeltaT_ij_x,
                      vector<vector<double>> &kdeltaT_ij_y,
                      vector<vector<double>> &kdeltaT_ij_z)
{
    int no_particle = x.size();

    vector<vector<double>> temp_x(no_particle);
    vector<vector<double>> temp_y(no_particle);
    vector<vector<double>> temp_z(no_particle);

    // #pragma omp parallel for
    for (int i = 0; i < no_particle; i++)
    {
        if (is_dummy[i] == 1)
        {
            continue;
        }
        
        int no_neighbor = neighbor[i].size();
        temp_x[i].resize(no_neighbor);
        temp_y[i].resize(no_neighbor);
        temp_z[i].resize(no_neighbor);

        double xi = x[i];
        double yi = y[i];
        double zi = z[i];

        #pragma omp parallel for
        for (int j = 0; j < no_neighbor; j++)
        {
            int idxj = neighbor[i][j];

            if(is_dummy[idxj] == 1)
            {
                temp_x[i][j] = (0);
                temp_y[i][j] = (0);
                temp_z[i][j] = (0);
            }

            double xj = x[idxj];
            double yj = y[idxj];
            double zj = z[idxj];

            double value = kdeltaT_ij[i][j];

            double dist = sqrt(pow(xj-xi,2)+pow(yj-yi,2)+pow(zj-zi,2));

            temp_x[i][j] = (value*(xj-xi)/dist);
            temp_y[i][j] = (value*(yj-yi)/dist);
            temp_z[i][j] = (value*(zj-zi)/dist);
        }
    }

    kdeltaT_ij_x = temp_x;
    kdeltaT_ij_y = temp_y;
    kdeltaT_ij_z = temp_z;
}

void calc_hybrid_3D(vector<double> x,
                    vector<double> y,
                    vector<double> z,
                    vector<vector<int>> neighbor,
                    vector<vector<double>> &alpha_ij,
                    vector<double> kdeltaT_x,
                    vector<double> kdeltaT_y,
                    vector<double> kdeltaT_z,
                    vector<int> is_dummy,
                    vector<vector<double>> MPS_like_value,
                    vector<vector<double>> &kdeltaT_ij_x,
                    vector<vector<double>> &kdeltaT_ij_y,
                    vector<vector<double>> &kdeltaT_ij_z)
{
    int no_particle = x.size();

    vector<vector<double>> temp(no_particle);
    vector<vector<double>> temp_x(no_particle);
    vector<vector<double>> temp_y(no_particle);
    vector<vector<double>> temp_z(no_particle);

    // #pragma omp parallel for
    for (int i = 0; i < no_particle; i++)
    {
        if(is_dummy[i] == 1)
        {
            continue;
        }
        
        int no_neighbor = neighbor[i].size();
        temp[i].resize(no_neighbor);
        temp_x[i].resize(no_neighbor);
        temp_y[i].resize(no_neighbor);
        temp_z[i].resize(no_neighbor);

        double xi = x[i];
        double yi = y[i];
        double zi = z[i];
        double kdti_x = kdeltaT_x[i];
        double kdti_y = kdeltaT_y[i];
        double kdti_z = kdeltaT_z[i];

        double alpha_star;

        #pragma omp parallel for
        for (int j = 0; j < no_neighbor; j++)
        {
            int idxj = neighbor[i][j];

            if(is_dummy[idxj] == 1)
            {
                temp[i][j] = (0.5);
                temp_x[i][j] = (0.5*kdeltaT_x[i]);
                temp_y[i][j] = (0.5*kdeltaT_y[i]);
                temp_z[i][j] = (0.5*kdeltaT_z[i]);
                continue;
            }

            double xj = x[idxj];
            double yj = y[idxj];
            double zj = z[idxj];
            double kdtj_x = kdeltaT_x[idxj];
            double kdtj_y = kdeltaT_y[idxj];
            double kdtj_z = kdeltaT_z[idxj];

            double RHS = MPS_like_value[i][j];

            double dist = sqrt(pow(xj-xi,2)+pow(yj-yi,2)+pow(zj-zi,2));
            double xdir = (xj-xi)/dist;
            double ydir = (yj-yi)/dist;
            double zdir = (zj-zi)/dist;
            

            RHS = RHS - (kdtj_x*xdir + kdtj_y*ydir + kdtj_z*zdir);

            double left_coef = (kdti_x-kdtj_x)*xdir + (kdti_y-kdtj_y)*ydir + (kdti_z-kdtj_z)*zdir;

            if (left_coef == 0)
            {
                alpha_star = 0.5;
            }
            else
            {
                alpha_star = RHS/left_coef;
            }

            if (alpha_star < 0)
            {
                temp[i][j] = (0);
            }
            else if (alpha_star > 1)
            {
                temp[i][j] = (1);
            }
            else
            {
                temp[i][j] = (alpha_star);
            }

            temp_x[i][j] = (temp[i][j]*kdeltaT_x[i]+(1-temp[i][j])*kdeltaT_x[idxj]);
            temp_y[i][j] = (temp[i][j]*kdeltaT_y[i]+(1-temp[i][j])*kdeltaT_y[idxj]);
            temp_z[i][j] = (temp[i][j]*kdeltaT_z[i]+(1-temp[i][j])*kdeltaT_z[idxj]);
            
        }
    }

    alpha_ij = temp;
    kdeltaT_ij_x = temp_x;
    kdeltaT_ij_y = temp_y;
    kdeltaT_ij_z = temp_z;
}