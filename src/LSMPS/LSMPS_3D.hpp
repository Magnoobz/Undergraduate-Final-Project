#include <vector>
#include "../eigen-3.3.7/Eigen/Dense"

using namespace std;
using namespace Eigen;

void calc_LSMPS_eta_3D(vector<vector<vector<double>>> &LSMPS_eta,
                       vector<double> x,
                       vector<double> y,
                       vector<double> z,
                       vector<double> hx,
                       vector<double> hy,
                       vector<double> hz,
                       vector<vector<int>> neighbor,
                       vector<vector<double>> weight_data)
{
    int no_particle = x.size();

    // #pragma omp parallel for
    for (int i = 0; i < no_particle; i++)
    {
        int no_neighbor = neighbor[i].size();
        
        double xi = x[i];
		double yi = y[i];
        double zi = z[i];
        double hxi = hx[i];
        double hyi = hy[i];
        double hzi = hz[i];

        MatrixXd Hrs = MatrixXd::Zero(9,9);
        
        Hrs(0,0) = pow(hxi,-1);
        Hrs(1,1) = pow(hyi,-1);
        Hrs(2,2) = pow(hzi,-1);
        Hrs(3,3) = pow(hxi*hyi,-1);
        Hrs(4,4) = pow(hxi*hzi,-1);
        Hrs(5,5) = pow(hyi*hzi,-1);
        Hrs(6,6) = pow(hxi,-2)*2.0;
        Hrs(7,7) = pow(hyi,-2)*2.0;
        Hrs(8,8) = pow(hzi,-2)*2.0;


        MatrixXd M = MatrixXd::Zero(9,9);
        MatrixXd P = MatrixXd::Zero(9,1);
        MatrixXd bi = MatrixXd::Zero(9,no_neighbor);

        #pragma omp parallel for
        for (int j = 0; j < no_neighbor; j++)
        {
            int idxj = neighbor[i][j];

            double xj = x[idxj];
            double yj = y[idxj];
            double zj = z[idxj];

            double dx_scaled = (xj-xi)/hxi;
            double dy_scaled = (yj-yi)/hyi;
            double dz_scaled = (zj-zi)/hzi;

            P(0,0) = dx_scaled;
            P(1,0) = dy_scaled;
            P(2,0) = dz_scaled;
            P(3,0) = dx_scaled*dy_scaled;
            P(4,0) = dx_scaled*dz_scaled;
            P(5,0) = dy_scaled*dz_scaled;
            P(6,0) = dx_scaled*dx_scaled;
            P(7,0) = dy_scaled*dy_scaled;
            P(8,0) = dz_scaled*dz_scaled;

            for (int k = 0; k < 9; k++)
            {
                for (int l = 0; l < 9; l++)
                {
                    M(k,l) = M(k,l) + weight_data[i][j]*P(k,0)*P(l,0);
                }
            }

            for (int k = 0; k < 9; k++)
            {
                bi(k,j) = weight_data[i][j]*P(k,0);
            }
        }

        MatrixXd M_Inv_Bi = M.bdcSvd(ComputeThinU | ComputeThinV).solve(bi);
        MatrixXd LSMPS = Hrs*M_Inv_Bi;

        // MatrixXd Minv = M.inverse();
        // MatrixXd LSMPS = Hrs*Minv*bi;

        #pragma omp parallel for
        for (int j = 0; j < 9; j++)
        {
            LSMPS_eta[i][j].resize(no_neighbor);
        }

        #pragma omp parallel for
        for (int j = 0; j < no_neighbor; j++)
        {
            LSMPS_eta[i][0][j] = (LSMPS(0,j)); 
            LSMPS_eta[i][1][j] = (LSMPS(1,j));
            LSMPS_eta[i][2][j] = (LSMPS(2,j));
            LSMPS_eta[i][3][j] = (LSMPS(3,j));
            LSMPS_eta[i][4][j] = (LSMPS(4,j));
            LSMPS_eta[i][5][j] = (LSMPS(5,j));
            LSMPS_eta[i][6][j] = (LSMPS(6,j));
            LSMPS_eta[i][7][j] = (LSMPS(7,j));
            LSMPS_eta[i][8][j] = (LSMPS(8,j));
        }
    }
}

void calc_LSMPS_eta_3D_2(vector<vector<vector<double>>> &LSMPS_eta,
                       vector<double> x,
                       vector<double> y,
                       vector<double> z,
                       vector<double> hx,
                       vector<double> hy,
                       vector<double> hz,
                       vector<vector<int>> neighbor,
                       vector<vector<double>> weight_data)
{
    int no_particle = x.size();

    double hxi = hx[0];
    double hyi = hy[0];
    double hzi = hz[0];

    MatrixXd Hrs = MatrixXd::Zero(9,9);
        
    Hrs(0,0) = pow(hxi,-1);
    Hrs(1,1) = pow(hyi,-1);
    Hrs(2,2) = pow(hzi,-1);
    Hrs(3,3) = pow(hxi*hyi,-1);
    Hrs(4,4) = pow(hxi*hzi,-1);
    Hrs(5,5) = pow(hyi*hzi,-1);
    Hrs(6,6) = pow(hxi,-2)*2.0;
    Hrs(7,7) = pow(hyi,-2)*2.0;
    Hrs(8,8) = pow(hzi,-2)*2.0;

    // #pragma omp parallel for
    for (int i = 0; i < no_particle; i++)
    {
        int no_neighbor = neighbor[i].size();
        
        double xi = x[i];
		double yi = y[i];
        double zi = z[i];
        // double hxi = hx[i];
        // double hyi = hy[i];
        // double hzi = hz[i];

        MatrixXd M = MatrixXd::Zero(9,9);
        MatrixXd P = MatrixXd::Zero(9,1);
        MatrixXd bi = MatrixXd::Zero(9,no_neighbor);

        #pragma omp parallel for
        for (int j = 0; j < no_neighbor; j++)
        {
            int idxj = neighbor[i][j];

            double xj = x[idxj];
            double yj = y[idxj];
            double zj = z[idxj];
            double hxj = hx[idxj];
            double hyj = hy[idxj];
            double hzj = hz[idxj];

            double v_ratio = hxj*hyj*hzj/(hxi*hyi*hzi);

            double dx_scaled = (xj-xi)/hxi;
            double dy_scaled = (yj-yi)/hyi;
            double dz_scaled = (zj-zi)/hzi;

            P(0,0) = dx_scaled;
            P(1,0) = dy_scaled;
            P(2,0) = dz_scaled;
            P(3,0) = dx_scaled*dy_scaled;
            P(4,0) = dx_scaled*dz_scaled;
            P(5,0) = dy_scaled*dz_scaled;
            P(6,0) = dx_scaled*dx_scaled;
            P(7,0) = dy_scaled*dy_scaled;
            P(8,0) = dz_scaled*dz_scaled;

            for (int k = 0; k < 9; k++)
            {
                for (int l = 0; l < 9; l++)
                {
                    M(k,l) = M(k,l) + v_ratio*weight_data[i][j]*P(k,0)*P(l,0);
                }
            }

            for (int k = 0; k < 9; k++)
            {
                bi(k,j) = v_ratio*weight_data[i][j]*P(k,0);
            }
        }

        MatrixXd M_Inv_Bi = M.bdcSvd(ComputeThinU | ComputeThinV).solve(bi);
        MatrixXd LSMPS = Hrs*M_Inv_Bi;

        // MatrixXd Minv = M.inverse();
        // MatrixXd LSMPS = Hrs*Minv*bi;

        #pragma omp parallel for
        for (int j = 0; j < 9; j++)
        {
            LSMPS_eta[i][j].resize(no_neighbor);
        }

        #pragma omp parallel for
        for (int j = 0; j < no_neighbor; j++)
        {
            LSMPS_eta[i][0][j] = (LSMPS(0,j)); 
            LSMPS_eta[i][1][j] = (LSMPS(1,j));
            LSMPS_eta[i][2][j] = (LSMPS(2,j));
            LSMPS_eta[i][3][j] = (LSMPS(3,j));
            LSMPS_eta[i][4][j] = (LSMPS(4,j));
            LSMPS_eta[i][5][j] = (LSMPS(5,j));
            LSMPS_eta[i][6][j] = (LSMPS(6,j));
            LSMPS_eta[i][7][j] = (LSMPS(7,j));
            LSMPS_eta[i][8][j] = (LSMPS(8,j));
        }
    }
}