#include <vector>
#include "../eigen-3.3.7/Eigen/Dense"

using namespace std;
using namespace Eigen;

void calc_LSMPS_eta(vector<vector<vector<double>>> &LSMPS_eta,
                    vector<double> x,
                    vector<double> y,
                    vector<double> hx,
                    vector<double> hy,
                    vector<vector<int>> neighbor,
                    vector<vector<double>> weight_data)
{    
    int no_particle = x.size();

    #pragma omp parallel for
    for (int i = 0; i < no_particle; i++)
    {
        int no_neighbor = neighbor[i].size();
        
        double xi = x[i];
		double yi = y[i];
        double hxi = hx[i];
        double hyi = hy[i];


        MatrixXd Hrs = MatrixXd::Zero(5,5);
        
        Hrs(0,0) = pow(hxi,-1);
        Hrs(1,1) = pow(hyi,-1);
        Hrs(2,2) = pow(hxi,-2)*2.0;
        Hrs(3,3) = pow(hxi*hyi,-1);
        Hrs(4,4) = pow(hyi,-2)*2.0;


        MatrixXd M = MatrixXd::Zero(5,5);
        MatrixXd P = MatrixXd::Zero(5,1);
        MatrixXd bi = MatrixXd::Zero(5,no_neighbor);

        for (int j = 0; j < no_neighbor; j++)
        {
            int idxj = neighbor[i][j];

            double xj = x[idxj];
            double yj = y[idxj];

            double dx_scaled = (xj-xi)/hxi;
            double dy_scaled = (yj-yi)/hyi;

            P(0,0) = dx_scaled;
            P(1,0) = dy_scaled;
            P(2,0) = dx_scaled*dx_scaled;
            P(3,0) = dx_scaled*dy_scaled;
            P(4,0) = dy_scaled*dy_scaled;

            for (int k = 0; k < 5; k++)
            {
                for (int l = 0; l < 5; l++)
                {
                    M(k,l) = M(k,l) + weight_data[i][j]*P(k,0)*P(l,0);
                }
            }

            for (int k = 0; k < 5; k++)
            {
                bi(k,j) = weight_data[i][j]*P(k,0);
            }
        }

        MatrixXd M_Inv_Bi = M.bdcSvd(ComputeThinU | ComputeThinV).solve(bi);
        MatrixXd LSMPS = Hrs*M_Inv_Bi;

        // MatrixXd Minv = M.inverse();
        // MatrixXd LSMPS = Hrs*Minv*bi;

        for (int j = 0; j < no_neighbor; j++)
        {
            LSMPS_eta[i][0].push_back(LSMPS(0,j)); 
            LSMPS_eta[i][1].push_back(LSMPS(1,j));
            LSMPS_eta[i][2].push_back(LSMPS(2,j));
            LSMPS_eta[i][3].push_back(LSMPS(3,j));
            LSMPS_eta[i][4].push_back(LSMPS(4,j));
        }
    }
}

void calc_LSMPS_eta_2(vector<vector<vector<double>>> &LSMPS_eta,
                    vector<double> x,
                    vector<double> y,
                    vector<double> hx,
                    vector<double> hy,
                    vector<vector<int>> neighbor,
                    vector<vector<double>> weight_data)
{
    int no_particle = x.size();

    double hxi = hx[0];
    double hyi = hy[0];

    MatrixXd Hrs = MatrixXd::Zero(5,5);
        
    Hrs(0,0) = pow(hxi,-1);
    Hrs(1,1) = pow(hyi,-1);
    Hrs(2,2) = pow(hxi,-2)*2.0;
    Hrs(3,3) = pow(hxi*hyi,-1);
    Hrs(4,4) = pow(hyi,-2)*2.0;

    #pragma omp parallel for
    for (int i = 0; i < no_particle; i++)
    {
        int no_neighbor = neighbor[i].size();
        
        double xi = x[i];
		double yi = y[i];
        // double hxi = hx[i];
        // double hyi = hy[i];

        MatrixXd M = MatrixXd::Zero(5,5);
        MatrixXd P = MatrixXd::Zero(5,1);
        MatrixXd bi = MatrixXd::Zero(5,no_neighbor);

        for (int j = 0; j < no_neighbor; j++)
        {
            int idxj = neighbor[i][j];

            double xj = x[idxj];
            double yj = y[idxj];
            double hxj = hx[idxj];
            double hyj = hy[idxj];

            double v_ratio = hxj*hyj/(hxi*hyi);

            double dx_scaled = (xj-xi)/hxi;
            double dy_scaled = (yj-yi)/hyi;

            P(0,0) = dx_scaled;
            P(1,0) = dy_scaled;
            P(2,0) = dx_scaled*dx_scaled;
            P(3,0) = dx_scaled*dy_scaled;
            P(4,0) = dy_scaled*dy_scaled;

            for (int k = 0; k < 5; k++)
            {
                for (int l = 0; l < 5; l++)
                {
                    M(k,l) = M(k,l) + v_ratio*weight_data[i][j]*P(k,0)*P(l,0);
                }
            }

            for (int k = 0; k < 5; k++)
            {
                bi(k,j) = v_ratio*weight_data[i][j]*P(k,0);
            }
        }

        MatrixXd M_Inv_Bi = M.bdcSvd(ComputeThinU | ComputeThinV).solve(bi);
        MatrixXd LSMPS = Hrs*M_Inv_Bi;

        // MatrixXd Minv = M.inverse();
        // MatrixXd LSMPS = Hrs*Minv*bi;

        for (int j = 0; j < no_neighbor; j++)
        {
            LSMPS_eta[i][0].push_back(LSMPS(0,j)); 
            LSMPS_eta[i][1].push_back(LSMPS(1,j));
            LSMPS_eta[i][2].push_back(LSMPS(2,j));
            LSMPS_eta[i][3].push_back(LSMPS(3,j));
            LSMPS_eta[i][4].push_back(LSMPS(4,j));
        }
    }
}