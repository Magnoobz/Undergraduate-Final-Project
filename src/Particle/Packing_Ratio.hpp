# define _USE_MATH_DEFINE
# include <vector>
# include <cmath>
# include "../eigen-3.3.7/Eigen/Dense"

using namespace std;
using namespace Eigen;

void calc_Ri_a(vector<double> h,
               double Re,
               int A,
               vector<vector<double>> &Ri_a)
{
    int no_particle = h.size();

    vector<vector<double>> temp(no_particle);

    #pragma omp parallel for
    for (int i = 0; i < no_particle; i++)
    {
        double Ri = h[i]*Re;

        for (int j = 0; j < A; j++)
        {
            temp[i].push_back(h[i]/2 + j/(A-1)*(Ri-h[i]));
        }
    }

    Ri_a = temp;
}

void calc_Ni(vector<double> x,
             vector<double> y,
             vector<double> h,
             vector<vector<int>> neighbor,
             vector<vector<double>> weight,
             double R_e,
             vector<vector<double>> Ri_a,
             vector<vector<double>> &Ni_a)
{
    int no_particle = x.size();
    int A = Ri_a[0].size();

    vector<vector<double>> temp(no_particle,vector<double>(A));

    #pragma omp parallel for
    for (int i = 0; i < no_particle; i++)
    {
        int no_neighbor = neighbor[i].size();

        double xi = x[i];
        double yi = y[i];
        double Li = h[i];

        for (int j = 0; j < A; j++)
        {
            for (int k = 0; k < no_neighbor; k++)
            {
                int idxk = neighbor[i][k];
                
                double xk = x[idxk];
                double yk = y[idxk];
                double Lk = h[idxk];

                double r_ik = sqrt(pow(xk-xi,2)+pow(yk-yi,2));

                if (Ri_a[i][j] >= r_ik - Lk/2 && Ri_a[i][j] <= r_ik + Lk/2)
                {
                    temp[i][j] += (Ri_a[i][j]-r_ik+Lk/2)*Lk;
                }
                else if (Ri_a[i][j] > r_ik + Lk/2)
                {
                    temp[i][j] += pow(Lk,2);
                }
                
            }        
        }
    }

    Ni_a = temp;    
}

void calc_ci (vector<double> x,
              vector<vector<double>> Ri_a,
              vector<vector<double>> Ni_a,
              vector<double> &ci)
{
    int no_particle = x.size();
    int A = Ri_a[0].size();

    double temp_num = 0;
    double temp_denum = 0;

    vector<double> temp(no_particle);

    #pragma omp parallel for
    for (int i = 0; i < no_particle; i++)
    {
        for (int j = 0; j < A; j++)
        {
            temp_num += Ni_a[i][j]*Ri_a[i][j];
            temp_denum += Ri_a[i][j];
        }

        temp[i] = temp_num/temp_denum;
    }

    ci = temp;
}