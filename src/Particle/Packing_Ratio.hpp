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

    #pragma omp parallel for num_threads(50)
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

    #pragma omp parallel for num_threads(50)
    for (int i = 0; i < A; i++)
    {
        for (int j = 0; j < no_particle; j++)
        {
            for (int k = 0; k < no_particle; k++)
            {
                if (k == j)
                {
                    continue;
                }
                else
                {
                    double xj = x[j];
                    double yj = y[j];
                    double Lj = h[j];

                    double xk = x[k];
                    double yk = y[k];
                    double Lk = h[k];

                    double r_jk = sqrt(pow(xk-xj,2)+pow(yk-yj,2));

                    if (Ri_a[j][i] < r_jk - Lk/2)
                    {
                        continue;
                    }
                    else if (Ri_a[j][i] >= r_jk - Lk/2 && Ri_a[j][i] <= r_jk + Lk/2)
                    {
                        temp[j][i] += (Ri_a[j][i]-r_jk+Lj/2)*Lj;
                    }
                    else if (Ri_a[j][i] > r_jk + Lk/2)
                    {
                        temp[j][i] += pow(Lj,2);
                    }


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

    #pragma omp parallel for num_threads(50)
    for (int i = 0; i < no_particle; i++)
    {
        for (int j = 0; j < A; j++)
        {
            temp_num += Ni_a[i][j]*Ri_a[i][j];
            temp_denum += Ri_a[i][j];
        }

        temp[i] = temp_num/temp_denum/M_PI;
    }

    ci = temp;
}