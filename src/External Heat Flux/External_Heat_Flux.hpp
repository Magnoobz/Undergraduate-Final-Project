#include <vector>
#include "../eigen-3.3.7/Eigen/Dense"

using namespace std;
using namespace Eigen;

void External_Heat_Flux_3D(vector<double> q_x,
                           vector<double> q_y,
                           vector<double> q_z,
                           vector<double> hx,
                           vector<double> hy,
                           vector<double> hz,
                           vector<double> &heat_flux)
{
    int no_particle = hx.size();

    vector<double> temp(no_particle);

    #pragma omp parallel for
    for (int i = 0; i < no_particle; i++)
    {
        temp[i] = hx[i]*hy[i]*q_z[i] + hx[i]*hz[i]*q_y[i]  + hy[i]*hz[i]*q_x[i];
    }

    heat_flux = temp;
}