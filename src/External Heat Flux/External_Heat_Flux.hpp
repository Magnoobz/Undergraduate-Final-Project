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

void Internal_Heat_flux_3D(vector<double> Q,
                           vector<double> hx,
                           vector<double> hy,
                           vector<double> hz,
                           vector<double> &heat_flux)
{
    int no_particle = hx.size();

    # pragma omp parallel for
    for (int i = 0; i < no_particle; i++)
    {
        heat_flux[i] += hx[i]*hy[i]*hz[i]*Q[i];
    }
}

void Convection_Heating_3D(vector<double> conv_x,
                           vector<double> conv_y,
                           vector<double> conv_z,
                           vector<double> T_conv,
                           vector<double> hx,
                           vector<double> hy,
                           vector<double> hz,
                           vector<double> &heat_flux)
{
    int no_particle = hx.size();

    #pragma omp parallel for
    for (int i = 0; i < no_particle; i++)
    {
        heat_flux[i] += (hx[i]*hy[i]*conv_z[i] + hx[i]*hz[i]*conv_y[i]  + hy[i]*hz[i]*conv_x[i])*T_conv[i];
    }
}