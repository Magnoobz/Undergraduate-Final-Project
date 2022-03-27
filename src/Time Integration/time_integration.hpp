#include <vector>

using namespace std;

void time_integration(double dt,
                      vector<double> &T,
                      vector<double> dTdt)
{
    #pragma omp parallel for
    
    int no_particle = T.size();

    for (int i = 0; i < no_particle; i++)
    {
        T[i] += dTdt[i]*dt;
    }
}