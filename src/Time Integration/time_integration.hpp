#include <vector>

using namespace std;

void time_integration(double dt,
                      vector<double> &T,
                      vector<double> dTdt)
{    
    int no_particle = T.size();

    #pragma omp parallel for num_threads(50)
    for (int i = 0; i < no_particle; i++)
    {
        T[i] += dTdt[i]*dt;
    }
}