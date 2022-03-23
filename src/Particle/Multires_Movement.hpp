# include <vector>
# include <cmath>

using namespace std;

void calc_DeltaX(vector<double> x,
                 vector<double> y,
                 vector<double> h,
                 double R_e,
                 double Cm,
                 vector<vector<int>> neighbor,
                 vector<vector<double>> weight,
                 vector<double> ci,
                 vector<double> &delta_x,
                 vector<double> &delta_y)
{
    #pragma omp parallel for
    
    int no_particle = x.size();

    vector<double> temp_x(no_particle);
    vector<double> temp_y(no_particle);

    for (int i = 0; i < no_particle; i++)
    {
        int no_neighbor = neighbor[i].size();

        double xi = x[i];
        double yi = y[i];
        double Li = h[i];
        double Pi = ci[i];

        for (int j = 0; j < no_neighbor; j++)
        {
            int idxj = neighbor[i][j];
            
            double xj = x[idxj];
            double yj = y[idxj];
            double Lj = h[idxj];
            double Pj = ci[idxj];
            
            double Rij = 0.5*(Li+Lj)*R_e;
            double wij = weight[i][j];

            double dx = xj-xi;
            double dy = yj-yi;

            double dist = sqrt(pow(dx,2)+pow(dy,2));

            temp_x[i] += wij*Cm*(dist-Rij)*Pj/Pi*dx/dist;
            temp_y[i] += wij*Cm*(dist-Rij)*Pj/Pi*dy/dist;
        }
    }

    delta_x = temp_x;
    delta_y = temp_y;
}