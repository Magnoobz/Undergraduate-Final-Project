#include <vector>
#include "../eigen-3.3.7/Eigen/Dense"

using namespace std;
using namespace Eigen;

void initialize_particle(double x_left,
                         double x_right,
                         double y_bottom,
                         double y_top,
                         double nx,
                         double ny,
                         vector<double> x,
                         vector<double> y,
                         vector<double> h)
{
    double x_length = x_right-x_left;
    double y_length = y_top-y_bottom;

    double dx = x_length/(nx-1);
    double dy = x_length/(ny-1);

    vector<double> x_temp, y_temp, h_temp;

    double x_par = x_left-dx;
    double y_par = y_bottom-dy;

    for (int i = 0; i < nx; i++)
    {
        x_par += dx;
        y_par = y_bottom-dy;

        for (int j = 0; j < ny; j++)
        {
            y_par += dy;

            x_temp.push_back(x_par);
            y_temp.push_back(y_par);
            h_temp.push_back(dx);          
        }
    }
}