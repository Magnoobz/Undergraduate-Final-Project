#include <vector>

using namespace std;

void initialize_particle(double x_left,
                         double x_right,
                         double y_bottom,
                         double y_top,
                         int nx,
                         int ny,
                         double dummy,
                         vector<int>& not_moving,
                         vector<double>& x,
                         vector<double>& y,
                         vector<double>& h)
{
    double x_length = x_right-x_left;
    double y_length = y_top-y_bottom;

    double dx = x_length/nx;
    double dy = y_length/ny;

    vector<double> x_temp, y_temp, h_temp;
    vector<int> not_moving_temp;

    double n_dummy = dummy+0.5;

    double x_par = x_left-n_dummy*dx;
    double y_par = y_bottom-n_dummy*dy;

    for (int i = 0; i < nx + 2*n_dummy-1; i++)
    {
        x_par += dx;
        y_par = y_bottom-n_dummy*dy;

        for (int j = 0; j < ny+2*n_dummy-1; j++)
        {
            y_par += dy;

            x_temp.push_back(x_par);
            y_temp.push_back(y_par);
            h_temp.push_back(dx);

            if ((x_par < x_left-(dummy-1)*dx) || (x_par > x_right+(dummy-1)*dx) || (y_par < y_bottom-(dummy-1)*dy) || (y_par > y_top+(dummy-1)*dy))
            {
                not_moving_temp.push_back(1);
            }
            else
            {
                not_moving_temp.push_back(0);
            }
        }
    }

    x = x_temp;
    y = y_temp;
    h = h_temp;

    not_moving = not_moving_temp;
}