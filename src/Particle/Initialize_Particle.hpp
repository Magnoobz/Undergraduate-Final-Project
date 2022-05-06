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


void initialize_particle_2(double x_left,
                         double x_right,
                         double y_bottom,
                         double y_top,
                         double x_cut1,
                         double x_cut2,
                         double y_cut1,
                         double y_cut2,
                         int nx,
                         int ny,
                         double dummy,
                         string str,
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

    string movement = str;

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

            if ((x_par>x_cut1) && (x_par<x_cut2) && (y_par>y_cut1) && (y_par<y_cut2))
            {
                continue;
            }
            
            x_temp.push_back(x_par);
            y_temp.push_back(y_par);
            h_temp.push_back(dx);

            if (movement == "all_move")
            {
                not_moving_temp.push_back(0);
            }
            else if (movement == "no_move")
            {
                not_moving_temp.push_back(1);
            }
            else
            {
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
    }

    x.insert(x.end(),x_temp.begin(),x_temp.end());
    y.insert(y.end(),y_temp.begin(),y_temp.end());
    h.insert(h.end(),h_temp.begin(),h_temp.end());

    not_moving.insert(not_moving.end(),not_moving_temp.begin(),not_moving_temp.end());
}