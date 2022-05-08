#include <vector>

using namespace std;

void initialize_particle_3D(double x_left,
                         double x_right,
                         double y_bottom,
                         double y_top,
                         double z_back,
                         double z_front,
                         int nx,
                         int ny,
                         int nz,
                         double dummy,
                         vector<int>& not_moving,
                         vector<double>& x,
                         vector<double>& y,
                         vector<double>& z,
                         vector<double>& h)
{
    double x_length = x_right-x_left;
    double y_length = y_top-y_bottom;
    double z_length = z_front-z_back;

    double dx = x_length/nx;
    double dy = y_length/ny;
    double dz = z_length/nz;

    vector<double> x_temp, y_temp, z_temp, h_temp;
    vector<int> not_moving_temp;

    double n_dummy = dummy+0.5;

    double x_par = x_left-n_dummy*dx;
    double y_par = y_bottom-n_dummy*dy;
    double z_par = z_back-n_dummy*dz;

    for (int i = 0; i < nx + 2*n_dummy-1; i++)
    {
        x_par += dx;
        y_par = y_bottom-n_dummy*dy;

        for (int j = 0; j < ny+2*n_dummy-1; j++)
        {
            y_par += dy;
            z_par = z_back - n_dummy*dz;
            
            for (int m = 0; m < nz+2*n_dummy-1; m++)
            {
                z_par += dz;

                x_temp.push_back(x_par);
                y_temp.push_back(y_par);
                z_temp.push_back(z_par);
                h_temp.push_back(dx);

                if ((x_par < x_left-(dummy-1)*dx) || (x_par > x_right+(dummy-1)*dx) || (y_par < y_bottom-(dummy-1)*dy) || (y_par > y_top+(dummy-1)*dy) || (z_par < z_back-(dummy-1)*dz) || (z_par > z_front+(dummy-1)*dz))
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

    x = x_temp;
    y = y_temp;
    z = z_temp;
    h = h_temp;

    not_moving = not_moving_temp;
}

void initialize_particle_3D_2(double x_left,
                         double x_right,
                         double y_bottom,
                         double y_top,
                         double z_back,
                         double z_front,
                         int nx,
                         int ny,
                         int nz,
                         double dummy_r,
                         double dummy_l,
                         double dummy_t,
                         double dummy_bo,
                         double dummy_f,
                         double dummy_ba,
                         vector<int>& not_moving,
                         vector<double>& x,
                         vector<double>& y,
                         vector<double>& z,
                         vector<double>& h)
{
    double x_length = x_right-x_left;
    double y_length = y_top-y_bottom;
    double z_length = z_front-z_back;

    double dx = x_length/nx;
    double dy = y_length/ny;
    double dz = z_length/nz;

    vector<double> x_temp, y_temp, z_temp, h_temp;
    vector<int> not_moving_temp;

    double n_dummy_r = dummy_r+0.5;
    double n_dummy_l = dummy_l+0.5;
    double n_dummy_t = dummy_t+0.5;
    double n_dummy_bo = dummy_bo+0.5;
    double n_dummy_f = dummy_t+0.5;
    double n_dummy_ba = dummy_ba+0.5;

    double x_par = x_left-n_dummy_l*dx;
    double y_par = y_bottom-n_dummy_bo*dy;
    double z_par = z_back-n_dummy_ba*dz;

    for (int i = 0; i < nx + n_dummy_r+n_dummy_l-1; i++)
    {
        x_par += dx;
        y_par = y_bottom-n_dummy_bo*dy;

        for (int j = 0; j < ny+n_dummy_t+n_dummy_bo-1; j++)
        {
            y_par += dy;
            z_par = z_back - n_dummy_ba*dz;
            
            for (int m = 0; m < nz+n_dummy_f+n_dummy_ba-1; m++)
            {
                z_par += dz;

                x_temp.push_back(x_par);
                y_temp.push_back(y_par);
                z_temp.push_back(z_par);
                h_temp.push_back(dx);

                // if ((x_par < x_left-(dummy-1)*dx) || (x_par > x_right+(dummy-1)*dx) || (y_par < y_bottom-(dummy-1)*dy) || (y_par > y_top+(dummy-1)*dy) || (z_par < z_back-(dummy-1)*dz) || (z_par > z_front+(dummy-1)*dz))
                if (true)
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
    z.insert(z.end(),z_temp.begin(),z_temp.end());
    h.insert(h.end(),h_temp.begin(),h_temp.end());

    not_moving.insert(not_moving.end(),not_moving_temp.begin(),not_moving_temp.end());
}