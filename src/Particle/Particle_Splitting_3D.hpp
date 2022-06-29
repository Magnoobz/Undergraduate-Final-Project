# include <vector>
# include <cmath>

using namespace std;

void Particle_Splitting_3D(vector<int> split_index,
                        vector<int>& not_moving,
                        vector<double> &x,
                        vector<double> &y,
                        vector<double> &z,
                        vector<double> &h)
{
    int no_particle = x.size();

    int split = 0;

    vector<double> x_temp, y_temp, z_temp, h_temp; 
    vector<int> not_moving_temp;

    double new_r = 1/cbrt(2);
    
    #pragma omp parallel for
    for (int i = 0; i < no_particle; i++)
    {
        if (not_moving[i] == 1)
        {
            x_temp.push_back(x[i]);
            y_temp.push_back(y[i]);
            z_temp.push_back(z[i]);
            h_temp.push_back(h[i]);

            not_moving_temp.push_back(1);
            continue;
        }

        split = 0;
        
        for (auto j : split_index)
        {
            if (j < i){}
            else if (j == i)
            {
                split = 1;
                break;
            }
            else if (j > i)
            {
                break;
            }
        }

        if (split == 1)
        {
            double r = new_r*h[i];
            double dx = 0.5*(h[i]-r);
                    
            x_temp.push_back(x[i]+dx);
            x_temp.push_back(x[i]-dx);

            y_temp.push_back(y[i]+dx);
            y_temp.push_back(y[i]-dx);

            z_temp.push_back(z[i]+dx);
            z_temp.push_back(z[i]-dx);            
            
            h_temp.push_back(r);
            h_temp.push_back(r);

            not_moving_temp.push_back(0);
            not_moving_temp.push_back(0);
        }
        else if (split == 0)
        {
            x_temp.push_back(x[i]);
            y_temp.push_back(y[i]);
            z_temp.push_back(z[i]);
            h_temp.push_back(h[i]);

            not_moving_temp.push_back(0);
        }
    }

    x.clear();
    y.clear();
    z.clear();
    h.clear();
    not_moving.clear();

    x = x_temp;
    y = y_temp;
    z = z_temp;
    h = h_temp;
    not_moving = not_moving_temp;
}

void Particle_Splitting_3D_2(vector<int> split_index,
                        vector<int>& not_moving,
                        vector<double> &x,
                        vector<double> &y,
                        vector<double> &z,
                        vector<double> &h)
{
    int no_particle = x.size();

    int split = 0;

    vector<double> x_temp, y_temp, z_temp, h_temp; 
    vector<int> not_moving_temp;

    double new_r = 1/cbrt(2);
    
    #pragma omp parallel for
    for (int i = 0; i < no_particle; i++)
    {
        if (not_moving[i] == 1)
        {
            x_temp.push_back(x[i]);
            y_temp.push_back(y[i]);
            z_temp.push_back(z[i]);
            h_temp.push_back(h[i]);

            not_moving_temp.push_back(1);
            continue;
        }

        split = 0;
        
        for (auto j : split_index)
        {
            if (j < i){}
            else if (j == i)
            {
                split = 1;
                break;
            }
            else if (j > i)
            {
                break;
            }
        }

        if (split == 1)
        {
            double r = new_r*h[i];
            double dx = 0.5*(h[i]-r);
                    
            x_temp.push_back(x[i]+dx);
            x_temp.push_back(x[i]-dx);

            y_temp.push_back(y[i]+dx);
            y_temp.push_back(y[i]-dx);

            z_temp.push_back(z[i]+dx);
            z_temp.push_back(z[i]-dx);            
            
            h_temp.push_back(r);
            h_temp.push_back(r);

            not_moving_temp.push_back(0);
            not_moving_temp.push_back(0);
        }
        else if (split == 0)
        {
            x_temp.push_back(x[i]);
            y_temp.push_back(y[i]);
            z_temp.push_back(z[i]);
            h_temp.push_back(h[i]);

            not_moving_temp.push_back(1);
        }
    }

    x.clear();
    y.clear();
    z.clear();
    h.clear();
    not_moving.clear();

    x = x_temp;
    y = y_temp;
    z = z_temp;
    h = h_temp;
    not_moving = not_moving_temp;
}