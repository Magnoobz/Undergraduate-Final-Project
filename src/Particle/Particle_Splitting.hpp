# include <vector>
# include <cmath>

using namespace std;

void Particle_Splitting(vector<int> split_index,
                        vector<int>& not_moving,
                        vector<double> &x,
                        vector<double> &y,
                        vector<double> &h)
{
    #pragma omp parallel

    int no_particle = x.size();

    int split = 0;

    vector<double> x_temp, y_temp, h_temp; 
    vector<int> not_moving_temp;

    for (int i = 0; i < no_particle; i++)
    {
        if (not_moving[i] == 1)
        {
            x_temp.push_back(x[i]);
            y_temp.push_back(y[i]);
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
            double dx = (0.5-0.25*sqrt(2))*h[i];
                    
            x_temp.push_back(x[i]+dx);
            x_temp.push_back(x[i]-dx);

            y_temp.push_back(y[i]+dx);
            y_temp.push_back(y[i]-dx);            
            
            h_temp.push_back(h[i]/sqrt(2));
            h_temp.push_back(h[i]/sqrt(2));

            not_moving_temp.push_back(0);
            not_moving_temp.push_back(0);
        }
        else if (split == 0)
        {
            x_temp.push_back(x[i]);
            y_temp.push_back(y[i]);
            h_temp.push_back(h[i]);

            not_moving_temp.push_back(0);
        }
    }

    x.clear();
    y.clear();
    h.clear();
    not_moving.clear();

    x = x_temp;
    y = y_temp;
    h = h_temp;
    not_moving = not_moving_temp;
}