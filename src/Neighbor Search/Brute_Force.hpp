# include <vector>
# include <cmath>

using namespace std;

void brute_force(vector<double> x_w,
                 vector<double> y_w,
                 vector<double> hx,
                 double eay,
                 vector<vector<int>> &neighbor,
                 vector<vector<double>> &weight,
                 double R_e)
{
    #pragma omp parallel for
    
    int no_particle = x_w.size();

    double xi, yi;
    double xj, yj;
    double hi, hj;

    vector<vector<int>> neighbor_temp(no_particle);
    vector<vector<double>> weight_temp(no_particle);

    for (int i = 0; i < no_particle-1; i++)
    {
        for (int j = i+1; j < no_particle; j++)
        {
            xi = x_w[i];
            xj = x_w[j];
            
            yi = y_w[i];
            yj = y_w[j];

            hi = hx[i]*eay;
            hj = hx[j]*eay;

            double h_avg = 0.5*(hi+hj);

            double dist = sqrt(pow(xi-xj,2)+pow(yi-yj,2));

            if (dist <= R_e*h_avg)
            {
                double weight = pow(1-dist/(R_e*h_avg),2);

                neighbor_temp[i].push_back(j);
                neighbor_temp[j].push_back(i);
                weight_temp[i].push_back(weight);
                weight_temp[j].push_back(weight);
            }
        }
    }

    neighbor = neighbor_temp;
    weight = weight_temp;
}

void brute_force_3D(vector<double> x_w,
                    vector<double> y_w,
                    vector<double> z_w,
                    vector<double> hx,
                    double eay,
                    double eaz,
                    vector<vector<int>> &neighbor,
                    vector<vector<double>> &weight,
                    double R_e)
{
    #pragma omp parallel for
    
    int no_particle = x_w.size();

    double xi, yi, zi;
    double xj, yj, zj;
    double hi, hj;

    vector<vector<int>> neighbor_temp(no_particle);
    vector<vector<double>> weight_temp(no_particle);

    for (int i = 0; i < no_particle-1; i++)
    {
        for (int j = i+1; j < no_particle; j++)
        {
            xi = x_w[i];
            xj = x_w[j];
            
            yi = y_w[i];
            yj = y_w[j];

            zi = z_w[i];
            zj = z_w[j];

            hi = hx[i]*eay*eaz;
            hj = hx[j]*eay*eaz;

            double h_avg = 0.5*(hi+hj);

            double dist = sqrt(pow(xi-xj,2)+pow(yi-yj,2)+pow(zi-zj,2));

            if (dist <= R_e*h_avg)
            {
                double weight = pow(1-dist/(R_e*h_avg),2);

                neighbor_temp[i].push_back(j);
                neighbor_temp[j].push_back(i);
                weight_temp[i].push_back(weight);
                weight_temp[j].push_back(weight);
            }
        }
    }

    neighbor = neighbor_temp;
    weight = weight_temp;
}

void brute_force_2(vector<double> x_w,
                 vector<double> y_w,
                 vector<double> hx,
                 double eay,
                 vector<vector<int>> &neighbor,
                 vector<vector<double>> &weight,
                 double R_e)
{
    #pragma omp parallel for
    
    int no_particle = x_w.size();

    double xi, yi;
    double xj, yj;
    double hi, hj;

    double h = hx[0]*eay;

    vector<vector<int>> neighbor_temp(no_particle);
    vector<vector<double>> weight_temp(no_particle);

    for (int i = 0; i < no_particle-1; i++)
    {
        for (int j = i+1; j < no_particle; j++)
        {
            xi = x_w[i];
            xj = x_w[j];
            
            yi = y_w[i];
            yj = y_w[j];

            // hi = hx[i]*eay;
            // hj = hx[j]*eay;

            // double h_avg = 0.5*(h+h);

            double dist = sqrt(pow(xi-xj,2)+pow(yi-yj,2));

            if (dist <= R_e*h)
            {
                double weight = pow(1-dist/(R_e*h),2);

                neighbor_temp[i].push_back(j);
                neighbor_temp[j].push_back(i);
                weight_temp[i].push_back(weight);
                weight_temp[j].push_back(weight);
            }
        }
    }

    neighbor = neighbor_temp;
    weight = weight_temp;
}

void brute_force_3D_2(vector<double> x_w,
                    vector<double> y_w,
                    vector<double> z_w,
                    vector<double> hx,
                    double eay,
                    double eaz,
                    vector<vector<int>> &neighbor,
                    vector<vector<double>> &weight,
                    double R_e)
{
    #pragma omp parallel for
    
    int no_particle = x_w.size();

    double xi, yi, zi;
    double xj, yj, zj;
    double hi, hj;

    double h = hx[0]*eay*eaz;

    vector<vector<int>> neighbor_temp(no_particle);
    vector<vector<double>> weight_temp(no_particle);

    for (int i = 0; i < no_particle-1; i++)
    {
        for (int j = i+1; j < no_particle; j++)
        {
            xi = x_w[i];
            xj = x_w[j];
            
            yi = y_w[i];
            yj = y_w[j];

            zi = z_w[i];
            zj = z_w[j];

            // hi = hx[i]*eay*eaz;
            // hj = hx[j]*eay*eaz;

            // double h_avg = 0.5*(hi+hj);

            double dist = sqrt(pow(xi-xj,2)+pow(yi-yj,2)+pow(zi-zj,2));

            if (dist <= R_e*h)
            {
                double weight = pow(1-dist/(R_e*h),2);

                neighbor_temp[i].push_back(j);
                neighbor_temp[j].push_back(i);
                weight_temp[i].push_back(weight);
                weight_temp[j].push_back(weight);
            }
        }
    }

    neighbor = neighbor_temp;
    weight = weight_temp;
}