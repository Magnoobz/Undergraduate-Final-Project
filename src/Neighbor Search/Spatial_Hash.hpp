# include <vector>
# include <cmath>
# include <algorithm>
# include <iostream>
# include <chrono>

using namespace std;

void hash_grid(vector<double> x,
               vector<double> y,
               double dx,
               int &ncell_x,
               int &ncell_y,
               int &ncell,
               vector<vector<vector<int>>> &hash_table,
               vector<int> &gridpos_x,
               vector<int> &gridpos_y)
{
    int no_particle = x.size();
    
    double x_min = *min_element(x.begin(),x.end());
    double x_max = *max_element(x.begin(),x.end());
    double y_min = *min_element(y.begin(),y.end());
    double y_max = *max_element(y.begin(),y.end());

    ncell_x = ceil((x_max-x_min)/dx);
    ncell_y = ceil((y_max-y_min)/dx);
    
    ncell = ncell_x*ncell_y;

    double lengrid_x = ncell_x*dx;
    double lengrid_y = ncell_y*dx;

    double x0 = x_min - (lengrid_x-(x_max-x_min))/2;
    double y0 = y_min - (lengrid_y-(y_max-y_min))/2;

    vector<vector<vector<int>>> hash_table_temp(ncell_x, vector<vector<int>>(ncell_y));
    vector<int> gridpos_x_temp(no_particle);
    vector<int> gridpos_y_temp(no_particle);

    # pragma omp parallel for
    for (int i = 0; i < no_particle; i++)
    {
        int pos_x = (x[i]-x0)/dx;
        int pos_y = (y[i]-y0)/dx;

        hash_table_temp[pos_x][pos_y].push_back(i);
        gridpos_x_temp[i] = pos_x;
        gridpos_y_temp[i] = pos_y;
    }

    hash_table = hash_table_temp;
    gridpos_x = gridpos_x_temp;
    gridpos_y = gridpos_y_temp;
}

void spatial_hash_neighbor(vector<double> x,
                           vector<double> y,
                           double dx,
                           int ncell_x,
                           int ncell_y,
                           vector<int> gridpos_x,
                           vector<int> gridpos_y,
                           vector<vector<vector<int>>> hash_table,
                           vector<vector<int>> &neighbor,
                           vector<vector<double>> &weight)
{
    int no_particle = x.size();

    vector<vector<int>> neighbor_temp(no_particle);
    vector<vector<double>> weight_temp(no_particle);

    // # pragma omp parallel for
    for (int i = 0; i < no_particle; i++)
    {
        int pos_x = gridpos_x[i];
        int pos_y = gridpos_y[i];

        double xi = x[i];
        double yi = y[i];

        for (int j = 0; j < 3; j++)
        {
            int gridX = pos_x + j - 1;
            if (gridX < 0 or gridX >= ncell_x){continue;}

            for (int k = 0; k < 3; k++)
            {
                int gridY = pos_y + k - 1;
                if (gridY < 0 or gridY >= ncell_y){continue;}

                # pragma omp parallel for
                for (int s : hash_table[gridX][gridY])
                {
                    if (s == i){continue;}

                    double xj = x[s];
                    double yj = y[s];

                    double dist = sqrt(pow(xj-xi,2)+pow(yj-yi,2));

                    if (dist <= dx)
                    {
                        neighbor_temp[i].push_back(s);
                    }
                }
            }
        }

        sort(neighbor_temp[i].begin(),neighbor_temp[i].end());
        weight_temp[i].resize(neighbor_temp[i].size());

        int count = 0;
        
        # pragma omp parallel for
        for (int t : neighbor_temp[i])
        {
            double xj = x[t];
            double yj = y[t];

            double dist = sqrt(pow(xj-xi,2)+pow(yj-yi,2));
            double weight = pow(1-dist/dx,2);
            weight_temp[i][count] = (weight);

            count++;
        }
    }

    neighbor = neighbor_temp;
    weight = weight_temp;
}

void spatial_hash_neighbor_2(vector<double> x,
                           vector<double> y,
                           vector<double> hx,
                           double eay,
                           double R_e,
                           int ncell_x,
                           int ncell_y,
                           vector<int> gridpos_x,
                           vector<int> gridpos_y,
                           vector<vector<vector<int>>> hash_table,
                           vector<vector<int>> &neighbor,
                           vector<vector<double>> &weight)
{
    int no_particle = x.size();

    vector<vector<int>> neighbor_temp(no_particle);
    vector<vector<double>> weight_temp(no_particle);

    // # pragma omp parallel for
    for (int i = 0; i < no_particle; i++)
    {
        int pos_x = gridpos_x[i];
        int pos_y = gridpos_y[i];

        double xi = x[i];
        double yi = y[i];
        double hi = hx[i]*eay*R_e;

        for (int j = 0; j < 3; j++)
        {
            int gridX = pos_x + j - 1;
            if (gridX < 0 or gridX >= ncell_x){continue;}

            for (int k = 0; k < 3; k++)
            {
                int gridY = pos_y + k - 1;
                if (gridY < 0 or gridY >= ncell_y){continue;}

                # pragma omp parallel for
                for (int s : hash_table[gridX][gridY])
                {
                    if (s == i){continue;}

                    double xj = x[s];
                    double yj = y[s];
                    double hj = hx[s]*eay*R_e;

                    double dist = sqrt(pow(xj-xi,2)+pow(yj-yi,2));

                    if (dist <= 0.5*(hi+hj))
                    {
                        neighbor_temp[i].push_back(s);
                    }
                }
            }
        }

        sort(neighbor_temp[i].begin(),neighbor_temp[i].end());
        weight_temp[i].resize(neighbor_temp[i].size());

        int count = 0;
        
        # pragma omp parallel for
        for (int t : neighbor_temp[i])
        {
            double xj = x[t];
            double yj = y[t];
            double hj = hx[t]*eay*R_e;

            double dist = sqrt(pow(xj-xi,2)+pow(yj-yi,2));
            double weight = pow(1-dist/(0.5*(hi+hj)),2);
            weight_temp[i][count] = (weight);

            count++;
        }
    }

    neighbor = neighbor_temp;
    weight = weight_temp;
}