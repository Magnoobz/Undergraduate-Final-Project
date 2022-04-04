# include <vector>
# include <cmath>
# include <algorithm>
# include <iostream>
# include <chrono>

using namespace std;

void hash_grid_3D(vector<double> x,
               vector<double> y,
               vector<double> z,
               double dx,
               int &ncell_x,
               int &ncell_y,
               int &ncell_z,
               int &ncell,
               vector<vector<vector<vector<int>>>> &hash_table,
               vector<int> &gridpos_x,
               vector<int> &gridpos_y,
               vector<int> &gridpos_z)
{
    int no_particle = x.size();
    
    double x_min = *min_element(x.begin(),x.end());
    double x_max = *max_element(x.begin(),x.end());
    double y_min = *min_element(y.begin(),y.end());
    double y_max = *max_element(y.begin(),y.end());
    double z_min = *min_element(z.begin(),z.end());
    double z_max = *max_element(z.begin(),z.end());

    ncell_x = ceil((x_max-x_min)/dx);
    ncell_y = ceil((y_max-y_min)/dx);
    ncell_z = ceil((z_max-z_min)/dx);
    
    ncell = ncell_x*ncell_y*ncell_z;

    double lengrid_x = ncell_x*dx;
    double lengrid_y = ncell_y*dx;
    double lengrid_z = ncell_z*dx;

    double x0 = x_min - (lengrid_x-(x_max-x_min))/2;
    double y0 = y_min - (lengrid_y-(y_max-y_min))/2;
    double z0 = z_min - (lengrid_z-(z_max-z_min))/2;

    vector<vector<vector<vector<int>>>> hash_table_temp(ncell_x, vector<vector<vector<int>>>(ncell_y, vector<vector<int>>(ncell_z)));
    vector<int> gridpos_x_temp(no_particle);
    vector<int> gridpos_y_temp(no_particle);
    vector<int> gridpos_z_temp(no_particle);

    # pragma omp parallel for
    for (int i = 0; i < no_particle; i++)
    {
        int pos_x = (x[i]-x0)/dx;
        int pos_y = (y[i]-y0)/dx;
        int pos_z = (z[i]-z0)/dx;

        hash_table_temp[pos_x][pos_y][pos_z].push_back(i);
        gridpos_x_temp[i] = pos_x;
        gridpos_y_temp[i] = pos_y;
        gridpos_z_temp[i] = pos_z;
    }

    hash_table = hash_table_temp;
    gridpos_x = gridpos_x_temp;
    gridpos_y = gridpos_y_temp;
    gridpos_z = gridpos_z_temp;
}

void spatial_hash_neighbor_3D(vector<double> x,
                           vector<double> y,
                           vector<double> z,
                           double dx,
                           int ncell_x,
                           int ncell_y,
                           int ncell_z,
                           vector<int> gridpos_x,
                           vector<int> gridpos_y,
                           vector<int> gridpos_z,
                           vector<vector<vector<vector<int>>>> hash_table,
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
        int pos_z = gridpos_z[i];

        double xi = x[i];
        double yi = y[i];
        double zi = z[i];

        for (int j = 0; j < 3; j++)
        {
            int gridX = pos_x + j - 1;
            if (gridX < 0 or gridX >= ncell_x){continue;}

            for (int k = 0; k < 3; k++)
            {
                int gridY = pos_y + k - 1;
                if (gridY < 0 or gridY >= ncell_y){continue;}

                for (int m = 0; m < 3; m++)
                {
                    int gridZ = pos_z +m -1;
                    if (gridZ < 0 or gridZ >= ncell_z){continue;}

                    # pragma omp parallel for
                    for (int s : hash_table[gridX][gridY][gridZ])
                    {
                        if (s == i){continue;}

                        double xj = x[s];
                        double yj = y[s];
                        double zj = z[s];

                        double dist = sqrt(pow(xj-xi,2)+pow(yj-yi,2)+pow(zj-zi,2));

                        if (dist <= dx)
                        {
                            neighbor_temp[i].push_back(s);
                        }
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
            double zj = z[t];

            double dist = sqrt(pow(xj-xi,2)+pow(yj-yi,2)+pow(zj-zi,2));
            double weight = pow(1-dist/dx,2);
            weight_temp[i][count] = (weight);

            count++;
        }
    }

    neighbor = neighbor_temp;
    weight = weight_temp;
}