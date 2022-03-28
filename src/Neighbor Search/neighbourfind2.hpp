#ifndef NEIGHBOURFIND_HPP
#define NEIGHBOURFIND_HPP

#include <vector>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <chrono>

using namespace std;

//////////////////////////////////////////////////////////////Basic Functions//////////////////////////////////////////////////////////

void zip   (vector<int> sorted,
            vector<double> sorter,
            vector<pair<int,double>> &zipped)
{
    for (int i = 0; i < sorted.size(); i++)
    {
        zipped.push_back(make_pair(sorted[i],sorter[i]));
    }
}

void unzip (vector<int> &sorted,
            vector<double> &sorter,
            vector<pair<int,double>> zipped)
{
    for (int i = 0; i < sorted.size(); i++)
        {
            sorted[i] = zipped[i].first;
            sorter[i] = zipped[i].second;
        }
}

void create_grid   (vector<vector<int>>& cellneighbouroutput,
                    vector<int>& cellneighbour_numoutput,
                    vector<double>& xstartoutput,
                    vector<double>& xendoutput,
                    vector<double>& ystartoutput,
                    vector<double>& yendoutput,
                    int& Ncell,
                    vector<double> x,
                    vector<double> y,
                    double rc,
                    double h)
{
    int size = x.size();

    double x_min = *min_element(x.begin(),x.end()) - 0.0000001;
    double x_max = *max_element(x.begin(),x.end()) + 0.0000001;
    double y_min = *min_element(y.begin(),y.end()) - 0.0000001;
    double y_max = *max_element(y.begin(),y.end()) + 0.0000001;
    
    double dY = 3.5*h; // Can be adjusted for efficiency
    double dX = 3.5*h; // Can be adjusted for efficiency

    // Checking (dX & dY must be larger than rc)
    if (dX < rc || dY < rc)
    {
        cout << "WARNING: dY or dX is set lower than rc" << endl;
    }

    int NcellX = ceil((x_max-x_min)/dX);
    int NcellY = ceil((y_max-y_min)/dY);

    Ncell = NcellX*NcellY;

    vector<vector<int>> cellneighbour(Ncell, vector<int> (9));
    vector<int> cellneighbour_num(Ncell);
    vector<double> xstart(Ncell);
    vector<double> xend(Ncell);
    vector<double> ystart(Ncell);
    vector<double> yend(Ncell);

    // Dividing Domain
    int counter1 = 0;
    int counter2 = 0;
    for (int i = 0; i < NcellX; i++)
    {
        for (int j = 0; j < NcellY; j++)
        {
            // Cell Neighbour Search
            counter2 = 0;
            if ((i != 0) && (j != 0))
            {
                cellneighbour[counter1][counter2] = (counter1-NcellY-1);
                counter2 = counter2+1;
            }
            if ((i != 0) && (j != NcellY-1))
            {
                cellneighbour[counter1][counter2] = (counter1-NcellY+1);
                counter2 = counter2+1;
            }
            if ((i != NcellX-1) && (j != 0))
            {
                cellneighbour[counter1][counter2] = (counter1+NcellY-1);
                counter2 = counter2+1;
            }
            if ((i != NcellX-1) && (j != NcellY-1))
            {
                cellneighbour[counter1][counter2] = (counter1+NcellY+1);
                counter2 = counter2+1;
            }
            if (i != 0)
            {
                cellneighbour[counter1][counter2] = (counter1-NcellY);
                counter2 = counter2+1;
            }
            if (i != NcellX-1)
            {
                cellneighbour[counter1][counter2] = (counter1+NcellY);
                counter2 = counter2+1;
            }
            if (j != 0)
            {
                cellneighbour[counter1][counter2] = (counter1-1);
                counter2 = counter2+1;
            }
            if (j != NcellY-1)
            {
                cellneighbour[counter1][counter2] = (counter1+1);
                counter2 = counter2+1;
            }

            cellneighbour[counter1][counter2] = counter1; // Self neighbour
            counter2 = counter2+1;

            cellneighbour_num[counter1] = counter2;

            // Cell Starting and Ending Coordinates
            xstart[counter1] = x_min + i*dX;
            xend[counter1] = x_min + (i+1)*dX;
            ystart[counter1] = y_min + j*dY;
            yend[counter1] = y_min + (j+1)*dY;

            counter1 = counter1 + 1;
        }
    }

    cellneighbouroutput = cellneighbour;
    cellneighbour_numoutput = cellneighbour_num;
    xstartoutput = xstart;
    xendoutput = xend;
    ystartoutput = ystart;
    yendoutput = yend;
}

void assign_grid   (vector<int>& cellidoutput,
                    vector<vector<int>>& cell_particleoutput,
                    vector<int>& cell_particle_numoutput,
                    vector<double> x,
                    vector<double> y,
                    vector<double> xstart,
                    vector<double> xend,
                    vector<double> ystart,
                    vector<double> yend,
                    int Ncell)
{
    int part_num = x.size();
    vector<int> cellid(part_num,-2);
    vector<vector<int>> cell_particle(Ncell);
    vector<int> cell_particle_num(Ncell);
    #pragma omp parallel for num_threads(30)
    for (int k = 0; k < part_num; k++)
    {
        for (int i = 0; i < Ncell; i++)
        {
            if ((x[k] >= xstart[i]) && (x[k] <= xend[i]) && (y[k] >= ystart[i]) && (y[k] <= yend[i]))
            {
                //#pragma omp critical
                cell_particle[i].push_back(k);
                cell_particle_num[i] = cell_particle_num[i] + 1;
                cellid[k] = i;
                break;
            }
        }
    }
    cellidoutput = cellid;
    cell_particleoutput = cell_particle;
    cell_particle_numoutput = cell_particle_num;
}

void find_neighbours_grid  (vector<vector<int>>& neighbourfull,
                            vector<vector<double>>& rangefull,
                            vector<vector<int>> cellneighbour,
                            vector<int> cellneighbour_num,
                            vector<double> x_i,
                            vector<double> y_i,
                            vector<double> x_j,
                            vector<double> y_j,
                            vector<int> cellid_i,
                            vector<vector<int>> cell_particle_j,
                            vector<int> cell_particle_num_j,
                            double rc)
{
    int i_size = x_i.size();

    double dx;
    double dy;
    double dr;
    int nbi;

    vector<vector<int>> neighbourdata(i_size);
    vector<vector<double>> rangedata(i_size); 

    vector<int> temp_neighbourdata;
    vector<double> temp_drdata;

    for (int k = 0; k < i_size; k++)
    {

        temp_neighbourdata.clear();
        temp_drdata.clear();

        for (int i = 0; i < cellneighbour_num[cellid_i[k]]; i++)
        {
            
            for (int j = 0; j < cell_particle_num_j[cellneighbour[cellid_i[k]][i]]; j++)
            {
                nbi = cell_particle_j[cellneighbour[cellid_i[k]][i]][j];
                dx = x_j[nbi] - x_i[k];
                dy = y_j[nbi] - y_i[k];
                dr = sqrt(dx*dx + dy*dy);

                if (dr < rc)
                {
                    temp_neighbourdata.push_back(nbi);
                    temp_drdata.push_back(dr);
                }
            }
        }

        neighbourdata[k] = temp_neighbourdata;
        rangedata[k] = temp_drdata;
    }

    neighbourfull = neighbourdata;
    rangefull = rangedata;
}

void neighbour_limiting(vector<vector<int>>& neighbour,
                        vector<vector<double>>& range,
                        vector<vector<int>> neighbourfull,
                        vector<vector<double>> rangefull,
                        int neighbourlimit)
{
    int i_size = neighbourfull.size();

    vector<vector<int>> neighbourdata(i_size);
    vector<vector<double>> rangedata(i_size); 

    vector<pair<int,double>> zipped;
    vector<int> temp_neighbourdata;
    vector<double> temp_drdata;
    
    for (int k = 0; k < i_size; k++)
    {
        temp_neighbourdata = neighbourfull[k];
        temp_drdata = rangefull[k];

        if (temp_neighbourdata.size() < neighbourlimit)
        {
            cout<<"WARNING: Neighbour Search Error: Not enough neighbour"<<endl;
        }

        zip(temp_neighbourdata,temp_drdata,zipped);

        sort(begin(zipped), end(zipped), 
        [&](const auto& a, const auto& b)
        {
            return a.second < b.second;
        });

        unzip(temp_neighbourdata,temp_drdata,zipped);
        zipped.clear();

        for (int i = 0; i < neighbourlimit; i++)
        {
            neighbourdata[k].push_back(temp_neighbourdata[i]);
            rangedata[k].push_back(temp_drdata[i]);
        }
    }

    neighbour = neighbourdata;
    range = rangedata;
}

void boundary_check(vector<int>& isboundary,
                    vector<vector<int>> neighbourfull,
                    int boundlimit)
{
    int i_size = neighbourfull.size();

    vector<int> isboundarytemp(i_size);

    for (int k = 0; k < i_size; k++)
    {
        if (neighbourfull[k].size() < boundlimit)
        {
            isboundarytemp[k] = 1;
        }
        else
        {
            isboundarytemp[k] = 0;
        }
    }

    isboundary = isboundarytemp;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void neighbourfind (vector<vector<int>>& neighbour,
                    vector<vector<double>>& range,
                    vector<double> x_i,
                    vector<double> y_i,
                    vector<double> x_j,
                    vector<double> y_j,
                    double rc,
                    double h)
{
    // Dividing domain into grids
    vector<vector<int>> cellneighbour;
    vector<int> cellneighbour_num;
    vector<double> xstart;
    vector<double> xend;
    vector<double> ystart;
    vector<double> yend;
    int Ncell;

    vector<double> x_all = x_i;
    x_all.insert(x_all.end(), x_j.begin(), x_j.end());
    vector<double> y_all = y_i;
    y_all.insert(y_all.end(), y_j.begin(), y_j.end());

    create_grid(cellneighbour,cellneighbour_num,xstart,xend,ystart,yend,Ncell,x_all,y_all,rc,h);

    // Assigning particles into grid
    vector<int> cellid_i;
    vector<vector<int>> cell_particle_i;
    vector<int> cell_particle_num_i;
    assign_grid(cellid_i,cell_particle_i,cell_particle_num_i,x_i,y_i,xstart,xend,ystart,yend,Ncell);

    vector<int> cellid_j;
    vector<vector<int>> cell_particle_j;
    vector<int> cell_particle_num_j;
    assign_grid(cellid_j,cell_particle_j,cell_particle_num_j,x_j,y_j,xstart,xend,ystart,yend,Ncell);

    // Finding particle_i's neighbouring particle_j
    find_neighbours_grid(neighbour,range,cellneighbour,cellneighbour_num,x_i,y_i,x_j,y_j,cellid_i,cell_particle_j,cell_particle_num_j,rc);
}

void neighbourfind_boundary(vector<vector<int>>& neighbour,
                                    vector<vector<double>>& range,
                                    vector<int>& isboundary,
                                    vector<double> x_i,
                                    vector<double> y_i,
                                    vector<double> x_j,
                                    vector<double> y_j,
                                    double rc,
                                    double h,
                                    int boundlimit)
{
    // Dividing domain into grids
    vector<vector<int>> cellneighbour;
    vector<int> cellneighbour_num;
    vector<double> xstart;
    vector<double> xend;
    vector<double> ystart;
    vector<double> yend;
    int Ncell;

    vector<double> x_all = x_i;
    x_all.insert(x_all.end(), x_j.begin(), x_j.end());
    vector<double> y_all = y_i;
    y_all.insert(y_all.end(), y_j.begin(), y_j.end());

    create_grid(cellneighbour,cellneighbour_num,xstart,xend,ystart,yend,Ncell,x_all,y_all,rc,h);

    // Assigning particles into grid
    vector<int> cellid_i;
    vector<vector<int>> cell_particle_i;
    vector<int> cell_particle_num_i;
    assign_grid(cellid_i,cell_particle_i,cell_particle_num_i,x_i,y_i,xstart,xend,ystart,yend,Ncell);

    vector<int> cellid_j;
    vector<vector<int>> cell_particle_j;
    vector<int> cell_particle_num_j;
    assign_grid(cellid_j,cell_particle_j,cell_particle_num_j,x_j,y_j,xstart,xend,ystart,yend,Ncell);

    // Finding particle_i's neighbouring particle_j
    find_neighbours_grid(neighbour,range,cellneighbour,cellneighbour_num,x_i,y_i,x_j,y_j,cellid_i,cell_particle_j,cell_particle_num_j,rc);
    
    // Boundary particles checking
    boundary_check(isboundary,neighbour,boundlimit);
}

void neighbourfind_limited_boundary(vector<vector<int>>& neighbour,
                                    vector<vector<double>>& range,
                                    vector<int>& isboundary,
                                    vector<double> x_i,
                                    vector<double> y_i,
                                    vector<double> x_j,
                                    vector<double> y_j,
                                    double rc,
                                    double h,
                                    int neighbourlimit,
                                    int boundlimit)
{
    // Dividing domain into grids
    vector<vector<int>> cellneighbour;
    vector<int> cellneighbour_num;
    vector<double> xstart;
    vector<double> xend;
    vector<double> ystart;
    vector<double> yend;
    int Ncell;

    vector<double> x_all = x_i;
    x_all.insert(x_all.end(), x_j.begin(), x_j.end());
    vector<double> y_all = y_i;
    y_all.insert(y_all.end(), y_j.begin(), y_j.end());

    create_grid(cellneighbour,cellneighbour_num,xstart,xend,ystart,yend,Ncell,x_all,y_all,rc,h);

    // Assigning particles into grid
    vector<int> cellid_i;
    vector<vector<int>> cell_particle_i;
    vector<int> cell_particle_num_i;
    assign_grid(cellid_i,cell_particle_i,cell_particle_num_i,x_i,y_i,xstart,xend,ystart,yend,Ncell);

    vector<int> cellid_j;
    vector<vector<int>> cell_particle_j;
    vector<int> cell_particle_num_j;
    assign_grid(cellid_j,cell_particle_j,cell_particle_num_j,x_j,y_j,xstart,xend,ystart,yend,Ncell);

    // Finding particle_i's neighbouring particle_j
    vector<vector<int>> neighbourfull;
    vector<vector<double>> rangefull;
    find_neighbours_grid(neighbourfull,rangefull,cellneighbour,cellneighbour_num,x_i,y_i,x_j,y_j,cellid_i,cell_particle_j,cell_particle_num_j,rc);

    // Limiting the maximum neighbouring particles
    neighbour_limiting(neighbour,range,neighbourfull,rangefull,neighbourlimit);
    
    // Boundary particles checking
    boundary_check(isboundary,neighbourfull,boundlimit);
}

#endif