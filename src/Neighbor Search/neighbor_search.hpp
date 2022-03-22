#include <bits/stdc++.h>
#include <fstream>

#ifndef INCLUDE_PARTICLE
#include "particle.hpp"
#endif

#ifndef INCLUDE_GLOBAL
#include "global.hpp"
#endif

using namespace std;

#ifndef INCLUDE_NEIGHBORSEARCH
#define INCLUDE_NEIGHBORSEARCH
class Grid
{
public:
    struct coordinate
    {
        double x, y;
    };

    vector<vector<int>> grid;
    Grid(): grid(table_size) 
    {
        cout << "Starting neighbor search \n";
    }
    ~Grid() 
    {
        cout << "Neighbor search finished \n";
    }
    const int hashkey(const double &x, const double &y)
    {
        return (int(x * conversion_factor) + int(y * conversion_factor) * rc);
    }
    void organize(const Particle &p)
    {
        for (int i = 0; i < N; ++i)
        {
            grid[hashkey(p.m_x[0][i], p.m_x[1][i])].push_back(p.m_index[i]);
        }   
    }
    const bool distance(const Particle &p, const int &i, const int &j)
    {
        //if (p.m_index[i] != p.m_index[j])
        //{
            if (pow(p.m_x[0][i] - p.m_x[0][j], 2) + pow(p.m_x[1][i] - p.m_x[1][j], 2) < pow(re, 2))
            {
                return true;
            }
        //}  

        else
        {
            return false;
        }
    }
    void range_query(Particle &particle)
    {
        vector<vector<int>> neighbors;
        vector<coordinate> coord(table_size);
        
        neighbors.clear();
        neighbors.resize(table_size);

        int k = 0;
        for (int i = 0; i < rc; ++i)
        {
            for (int j = 0; j < rc; ++j)
            {
                coord[k].x = xmin + h * (0.5 + j);
                coord[k].y = ymin + h * (0.5 + i);
                k++;
            }
        }

        for (int p = 0; p < table_size; ++p)
        {
            for (int i = -1; i <= 1; ++i)
            {
                for (int j = -1; j <= 1; ++j)
                {
                    double x = coord[p].x + j * h;
                    double y = coord[p].y + i * h;
                    if (x >= xmin && x <= xmax && y >= ymin && y <= ymax)
                    {
                        int key = hashkey(x, y);
                        neighbors[p].push_back(key);
                    }
                }
            }
        }

        /*ofstream ss("Results/buckets.dat");
        ss << "Neighbors of grid: \n";
        for (int i = 0; i < neighbors.size(); ++i)
        {
            ss << "bucket[" << i << "] neighbor : " << neighbors[i].size() << "\n";
            for (int j = 0; j < neighbors[i].size(); ++j)
            {
                ss << neighbors[i][j] << '\n';
            }
        }
        ss.close();*/


        int index_neighbor, index_main, size_neighbor, grid_neighbor, grid_main;
        for (int i = 0; i < grid.size(); ++i)
        {
            for (int j = 0; j < neighbors[i].size(); ++j)
            {
                grid_neighbor = neighbors[i][j];
                size_neighbor = grid[grid_neighbor].size();
                //cout << size_neighbor << "\n";
                for (int ii = 0; ii < grid[i].size(); ++ii)
                {
                    index_main = grid[i][ii];
                    for (int jj = 0; jj < size_neighbor; ++jj)
                    {
                        index_neighbor = grid[grid_neighbor][jj];
                        if (distance(particle, index_main, index_neighbor))
                        {
                            particle.m_neighbors[index_main].push_back(index_neighbor);
                        }
                    }
                }
            }
        }
        
    }

    void bruteforce(Particle &particle)
    {
        for (int i = 0; i < N; ++i)
        {
            for (int j = 0; j < N; ++j)
            {
                //if (i != j)
                //{
                    if (pow(particle.m_x[0][i] - particle.m_x[0][j], 2) + pow(particle.m_x[1][i] - particle.m_x[1][j], 2) < pow(re, 2))
                    {
                        particle.m_neighbors[i].push_back(particle.m_index[j]);
                    }
                //}
            }
        }
    }

    void display(const Particle &p)
    {
        ofstream ss("Results/grids.dat");
        ss << "Particles inside grid: \n";
        for (int i = 0; i < grid.size(); ++i)
        {
            ss << "grid[" << i << "] : \n";
            for (int j = 0; j < grid[i].size(); ++j)
            {
                ss << grid[i][j] << '\n';
            }
            ss << '\n';
        }
        ss.close();
    }

    void show_particles(const Particle &p)
    {
        ofstream ss("Results/hash.dat");
        ss << "Neighbors of particle: \n";
        for (int i = 0; i < p.m_neighbors.size(); ++i)
        {
            ss << "particle[" << p.m_index[i] << "] neighbor :" << p.m_neighbors[i].size() << "\n";
            vector<int> cop = p.m_neighbors[i];
            sort(cop.begin(), cop.end());
            for (int j = 0; j < cop.size(); ++j)
            {
                ss << cop[j] << '\n';
            }
        }
        ss.close();
    }
};
#endif