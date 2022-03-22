#include <vector>

class Particle{
    public:

    // Particle Coordinates
    double x;
    double y;
    double z;

    // Particle Radius
    double r;

    // Particle Temperature
    double T;

    // Particle Neighbor
    std::vector<int> neighbor;

    Particle(double _x, double _y, double _z, double _r, double _T){

        x = _x;
        y = _y;
        z = _z;
        r = _r;
        T = _T;

    }

};