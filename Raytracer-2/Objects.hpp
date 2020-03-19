
#include <stdio.h>
#include <vector>




class Ray {
    public :
    Ray(const Vector& o, const Vector& d ) : origin(o), direction(d) {};
    Vector origin, direction;
};

class Sphere {
    public :
    Sphere(const Vector& o, const double& r ) : origin(o), rayon(r) {};
    bool intersection (const Ray& rayCam,  Vector& P, Vector& N) const ; 
    Vector origin;
    double rayon;
};

