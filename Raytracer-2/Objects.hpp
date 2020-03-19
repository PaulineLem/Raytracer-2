
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
    bool intersection (const Ray& rayCam,  Vector& P, Vector& N, double &t) const ; 
    Vector origin;
    double rayon;
};

class Scene {
    public :
    Scene() {};
    void addSphere(const Sphere& s) {spheres.push_back(s);}
    bool intersection (const Ray& r,  Vector& P, Vector& N, int& sphere_id, double& min_t) const;


    std::vector<Sphere> spheres;

};
