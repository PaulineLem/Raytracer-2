
#include <stdio.h>
#include <vector>




class Ray {
    public :
    Ray(const Vector& o, const Vector& d ) : origin(o), direction(d) {};
    Vector origin, direction;
};

class Objet {
    public :
    Objet(){};
    virtual bool intersection (const Ray& rayCam,  Vector& P, Vector& N, double &t) const=0 ;
    Vector albedo;
    bool mirror;
    bool transparent;

};

class Sphere :public Objet{
    public :
    Sphere(const Vector& o, const double& r, const Vector &color, bool is_mirror = false, bool is_transp = false) : origin(o), rayon(r){
        albedo = color;
        mirror = is_mirror;
        transparent = is_transp;
    };
    bool intersection (const Ray& rayCam,  Vector& P, Vector& N, double &t) const ;
    Vector origin;
    double rayon;


};

class Triangle : public Objet {
    public :
    Triangle(const Vector& A, const Vector& B,const Vector& C, const Vector &color, bool is_mirror = false, bool is_transp = false) : A(A), B(B), C(C){
        albedo = color;
        mirror = is_mirror;
        transparent = is_transp;
    };
    bool intersection (const Ray& rayCam,  Vector& P, Vector& N, double &t) const ;
    const Vector &A;
    const Vector &B;
    const Vector &C;

};

class TriangleIndices {
public:
    TriangleIndices(int vtxi = -1, int vtxj = -1, int vtxk = -1, int ni = -1, int nj = -1, int nk = -1, int uvi = -1, int uvj = -1, int uvk = -1) : vtxi(vtxi), vtxj(vtxj), vtxk(vtxk), uvi(uvi), uvj(uvj), uvk(uvk), ni(ni), nj(nj), nk(nk) {
    };
    int vtxi, vtxj, vtxk;
    int uvi, uvj, uvk;
    int ni, nj, nk;
    int faceGroup;
};

class Bbox{
public:
    Bbox(){};
    Bbox (const Vector& bmin, const Vector& bmax): bmin(bmin), bmax(bmax){}
    
    bool intersection(const Ray& r) const;
    
    Vector bmin, bmax;
};



class Geometry : public Objet {
public:
    Geometry() {};
    Geometry(const char* obj, double scaling,const Vector& offset, const Vector &color, bool is_mirror = false, bool is_transp = false) {
        
        albedo = color;
        mirror =  is_mirror;
        transparent = is_transp;
        
        readOBJ(obj);
        
        for (int i = 0; i < vertices.size(); i++) {
            vertices[i] = scaling * vertices[i] + offset;
        }
        
        build_bb(0, indices.size());
        
    };
    
    
    
    void readOBJ(const char* obj);
    Bbox build_bb( int i0, int i1);

    bool intersection(const Ray& rayCam, Vector& P, Vector &N, double &t) const;
    
    std::vector<TriangleIndices> indices;
    std::vector<Vector> vertices;
    std::vector<Vector> normals;
    std::vector<Vector> uvs; // Vector en 3D mais on n'utilise que 2 composantes
    std::vector<Vector> vertexcolors;

    
    
    private :
    Bbox bbox;
    std::vector<int> w,h;
};



class Scene {
    public :
    Scene() {};

    void addSphere(const Sphere& s) {objets.push_back((Objet*)&s);}
    void addTriangle(const Triangle& t) {objets.push_back((Objet*)&t);}

    bool intersection (const Ray& r,  Vector& P, Vector& N, int& sphere_id, double& min_t) const;
    double lumIntensite;
    Sphere *lumiere;
    std::vector<Objet*> objets;

};
