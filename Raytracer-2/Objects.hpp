
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
    virtual bool intersection (const Ray& rayCam,  Vector& P, Vector& N, double &t, Vector &color) const=0 ;
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
    bool intersection (const Ray& rayCam,  Vector& P, Vector& N, double &t, Vector &color) const ;
    Vector origin;
    double rayon;


};

class Triangle : public Objet {
    public :
    Triangle(const Vector& A, const Vector& B,const Vector& C, const Vector &color, bool is_mirror = false, bool is_transp = false) : A(A), B(B), C(C){
        albedo = color;
        mirror = is_mirror;
        transparent = is_transp;
    }
    bool intersection(const Ray& rayCam, Vector& P, Vector &N, double &t, Vector &color) const;
    bool intersection(const Ray& rayCam, Vector& P, Vector &N, double &t, double &alpha, double &beta, double &gamma) const ;

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

class BVH{
public:

    int i0, i1;
    Bbox bbox;
    
    BVH *arbreg, *arbred;

    
};


class Geometry : public Objet {
public:
    Geometry() {};
    Geometry(const char* obj, double scaling,const Vector& offset, const Vector &color, bool is_mirror = false, bool is_transp = false) {
        
        albedo = color;
        mirror =  is_mirror;
        transparent = is_transp;
        
        readOBJ(obj);
        
        //on applique le scaling et l'offset
        for (int i = 0; i < vertices.size(); i++) {
            vertices[i] = scaling * vertices[i] + offset;
        }
        
        //on construit la bvh de l'objet
        constr_bvh(&bvh, 0, indices.size());
        

    };
    
    
    
    void readOBJ(const char* obj);
    bool intersection(const Ray& rayCam, Vector& P, Vector &N, double &t, Vector &color) const;
    Bbox constr_bbox(int i0, int i1);
    void constr_bvh(BVH *noeud, int i0, int i1);
    void add_texture(const char* filename);

    
    std::vector<TriangleIndices> indices;
    std::vector<Vector> vertices;
    std::vector<Vector> normals;
    std::vector<Vector> uvs; // Vector en 3D mais on n'utilise que 2 composantes
    std::vector<Vector> vertexcolors;

    
    
    private :
    BVH bvh;
    std::vector<int> w,h;
    std::vector<std::vector<unsigned char> > textures;

};



class Scene {
    public :
    Scene() {};

    void addSphere(const Sphere& s) {objets.push_back((Objet*)&s);}
    void addTriangle(const Triangle& t) {objets.push_back((Objet*)&t);}
    void addGeometry(const Geometry& g) {objets.push_back((Objet*)&g);}


    bool intersection (const Ray& r,  Vector& P, Vector& N, int& sphere_id, double& min_t, Vector &color) const;
    double lumIntensite;
    Sphere *lumiere;
    std::vector<Objet*> objets;

};
