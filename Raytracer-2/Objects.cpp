//
//  Objets.cpp
//  Raytracer-2
//
//  Created by Pauline Lemeille on 19/03/2020.
//  Copyright © 2020 Pauline Lemeille. All rights reserved.
//
#include "Vector.hpp"
#include "Objects.hpp"
#include <vector>


bool Sphere::intersection(const Ray& rayCam,  Vector& P, Vector& N, double &t) const {
    
        double a = 1;
       double b = 2*dot(rayCam.direction, rayCam.origin - origin);
       double c = ( rayCam.origin - origin).getNorm2() -rayon*rayon;
       
       double delta = b*b - 4*a*c;
       
       if (delta<0) return false ;
       double t1 = (-b - sqrt(delta))/ (2*a);
       double t2 = (-b + sqrt(delta))/ (2*a);
    if(t2<0){
        return false;

    }
       if (t1>0)
           t=t1;
       else t=t2;
    
    P = rayCam.origin + t * rayCam.direction;
    N = (P-origin).getNormalized();
    return true;
}

bool Scene::intersection (const Ray& r,  Vector& P, Vector& N, int& sphere_id, double &t_min) const{

    bool has_inter = false;
    t_min = 1E99;
    
    for (int i=0; i<spheres.size(); i++) {
        Vector localP, localN;
        double localt;
        bool local_has_inter = spheres[i].intersection(r, localP, localN, localt);
        
        if (local_has_inter){
            has_inter = true;
            if (localt<t_min){
                t_min =localt;
                P= localP;
                N=  localN;
                sphere_id = i;
            }
        }
    }
    return has_inter;
};
