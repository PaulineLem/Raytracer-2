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
#include <omp.h>
#include <random>
#include <map>




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

bool Triangle::intersection(const Ray& rayCam,  Vector& P, Vector& N, double &t) const {
    
            N = cross(A-B, C-A);
            N.normalize();
            
            double denom = dot(rayCam.direction, N);
            if (std::abs(denom)< 0) return false; // ray parallele
            
            t = dot(C-rayCam.origin, N)/denom;
            if (t<0) return false ; //intersect derriere
            
            P = rayCam.origin + t*rayCam.direction; // Point d'intersection
            
            double APAB = dot(P-A, B-A);
            double ACAB = dot(C-A, B-A);
            double ABAB = dot(B-A, B-A);
            double APAC = dot(P-A, C-A);
            double ACAC = dot(C-A, C-A);
            double det = ABAB*ACAC-ACAB*ACAB;
    
            double alpha, beta, gamma;
            
            //Systeme de Crammer
            beta = (APAB*ACAC-APAC*ACAB)/det;
            gamma = (ABAB*APAC -ACAB*APAB)/det;
            alpha = 1-beta-gamma;
            if (beta<0 || beta>1) return false;
            if (gamma<0 || gamma>1) return false;
            if (alpha<0 || alpha>1) return false;
            return true;
        
}

bool Bbox::intersection(const Ray& r) const {
    double t_1_x = (bmin[0]-r.origin[0])/r.direction[0];
    double t_2_x = (bmax[0]-r.origin[0])/r.direction[0];
    double t_min_x = std::min(t_1_x, t_2_x);
    double t_max_x = std::max(t_1_x, t_2_x);
    
    double t_1_y = (bmin[1]-r.origin[1])/r.direction[1];
    double t_2_y = (bmax[1]-r.origin[1])/r.direction[1];
    double t_min_y = std::min(t_1_y, t_2_y);
    double t_max_y = std::max(t_1_y, t_2_y);
    
    double t_1_z = (bmin[2]-r.origin[2])/r.direction[2];
    double t_2_z = (bmax[2]-r.origin[2])/r.direction[2];
    double t_min_z = std::min(t_1_z, t_2_z);
    double t_max_z = std::max(t_1_z, t_2_z);
    
    double t_max = std::min(std::min(t_max_x,t_max_y),t_max_z) ;
    double t_min = std::max(std::max(t_min_x, t_min_y), t_min_z) ;

    if(t_max<0) return false;
    
    if (t_max - t_min > 0) return true ;
    return false;
}

void Geometry::readOBJ(const char* obj) {
    
    char matfile[255];
    char grp[255];
    
    FILE* f;
    f = fopen(obj, "r");
    
    std::map<std::string, int> groupNames;
    int curGroup = -1;
    while (!feof(f)) {
        char line[255];
        if (!fgets(line, 255, f)) break;
        
        std::string linetrim(line);
        linetrim.erase(linetrim.find_last_not_of(" \r\t") + 1);
        strcpy(line, linetrim.c_str());
        
        if (line[0] == 'u' && line[1] == 's') {
            sscanf(line, "usemtl %[^\n]\n", grp);
            if (groupNames.find(std::string(grp)) != groupNames.end()) {
                curGroup = groupNames[std::string(grp)];
            }
            else {
                curGroup = groupNames.size();
                groupNames[std::string(grp)] = curGroup;
            }
        }
        if (line[0] == 'm' && line[1] == 't' && line[2] == 'l') {
            sscanf(line, "mtllib %[^\n]\n", matfile);
        }
        if (line[0] == 'v' && line[1] == ' ') {
            Vector vec;
            Vector col;
            if (sscanf(line, "v %lf %lf %lf %lf %lf %lf\n", &vec[0], &vec[2], &vec[1], &col[0], &col[1], &col[2]) == 6) {
                vertices.push_back(vec);
                vertexcolors.push_back(col);
            }
            else {
                sscanf(line, "v %lf %lf %lf\n", &vec[0], &vec[2], &vec[1]);  // helmet
                //vec[2] = -vec[2]; //car2
                vertices.push_back(vec);
            }
        }
        if (line[0] == 'v' && line[1] == 'n') {
            Vector vec;
            sscanf(line, "vn %lf %lf %lf\n", &vec[0], &vec[2], &vec[1]); //girl
            normals.push_back(vec);
        }
        if (line[0] == 'v' && line[1] == 't') {
            Vector vec;
            sscanf(line, "vt %lf %lf\n", &vec[0], &vec[1]);
            uvs.push_back(vec);
        }
        if (line[0] == 'f') {
            TriangleIndices t;
            int i0, i1, i2, i3;
            int j0, j1, j2, j3;
            int k0, k1, k2, k3;
            int nn;
            
            char* consumedline = line + 1;
            int offset;
            t.faceGroup = curGroup;
            nn = sscanf(consumedline, "%u/%u/%u %u/%u/%u %u/%u/%u%n", &i0, &j0, &k0, &i1, &j1, &k1, &i2, &j2, &k2, &offset);
            if (nn == 9) {
                if (i0 < 0) t.vtxi = vertices.size() + i0; else    t.vtxi = i0 - 1;
                if (i1 < 0) t.vtxj = vertices.size() + i1; else    t.vtxj = i1 - 1;
                if (i2 < 0) t.vtxk = vertices.size() + i2; else    t.vtxk = i2 - 1;
                if (j0 < 0) t.uvi = uvs.size() + j0; else    t.uvi = j0 - 1;
                if (j1 < 0) t.uvj = uvs.size() + j1; else    t.uvj = j1 - 1;
                if (j2 < 0) t.uvk = uvs.size() + j2; else    t.uvk = j2 - 1;
                if (k0 < 0) t.ni = normals.size() + k0; else    t.ni = k0 - 1;
                if (k1 < 0) t.nj = normals.size() + k1; else    t.nj = k1 - 1;
                if (k2 < 0) t.nk = normals.size() + k2; else    t.nk = k2 - 1;
                
                indices.push_back(t);
            }
            else {
                nn = sscanf(consumedline, "%u/%u %u/%u %u/%u%n", &i0, &j0, &i1, &j1, &i2, &j2, &offset);
                if (nn == 6) {
                    if (i0 < 0) t.vtxi = vertices.size() + i0; else    t.vtxi = i0 - 1;
                    if (i1 < 0) t.vtxj = vertices.size() + i1; else    t.vtxj = i1 - 1;
                    if (i2 < 0) t.vtxk = vertices.size() + i2; else    t.vtxk = i2 - 1;
                    if (j0 < 0) t.uvi = uvs.size() + j0; else    t.uvi = j0 - 1;
                    if (j1 < 0) t.uvj = uvs.size() + j1; else    t.uvj = j1 - 1;
                    if (j2 < 0) t.uvk = uvs.size() + j2; else    t.uvk = j2 - 1;
                    indices.push_back(t);
                }
                else {
                    nn = sscanf(consumedline, "%u %u %u%n", &i0, &i1, &i2, &offset);
                    if (nn == 3) {
                        if (i0 < 0) t.vtxi = vertices.size() + i0; else    t.vtxi = i0 - 1;
                        if (i1 < 0) t.vtxj = vertices.size() + i1; else    t.vtxj = i1 - 1;
                        if (i2 < 0) t.vtxk = vertices.size() + i2; else    t.vtxk = i2 - 1;
                        indices.push_back(t);
                    }
                    else {
                        nn = sscanf(consumedline, "%u//%u %u//%u %u//%u%n", &i0, &k0, &i1, &k1, &i2, &k2, &offset);
                        if (i0 < 0) t.vtxi = vertices.size() + i0; else    t.vtxi = i0 - 1;
                        if (i1 < 0) t.vtxj = vertices.size() + i1; else    t.vtxj = i1 - 1;
                        if (i2 < 0) t.vtxk = vertices.size() + i2; else    t.vtxk = i2 - 1;
                        if (k0 < 0) t.ni = normals.size() + k0; else    t.ni = k0 - 1;
                        if (k1 < 0) t.nj = normals.size() + k1; else    t.nj = k1 - 1;
                        if (k2 < 0) t.nk = normals.size() + k2; else    t.nk = k2 - 1;
                        indices.push_back(t);
                    }
                }
            }
            
            
            consumedline = consumedline + offset;
            
            while (true) {
                if (consumedline[0] == '\n') break;
                if (consumedline[0] == '\0') break;
                nn = sscanf(consumedline, "%u/%u/%u%n", &i3, &j3, &k3, &offset);
                TriangleIndices t2;
                t2.faceGroup = curGroup;
                if (nn == 3) {
                    if (i0 < 0) t2.vtxi = vertices.size() + i0; else    t2.vtxi = i0 - 1;
                    if (i2 < 0) t2.vtxj = vertices.size() + i2; else    t2.vtxj = i2 - 1;
                    if (i3 < 0) t2.vtxk = vertices.size() + i3; else    t2.vtxk = i3 - 1;
                    if (j0 < 0) t2.uvi = uvs.size() + j0; else    t2.uvi = j0 - 1;
                    if (j2 < 0) t2.uvj = uvs.size() + j2; else    t2.uvj = j2 - 1;
                    if (j3 < 0) t2.uvk = uvs.size() + j3; else    t2.uvk = j3 - 1;
                    if (k0 < 0) t2.ni = normals.size() + k0; else    t2.ni = k0 - 1;
                    if (k2 < 0) t2.nj = normals.size() + k2; else    t2.nj = k2 - 1;
                    if (k3 < 0) t2.nk = normals.size() + k3; else    t2.nk = k3 - 1;
                    indices.push_back(t2);
                    consumedline = consumedline + offset;
                    i2 = i3;
                    j2 = j3;
                    k2 = k3;
                }
                else {
                    nn = sscanf(consumedline, "%u/%u%n", &i3, &j3, &offset);
                    if (nn == 2) {
                        if (i0 < 0) t2.vtxi = vertices.size() + i0; else    t2.vtxi = i0 - 1;
                        if (i2 < 0) t2.vtxj = vertices.size() + i2; else    t2.vtxj = i2 - 1;
                        if (i3 < 0) t2.vtxk = vertices.size() + i3; else    t2.vtxk = i3 - 1;
                        if (j0 < 0) t2.uvi = uvs.size() + j0; else    t2.uvi = j0 - 1;
                        if (j2 < 0) t2.uvj = uvs.size() + j2; else    t2.uvj = j2 - 1;
                        if (j3 < 0) t2.uvk = uvs.size() + j3; else    t2.uvk = j3 - 1;
                        consumedline = consumedline + offset;
                        i2 = i3;
                        j2 = j3;
                        indices.push_back(t2);
                    }
                    else {
                        nn = sscanf(consumedline, "%u//%u%n", &i3, &k3, &offset);
                        if (nn == 2) {
                            if (i0 < 0) t2.vtxi = vertices.size() + i0; else    t2.vtxi = i0 - 1;
                            if (i2 < 0) t2.vtxj = vertices.size() + i2; else    t2.vtxj = i2 - 1;
                            if (i3 < 0) t2.vtxk = vertices.size() + i3; else    t2.vtxk = i3 - 1;
                            if (k0 < 0) t2.ni = normals.size() + k0; else    t2.ni = k0 - 1;
                            if (k2 < 0) t2.nj = normals.size() + k2; else    t2.nj = k2 - 1;
                            if (k3 < 0) t2.nk = normals.size() + k3; else    t2.nk = k3 - 1;
                            consumedline = consumedline + offset;
                            i2 = i3;
                            k2 = k3;
                            indices.push_back(t2);
                        }
                        else {
                            nn = sscanf(consumedline, "%u%n", &i3, &offset);
                            if (nn == 1) {
                                if (i0 < 0) t2.vtxi = vertices.size() + i0; else    t2.vtxi = i0 - 1;
                                if (i2 < 0) t2.vtxj = vertices.size() + i2; else    t2.vtxj = i2 - 1;
                                if (i3 < 0) t2.vtxk = vertices.size() + i3; else    t2.vtxk = i3 - 1;
                                consumedline = consumedline + offset;
                                i2 = i3;
                                indices.push_back(t2);
                            }
                            else {
                                consumedline = consumedline + 1;
                            }
                        }
                    }
                }
            }
            
        }
        
        
    }
    fclose(f);
    
}


Bbox Geometry::build_bb(int i0, int i1)
{
    Bbox bb;
    bb.bmin =vertices[indices[i0].vtxi];
    bb.bmax =vertices[indices[i0].vtxi];
    for (int i=i0; i<i1; i++) {// triangle
        for( int k=0; k<3; k++){//dimension
            bb.bmin[k] = std::min(bb.bmin[k],vertices[indices[i].vtxi][k]);
            bb.bmax[k] = std::max(bb.bmax[k],vertices[indices[i].vtxi][k]);
            
            bb.bmin[k] = std::min(bb.bmin[k],vertices[indices[i].vtxj][k]);
            bb.bmax[k] = std::max(bb.bmax[k],vertices[indices[i].vtxj][k]);
            
            bb.bmin[k] = std::min(bb.bmin[k],vertices[indices[i].vtxk][k]);
            bb.bmax[k] = std::max(bb.bmax[k],vertices[indices[i].vtxk][k]);
            
        }
    }
    
    return bb;
}

    bool Geometry::intersection(const Ray& rayCam, Vector& P, Vector &N, double &t) const {
        
        if(!bbox.intersection(rayCam)) return false;

        t=1E99;
        bool intersection = false;

        for (int i=0; i<indices.size(); i++){
            
            int a = indices[i].vtxi ;
             int b = indices[i].vtxj;
             int c = indices[i].vtxk;
             Triangle tri(vertices[a], vertices[b], vertices[c], albedo, mirror, transparent);
             Vector localP, localN;
             double localt;
             if (tri.intersection(rayCam, localP, localN, localt)) {
                
                 intersection = true;
                 if(localt<t){
                     t=localt;
                     P=localP;
                     N=localN;
            
        }
             }
        }
            return intersection;
        }
                 
    



bool Scene::intersection (const Ray& r,  Vector& P, Vector& N, int& sphere_id, double &t_min) const{

    bool has_inter = false;
    t_min = 1E99;
    
    for (int i=0; i<objets.size(); i++) {
        Vector localP, localN;
        double localt;
        bool local_has_inter = objets[i]->intersection(r, localP, localN, localt);
        
        if (local_has_inter){
            has_inter = true;
            if (localt<t_min){
                t_min =localt;
                P= localP;
                N= localN;
                sphere_id = i;
            }
        }
    }
    
    return has_inter;
};
