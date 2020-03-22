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
#include <list>






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
    
            N = -cross(B-A, C-A);
            
            double denom = dot(rayCam.direction, N);
            if (std::abs(denom)< 0) return false; // ray parallele
            
            t = dot(C-rayCam.origin, N)/denom;
            if (t<0) return false ; //intersect derriere
            
            P = rayCam.origin + t*rayCam.direction; // Point d'intersection
            N.normalize();

            
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

    if (t_max<0)return false;
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




bool Geometry::intersection(const Ray& rayCam, Vector& P, Vector &N, double &t) const {
    t=1E99;
    bool intersection = false;
        
    //si pas d'intersection avec la racine, pas la peine d'aller plus loin
    if(!bvh.bbox.intersection(rayCam)) return false;
    
    // on crée une liste des fils
    std::list<const BVH*> l;
    //on met le bvh racine dans la liste
    l.push_front(&bvh);
    
    //Tant que notre liste n'est pas vide, on teste le BVH en tête de cette liste et on regarde si un de ses fils a une intersection avec le rayon, qu'on ajoute à la liste si c'est le cas. Si notre test de liste est une feuille, alors c'est qu'il y a intersection et on regarde, pour chaque triangle de cette bbox, s'il y a intersection
    while(!l.empty()){
        const BVH* bvh_test = l.front();
        l.pop_front();
        
        //si mon bvh en cours de test a un fils dont la bbox a aussi une intersection avec mon rayon, je l'ajoute à la liste
        if (bvh_test->arbreg && bvh_test->arbreg->bbox.intersection(rayCam)){
            l.push_back(bvh_test->arbreg);
        }
        //même chose avec le fils droit
        if (bvh_test->arbred && bvh_test->arbred->bbox.intersection(rayCam)){
            l.push_back(bvh_test->arbred);
        }
        
        
        if(!bvh_test->arbreg){// pas besoin de tester si arbre droit, car si pas d'arbre gauche, pas d'arbre droit
            
            // on est sur une feuille, on fait donc le même test qu'avant, sur l'ensemble des éléments (triangles) de la feuille
            for (int i=bvh_test->i0; i<bvh_test->i1; i++){
                
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
            
        }
        
    }




            return intersection;
        }


Bbox Geometry::constr_bbox(int i0, int i1){
    Bbox bbox;
    bbox.bmin =vertices[indices[i0].vtxi];
    bbox.bmax =vertices[indices[i0].vtxi];
    for (int i=i0; i<i1; i++) {// triangle
        for( int k=0; k<3; k++){//dimension
            bbox.bmin[k] = std::min(bbox.bmin[k],vertices[indices[i].vtxi][k]);
            bbox.bmax[k] = std::max(bbox.bmax[k],vertices[indices[i].vtxi][k]);
            
            bbox.bmin[k] = std::min(bbox.bmin[k],vertices[indices[i].vtxj][k]);
            bbox.bmax[k] = std::max(bbox.bmax[k],vertices[indices[i].vtxj][k]);
            
            bbox.bmin[k] = std::min(bbox.bmin[k],vertices[indices[i].vtxk][k]);
            bbox.bmax[k] = std::max(bbox.bmax[k],vertices[indices[i].vtxk][k]);
            
        }
    }
    
    return bbox;
}

void Geometry::constr_bvh(BVH *noeud, int i0, int i1){
    noeud->i0 = i0;
    noeud->i1 = i1;
    noeud->bbox = constr_bbox(i0, i1);
    noeud -> arbreg = NULL;
    noeud -> arbred = NULL;
    
    Vector diag = noeud->bbox.bmax-noeud->bbox.bmin;
    int dimension;
//    choix de la plus grande dimension de la bbox
    if((diag[0]>diag[1]) && (diag[0]>diag[2])){
        dimension = 0;
    }
    else {
        if ((diag[1]>diag[0]) && (diag[1]>diag[2]))
        {
            dimension = 1;
        }
        
        else {
            dimension =2 ;
        }
        
    }
    
    double valeur_coupe = noeud -> bbox.bmin[dimension] + diag[dimension]/2; //la moitier de la bbox sur la dimension de coupe
    
//    pivot du quisk sort
    int pivot = i0-1;

    for(int i=i0; i<i1; i++){
        
        double bary_dim = (vertices[indices[i].vtxi][dimension] + vertices[indices[i].vtxj][dimension] + vertices[indices[i].vtxk][dimension])/3; //barycentre projeté sur la dimension de coupe de la bbox
        //algo du quick sort
        if(bary_dim <valeur_coupe){
            pivot ++;
            
            //on remet directement en ordre le vecteur (TriangleIndices) indices, ainsi, on a juste besoin de travailler sur les indices et pas crée de nouveau vecteur trié
            
            std::swap(indices[i].vtxi, indices[pivot].vtxi);
            std::swap(indices[i].vtxj, indices[pivot].vtxj);
            std::swap(indices[i].vtxk, indices[pivot].vtxk);
            
            std::swap(indices[i].ni, indices[pivot].ni);
            std::swap(indices[i].nj, indices[pivot].nj);
            std::swap(indices[i].nk, indices[pivot].nk);
            
            std::swap(indices[i].uvi, indices[pivot].uvi);
            std::swap(indices[i].uvj,indices[pivot].uvj);
            std::swap(indices[i].uvk, indices[pivot].uvk);
            
            std::swap(indices[i].faceGroup, indices[pivot].faceGroup);
            
        }
    }
    
    // si tous les éléments sont bien ordonnés ou si on a déja une feuille on ne faut rien
    if(pivot < i0|| pivot>=i1-1 || i1==i0+1){
        return;
    }
    // sinon on  crée les arbres fils gauche et droit de notre bvh en appelant recursivement la fonction
    noeud->arbreg= new BVH();
    constr_bvh(noeud->arbreg, i0, pivot+1);
    
    noeud->arbred= new BVH();
    constr_bvh(noeud->arbred, pivot+1, i1);
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
