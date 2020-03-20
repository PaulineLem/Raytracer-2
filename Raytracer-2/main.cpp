//
//  main.cpp
//  Raytracer-2
//
//  Created by Pauline Lemeille on 19/03/2020.
//  Copyright © 2020 Pauline Lemeille. All rights reserved.
//


#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

// Raytracer.cpp : Defines the entry point for the console application.
#define _CRT_SECURE_NO_WARNINGS // for Visual Studio 2017 (maybe 2015 as well)

#include <iostream>
#include <random>

#include <vector>
#include "Vector.hpp"
#include "Objects.hpp"




void save_image(const char* filename, const unsigned char* tableau, int w, int h) { // (0,0) is top-left corner
    
    FILE *f;
    
    int filesize = 54 + 3 * w*h;
    
    unsigned char bmpfileheader[14] = { 'B','M', 0,0,0,0, 0,0, 0,0, 54,0,0,0 };
    unsigned char bmpinfoheader[40] = { 40,0,0,0, 0,0,0,0, 0,0,0,0, 1,0, 24,0 };
    unsigned char bmppad[3] = { 0,0,0 };
    
    bmpfileheader[2] = (unsigned char)(filesize);
    bmpfileheader[3] = (unsigned char)(filesize >> 8);
    bmpfileheader[4] = (unsigned char)(filesize >> 16);
    bmpfileheader[5] = (unsigned char)(filesize >> 24);
    
    bmpinfoheader[4] = (unsigned char)(w);
    bmpinfoheader[5] = (unsigned char)(w >> 8);
    bmpinfoheader[6] = (unsigned char)(w >> 16);
    bmpinfoheader[7] = (unsigned char)(w >> 24);
    bmpinfoheader[8] = (unsigned char)(h);
    bmpinfoheader[9] = (unsigned char)(h >> 8);
    bmpinfoheader[10] = (unsigned char)(h >> 16);
    bmpinfoheader[11] = (unsigned char)(h >> 24);
    
    f = fopen(filename, "wb");
    fwrite(bmpfileheader, 1, 14, f);
    fwrite(bmpinfoheader, 1, 40, f);
    unsigned char *row = new unsigned char[w * 3];
    for (int i = 0; i < h; i++)
    {
        for (int j = 0; j < w; j++) {
            row[j * 3] = tableau[(w*(h - i - 1) * 3) + j * 3+2];
            row[j * 3+1] = tableau[(w*(h - i - 1) * 3) + j * 3+1];
            row[j * 3+2] = tableau[(w*(h - i - 1) * 3) + j * 3];
        }
        fwrite(row, 3, w, f);
        fwrite(bmppad, 1, (4 - (w * 3) % 4) % 4, f);
    }
    fclose(f);
    delete[] row;
}

Vector getColor(const Ray rayCam, const Scene s, Vector lumOrigin, int nb_rebond){
    double eps(0.001);
    Vector pixelColor(0.,0.,0.);
    Vector P, N;
    double t;
    int sphere_id;
    bool intersect = s.intersection(rayCam, P, N, sphere_id, t);
    if (!intersect) return pixelColor;
    
//    if (s.spheres[sphere_id].mirror && nb_rebond>0) {
//        pixelColor = getColor(ray_refl, s, lumOrigin, nb_rebon)
//        
//    }
    Vector PL =(lumOrigin - P);
    double d2 = PL.getNorm2();

    
//Gestion des ombres portées
    Ray new_ray(P+eps*N, PL.getNormalized());
    Vector newP, newN;
    double new_t;
    int new_sphere_id;
    
    bool new_intersect = s.intersection(new_ray, newP, newN, new_sphere_id, new_t);
    if (new_intersect) {
        Vector newPP =(P - newP);
        double newd2 = newPP.getNorm2();
        if (newd2<d2){
            return pixelColor;
        }
        
    }
    
    Vector intensite;
    intensite = s.spheres[sphere_id].albedo * 1000 * std::max(0., dot(PL.getNormalized(), N))/d2;
        
    pixelColor = Vector(std::min(255., std::max(0., intensite[0])), std::min(255., std::max(0., intensite[1])),std::min(255., std::max(0., intensite[2])));
    
       return pixelColor;
}


int main() {
    
    int W = 1024;
    int H = 1024;
    double fov = 60 * M_PI / 180;
    std::vector<unsigned char> image(W*H * 3);
    
    Vector cameraPos(0., 0., 55);
    Sphere sphere_1(Vector(0,0, 0),10, Vector(255., 255., 255));
    Sphere sphere_2(Vector(0,-1000, 0),990, Vector(0., 0., 255));

    Scene s;
    Vector lumOrigin(-10, 20, 40);
    
    s.addSphere(sphere_1);
    s.addSphere(sphere_2);


  
    // Lance des thread pour la boucle for suivante
    for (int i = 0; i < H; i++) {
        
        for (int j = 0; j < W; j++) {
            Vector pixColor=(0.,0.,0.);
            
            Vector direction(j-W/2 +0.5 , i-H/2+0.5, -H/ (2*tan(fov/2)));
            direction.normalize();
            Ray rayCam(cameraPos, direction);
            
            pixColor=getColor(rayCam, s, lumOrigin, 5);

            image[((H-i-1)*W + j) * 3 + 0] = pixColor[0];
            image[((H-i-1)*W + j) * 3 + 1] = pixColor[1];
            image[((H-i-1)*W + j) * 3 + 2] = pixColor[2];
        }
    }
    save_image("seance2-ombres-portées-sans vruit.bmp",&image[0], W, H);

    return 0;
}
