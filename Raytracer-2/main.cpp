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
#include "random.hpp"






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


Vector getColor(const Ray rayCam, const Scene s,  int nb_rebond){
    double eps(0.001);
    Vector pixelColor(0.,0.,0.);
    Vector P, N, albedo;
    double t;
    int sphere_id;
    
    if (nb_rebond ==0) return pixelColor;

    bool intersect = s.intersection(rayCam, P, N, sphere_id, t, albedo);
    if (!intersect) return pixelColor;

    
    if (s.objets[sphere_id]->mirror ) {
            Vector dir_refl(rayCam.direction - 2*dot(N, rayCam.direction)*N);
            Ray ray_refl(P+eps*N, dir_refl);
            pixelColor = getColor(ray_refl, s, nb_rebond-1);
        }
    else {
        if (s.objets[sphere_id]->transparent ) {
            double n1=1;
            double n2=1.5; //verre
            Vector Ntransp(N);
            if(dot(rayCam.direction,N)>0) {
                //si le rayon sort de la sphere
                n1=1.5;
                n2=1;
                Ntransp = -N;
            }
            double radical = 1-pow((n1/n2),2)*(1-pow((dot(Ntransp , rayCam.direction)),2));
            if (radical > 0) {
                Vector dir_refr = (n1/n2)*(rayCam.direction - dot(rayCam.direction, Ntransp)*Ntransp) - Ntransp * sqrt (radical);
                Ray refrRay(P - eps*Ntransp, dir_refr);
                pixelColor = getColor(refrRay, s, nb_rebond-1);
            }
            
        }
        else {
            
                //                    Eclairage direct
            //On ne génère les rayons que vers la source de lumière
                Vector LumOr = s.lumiere->origin;
                //P point d'intersection
                Vector LP = P-LumOr;
                LP.normalize();
                Vector randSdir = randomcos(LP);
                Vector Pi = LumOr + randSdir*s.lumiere->rayon; //point aleatoire sur la sphere (changement de variable dans l'equation du rendu)
                Vector wi = Pi-P; // dir aleatoire
                wi.normalize();
        
                double d2 = (Pi-P).getNorm2();
                double costheta, costhetaprime, costhetasecond;
                 costheta = std::max(0., dot(N, wi));
                 costhetaprime = dot(randSdir,-wi );
                 costhetasecond = dot(LP,randSdir);
            
                Ray ray_lum(P + eps * N, wi);
                Vector P_lum, N_lum, albedo_lum;
                int sphere_id_lum;
                double t_lum;
            
                bool intersect_lum= s.intersection(ray_lum, P_lum, N_lum, sphere_id_lum, t_lum, albedo_lum);

                if ( intersect_lum && t_lum*t_lum < d2*0.99 ){
                    pixelColor = Vector(0,0,0);
                }
                
                else {
                    pixelColor = (s.lumIntensite/(4*M_PI*d2)*costheta*costhetaprime/costhetasecond)* albedo;
                }

            
//             contribution indirect (meme chose que miroir mais rayons aleatoires

            Vector dir_alea= randomcos(N);
            Ray ray_alea(P+eps*N, dir_alea);
            pixelColor += getColor(ray_alea, s, nb_rebond-1)*albedo ;
        }
        
    }
//    if (sphere_id ==6 ){ std::cout<<"jesuisla"<< pixelColor[0]<<"\n" ;}

    return pixelColor;

}


int main() {
    
    int W = 1024;
    int H = 1024;
    double fov = 60 * M_PI / 180;
    std::vector<unsigned char> image(W*H * 3);
    int nb_rayon = 500;
    int focus_cam = 25;
    
    
    Vector cameraPos(0., 0., 0. );
    
    Sphere sphere_lum(Vector(15, 70, 30), 15,Vector (1,1,1));

    Sphere sphere_1(Vector(10, 10, -focus_cam-5 ),5, Vector(1., 1.,1.), true);
    Sphere sphere_7(Vector(10, -5, -focus_cam-5 ),5, Vector(1., 1.,1.), false, true);

    Sphere sphere_2(Vector(0,-2000-20, 0),2000, Vector (0.7,0.7,0.7)); //ground
    Sphere sphere_3(Vector(0,2000+100, 0),2000, Vector (0.7,0.7,0.7)); //ceiling
    Sphere sphere_4(Vector(-2000-50,0, 0),2000, Vector (0.,1.,0)); // left wall
    Sphere sphere_5(Vector(2000+50,0, 0),2000,  Vector (0.,0,1.)); // right wall
    Sphere sphere_6(Vector(0, 0, -2000-100),2000, Vector (0.,0.8,0.8)); // back wall
    
//    Triangle triangle_1(Vector(-10, -10, -focus_cam), Vector(10, -10, -focus_cam), Vector(0, 10, -focus_cam), Vector(1, 0, 0));
    
    Geometry geometry_1("Beautiful Girl.obj", 10,Vector(-10, -10, -focus_cam), Vector (1,1,1));


    




    Scene s;
    s.addSphere(sphere_lum);
    s.addSphere(sphere_1);
    s.addSphere(sphere_7);
    s.addSphere(sphere_2);
    s.addSphere(sphere_3);
    s.addSphere(sphere_4);
    s.addSphere(sphere_5);
    s.addSphere(sphere_6);
//    s.addTriangle(triangle_1);
    s.addGeometry(geometry_1);
    
    s.lumiere = &sphere_lum;
    s.lumIntensite = 8000000000;




  #pragma omp parallel for

    // Lance des thread pour la boucle for suivante
    for (int i = 0; i < H; i++) {
        
        for (int j = 0; j < W; j++) {
            Vector pixColor=(0.,0.,0.);
            
            for (int k=0; k<nb_rayon; k++) {
                double  dx, dy;
                
                //aiti-aisling
                Vector rand = random_vect();

                dx = cos(2 * M_PI * rand[1])* sqrt(-2*log(rand[0])) * 0.5;
                dy = sin(2 * M_PI * rand[1])* sqrt(-2*log(rand[0])) *0.5;
                
                

                Vector direction(j-W/2 +0.5 +dx , i-H/2+0.5 +dy, -H/ (2*tan(fov/2)));
                direction.normalize();
                
                Vector rand2 = random_vect();
                Vector destination = cameraPos + focus_cam * direction;
                Vector origine = cameraPos +Vector((rand2[0] -0.5) *0.5 , (rand2[1] -0.5) *0.5 , 0);
                Ray rayCam(origine, (destination-origine).getNormalized());
        
                pixColor+=getColor(rayCam, s, 5) ;
            }
            pixColor = pixColor/nb_rayon;
            
            image[((H-i-1)*W + j) * 3 + 0] = std::min(255., std::max(0.,pow(pixColor[0], 1/2.2)));
            image[((H-i-1)*W + j) * 3 + 1] = std::min(255., std::max(0.,pow(pixColor[1], 1/2.2)));
            image[((H-i-1)*W + j) * 3 + 2] = std::min(255., std::max(0.,pow(pixColor[2], 1/2.2)));
        }
    }
    save_image("seance5-beautiful-girl-texture-500r.bmp",&image[0], W, H);

    return 0;
}
