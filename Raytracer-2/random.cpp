#include "Vector.hpp"
#include <vector>
#include <omp.h>

#include <random>

std::default_random_engine engine[8];
std::uniform_real_distribution<double> distrib(0,1);

Vector random_vect() {
    double r1, r2;
    r1 = distrib(engine[omp_get_thread_num()]);
    r2 = distrib(engine[omp_get_thread_num()]);
    return Vector(r1,r2);
}

Vector randomcos(const Vector &N) {
    
    Vector Randvec = random_vect();

    Vector V; //repere local
    V[0] = cos(2 * M_PI * Randvec[0])* sqrt(1-Randvec[1]);
    V[1] = sin(2 * M_PI * Randvec[0])* sqrt(1-Randvec[1]);
    V[2] = sqrt(Randvec[1]);

    
    Vector aleatoire (distrib(engine[omp_get_thread_num()]) - 0.5, distrib(engine[omp_get_thread_num()])-0.5, distrib(engine[omp_get_thread_num()])-0.5);
    Vector tangent1= cross(N, aleatoire); // perpendiculaire a N et aleatoire donc tangent a la surface
    tangent1.normalize();
    Vector tangent2 = cross(tangent1, N); // perpendiculaire à t1 et N
    
    return V[0]*tangent1 + V[1] * tangent2+ V[2] *N;
}




