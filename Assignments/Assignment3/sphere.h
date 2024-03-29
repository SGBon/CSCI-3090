//
//  Framework for a raytracer
//  File: sphere.h
//
//  Created for the Computer Science course "Introduction Computer Graphics"
//  taught at the University of Groningen by Tobias Isenberg.
//
//  Authors:
//    Maarten Everts
//    Jasper van de Gronde
//
//  This framework is inspired by and uses code of the raytracer framework of
//  Bert Freudenberg that can be found at
//  http://isgwww.cs.uni-magdeburg.de/graphik/lehre/cg2/projekt/rtprojekt.html
//

#ifndef SPHERE_H_115209AE
#define SPHERE_H_115209AE

#include "object.h"

class Sphere : public Object
{
public:
    Sphere(Point position,double r) : position(position), r(r) { }

    // calculates roots for quadratic formula
    // coefficents of quadratic formula A,B,C
    // neg controls the +/-
    // true for -, false for +
    double quadf(double A, double B, double C, bool neg);

    virtual Hit intersect(const Ray &ray);

    const Point position;
    const double r;
};

#endif /* end of include guard: SPHERE_H_115209AE */
