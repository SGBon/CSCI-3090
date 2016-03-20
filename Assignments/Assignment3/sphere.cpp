//
//  Framework for a raytracer
//  File: sphere.cpp
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

#include "sphere.h"
#include <iostream>
#include <math.h>

/************************** Sphere **********************************/

Hit Sphere::intersect(const Ray &ray)
{
    /****************************************************
    * RT1.1: INTERSECTION CALCULATION
    *
    * Given: ray, position, r
    * Sought: intersects? if true: *t
    *
    * Insert calculation of ray/sphere intersection here.
    *
    * You have the sphere's center (C) and radius (r) as well as
    * the ray's origin (ray.O) and direction (ray.D).
    *
    * If the ray does not intersect the sphere, return false.
    * Otherwise, return true and place the distance of the
    * intersection point from the ray origin in *t (see example).
    ****************************************************/

    double t; // intersect parameter
    Point p; // intersect point

    // get quadratic formula coefficients
    double A = (ray.D.dot(ray.D));
    double B = 2.0*(ray.D.dot(ray.O - position));
    double C = (ray.O - position).dot(ray.O - position) - pow(r,2.0);
    double root0,root1;

    // test if roots are imaginary/complex
    // if they are, return false
    // otherwise get roots
    double test = pow(B,2.0) - (4.0*A*C);
    if(test < 0){
      return Hit::NO_HIT();
    }else{
      root0 = quadf(A,B,C,false);
      root1 = quadf(A,B,C,true);
    }

    // check for negative roots
    // get smallest root
    if(root0 < 0){
      if(root1 < 0){
        // if both roots are negative
        return Hit::NO_HIT();
      }else{
        // if only root0 is negative
        t = root1;
      }
      // if root0 is positive
    }else{
      if(root1<0){
        // if only root0 is positive
        t = root0;
      }else{
        // if both are positive, get smallest root
        t = fmin(root0,root1);
      }
    }

    /****************************************************
    * RT1.2: NORMAL CALCULATION
    *
    * Given: t, C, r
    * Sought: N
    *
    * Insert calculation of the sphere's normal at the intersection point.
    ****************************************************/

    p = ray.O + t*ray.D;
    Vector N = position - p;
    N.normalize();

    return Hit(t,N);
}

double Sphere::quadf(double A, double B, double C, bool neg){
  if(neg){
    return ((-B - sqrt(pow(B,2.0)-(4.0*A*C)))/(2.0*A));
  }else{
    return ((-B + sqrt(pow(B,2.0)-(4.0*A*C)))/(2.0*A));
  }
}
