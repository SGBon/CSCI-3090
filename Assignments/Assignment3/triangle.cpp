#include "triangle.h"

Hit Triangle::intersect(const Ray &ray){
  double t; // parameter
  Vector N; // hit point

  // get 2 vectors on plane
  Vector U = V2 - V1;
  Vector V = V3 - V1;

  N = U.cross(V);

  double dndot = ray.D.dot(N);
  // check if ray is parallel to plane
  if(dndot == 0.0){
    return Hit::NO_HIT();
  }


  t = ((V1 - ray.O).dot(N))/dndot;
  // check if intersection is behind origin
  if(t < 0.0){
    return Hit::NO_HIT();
  }

  Point sect = ray.O + t*ray.D; //intersect point
  Vector w = sect - V1; // vector in plane

  //v(s,t) = v0 + su + tv
  // method found on geomalgorithms.com
  double uu,uv,vv,wu,wv,D;
  uu = U.dot(U);
  uv = U.dot(V);
  vv = V.dot(V);
  wu = w.dot(U);
  wv = w.dot(V);
  D = (uv * uv) - (uu * vv);

  double s = (uv * wv - vv * wu) / D;

  double r = (uv * wu - uu * wv) / D;

  // the point T is within the triangle when both s and t are >= 0
  // and when s+t is <= 1
  // s = 0 or t = 0 or s+t = 1 signify that the point is on the edge
  // I used r instead of t for the variable since I already used t
  if(s >= 0.0 && r >= 0.0 && (s+r) <= 1.0){
    return Hit(t,N);
  }else{
    return Hit::NO_HIT();
  }

}
