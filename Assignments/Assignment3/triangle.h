#ifndef TRIANGLE_H
#define TRIANGLE_H

#include "object.h"

class Triangle : public Object{
public:
  Triangle(Point V1, Point V2, Point V3)
  : V1(V1), V2(V2), V3(V3){}

  virtual Hit intersect(const Ray &ray);

  // vertices
  const Point V1;
  const Point V2;
  const Point V3;
};

#endif
