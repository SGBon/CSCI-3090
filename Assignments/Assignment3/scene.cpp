//
//  Framework for a raytracer
//  File: scene.cpp
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

#include "scene.h"
#include "material.h"

Color Scene::trace(const Ray &ray, int steps,double eta1)
{
    // Find hit object and distance
    Hit min_hit(std::numeric_limits<double>::infinity(),Vector());
    Object *obj = NULL;
    for (unsigned int i = 0; i < objects.size(); ++i) {
        Hit hit(objects[i]->intersect(ray));
        if (hit.t<min_hit.t) {
            min_hit = hit;
            obj = objects[i];
        }
    }

    // No hit? Return background color.
    if (!obj) return Color(0.0, 0.0, 0.0);

    Material *material = obj->material;            //the hit objects material
    Point hit = ray.at(min_hit.t);                 //the hit point
    Vector N = min_hit.N;                          //the normal at hit point
    Vector V = -ray.D;                             //the view vector

    /****************************************************
    * This is where you should insert the color
    * calculation (Phong model).
    *
    * Given: material, hit, N, V, lights[]
    * Sought: color
    *
    * Hints: (see triple.h)
    *        Triple.dot(Vector) dot product
    *        Vector+Vector      vector sum
    *        Vector-Vector      vector difference
    *        Point-Point        yields vector
    *        Vector.normalize() normalizes vector, returns length
    *        double*Color        scales each color component (r,g,b)
    *        Color*Color        dito
    *        pow(a,b)           a to the power of b
    ****************************************************/

    // extract and/or declare variables
    // for lighting calculations
    Color ambient = Color(1.0,1.0,1.0); // ambient colour
    Color base = material->color; // material colour
    double ka = material->ka; // ambient intensity;
    double kd = material->kd; // diffuse intensity
    double ks = material->ks; // specular intensity
    double e = material->n; // exponent of specular highlight size
    double reflect = material->reflect; // reflect coefficient
    double refract = material->refract; // refraction coefficient

    double eta2 = material->eta; // refraction index
    if(eta1 == eta2){
      eta2 = 1.0;
    }

    // get reflected ray
    Vector vec_ref = ray.D - (2.0 *(ray.D.dot(N))*N); // reflect ray direction
    Ray ray_ref(hit,vec_ref); //reflect ray
    // jiggle the ray
    jiggle(ray_ref);

    // hack
    Vector frac_n;
    if(ray.D.dot(N) < 0.0){
      frac_n = N;
    }else{
      frac_n = -N;
    }

    // get refracted ray
    bool frac_flag;
    Vector frac_dir = fractf(eta1,eta2,ray.D,frac_n); // direction of refraction

    Ray ray_frac(hit,frac_dir); // ray going out of the material
    if(frac_dir.length_2() > 0.0 && refract > 0.0){
      frac_flag = true;
    }else{
      frac_flag = false;
    }

    // jiggle the ray
    jiggle(ray_frac);

    Color c_ref; // colour of reflected ray
    Color c_frac; // colour of refracted ray
    // recursively trace reflected/refracted rays up to steps times
    if(steps > 0){
      if(reflect > 0.0) c_ref = trace(ray_ref,steps-1,eta1);
      if(frac_flag) c_frac = trace(ray_frac, steps-1,eta2);
    }

    Color color = ka * base * ambient; // set ambient colour
    for(unsigned int i = 0;i<lights.size();i++){
      bool shadows = false; // flag if the current light cast a shadow
      Vector L = hit - lights[i]->position; // vector of light direction
      Vector SL = lights[i]->position - hit; // vector of shadow feeler
      L.normalize();
      SL.normalize();

      // get shadow feelers
      Ray feeler = Ray(hit,SL);
      //jiggle Ray
      jiggle(feeler);
      // test to see if object is in shadow
      for(unsigned int i = 0;i<objects.size();i++){
        Hit test(objects[i]->intersect(feeler));
        if(test.t >= 0.0)  {
          shadows = true;
          break;
        }
      }

      if(!shadows){
        Color lc = lights[i]->color; // colour of light

        double lnDot = L.dot(N); // dot product of light

        Vector R = L - (2.0*N*lnDot); // reflection vector
        R.normalize();
        double rvDot = R.dot(V) ;

        color += (kd*base*lc*max(0.0,lnDot)) + (ks*lc*pow(max(0.0,rvDot),e));
      }
    }

    color += reflect*c_ref + refract*c_frac;

    return color;
}

void Scene::render(Image &img)
{
    int w = img.width();
    int h = img.height();
    for (int y = 0; y < h; y++) {
        for (int x = 0; x < w; x++) {
            Point pixel(x, h-1-y, 0);
            Ray ray(eye, (pixel-eye).normalized());
            Color col = trace(ray,10,1.0);
            col.clamp();
            img(x,y) = col;
        }
    }
}

void Scene::addObject(Object *o)
{
    objects.push_back(o);
}

void Scene::addLight(Light *l)
{
    lights.push_back(l);
}

void Scene::setEye(Triple e)
{
    eye = e;
}

void Scene::jiggle(Ray& ray){
  Point pjig = ray.at(pow(2,-32));
  ray = Ray(pjig,ray.D);
}

Vector Scene::fractf(double eta1, double eta2, Vector dir, Vector norm){
  Vector t;

  // check if root is negative
  double root = 1.0 - ((pow(eta1,2.0)*(1.0 - (pow(dir.dot(norm),2.0))))/pow(eta2,2.0));

  if(root < 0.0){
    return t;
  }

  // first half of equation
  t = ((eta1*(dir - norm*(dir.dot(norm))))/eta2) - (norm*sqrt(root));

  return t;
}
