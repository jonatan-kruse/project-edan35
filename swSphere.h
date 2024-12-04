#pragma once

#include "swIntersection.h"
#include "swPrimitive.h"
#include "swRay.h"

namespace sw {

class Sphere : public Primitive {
  public:
    Sphere() = default;
    Sphere(const Vec3 &c, const float &r, const Material &m) : center(c), radius(r), material(m) {velocity = Vec3(0.f, 0.f, 0.f);}
    Sphere(const Sphere &s) = default;
    Sphere(Sphere &&) = default;
    Sphere &operator=(Sphere &&) = default;

    bool intersect(const Ray &r, Intersection &isect) const;
    void tick(const float dt) {center = center + velocity * dt;}
    void set_velocity(const Vec3 &v) {velocity = v;}

  public:
    Vec3 center;
    float radius{0.0f};
    Material material;
    Vec3 velocity;
};

} // namespace sw
