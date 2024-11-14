#include "swTriangle.h"

namespace sw {

bool Triangle::intersect(const Ray &ray, Intersection &isect) const {
    Vec3 edge1 = vertices[1] - vertices[0];
    Vec3 edge2 = vertices[2] - vertices[0];
    Vec3 v0 = vertices[0];
    Vec3 n = edge1 % edge2;
    float m = -n * v0;
    float t = (n * ray.orig + m) / (-n * ray.dir);
    if (t < ray.minT || t > ray.maxT) return false;
    Vec3 Q = ray.orig + t * ray.dir;
    Vec3 r = Q - v0;
    float v = (edge1 % r).norm() / n.norm();
    float w = (r % edge2).norm() / n.norm();

    if (v + w < 1.0f && (edge1 % r) * n >= 0.0f && (r % edge2) * n >= 0.0f) {
        isect.hitT = t;
        isect.position = Q;
        isect.normal = n.normalize();
        isect.material = material;
        isect.frontFacing = (-ray.dir * isect.normal) > 0.0f;
        if (!isect.frontFacing) isect.normal = -isect.normal;
        isect.ray = ray;
        return true;
    }

    return false;
}

} 
