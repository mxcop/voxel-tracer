#include "precomp.h"
#include "sphere.h"

Sphere::Sphere(const float3 pos, const f32 radius) : pos(pos), radius(radius) {}

AABB Sphere::get_aabb() const { return AABB(pos - float3(radius), pos + float3(radius)); }

float3 Sphere::center() const { return pos; }

HitInfo Sphere::intersect(const Ray& ray) const {
    HitInfo hit;

    const float3 oc = ray.origin - pos;
    const f32 a = dot(ray.dir, ray.dir);
    const f32 b = 2.0 * dot(oc, ray.dir);
    const f32 c = dot(oc, oc) - radius * radius;
    const f32 discriminant = b * b - 4 * a * c;
    if (discriminant < 0) {
        return hit;
    }
    const f32 tmin = (-b - sqrt(discriminant)) / (2.0 * a);

    hit.depth = tmin;
    return hit;
}
