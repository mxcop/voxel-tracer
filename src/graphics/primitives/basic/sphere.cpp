#include "precomp.h"
#include "sphere.h"

Sphere::Sphere(const float3 pos, const f32 radius) : pos(pos), radius(radius) {}

AABB Sphere::get_aabb() const { return AABB(pos - float3(radius), pos + float3(radius)); }

float3 Sphere::center() const { return pos; }

HitInfo Sphere::intersect(const Ray& ray) const {
    HitInfo hit;

    const float3 oc = ray.origin - pos;
    const f32 a = sqrLength(ray.dir);
    const f32 half_b = dot(oc, ray.dir);
    const f32 c = sqrLength(oc) - radius * radius;
    const f32 discriminant = half_b * half_b - a * c;
    if (discriminant < 0) {
        return hit;
    }
    const f32 tmin = (-half_b - sqrt(discriminant)) / a;
    if (tmin < 0) {
        return hit;
    }

    hit.depth = tmin;
    hit.albedo = float4(1, 1, 1, 0.01f);
    const float3 p = (ray.origin + ray.dir * tmin);
    hit.normal = normalize(p - pos);
    // TEMP: for fun.
    // hit.albedo = hit.normal + 1.0f * 0.5f; 
    return hit;
}
