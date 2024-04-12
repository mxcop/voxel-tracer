#include "capsule.h"

Capsule::Capsule(const float3 a, const float3 b, const f32 radius) : a(a), b(b), radius(radius) {}

AABB Capsule::get_aabb() const {
    const float3 min = fminf(a, b) - radius;
    const float3 max = fmaxf(a, b) + radius;
    return AABB(min, max);
}

float3 Capsule::center() const { return (a + b) * 0.5f; }

static inline f32 cap_intersect(const float3& ro, const float3& rd, const float3& pa,
                               const float3& pb, f32 r) {
    const float3 ba = pb - pa;
    const float3 oa = ro - pa;

    const f32 baba = dot(ba, ba);
    const f32 bard = dot(ba, rd);
    const f32 baoa = dot(ba, oa);
    const f32 rdoa = dot(rd, oa);
    const f32 oaoa = dot(oa, oa);

    f32 a = baba - bard * bard;
    f32 b = baba * rdoa - baoa * bard;
    f32 c = baba * oaoa - baoa * baoa - r * r * baba;
    f32 h = b * b - a * c;
    if (h >= 0.0f && b < 1.0f) {
        const f32 t = (-b - sqrt(h)) / a;
        const f32 y = baoa + t * bard;
        // body
        if (y > 0.0f && y < baba) return t;
        // caps
        const float3 oc = (y <= 0.0f) ? oa : ro - pb;
        const f32 d = dot(rd, rd);
        b = dot(rd, oc);
        c = dot(oc, oc) - r * r;
        h = b * b - d * c;
        if (h < 0) return BIG_F32;

        const f32 tmin = (- b - sqrt(h));
        if (tmin < 0) return BIG_F32;

        return tmin;
    }
    return BIG_F32;
}

static inline float3 cap_normal(const float3& pos, const float3& a, const float3& b, f32 r) {
    const float3 ba = b - a;
    const float3 pa = pos - a;
    const f32 h = clamp(dot(pa, ba) / dot(ba, ba), 0.0f, 1.0f);
    return (pa - h * ba) / r;
}

HitInfo Capsule::intersect(const Ray& ray) const {
    HitInfo hit;

    const f32 t = cap_intersect(ray.origin, ray.dir, a, b, radius);
    if (t < BIG_F32) {
        hit.depth = t - 0.01f;

        const float3 point = ray.origin + ray.dir * (t - 0.01f);
        hit.normal = cap_normal(point, a, b, radius);
        hit.albedo = {50, 0, 0};
        hit.material = 0xFF;
    }

    return hit;
}
