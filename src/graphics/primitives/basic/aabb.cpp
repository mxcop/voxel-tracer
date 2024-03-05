#include "aabb.h"

AABB::AABB(const f32 min, const f32 max) : min(min), max(max) {}

AABB::AABB(const float3 min, const float3 max) : min(min), max(max) {}

void AABB::grow(const float3 p) {
    min = fminf(min, p);
    max = fmaxf(max, p);
}

void AABB::grow(const AABB& aabb) {
    if (aabb.min.x != 1e30f) {
        this->grow(aabb.min);
        this->grow(aabb.max);
    }
}

float3 AABB::center() const { return min + (max - min) * 0.5f; }

float AABB::area() const {
    const float3 e = max - min;
    return e.x * e.x + e.y * e.y + e.z * e.z;
}

AABB AABB::get_aabb() const { return *this; }

HitInfo AABB::intersect(const Ray& ray) const {
    HitInfo hit;
    /* Idea to use fmsub to save 1 instruction came from
     * <http://www.joshbarczak.com/blog/?p=787> */
    const f128 ord = _mm_mul_ps(ray.o, ray.rd);
    const f128 t1 = _mm_fmsub_ps(bmin, ray.rd, ord);
    const f128 t2 = _mm_fmsub_ps(bmax, ray.rd, ord);

    /* Find the near and far intersection point */
    const f128 vmax4 = _mm_max_ps(t1, t2), vmin4 = _mm_min_ps(t1, t2);

    /* Set the 4th element to 0 */
    const f128 tmin4 = (f128&)_mm_slli_si128((i128&)vmin4, 4);

    /* Get the horizontal minimum and maximum "t" */
    const f32 tmax = _mm_hmin3_ps(vmax4), tmin = _mm_hmax_ps(tmin4);

    if (tmax >= tmin) hit.depth = tmin;
    return hit;
}
