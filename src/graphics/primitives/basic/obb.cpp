#include "precomp.h"
#include "obb.h"

OBB::OBB(const float3 pos, const float3 size, const float3 axis, const f32 angle) : pos(pos), size(size) {
    /* Rotate around the center of the box */
    model = mat4::Translate(center()) * mat4::Rotate(axis, angle) * mat4::Translate(-center());
}

float3 OBB::center() const { return pos + size * 0.5f; }

float OBB::area() const {
    const float3 e = size;
    return e.x * e.x + e.y * e.y + e.z * e.z;
}

AABB OBB::get_aabb() const {
    const float3 extent = size * 0.5f;

    /* Inspired by : <https://zeux.io/2010/10/17/aabb-from-obb-with-component-wise-abs/> */
    /* Get the transformed center and extend */
    const float3 t_center = TransformPosition(center(), model);
    const float3 t_extent = TransformVector(extent, fabs(model));

    return AABB(t_center - t_extent, t_center + t_extent);
}

HitInfo OBB::intersect(const Ray& ray) const {
    HitInfo hit;
    f32 tmin = 0.0f, tmax = BIG_F32;

    /* "model[3]" holds the world position of the box */
    const float3 delta = model.GetTranslation() - ray.origin;

    /* Loop to be unrolled by the compiler */
    for (u32 d = 0; d < 3; ++d) {
        const float3 axis = float3(model.GetColumn(d));
        const f32 e = dot(axis, delta), f_inv = 1.0f / dot(ray.dir, axis);

        f32 t1 = (e + pos[d]) * f_inv;
        f32 t2 = (e + pos[d] + size[d]) * f_inv;

        /* Swap t1 & t2 so t1 is always the smallest */
        if (t1 > t2) {
            const f32 temp = t1;
            t1 = t2, t2 = temp;
        }

        tmin = fmaxf(t1, tmin);
        tmax = fminf(t2, tmax);

        /* Early out check */
        if (tmax < tmin) return hit;
    }

    hit.depth = tmin;
    return hit;
}
