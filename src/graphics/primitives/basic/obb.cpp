#include "precomp.h"
#include "obb.h"

OBB::OBB(const float3 pos, const float3 size, const float3 pivot, const quat rot)
    : pos(pos), pivot(pivot), rot(rot), size(size) {
    set_rotation_pivot(pivot, rot);
}

float3 OBB::center() const { return pos; }

float OBB::area() const {
    const float3 e = size;
    return e.x * e.x + e.y * e.y + e.z * e.z;
}

void OBB::set_position(const float3 position) { 
    pos = position;
    set_rotation_pivot(pivot, rot);
}

void OBB::set_rotation(const quat& rotation) {
    rot = rotation;
    set_rotation_pivot(pivot, rot);
}

void OBB::set_rotation_pivot(const float3 pivot, const quat& rotation) {
    rot = rotation;
    /* Rotate around the center of the box */
    model = mat4::Identity();
    model = model * mat4::Translate(pos);
    model = model * rotation.toMatrix();
    model = model * mat4::Translate(-pivot);

    imodel = model.Inverted();
}

AABB OBB::get_aabb() const {
    const float3 extent = size * 0.5f;

    /* Inspired by : <https://zeux.io/2010/10/17/aabb-from-obb-with-component-wise-abs/> */
    /* Get the transformed center and extend */
    const float3 t_center = TransformPosition(size * 0.5f, model);
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

        f32 t1 = (e + 0) * f_inv;
        f32 t2 = (e + size[d]) * f_inv;

        /* Swap t1 & t2 so t1 is always the smallest */
        if (t1 > t2) {
            const f32 temp = t1;
            t1 = t2, t2 = temp;
        }

        tmin = fmaxf(t1, tmin);
        tmax = fminf(t2, tmax);

        /* Early out check (small value fixes strange black edges when shading) */
        if (tmax - 0.0001f < tmin) return hit;
    }

    hit.depth = tmin;
    hit.albedo = float4(1, 1, 1, 0);
    hit.normal = intersection_normal(ray, tmin);
    return hit;
}

float3 OBB::intersection_normal(const Ray& ray, const f32 tmin) const {
    const float3 p = TransformPosition(ray.origin + ray.dir * tmin, imodel);
    const float3 min = 0, max = size;
    float3 normal = 0;
    if (fabs(p.x - min.x) < 0.001f) {
        normal = float3(-1, 0, 0);
    } else if (fabs(p.x - max.x) < 0.001f) {
        normal = float3(1, 0, 0);
    } else if (fabs(p.y - min.y) < 0.001f) {
        normal = float3(0, -1, 0);
    } else if (fabs(p.y - max.y) < 0.001f) {
        normal = float3(0, 1, 0);
    } else if (fabs(p.z - min.z) < 0.001f) {
        normal = float3(0, 0, -1);
    } else if (fabs(p.z - max.z) < 0.001f) {
        normal = float3(0, 0, 1);
    }
    return TransformVector(normal, model);
}

Ray OBB::world_to_local(const Ray& ray) const {
    return Ray(TransformPosition(ray.origin, imodel), TransformVector(ray.dir, imodel));
}
