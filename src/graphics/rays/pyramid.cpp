#include "pyramid.h"

#include "dev/debug.h"

Pyramid::Pyramid(const float3& o, const float3 f, const float3 tl, const float3 tr,
                 const float3 bl) {
    const float4 origin4 = float4(o);
    const float3 br = tr - (tl - bl);
    origin = o;

    /* Corners */
    far_corners[0] = o + tl * 10.0f, rays[0] = tl;
    far_corners[1] = o + tr * 10.0f, rays[1] = tr;
    far_corners[2] = o + bl * 10.0f, rays[2] = bl;
    far_corners[3] = o + br * 10.0f, rays[3] = br;

    /* Left pyramid plane */
    planes[0].normal = normalize(cross(bl, tl));
    planes[0].normal.w = -dot(planes[0].normal, origin4);

    /* Right pyramid plane */
    planes[1].normal = normalize(cross(tr, br));
    planes[1].normal.w = -dot(planes[1].normal, origin4);
    rays[4] = planes[1].normal;

    /* Top pyramid plane */
    planes[2].normal = normalize(cross(tl, tr));
    planes[2].normal.w = -dot(planes[2].normal, origin4);
    rays[5] = planes[2].normal;

    /* Bottom pyramid plane */
    planes[3].normal = normalize(cross(br, bl));
    planes[3].normal.w = -dot(planes[3].normal, origin4);

    /* Forward plane */
    forward = f;
    forward.w = -dot(forward, origin4);
}

void Pyramid::db_draw() const {
    for (u32 i = 0; i < 4; i++) {
        db::draw_line(origin, far_corners[i], 0xFF0000FF);
    }
    for (u32 i = 0; i < 4; i++) {
        db::draw_normal(origin, planes[i].normal, 0xFFFF0000);
    }
}

float2 Pyramid::project(const float3& point) const {
    const float4 p = float4(point, 1);

    /* Left & Right */
    const f32 d1 = dot(planes[0].normal, p);
    const f32 d2 = dot(planes[1].normal, p);
    const f32 u = d1 / (d1 + d2);

    /* Top & Bottom */
    const f32 d3 = dot(planes[2].normal, p);
    const f32 d4 = dot(planes[3].normal, p);
    const f32 v = d3 / (d3 + d4);

    return float2(u, v);
}

float2 Pyramid::safe_project(const float3& point) const {
    const float4 p = float4(point, 1);

    /* Exit if the point is behind the camera */
    if (dot(forward, p) < 0) return 100'000.0f;

    return project(point);
}

/**
 * @brief Project a 3D point onto a vector. (good enough for comparisons)
 *
 * @param p The point to project.
 * @param n The unit vector to project onto.
 * @return The projection of P onto N.
 */
inline f32 project_onto(const float3& p, const float3& n) { return dot(p, n); }

void Pyramid::projected_minmax(const float3& n, f32& min, f32& max) const {
    min = BIG_F32, max = -BIG_F32;

    for (u32 c = 0; c < 4; ++c) {
        /* Save min and max projection */
        const f32 proj = project_onto(far_corners[c], n);
        min = fminf(min, proj), max = fmaxf(max, proj);
    }

    const f32 proj_o = project_onto(origin, n);
    min = fminf(min, proj_o), max = fmaxf(max, proj_o);
}
