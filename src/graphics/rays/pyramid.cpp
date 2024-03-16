#include "pyramid.h"

Pyramid::Pyramid(const float3& o, const float3 f, const float3 tl, const float3 tr,
                 const float3 bl) {
    const float4 origin = float4(o);

    /* Left frustum plane */
    planes[0].normal = normalize(cross(bl, tl));
    planes[0].normal.w = -dot(planes[0].normal, origin);

    /* Right frustum plane */
    const float3 br = tr - (tl - bl);
    planes[1].normal = normalize(cross(tr, br));
    planes[1].normal.w = -dot(planes[1].normal, origin);

    /* Top frustum plane */
    planes[2].normal = normalize(cross(tl, tr));
    planes[2].normal.w = -dot(planes[2].normal, origin);

    /* Bottom frustum plane */
    planes[3].normal = normalize(cross(br, bl));
    planes[3].normal.w = -dot(planes[3].normal, origin);

    forward = f;
    forward.w = -dot(forward, origin);
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
