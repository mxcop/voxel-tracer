#include "precomp.h"
#include "pyramid.h"

Pyramid::Pyramid(const float3& o, float3 tl, float3 tr, float3 bl) {
    origin = o;

    /* Left frustum plane */
    planes[0].normal = normalize(cross(bl, tl));
    planes[0].d = dot(planes[0].normal, origin);

    /* Right frustum plane */
    const float3 br = normalize(tr - (tl - bl));
    planes[1].normal = normalize(cross(tr, br));
    planes[1].d = dot(planes[1].normal, origin);

    /* Top frustum plane */
    planes[2].normal = normalize(cross(tl, tr));
    planes[2].d = dot(planes[2].normal, origin);

    /* Bottom frustum plane */
    planes[3].normal = normalize(cross(br, bl));
    planes[3].d = dot(planes[3].normal, origin);
}

float2 Pyramid::project(const float3& point) const {
    /* Left & Right */
    const f32 d1 = dot(planes[0].normal, point) - planes[0].d;
    const f32 d2 = dot(planes[1].normal, point) - planes[1].d;

    /* Top & Bottom */
    const f32 d3 = dot(planes[2].normal, point) - planes[2].d;
    const f32 d4 = dot(planes[3].normal, point) - planes[3].d;

    const f32 u = d1 / (d1 + d2);
    const f32 v = d3 / (d3 + d4);

    return float2(u, v);
}
