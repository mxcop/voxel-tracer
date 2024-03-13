#include "frustum.h"

/**
 * @brief Create a frustum from an origin and 4 corner rays.
 *
 * @param fo Frustum origin point.
 * @param fd_tl Frustum direction top left. (normalized)
 * @param fd_tr Frustum direction top right. (normalized)
 * @param fd_bl Frustum direction bottom left. (normalized)
 * @param fd_br Frustum direction bottom right. (normalized)
 * @param extend How far away the far plane should be from the origin.
 */
Frustum::Frustum(const float3 fo, const float3 fd_tl, const float3 fd_tr, const float3 fd_bl,
                 const float3 fd_br, const f32 extend) {
    /* Far points */
    float3 ftl = fo + fd_tl * extend;
    float3 ftr = fo + fd_tr * extend;
    float3 fbl = fo + fd_bl * extend;
    float3 fbr = fo + fd_br * extend;
    /* Near points */
    float3 ntl = fo + fd_tl;
    float3 ntr = fo + fd_tr;
    float3 nbl = fo + fd_bl;
    float3 nbr = fo + fd_br;

    float3 p0, p1, p2;

    /* Left frustum plane */
    p0 = nbl, p1 = fbl, p2 = ftl;
    planes[0] = normalize(cross(p1 - p0, p2 - p1));
    planes[0].w = dot(float3(planes[0]), p0);

    /* Top frustum plane */
    p0 = ntl, p1 = ftl, p2 = ftr;
    planes[1] = normalize(cross(p1 - p0, p2 - p1));
    planes[1].w = dot(float3(planes[1]), p0);

    /* Right frustum plane */
    p0 = ntr, p1 = ftr, p2 = fbr;
    planes[2] = normalize(cross(p1 - p0, p2 - p1));
    planes[2].w = dot(float3(planes[2]), p0);

    /* Bottom frustum plane */
    p0 = nbr, p1 = fbr, p2 = fbl;
    planes[3] = normalize(cross(p1 - p0, p2 - p1));
    planes[3].w = dot(float3(planes[3]), p0);

    /* Far corners */
    corners[0] = ftl;
    corners[1] = ftr;
    corners[2] = fbl;
    corners[3] = fbr;

    /* Near corners */
    corners[4] = ntl;
    corners[5] = ntr;
    corners[6] = nbl;
    corners[7] = nbr;
}

bool Frustum::intersect_unitcube() const {
    const float3 cube_min = float3(0);
    const float3 cube_max = float3(8);

    /* Cube vertices */
    float3 cube_verts[8] = {};
    for (u32 z = 0; z < 2; z++) {
        for (u32 y = 0; y < 2; y++) {
            for (u32 x = 0; x < 2; x++) {
                cube_verts[(z * 2 * 2) + (y * 2) + x] = -float3(x, y, z) * cube_max;
            }
        }
    }
    
    /* Cube planes */
    float4 cube_planes[6] = {};
    cube_planes[0] = float4 (-1.0f, 0.0f, 0.0f,cube_max.x);
    cube_planes[1] = float4(0.0f, -1.0f, 0.0f, cube_max.y);
    cube_planes[2] = float4(0.0f, 0.0f, -1.0f, cube_max.z);
    cube_planes[3] = float4(1.0f, 0.0f, 0.0f, cube_min.x);
    cube_planes[4] = float4(0.0f, 1.0f, 0.0f, cube_min.y);
    cube_planes[5] = float4(0.0f, 0.0f, 1.0f, cube_min.z);

    bool intersects = true;

    /* Test cube vertices vs frustum planes */
    for (int i = 0; i < 4; ++i) {
        bool isAnyVertexInPositiveSide = false;
        for (int j = 0; j < 8; ++j) {
            float dotResult = dot(float3(planes[i]), cube_verts[j]) + planes[i].w;
            isAnyVertexInPositiveSide |= dotResult > 0;
        }

        intersects &= isAnyVertexInPositiveSide;
    }

    /* Test frustum vertices vs cube planes */
    /* NOTE: this isn't very effective, only cuts out corners IF angles are grazing. */
    /* TODO: improve accuracy of this check... */
    for (int i = 0; i < 6; ++i) {
        bool isAnyVertexInPositiveSide = false;
        for (int j = 0; j < 8; ++j) {
            float dotResult = dot(float3(cube_planes[i]), corners[j]) + cube_planes[i].w;
            isAnyVertexInPositiveSide |= dotResult > 0;
        }

        intersects &= isAnyVertexInPositiveSide;
    }
    return intersects;
}

float2 Frustum::project(const float3& point) const {
    /* Left & Right */
    const f32 d1 = dot(float3(planes[0]), point - corners[4]);
    const f32 d2 = dot(float3(planes[2]), point - corners[4]);

    /* Top & Bottom */
    const f32 d3 = dot(float3(planes[1]), point - corners[4]);
    const f32 d4 = dot(float3(planes[3]), point - corners[4]);

    const f32 u = d1 / (d1 + d2);
    const f32 v = d3 / (d3 + d4);
    
    return make_float2(u, v);
}
