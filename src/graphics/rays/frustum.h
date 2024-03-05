#pragma once

struct FrustumPlane {
    float3 normal, point;
    float d;
};

struct Frustum {
    float4 planes[4] = {};

    Frustum() = default;
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
    Frustum(float3 fo, float3 fd_tl, float3 fd_tr, float3 fd_bl, float3 fd_br,
            f32 extend = 1'000.0f) {
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
    }
};
