#pragma once

struct FrustumPlane {
    float3 normal, point;
    float d;
};

class Frustum {
    float4 planes[4] = {};
    float3 corners[8] = {};

   public:
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
    Frustum(const float3 fo, const float3 fd_tl, const float3 fd_tr, const float3 fd_bl,
            const float3 fd_br, const f32 extend = 1'000.0f);

    bool intersect_unitcube() const;

    float2 project(const float3& point) const;
};
