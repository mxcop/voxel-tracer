#pragma once

/**
 * @brief View pyramid. (used for reprojection)
 */
class Pyramid {
    /* Pyramid plane */
    struct Plane {
        float4 normal;
    };

    Plane planes[4] = {};
    float4 forward = {};

   public:
    Pyramid() = default;
    /**
     * @brief Create a pyramid from an origin and 3 corner rays.
     *
     * @param o Frustum origin point.
     * @param f Frustum forward direction.
     * @param tl Frustum direction top left.
     * @param tr Frustum direction top right.
     * @param bl Frustum direction bottom left.
     */
    Pyramid(const float3& o, const float3 f, const float3 tl, const float3 tr, const float3 bl);

    /**
     * @brief Project a world point onto the pyramid view.
     *
     * @param point The point to project.
     * @return UV coordinates of the point. (outside view if not between 0 and 1)
     */
    float2 project(const float3& point) const;

    /**
     * @brief Project a world point onto the pyramid view. (only if its in front of the pyramid)
     *
     * @param point The point to project.
     * @return UV coordinates of the point. (outside view if not between 0 and 1)
     */
    float2 safe_project(const float3& point) const;
};
