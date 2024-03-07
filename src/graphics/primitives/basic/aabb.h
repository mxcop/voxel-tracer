#pragma once

#include "graphics/primitives/traceable.h"

/**
 * @brief Axis aligned bounding box. (Box without rotation)
 */
struct AABB : public Traceable {
    union {
        float4 corners[2] = {float4(BIG_F32), float4(-BIG_F32)};
        struct {
            float3 min;
            f32 d0;
            float3 max;
            f32 d1;
        };
        struct {
            f128 bmin, bmax;
        };
    };
    float3 color = 0;

    AABB() = default;
    AABB(const f32 min, const f32 max);
    AABB(const float3 min, const float3 max, const float3 color = 0);

    /* Trace-able functions */
    AABB get_aabb() const override;
    float3 center() const override;
    HitInfo intersect(const Ray& ray) const override;

    /**
     * @brief Grow AABB to include given point.
     */
    void grow(const float3 p);
    void grow(const AABB& aabb);
    f32 area() const;

   private:
    float3 intersection_normal(const Ray& ray, const f32 tmin) const;
};
