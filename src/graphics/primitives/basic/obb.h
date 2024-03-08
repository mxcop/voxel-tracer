#pragma once

#include "graphics/primitives/traceable.h"

/**
 * @brief Oriented bounding box. (Box with rotation)
 */
struct OBB : public Traceable {
    /* Model matrix & inverted model matrix */
    mat4 model, imodel;
    float3 pos = 0, size = 0;

    OBB() = default;
    OBB(const float3 pos, const float3 size, const float3 axis = 0, const f32 angle = 0);

    /* Trace-able functions */
    AABB get_aabb() const override;
    float3 center() const override;
    HitInfo intersect(const Ray& ray) const override;

    f32 area() const;
    void set_rotation(const float3& axis, const f32 angle);
    void set_rotation_pivot(const float3& pivot, const float3& axis, const f32 angle);

    /* Transform a ray from world space to the OBB local space */
    Ray world_to_local(const Ray& ray) const;

   private:
    float3 intersection_normal(const Ray& ray, const f32 tmin) const;
};
