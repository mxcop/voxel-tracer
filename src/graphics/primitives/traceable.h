#pragma once

#include "graphics/rays/hit.h"
#include "graphics/rays/ray.h"

struct AABB;

/**
 * @brief A ray trace-able object interface.
 */
class Traceable {
   public:
    /**
     * @brief Get the axis aligned bounding box of this trace-able.
     */
    virtual AABB get_aabb() const = 0;

    /**
     * @brief Get the center point of the trace-able.
     */
    virtual float3 center() const = 0;

    /**
     * @brief Perform ray intersection test with this trace-able.
     * @return Information about the intersection. (as HitInfo)
     */
    virtual HitInfo intersect(const Ray& ray) const = 0;

    virtual PacketHit8x8 intersect(const RayPacket8x8& packet, const bool debug = false) const {
        return PacketHit8x8();
    };
};
