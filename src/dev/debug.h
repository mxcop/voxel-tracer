#pragma once

namespace db {

constexpr u32 DEFAULT_COLOR = 0xFF00FF00; /* Green */

#ifdef DEV

/**
 * @brief Draw a line from world point A to B.
 */
extern void draw_line(const float3& a, const float3& b, const u32 c = DEFAULT_COLOR);

/**
 * @brief Draw a normal at world point.
 */
extern void draw_normal(const float3& point, const float3& n, const u32 c = DEFAULT_COLOR);

/**
 * @brief Draw an AABB from world points min and max.
 */
extern void draw_aabb(const float3& min, const float3& max, const u32 c = DEFAULT_COLOR);

/**
 * @brief Draw an OBB from a world point, size, and a rotation quaternion.
 */
extern void draw_obb(const float3& point, const float3& size, const quat& rot,
                     const u32 c = DEFAULT_COLOR);

#else

/**
 * @brief DOES NOTHING IF "DEV" IS UNDEFINED.
 */
inline void draw_line(const float3& a, const float3& b, const u32 c = DEFAULT_COLOR){};

/**
 * @brief DOES NOTHING IF "DEV" IS UNDEFINED.
 */
inline void draw_normal(const float3& point, const float3& n, const u32 c = DEFAULT_COLOR){};

/**
 * @brief DOES NOTHING IF "DEV" IS UNDEFINED.
 */
inline void draw_aabb(const float3& min, const float3& max, const u32 c = DEFAULT_COLOR){};

/**
 * @brief DOES NOTHING IF "DEV" IS UNDEFINED.
 */
inline void draw_obb(const float3& point, const float3& size, const quat& rot,
                     const u32 c = DEFAULT_COLOR){};

#endif

}  // namespace db
