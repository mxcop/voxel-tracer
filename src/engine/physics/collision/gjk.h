#pragma once

/*
 * Huge props to this article for explaining GJK with code snippets!
 * Source : <https://winter.dev/articles/gjk-algorithm>
 */

/* Get point on the Minkowski difference. */
inline float3 gjk_support(const Collider& a, const Transform& ta, const Collider& b,
                          const Transform& tb, const float3& dir) {
    return a.furthest_point(ta, dir) - b.furthest_point(tb, -dir);
}

/**
 * @brief Check if 2 colliders are intersecting.
 */
bool gjk(const Collider& a, const Transform& ta, const Collider& b, const Transform& tb);
