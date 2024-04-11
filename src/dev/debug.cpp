#include "debug.h"

#ifdef DEV

/* Internal development variables */
namespace inter {}  // namespace inter

namespace db {

/**
 * @brief Project a world point onto the main camera.
 */
inline static int2 world_to_screen(const float3& p) {
    const Pyramid& view_pyramid = dev::main_camera->pyramid;
    return floori(view_pyramid.safe_project(p) * WIN_SIZE);
}

/**
 * @brief Draw a line from pixel A to B.
 */
inline static void draw_line_screen(const int2& a, const int2& b, const u32 c) {
    dev::db_screen->Line((f32)a.x, (f32)a.y, (f32)b.x, (f32)b.y, c);
}

/**
 * @brief Draw a line from world point A to B.
 */
void draw_line(const float3& a, const float3& b, const u32 c) {
    const int2 a_pix = world_to_screen(a);
    if (a_pix.x == WIN_WIDTH * 100'000) return;
    const int2 b_pix = world_to_screen(b);
    if (b_pix.x == WIN_WIDTH * 100'000) return;

    draw_line_screen(a_pix, b_pix, c);
}

/**
 * @brief Draw a normal at world point.
 */
void draw_normal(const float3& point, const float3& n, const u32 c) { 
    const float3 end_point = point + n;
    draw_line(point, end_point, c);

    const float3 side = normalize(cross(n, UP)) * 0.2f;
    draw_line(end_point, point + n * 0.8f + side, c);
    draw_line(end_point, point + n * 0.8f - side, c);
}

/**
 * @brief Draw an AABB from world points min and max.
 */
void draw_aabb(const float3& min, const float3& max, const u32 c) {
    /* Bottom plane */
    const int2 p000 = world_to_screen({min.x, min.y, min.z});
    const int2 p100 = world_to_screen({max.x, min.y, min.z});
    const int2 p101 = world_to_screen({max.x, min.y, max.z});
    const int2 p001 = world_to_screen({min.x, min.y, max.z});
    draw_line_screen(p000, p100, c);
    draw_line_screen(p000, p001, c);
    draw_line_screen(p101, p001, c);
    draw_line_screen(p101, p100, c);

    /* Top plane */
    const int2 p111 = world_to_screen({max.x, max.y, max.z});
    const int2 p110 = world_to_screen({max.x, max.y, min.z});
    const int2 p010 = world_to_screen({min.x, max.y, min.z});
    const int2 p011 = world_to_screen({min.x, max.y, max.z});
    draw_line_screen(p010, p110, c);
    draw_line_screen(p010, p011, c);
    draw_line_screen(p111, p011, c);
    draw_line_screen(p111, p110, c);

    /* Side edges */
    draw_line_screen(p000, p010, c);
    draw_line_screen(p100, p110, c);
    draw_line_screen(p001, p011, c);
    draw_line_screen(p101, p111, c);
}

/**
 * @brief Draw an OBB from a world point, size, and a rotation quaternion.
 */
void draw_obb(const float3& point, const float3& size, const quat& rot, const u32 c) {
    const float3 max = size * 0.5f;
    const float3 min = -max;

    /* Bottom plane */
    const int2 p000 = world_to_screen(point + rot.rotateVector({min.x, min.y, min.z}));
    const int2 p100 = world_to_screen(point + rot.rotateVector({max.x, min.y, min.z}));
    const int2 p101 = world_to_screen(point + rot.rotateVector({max.x, min.y, max.z}));
    const int2 p001 = world_to_screen(point + rot.rotateVector({min.x, min.y, max.z}));
    draw_line_screen(p000, p100, c);
    draw_line_screen(p000, p001, c);
    draw_line_screen(p101, p001, c);
    draw_line_screen(p101, p100, c);

    /* Top plane */
    const int2 p111 = world_to_screen(point + rot.rotateVector({max.x, max.y, max.z}));
    const int2 p110 = world_to_screen(point + rot.rotateVector({max.x, max.y, min.z}));
    const int2 p010 = world_to_screen(point + rot.rotateVector({min.x, max.y, min.z}));
    const int2 p011 = world_to_screen(point + rot.rotateVector({min.x, max.y, max.z}));
    draw_line_screen(p010, p110, c);
    draw_line_screen(p010, p011, c);
    draw_line_screen(p111, p011, c);
    draw_line_screen(p111, p110, c);

    /* Side edges */
    draw_line_screen(p000, p010, c);
    draw_line_screen(p100, p110, c);
    draw_line_screen(p001, p011, c);
    draw_line_screen(p101, p111, c);
}

}  // namespace db

#endif
