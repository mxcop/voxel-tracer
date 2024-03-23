#include "precomp.h"

#include "dev/debug.h"

bool CoherentPacked8x8::setup_slice(const float3& min, const float3& max, const f32 vpu) {
    // TODO: setup slice for traversal.

    /* Grab the 2 corner rays */
    const float3 ray_tl = rays[0], ray_br = rays[15];

    /* Packet major axis */
    k = major_axis(ray_tl);
    /* Remaining axis */
    u = (k + 1) < 3 ? (k + 1) : 0;
    v = ~(k | u) & 0b11;

    /* Find intersection time of U and V rays */
    const f32 entry_u = entry(origin[k], ray_br[k], min[k], max[k]);
    const f32 entry_v = entry(origin[k], ray_tl[k], min[k], max[k]);
    k_t = origin[k] + ray_br[k] * entry_u;

    const f32 upv = 1.0f / vpu;
    const f32 next_u = entry(origin[k], ray_br[k], min[k] + upv, max[k] - upv);
    const f32 next_v = entry(origin[k], ray_tl[k], min[k] + upv, max[k] - upv);

    /* Slice U and V minimum and maximum */
    const f32 u_min = origin[v] + ray_br[v] * entry_u;
    const f32 u_max = origin[u] + ray_br[u] * entry_u;
    const f32 v_min = origin[u] + ray_tl[u] * entry_v;
    const f32 v_max = origin[v] + ray_tl[v] * entry_v;

    /* Slice U and V delta */
    const f32 du_min = (origin[v] + ray_br[v] * next_u) - u_min;
    const f32 du_max = (origin[u] + ray_br[u] * next_u) - u_max;
    const f32 dv_min = (origin[u] + ray_tl[u] * next_v) - v_min;
    const f32 dv_max = (origin[v] + ray_tl[v] * next_v) - v_max;

    /* Slice U and V point */
    float3 u_pos;
    u_pos[k] = k_t;
    u_pos[u] = u_max, u_pos[v] = u_min;
    float3 v_pos;
    v_pos[k] = k_t;
    v_pos[u] = v_min, v_pos[v] = v_max;

    slice.xmm = _mm_setr_ps(u_min, u_max, v_min, v_max);
    delta_slice = _mm_setr_ps(du_min, du_max, dv_min, dv_max);

    // draw_slice();

    return true;
}

void CoherentPacked8x8::traverse(const float3& min, const float3& max, const f32 vpu) {
    const f32 sign = -getsign(rays[0][k]);
    const f32 upv = (1.0f / vpu) * sign;
    const f32 min_t = sign ? min[k] : max[k];
    const f32 max_t = sign ? max[k] : min[k];
    k_t += upv * 0.01f;
    for (k_t; k_t < max_t && k_t > min_t; k_t += upv) {
        slice.xmm = _mm_add_ps(slice.xmm, delta_slice);
        draw_slice(vpu);
    }
}

void CoherentPacked8x8::draw_slice(const f32 vpu) const {
    float3 a_pos;
    a_pos[k] = k_t, a_pos[u] = slice.u_max, a_pos[v] = slice.u_min;
    float3 b_pos;
    b_pos[k] = k_t, b_pos[u] = slice.v_min, b_pos[v] = slice.v_max;

    const f32 sign = -getsign(rays[0][k]);
    const f32 upv = 1.0f / vpu;
    b_pos[k] += upv * sign + upv * 0.01f;
    db::draw_aabb(a_pos, b_pos);

    db::draw_aabb(floorf(a_pos * vpu) * upv, floorf(b_pos * vpu) * upv + upv, 0xFFFF0000);
}

f32 CoherentPacked8x8::entry(const f32 ro, const f32 rd, const f32 min, const f32 max) const {
    bool sign = rd < 0.0f;
    const f32 bmin = sign ? max : min;
    // TODO: maybe try getting rid of division here?
    const f32 dmin = (bmin - ro) / rd;
    return fmaxf(dmin, 0.0f);
}

// ext CoherentPacked8x8::intersection(const float3& rd, const float3& min, const float3& max) const
// {
//     ext t = {0.0f, BIG_F32};
//     const float3 corners[2] = {min, max};
//
//     for (u32 d = 0; d < 3; ++d) {
//         bool sign = signs[d];
//         const f32 bmin = corners[sign][d];
//         const f32 bmax = corners[!sign][d];
//
//         // TODO: maybe try getting rid of division here?
//         const f32 dmin = (bmin - origin[d]) / rd[d];
//         const f32 dmax = (bmax - origin[d]) / rd[d];
//
//         t.min = fmaxf(dmin, t.min);
//         t.max = fminf(dmax, t.max);
//     }
//     return t;
// }
