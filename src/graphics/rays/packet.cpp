#include "precomp.h"

#include "dev/debug.h"

bool CoherentPacked8x8::setup_slice(const float3& min, const float3& max, const f32 vpu) {
    // TODO: setup slice for traversal.

    /* Grab the 4 corner rays */
    const float3 ray_tl = rays[0], ray_tr = rays[3], ray_bl = rays[12], ray_br = rays[15];

    /* Packet major axis */
    const u32 k = major_axis(ray_tl);

    /* Remaining axis */
    const u32 u = (k + 1) < 3 ? (k + 1) : 0;
    const u32 v = (~k & ~u) & 0b11;

    /* Find intersection extend of each corner ray */
    const f32 entry_tl = entry(origin[k], ray_tl[k], min[k], max[k]);
    const f32 entry_br = entry(origin[k], ray_br[k], min[k], max[k]);

    /* Slice U and V */
    const f32 u_min = origin[v] + ray_br[v] * entry_br;
    const f32 u_max = origin[u] + ray_br[u] * entry_br;
    const f32 v_min = origin[u] + ray_tl[u] * entry_tl;
    const f32 v_max = origin[v] + ray_tl[v] * entry_tl;

    float3 u_pos = 0;
    u_pos[k] = origin[k] + ray_tl[k] * entry_tl;
    u_pos[u] = v_min;
    u_pos[v] = v_max;
    float3 v_pos = 0;
    v_pos[k] = origin[k] + ray_br[k] * entry_br;
    v_pos[u] = u_max;
    v_pos[v] = u_min;
    
    { /* Draw the entry slice */
        float3 uo_pos = 0;
        uo_pos[k] = origin[k] + ray_tl[k] * entry_tl;
        uo_pos[u] = u_max;
        uo_pos[v] = v_max;
        float3 vo_pos = 0;
        vo_pos[k] = origin[k] + ray_br[k] * entry_br;
        vo_pos[u] = v_min;
        vo_pos[v] = u_min;
        db::draw_line(u_pos, uo_pos, 0xFF00FF00);
        db::draw_line(v_pos, vo_pos, 0xFF00FF00);
        db::draw_line(vo_pos, u_pos, 0xFF00FF00);
        db::draw_line(uo_pos, v_pos, 0xFF00FF00);
    }

    slice = _mm_setr_ps(u_min, u_max, v_min, v_max);
    return true;
}

f32 CoherentPacked8x8::entry(const f32 ro, const f32 rd, const f32 min, const f32 max) const {
    bool sign = rd;
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
