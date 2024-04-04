#include "precomp.h"

#include "dev/debug.h"

void CoherentPacket4x4::setup_slice(const float3& min, const float3& max, const f32 vpu) {
    /* Grab the 2 corner rays */
    const float3 ray_tl = rays[0], ray_br = rays[15];
    const f32 upv = 1.0f / vpu;

    /* Packet major axis */
    k = major_axis(ray_tl);
    /* Remaining axis */
    u = (k + 1) < 3 ? (k + 1) : 0;
    v = ~(k | u) & 0b11;

    /* Compute entry time of top left & bottom left rays along major axis K */
    const f32 entry_br = entry(origin[k], ray_br[k], min[k], max[k]);
    const f32 entry_tl = entry(origin[k], ray_tl[k], min[k], max[k]);

    /* Compute next entry time of top left & bottom left rays along major axis K */
    const f32 next_br = entry(origin[k], ray_br[k], min[k] + upv, max[k] - upv);
    const f32 next_tl = entry(origin[k], ray_tl[k], min[k] + upv, max[k] - upv);

    /* Entry time along major axis K */
    k_t = (origin[k] + ray_br[k] * entry_br) * vpu;

    /* Top left & bottom right U,V entry points */
    const f32 u_tl = (origin[u] + ray_tl[u] * entry_tl) * vpu;
    const f32 u_br = (origin[u] + ray_br[u] * entry_br) * vpu;
    const f32 v_tl = (origin[v] + ray_tl[v] * entry_tl) * vpu;
    const f32 v_br = (origin[v] + ray_br[v] * entry_br) * vpu;

    /* Next top left & bottom right U,V bounding points */
    const f32 nu_tl = (origin[u] + ray_tl[u] * next_tl) * vpu;
    const f32 nu_br = (origin[u] + ray_br[u] * next_br) * vpu;
    const f32 nv_tl = (origin[v] + ray_tl[v] * next_tl) * vpu;
    const f32 nv_br = (origin[v] + ray_br[v] * next_br) * vpu;

    /* Entry min and max U,V */
    u_min = fminf(u_tl, u_br), u_max = fmaxf(u_tl, u_br);
    v_min = fminf(v_tl, v_br), v_max = fmaxf(v_tl, v_br);

    /* Next min and max U,V */
    const f32 nu_min = fminf(nu_tl, nu_br), nu_max = fmaxf(nu_tl, nu_br);
    const f32 nv_min = fminf(nv_tl, nv_br), nv_max = fmaxf(nv_tl, nv_br);

    /* Slice delta U,V */
    du_min = nu_min - u_min, du_max = nu_max - u_max;
    dv_min = nv_min - v_min, dv_max = nv_max - v_max;

    /* Slice entry U,V */
    u_min = fminf(u_min, nu_min), u_max = fmaxf(u_max, nu_max);
    v_min = fminf(v_min, nv_min), v_max = fmaxf(v_max, nv_max);
}

void CoherentPacket4x4::traverse(const float3& min, const float3& max, const f32 vpu) {
    const f32 sign = -getsign(rays[0][k]);
    // const f32 upv = (1.0f / vpu) * sign;
    const f32 min_t = (sign ? min[k] : max[k]) * vpu;
    const f32 max_t = (sign ? max[k] : min[k]) * vpu;
    k_t += sign * 0.01f;
    for (k_t; k_t >= min_t && k_t < max_t; k_t += sign) {
        slice = _mm_add_ps(slice, delta_slice);
        draw_slice(vpu);
    }
}

void CoherentPacket4x4::draw_slice(const f32 vpu) const {
    float3 a_pos = 0;
    a_pos[k] = k_t, a_pos[u] = u_min, a_pos[v] = v_min;
    float3 b_pos = 0;
    b_pos[k] = k_t, b_pos[u] = u_max, b_pos[v] = v_max;

    const f32 sign = -getsign(rays[0][k]);
    const f32 upv = 1.0f / vpu;
    b_pos[k] += sign;
    /* Draw the floating point slice */
    db::draw_aabb(a_pos * upv, b_pos * upv);

    const i128 islice = _mm_cvttps_epi32(slice);

    a_pos[k] = k_t, a_pos[u] = islice.m128i_u32[0], a_pos[v] = islice.m128i_u32[2];
    b_pos[k] = k_t, b_pos[u] = islice.m128i_u32[1], b_pos[v] = islice.m128i_u32[3];
    b_pos[k] += 1, b_pos[u] += 1, b_pos[v] += 1;

    /* The minimum and maximum cell of the slice in the grid */
    const float3 cell_min = a_pos * upv;
    const float3 cell_max = b_pos * upv;

    /* Draw the grid slice */
    db::draw_aabb(cell_min, cell_max, 0xFFFF0000);
}

f32 CoherentPacket4x4::entry(const f32 ro, const f32 rd, const f32 min, const f32 max) const {
    const bool sign = (rd < 0);
    const f32 bmin = sign ? max : min;
    return (bmin - ro) / rd;  // TODO: maybe try getting rid of division here?
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
