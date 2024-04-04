#include "cvv.h"

#include "dev/debug.h"

void CoherentVoxelVolume::set_voxel(const u32 x, const u32 y, const u32 z, const u8 mat) {
    const u32 i = (z * grid_size.y * grid_size.x) + (y * grid_size.x) + x;
    voxels[i] = mat;
}

CoherentVoxelVolume::CoherentVoxelVolume(const float3& pos, const int3& grid_size, const f32 vpu)
    : bb(OBB(pos, float3(grid_size) / vpu)), grid_size(grid_size), vpu(vpu) {
    voxels = new u8[grid_size.x * grid_size.y * grid_size.z]{};

    /* Fill the grid with some noise */
    for (u32 z = 0; z < grid_size.z; z++) {
        for (u32 y = 0; y < grid_size.y; y++) {
            for (u32 x = 0; x < grid_size.x; x++) {
                const u32 i = (z * grid_size.y * grid_size.x) + (y * grid_size.x) + x;
#if 1
                const f32 noise =
                    noise3D((f32)x / grid_size.x, (f32)y / grid_size.y, (f32)z / grid_size.z);
                if (noise > 0.09f) {
                    set_voxel(x, y, z, 0xFF);
                } else {
                    set_voxel(x, y, z, 0x00);
                }
#else
                if (z % 8 == 0) {
                    set_voxel(x, y, z, 0xFF);
                } else {
                    set_voxel(x, y, z, 0x00);
                }
#endif
            }
        }
    }
    pivot = bb.size * 0.5f; /* center pivot */

    /* Initialize palette to 0,1,1,1 */
    palette = new u32[256];
    memset(palette, 0x00FFFFFF, sizeof(u32) * 256);
}

static f32 entry(const f32 ro, const f32 rd, const f32 min, const f32 max) {
    const bool sign = (rd < 0);
    const f32 bmin = sign ? max : min;
    return (bmin - ro) / rd;  // TODO: maybe try getting rid of division here?
}

static f32 axis_delta(const f32 axis_dir) { return fabs(1.0f / axis_dir); }

static f32 safe_entry(const f32 ro, const f32 rd, const f32 min, const f32 max) {
    const bool sign = (rd < 0);
    const f32 bmin = sign ? max : min;
    const f32 tmin = (bmin - ro) / rd;
    return fmaxf(0, tmin);  // TODO: maybe try getting rid of division here?
}

PacketHit8x8 CoherentVoxelVolume::intersect(const RayPacket8x8& packet,
                                              const bool debug) const {
    PacketHit8x8 hit;

    /* Grab the 4 corner rays */
    const float3 origin = TransformPosition(packet.rays[0].origin, bb.imodel);
    const float3 ray_tl = TransformVector(packet.rays[0].dir, bb.imodel);
    const float3 ray_tr = TransformVector(packet.rays[7].dir, bb.imodel);
    const float3 ray_bl = TransformVector(packet.rays[56].dir, bb.imodel);
    const float3 ray_br = TransformVector(packet.rays[63].dir, bb.imodel);
    const f32 upv = 1.0f / vpu;

    /* Packet major axis */
    const u32 k = major_axis(ray_tl);
    /* Remaining axis */
    const u32 u = (k + 1) < 3 ? (k + 1) : 0;
    const u32 v = ~(k | u) & 0b11;

    /* Bounding box min & max on major axis K */
    const f32 k_sign = -getsign(ray_tl[k]);
    const f32 k_min = bb.pos[k], k_max = bb.pos[k] + bb.size[k] - upv;

    /* Compute entry time of the corner rays along major axis K */
    const f32 entry_t = safe_entry(origin[k], ray_tl[k], k_min, k_max);

    /* Entry time along major axis K (floor + 0.5 is to align k_t to the grid) */
    f32 k_t = floorf((origin[k] + ray_tl[k] * entry_t) * vpu + 0.0001f);
    const f32 k_o = (k_t - origin[k] * vpu) * k_sign;

    //if (debug) {
    //    db::draw_line(origin, origin + packet.rays[0] * 1.5f, 0xFFFF0000);
    //    db::draw_line(origin, origin + packet.rays[63] * 1.5f, 0xFFFF0000);
    //    db::draw_line(origin, origin + packet.rays[7] * 1.5f, 0xFFFF0000);
    //    db::draw_line(origin, origin + packet.rays[56] * 1.5f, 0xFFFF0000);
    //}

    /* [du_min, du_max, dv_min, dv_max] */
    /* Factor by which the slice grows each step. */
    union {
        f128 delta_slice;
        struct {
            f32 du_min, du_max, dv_min, dv_max;
        };
    };
    /* [u_min, u_max, v_min, v_max] */
    /* Extend of the current slice. */
    union {
        f128 slice;
        struct {
            f32 u_min, u_max, v_min, v_max;
        };
    };

    /* Calculate the slice deltas */
    const f32 tl_dk = axis_delta(ray_tl[k]), br_dk = axis_delta(ray_br[k]);
    const f32 tl_du = ray_tl[u] * tl_dk, br_du = ray_br[u] * br_dk;
    const f32 tl_dv = ray_tl[v] * tl_dk, br_dv = ray_br[v] * br_dk;
    const f32 tr_dk = axis_delta(ray_tr[k]), bl_dk = axis_delta(ray_bl[k]);
    const f32 tr_du = ray_tr[u] * tr_dk, bl_du = ray_bl[u] * bl_dk;
    const f32 tr_dv = ray_tr[v] * tr_dk, bl_dv = ray_bl[v] * bl_dk;
    du_min = fminf(fminf(tl_du, br_du), fminf(tr_du, bl_du));
    du_max = fmaxf(fmaxf(tl_du, br_du), fmaxf(tr_du, bl_du));
    dv_min = fminf(fminf(tl_dv, br_dv), fminf(tr_dv, bl_dv));
    dv_max = fmaxf(fmaxf(tl_dv, br_dv), fmaxf(tr_dv, bl_dv));

    const f32 origin_u = origin[u] * vpu;
    const f32 origin_v = origin[v] * vpu;

    const f32 k_pos = -fminf(k_sign, 0) + 1;
    const f32 t = k_o - k_pos;

    /* Point 0 (adjust for k_max - 1) */
    const f32 u_a = origin_u + du_min * t;
    const f32 u_b = origin_u + du_max * t;
    const f32 v_a = origin_v + dv_max * t;
    const f32 v_b = origin_v + dv_min * t;

    /* Point 1 */
    const f32 u_c = origin_u + du_min * (t + 1);
    const f32 u_d = origin_u + du_max * (t + 1);
    const f32 v_c = origin_v + dv_max * (t + 1);
    const f32 v_d = origin_v + dv_min * (t + 1);

    /* Merge point 0 and 1 for our starting position */
    u_min = fminf(fminf(u_a, u_b), fminf(u_c, u_d));
    v_min = fminf(fminf(v_a, v_b), fminf(v_c, v_d));
    u_max = fmaxf(fmaxf(u_a, u_b), fmaxf(u_c, u_d));
    v_max = fmaxf(fmaxf(v_a, v_b), fmaxf(v_c, v_d));

    if (debug) {
        float3 min_p, max_p;
        min_p[k] = k_t, min_p[u] = u_min, min_p[v] = v_min;
        max_p[k] = k_t + 1, max_p[u] = u_max, max_p[v] = v_max;

        /* Draw floating point grid slice */
        db::draw_aabb(min_p * upv, max_p * upv, 0xFF00FF00);
    }

    /* Time to trace! */
    const f32 min_t = k_min * vpu;
    const f32 max_t = k_max * vpu;

    for (k_t; k_t >= min_t && k_t <= max_t; k_t += k_sign) {
        /* Move to the next slice */
        slice = _mm_add_ps(slice, delta_slice);

        /* Skip any slice 100% outside the grid */
        if (u_max < 0.0f || v_max < 0.0f) continue;
        if (u_min >= grid_size[u] || v_min >= grid_size[v]) continue;

        /* Truncate the floating point slice */
        const f128 fslice = _mm_floor_ps(slice);
        const i128 islice = _mm_cvttps_epi32(fslice);

        /* Grid slice extend */
        const i32 y_min = max(0, islice.m128i_i32[2]);
        const i32 y_max = min(grid_size[v] - 1, islice.m128i_i32[3]);
        const i32 x_min = max(0, islice.m128i_i32[0]);
        const i32 x_max = min(grid_size[u] - 1, islice.m128i_i32[1]);

        /* DEBUG: draw the grid slice */
        if (debug) {
            float3 min_p, max_p;
            min_p[k] = k_t, min_p[u] = u_min, min_p[v] = v_min;
            max_p[k] = k_t + 1, max_p[u] = u_max, max_p[v] = v_max;

            /* Draw floating point grid slice */
            db::draw_aabb(min_p * upv, max_p * upv, 0xFF0000FF);

            min_p[k] = k_t, min_p[u] = x_min, min_p[v] = y_min;
            max_p[k] = k_t + 1, max_p[u] = x_max + 1, max_p[v] = y_max + 1;

            /* Draw grid slice */
            db::draw_aabb(min_p * upv, max_p * upv, 0xFF00FF00);
        }

        /* Skip any slice 100% outside the grid */
        //if (u_max < 0.0f || v_max < 0.0f) continue;
        //if (u_min >= grid_size[u] || v_min >= grid_size[v]) continue;

        /* Iterate over cells in slice */
        for (i32 y = y_min; y <= y_max; y++) {
            for (i32 x = x_min; x <= x_max; x++) {
                int3 i; /* Cell coordinate */
                i[k] = k_t, i[u] = x, i[v] = y;

                // if (debug) {
                //     float3 min_p = i, max_p = i + 1;

                //    /* Draw floating point grid slice */
                //    db::draw_aabb(min_p * upv, max_p * upv, 0xFFFF0000);
                //}

                /* If the cell is a solid voxel */
                if (voxels[i.z * grid_size.y * grid_size.x + i.y * grid_size.x + i.x]) {
                    if (debug) {
                        float3 min_p = i, max_p = i + 1;

                        /* Draw floating point grid slice */
                        db::draw_aabb(min_p * upv, max_p * upv, 0xFFFF00FF);
                    }

                    /* Find which rays intersect this voxel */
                    for (u32 r = 0; r < 8 * 8; r++) {
                        // TODO: maybe don't do this transform every time?
                        const float3 rd = TransformVector(packet.rays[r].dir, bb.imodel);
                        const float3 ird = 1.0f / rd;
                        const float3 grid_o = origin * vpu;

                        /* Entry point */
                        const f32 entry_k = k_t - fminf(k_sign, 0.0f);
                        const f32 entry_t = entry(grid_o[k], rd[k], entry_k, entry_k) + 0.00001f;
                        float3 entry_p = grid_o + rd * entry_t;
                        int3 entry_c = floori(entry_p);

                        if (debug) {
                            float3 min_p = floorf(entry_p), max_p = min_p + 1;

                            /* Draw floating point grid slice */
                            db::draw_aabb(min_p * upv, max_p * upv, 0xFFFF0000);
                        }

                        /* DDA */
                        const int3 step = make_int3(-getsign(rd.x), -getsign(rd.y), -getsign(rd.z));
                        const float3 pos_step = fmaxf(step, 0);
                        const float3 delta = fabs(ird);
                        float3 tmax = ((float3(entry_c) - entry_p) + pos_step) * ird;
                        i32 axis = k;
                        f32 inner_t = 0;

                        /* Up to 5 DDA steps */
                        for (u32 q = 0; q < 5; q++) {
                            /* Stop if depth is higher or equal to previously found depth */
                            if (hit.hits[r].depth <= (entry_t + inner_t) * upv) break;

                            /* Hit! */
                            if (entry_c[u] == i[u] && entry_c[v] == i[v]) {
                                hit.hits[r].depth = (entry_t + inner_t) * upv;
                                /* Normal */
                                hit.hits[r].normal = 0, hit.hits[r].normal[axis] = -step[axis];
                                break;
                            }

                            /* Amanatides & Woo */
                            /* <http://www.cse.yorku.ca/~amana/research/grid.pdf> */
                            if (tmax[u] < tmax[v]) {
                                if (tmax[u] < tmax[k]) {
                                    entry_c[u] += step[u];
                                    inner_t = tmax[u];
                                    tmax[u] += delta[u];
                                    axis = u;
                                } else {
                                    axis = k;
                                    break;
                                }
                            } else {
                                if (tmax[v] < tmax[k]) {
                                    entry_c[v] += step[v];
                                    inner_t = tmax[v];
                                    tmax[v] += delta[v];
                                    axis = v;
                                } else {
                                    axis = k;
                                    break;
                                }
                            }
                        }
                    }
                }
            }
        }

        u32 inactive_rays = 0;
        for (u32 r = 0; r < 8 * 8; r++) {
            /* Only check active rays */
            if (hit.hits[r].depth != BIG_F32) {
                inactive_rays++;
            }
        }

        /* Early exit if all rays are done */
        if (inactive_rays == (8 * 8)) {
            break;
        }
    }

    /* Finalize the normals */
    for (u32 r = 0; r < 8 * 8; r++) {
        if (hit.hits[r].depth != BIG_F32) {
            hit.hits[r].normal = normalize(TransformVector(hit.hits[r].normal, bb.model));
        }
    }

    // TODO: whenever this happens, artifacts apear because of a misalingment.
    // if (entry_br == 0 or entry_tl == 0 or entry_tr == 0 or entry_bl == 0) {
    //    // if (debug) {
    //    //     db::draw_normal(0, 0);
    //    // }
    //    for (u32 r = 0; r < 8 * 8; r++) {
    //        hit.normal[r] = 0;
    //    }
    //    return hit;
    //}

    return hit;
}
