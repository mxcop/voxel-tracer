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

CoherentHit4x4 CoherentVoxelVolume::intersect(const CoherentPacket4x4& packet,
                                              const bool debug) const {
    CoherentHit4x4 hit;

    /* Grab the 2 corner rays */
    const float3 origin = TransformPosition(packet.origin, bb.imodel);
    const float3 ray_tl = TransformVector(packet.rays[0], bb.imodel);
    const float3 ray_br = TransformVector(packet.rays[15], bb.imodel);
    const f32 upv = 1.0f / vpu;

    /* Packet major axis */
    const u32 k = major_axis(ray_tl);
    /* Remaining axis */
    const u32 u = (k + 1) < 3 ? (k + 1) : 0;
    const u32 v = ~(k | u) & 0b11;

    /* Bounding box min & max on major axis K */
    const f32 k_min = bb.pos[k], k_max = bb.pos[k] + bb.size[k];

    /* Compute entry time of top left & bottom left rays along major axis K */
    const f32 entry_br = entry(origin[k], ray_br[k], k_min, k_max);
    const f32 entry_tl = entry(origin[k], ray_tl[k], k_min, k_max);

    /* Compute next entry time of top left & bottom left rays along major axis K */
    const f32 next_br = entry(origin[k], ray_br[k], k_min + upv, k_max - upv);
    const f32 next_tl = entry(origin[k], ray_tl[k], k_min + upv, k_max - upv);

    /* Entry time along major axis K */
    f32 k_t = (origin[k] + ray_tl[k] * entry_tl) * vpu;

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

    if (debug) {
        db::draw_line(packet.origin, packet.origin + packet.rays[0] * 2.0f, 0xFFFF0000);
        db::draw_line(packet.origin, packet.origin + packet.rays[15] * 2.0f, 0xFFFF0000);
        db::draw_line(packet.origin, packet.origin + packet.rays[3] * 2.0f, 0xFFFF0000);
        db::draw_line(packet.origin, packet.origin + packet.rays[12] * 2.0f, 0xFFFF0000);
    }

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

    /* Entry min and max U,V */
    u_min = fminf(u_tl, u_br), u_max = fmaxf(u_tl, u_br);
    v_min = fminf(v_tl, v_br), v_max = fmaxf(v_tl, v_br);

    if (debug) {
        float3 min_p, max_p;
        min_p[k] = k_t, min_p[u] = u_min, min_p[v] = v_min;
        max_p[k] = k_t, max_p[u] = u_max, max_p[v] = v_max;

        /* Draw floating point grid slice */
        db::draw_aabb(min_p * upv, max_p * upv, 0xFFFF00FF);
    }

    /* Next min and max U,V */
    const f32 nu_min = fminf(nu_tl, nu_br), nu_max = fmaxf(nu_tl, nu_br);
    const f32 nv_min = fminf(nv_tl, nv_br), nv_max = fmaxf(nv_tl, nv_br);

    if (debug) {
        float3 min_p, max_p;
        min_p[k] = k_t, min_p[u] = nu_min, min_p[v] = nv_min;
        max_p[k] = k_t, max_p[u] = nu_max, max_p[v] = nv_max;

        /* Draw floating point grid slice */
        db::draw_aabb(min_p * upv, max_p * upv, 0xFFFF00FF);
    }

    /* Slice delta U,V */
    du_min = nu_min - u_min, du_max = nu_max - u_max;
    dv_min = nv_min - v_min, dv_max = nv_max - v_max;

    /* Slice entry U,V */
    u_min = fminf(u_min, nu_min), u_max = fmaxf(u_max, nu_max);
    v_min = fminf(v_min, nv_min), v_max = fmaxf(v_max, nv_max);

    /* Time to trace! */
    const f32 sign = -getsign(ray_tl[k]);
    const f32 min_t = k_min * vpu;
    const f32 max_t = k_max * vpu;
    k_t += sign * 0.01f;

    /* Move back by 1 slice, because we first move then check! */
    slice = _mm_sub_ps(slice, delta_slice);

    for (k_t; k_t > min_t && k_t < max_t; k_t += sign) {
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
            max_p[k] = k_t, max_p[u] = u_max, max_p[v] = v_max;

            /* Draw floating point grid slice */
            db::draw_aabb(min_p * upv, max_p * upv, 0xFF0000FF);

            min_p[k] = k_t, min_p[u] = x_min, min_p[v] = y_min;
            max_p[k] = k_t + sign, max_p[u] = x_max + 1, max_p[v] = y_max + 1;

            /* Draw grid slice */
            db::draw_aabb(min_p * upv, max_p * upv, 0xFF00FF00);
        }

        const int3 step = make_int3(-getsign(ray_tl.x), -getsign(ray_tl.y), -getsign(ray_tl.z));
        const float3 pos_step = fmaxf(step, 0);

        /* Iterate over cells in slice */
        u32 inactive_rays = 0;
        for (i32 y = y_min; y <= y_max; y++) {
            for (i32 x = x_min; x <= x_max; x++) {
                int3 i; /* Cell coordinate */
                i[k] = k_t, i[u] = x, i[v] = y;

                /* If the cell is a solid voxel */
                if (voxels[i.z * grid_size.y * grid_size.x + i.y * grid_size.x + i.x]) {
                    /* Find which rays intersect this voxel */
                    inactive_rays = 0;
                    for (u32 r = 0; r < 4 * 4; r++) {
                        /* Only check active rays */
                        if (hit.depth[r] != BIG_F32) {
                            inactive_rays++;
                            // continue;
                        }

                        // TODO: maybe don't do this transform every time?
                        const float3 rd = TransformVector(packet.rays[r], bb.imodel);
                        const float3 ird = 1.0f / rd;
                        const float3 grid_o = origin * vpu;

                        /* Entry point */
                        const f32 entry_t = entry(grid_o[k], rd[k], k_t, k_t);
                        const float3 entry_p = grid_o + rd * entry_t;
                        int3 entry_c = floori(entry_p);
                        entry_c[k] = k_t;

                        /* DDA */
                        const float3 delta = fabs(ird);
                        float3 tmax = ((float3(entry_c) - entry_p) + pos_step) * ird;
                        u32 axis = k;
                        f32 inner_t = 0;

                        /* Up to 5 DDA steps */
                        for (u32 q = 0; q < 5; q++) {
                            /* Stop if depth is higher or equal to previously found depth */
                            if (hit.depth[r] <= (entry_t + inner_t) * upv) break;

                            /* Hit! */
                            if (entry_c[u] == i[u] && entry_c[v] == i[v]) {
                                inactive_rays++;
                                hit.depth[r] = (entry_t + inner_t) * upv;
                                /* Normal */
                                hit.normal[r] = 0, hit.normal[r][axis] = -step[axis];
                                break;
                            }

                            /* Amanatides & Woo */
                            /* <http://www.cse.yorku.ca/~amana/research/grid.pdf> */
                            if (tmax[u] < tmax[v]) {
                                if (tmax[u] < tmax[k]) {
                                    entry_c[u] += step[u];
                                    axis = u;
                                    inner_t = tmax[u];
                                    tmax[u] += delta[u];
                                } else {
                                    break;
                                }
                            } else {
                                if (tmax[v] < tmax[k]) {
                                    entry_c[v] += step[v];
                                    axis = v;
                                    inner_t = tmax[v];
                                    tmax[v] += delta[v];
                                } else {
                                    break;
                                }
                            }
                        }
                    }
                }
            }
        }

        /* Early exit if all rays are done */
        if (inactive_rays == (4 * 4)) {
            return hit;
        }
    }

    /* Finalize the normals */
    for (u32 r = 0; r < 4 * 4; r++) {
        if (hit.depth[r] != BIG_F32) {
            hit.normal[r] = normalize(TransformVector(hit.normal[r], bb.model));
        }
    }

    return hit;
}

CoherentHit8x8 CoherentVoxelVolume::intersect(const CoherentPacket8x8& packet,
                                              const bool debug) const {
    CoherentHit8x8 hit;

    /* Grab the 4 corner rays */
    const float3 origin = TransformPosition(packet.origin, bb.imodel);
    const float3 ray_tl = TransformVector(packet.rays[0], bb.imodel);
    const float3 ray_tr = TransformVector(packet.rays[7], bb.imodel);
    const float3 ray_bl = TransformVector(packet.rays[56], bb.imodel);
    const float3 ray_br = TransformVector(packet.rays[63], bb.imodel);
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

    if (debug) {
        db::draw_line(packet.origin, packet.origin + packet.rays[0] * 1.5f, 0xFFFF0000);
        db::draw_line(packet.origin, packet.origin + packet.rays[63] * 1.5f, 0xFFFF0000);
        db::draw_line(packet.origin, packet.origin + packet.rays[7] * 1.5f, 0xFFFF0000);
        db::draw_line(packet.origin, packet.origin + packet.rays[56] * 1.5f, 0xFFFF0000);
    }

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
                        const float3 rd = TransformVector(packet.rays[r], bb.imodel);
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
                            if (hit.depth[r] <= (entry_t + inner_t) * upv) break;

                            /* Hit! */
                            if (entry_c[u] == i[u] && entry_c[v] == i[v]) {
                                hit.depth[r] = (entry_t + inner_t) * upv;
                                /* Normal */
                                hit.normal[r] = 0, hit.normal[r][axis] = -step[axis];
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
            if (hit.depth[r] != BIG_F32) {
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
        if (hit.depth[r] != BIG_F32) {
            hit.normal[r] = normalize(TransformVector(hit.normal[r], bb.model));
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
