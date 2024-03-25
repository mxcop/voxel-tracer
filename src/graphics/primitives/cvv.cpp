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
                const f32 noise =
                    noise3D((f32)x / grid_size.x, (f32)y / grid_size.y, (f32)z / grid_size.z);
                if (noise > 0.09f) {
                    set_voxel(x, y, z, 0xFF);
                } else {
                    set_voxel(x, y, z, 0x00);
                }
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

    /* Next min and max U,V */
    const f32 nu_min = fminf(nu_tl, nu_br), nu_max = fmaxf(nu_tl, nu_br);
    const f32 nv_min = fminf(nv_tl, nv_br), nv_max = fmaxf(nv_tl, nv_br);

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
        const i32 y_min = islice.m128i_i32[2];
        const i32 y_max = islice.m128i_i32[3];
        const i32 x_min = islice.m128i_i32[0];
        const i32 x_max = islice.m128i_i32[1];

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

        /* Iterate over cells in slice */
        for (i32 y = y_min; y <= y_max; y++) {
            for (i32 x = x_min; x <= x_max; x++) {
                uint3 i; /* Cell coordinate */
                i[k] = k_t, i[u] = x, i[v] = y;

                // TODO: maybe avoid these checks? using clamp?
                /* Safety bounds check */
                if (x < 0 || y < 0) continue;
                if (x >= grid_size[u] || y >= grid_size[v]) continue;

                /* If the cell is a solid voxel */
                if (voxels[i.z * grid_size.y * grid_size.x + i.y * grid_size.x + i.x]) {
                    // TESTING
                    // for (u32 r = 0; r < 4 * 4; r++) {
                    //    hit.depth[r] = k_t;
                    //}
                    // return hit;
                    // TESTING

                    /* Find which rays intersect this voxel */
                    u32 inactive_rays = 0;
                    for (u32 r = 0; r < 4 * 4; r++) {
                        /* Only check active rays */
                        if (hit.depth[r] != BIG_F32) {
                            inactive_rays++;
                            continue;
                        }

                        // TODO: maybe don't do this transform every time?
                        const float3 rd = TransformVector(packet.rays[r], bb.imodel);
                        const float3 grid_o = origin * vpu;
                        const f32 entry_t = entry(grid_o[k], rd[k], k_t, k_t);
                        const f32 exit_t = entry(grid_o[k], rd[k], k_t + sign, k_t + sign);

                        /* Ray entry and exit point in the grid slice */
                        const int3 entry_p = floori(grid_o + rd * entry_t);
                        const int3 exit_p = floori(grid_o + rd * exit_t);

                        /* DEBUG DRAW */
                        if (debug) {
                            db::draw_line((grid_o + rd * entry_t) * upv,
                                          (grid_o + rd * exit_t) * upv,
                                          0xFFFF00FF);
                            db::draw_aabb(float3(entry_p) * upv, float3(entry_p + 1) * upv,
                                          0xFFFFFF00);
                            db::draw_aabb(float3(exit_p) * upv, float3(exit_p + 1) * upv,
                                          0xFFFFFF00);
                        }

                        /* THIS METHOD IS ALSO NOT CORRECT!!! */
                        const i32 pu_min = fminf(entry_p[u], exit_p[u]);
                        const i32 pu_max = fmaxf(entry_p[u], exit_p[u]);
                        const i32 pv_min = fminf(entry_p[v], exit_p[v]);
                        const i32 pv_max = fmaxf(entry_p[v], exit_p[v]);

                        if (x >= pu_min && x <= pu_max) {
                            hit.depth[r] = entry_t;
                            inactive_rays++;
                            continue;
                        }
                        if (y >= pv_min && y <= pv_max) {
                            hit.depth[r] = entry_t;
                            inactive_rays++;
                            continue;
                        }

                        ///* If the ray (entry) intersects this voxel */
                        //if (x == entry_p[u] && y == entry_p[v]) {
                        //    hit.depth[r] = entry_t;
                        //    inactive_rays++;
                        //    // TODO: store normal (which in this case is major axis K)
                        //    continue;
                        //}

                        ///* If the ray (exit) intersects this voxel */
                        //if (x == exit_p[u] && y == exit_p[v]) {
                        //    /* TEMP: is incorrect! */
                        //    hit.depth[r] = entry_t;
                        //    inactive_rays++;
                        //    // TODO: store normal (diff between entry and exit point)
                        //    continue;
                        //}
                    }

                    /* Early exit if all rays are done */
                    if (inactive_rays == (4 * 4)) {
                        return hit;
                    }
                }
            }
        }
    }

    return hit;
}
