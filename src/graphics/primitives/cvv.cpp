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
    const f32 sign_u = fmaxf(getsign(ray_tl[u]), 0);
    const f32 sign_v = fmaxf(getsign(ray_tl[v]), 0);
    // const f32 min_t = (sign ? k_min : k_max) * vpu;
    // const f32 max_t = (sign ? k_max : k_min) * vpu;
    const f32 min_t = k_min * vpu;
    const f32 max_t = k_max * vpu;
    k_t += sign * 0.01f;

    // if (debug) db::draw_aabb(cell_min, cell_max, 0xFFFF0000);

    for (k_t; k_t > min_t && k_t < max_t; k_t += sign) {
        slice = _mm_add_ps(slice, delta_slice);
        const i128 islice = _mm_cvtps_epi32(slice);

        float3 min_p, max_p;
        min_p[k] = k_t, min_p[u] = islice.m128i_i32[0], min_p[v] = islice.m128i_i32[2];
        max_p[k] = k_t, max_p[u] = islice.m128i_i32[1], max_p[v] = islice.m128i_i32[3];

        /* Skip any slice 100% outside the grid */
        // if (max_p[u] < 0 || max_p[v] < 0) continue;
        // if (min_p[u] > grid_size[u] || min_p[v] > grid_size[v]) continue;

        /* Draw the grid slice */
        // if (debug) {
        //     /* DEBUG DRAW */
        //     float3 a_pos, b_pos;
        //     a_pos[k] = k_t, a_pos[u] = islice.m128i_i32[0], a_pos[v] = islice.m128i_i32[2];
        //     b_pos[k] = k_t, b_pos[u] = islice.m128i_i32[1], b_pos[v] = islice.m128i_i32[3];
        //     b_pos[k] += sign, b_pos[u] += 1, b_pos[v] += 1;

        //    /* The minimum and maximum cell of the slice in the grid */
        //    const float3 cell_min = min_p * upv;
        //    const float3 cell_max = max_p * upv;

        //    db::draw_aabb(cell_min, cell_max, 0xFFFF0000);
        //}

        /* Clamp the slice inside the grid extend */
        const f32 y_min = islice.m128i_i32[2];  // fmaxf(0, min_p[v]);
        const f32 y_max = islice.m128i_i32[3];  // fminf(grid_size[v], max_p[v]);
        const f32 x_min = islice.m128i_i32[0];  // fmaxf(0, min_p[u]);
        const f32 x_max = islice.m128i_i32[1];  // fminf(grid_size[u], max_p[u]);

        // if (debug) {
        //     float3 cell_min, cell_max;
        //     cell_min[k] = k_t, cell_min[u] = x_min, cell_min[v] = y_min;
        //     cell_max[k] = k_t + sign, cell_max[u] = x_max + 1, cell_max[v] = y_max + 1;

        //    db::draw_aabb(cell_min * upv, cell_max * upv, 0xFFFF0000);
        //}

        /* Iterate over cells in slice */
        for (u32 y = y_min; y <= y_max; y++) {
            for (u32 x = x_min; x <= x_max; x++) {
                uint3 i; /* Cell coordinate */
                i[k] = k_t, i[u] = x, i[v] = y;

                /* DEBUG DRAW */
                if (debug) {
                    float3 c_pos, d_pos;
                    c_pos[k] = k_t, c_pos[u] = x, c_pos[v] = y;
                    d_pos[k] = k_t, d_pos[u] = x, d_pos[v] = y;
                    d_pos[k] += sign, d_pos[u] += 1, d_pos[v] += 1;
                    db::draw_aabb(c_pos * upv, d_pos * upv, 0xFF00FF00);
                }

                if (x < 0 || y < 0) continue;
                if (x >= grid_size[u] || y >= grid_size[v]) continue;

                /* If the cell is a solid voxel */
                if (voxels[i.z * grid_size.y * grid_size.x + i.y * grid_size.x + i.x]) {
                    // TESTING
                    for (u32 r = 0; r < 4 * 4; r++) {
                        hit.depth[r] = 100.0f;
                    }
                    return hit;
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

                        /* If the ray (entry) intersects this voxel */
                        if (x == entry_p[u] || y == entry_p[v]) {
                            hit.depth[r] = entry_t;
                            // TODO: store normal (which in this case is major axis K)
                            continue;
                        }

                        /* If the ray (exit) intersects this voxel */
                        if (x == exit_p[u] || y == exit_p[v]) {
                            /* TEMP: is incorrect! */
                            hit.depth[r] = entry_t;
                            // TODO: store normal (diff between entry and exit point)
                        }
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
