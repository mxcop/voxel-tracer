#include "vv.h"

#define OGT_VOX_IMPLEMENTATION
#include <ogt/ogt_vox.h>
#include <dev/debug.h>

constexpr u32 MAX_STEPS = 256;

/* Convert RGBA to u32 */
static inline u32 merge_u8_u32(ogt_vox_rgba c) { return (c.r << 16) | (c.g << 8) | c.b; }

OVoxelVolume::OVoxelVolume(const float3& pos, const char* vox_path, const i32 model_id,
                           const f32 vpu)
    : vpu(vpu) {
    /* Load the model file */
    FILE* fp = fopen(vox_path, "rb");
    uint32_t buffer_size = _filelength(_fileno(fp));
    uint8_t* buffer = new uint8_t[buffer_size];
    fread(buffer, buffer_size, 1, fp);
    fclose(fp);

    /* Parse the model file */
    const ogt_vox_scene* scene = ogt_vox_read_scene(buffer, buffer_size);
    delete[] buffer; /* Cleanup */

    /* Ignore everything except the first model */
    const ogt_vox_model* model = scene->models[model_id];

    /* Initialize the volume */
    grid_size = int3(model->size_y, model->size_z, model->size_x);
    brickmap = Brickmap(grid_size);
    // #if USE_BITPACKING
    voxels = new u8[grid_size.x * grid_size.y * grid_size.z]{};
    // #endif
    bb = OBB(pos, float3(grid_size) / vpu);
    pivot = bb.size * 0.5f; /* center pivot */

    /* Format the grid */
    for (u32 z = 0; z < model->size_z; z++) {
        for (u32 y = 0; y < model->size_y; y++) {
            for (u32 x = 0; x < model->size_x; x++) {
                /* .vox model axis are weird... */
                const u32 mi = (z * model->size_y * model->size_x) +
                               ((model->size_y - y - 1) * model->size_x) + x;

                set_voxel(int3(y, z, x), model->voxel_data[mi]);
            }
        }
    }

    /* Initialize the grid palette */
    palette = new u32[256];
    for (u32 c = 0; c < 256; ++c) palette[c] = merge_u8_u32(scene->palette.color[c]);
}

OVoxelVolume::OVoxelVolume(const float3& pos, const int3& grid_size, const f32 vpu)
    : bb(OBB(pos, float3(grid_size) / vpu)),
      brickmap(Brickmap(grid_size)),
      grid_size(grid_size),
      vpu(vpu) {
    // #if USE_BITPACKING
    voxels = new u8[grid_size.x * grid_size.y * grid_size.z]{};
    // #endif

    /* Fill the grid with some noise */
    for (u32 z = 0; z < grid_size.z; z++) {
        for (u32 y = 0; y < grid_size.y; y++) {
            for (u32 x = 0; x < grid_size.x; x++) {
                const u32 i = (z * grid_size.y * grid_size.x) + (y * grid_size.x) + x;
                const f32 noise =
                    noise3D((f32)x / grid_size.x, (f32)y / grid_size.y, (f32)z / grid_size.z);
                if (noise > 0.09f) {
                    set_voxel(int3(x, y, z), 16);
                } else {
                    set_voxel(int3(x, y, z), 0);
                }
            }
        }
    }
    pivot = bb.size * 0.5f; /* center pivot */

    /* Initialize palette to 0,1,1,1 */
    palette = new u32[256];
    memset(palette, 0x00FFFFFF, sizeof(u32) * 256);
}

/* Trace-able functions */
AABB OVoxelVolume::get_aabb() const { return bb.get_aabb(); }
float3 OVoxelVolume::center() const { return bb.center(); }

void OVoxelVolume::set_position(const float3& pos) { bb.set_position(pos); }

void OVoxelVolume::set_rotation(const quat& rot) { bb.set_rotation_pivot(pivot, rot); }

HitInfo OVoxelVolume::intersect(const Ray& ray) const {
    HitInfo hit = bb.intersect(ray);

    /* Traverse the voxel grid if we hit the bounding box */
    if (hit.depth != BIG_F32) {
        const Ray grid_ray = bb.world_to_local(ray);

        /* Reciprocal bricks per unit */
        const f32 bpu = vpu / 8, rbpu = 1.0f / bpu;
        const float3 entry = (grid_ray.origin + grid_ray.dir * hit.depth) * bpu;

        /* Clamp the entry point inside the volume grid */
        int3 cell = clamp(floori(entry), 0, brickmap.size - 1);

        /* Which direction each axis will step in -1 or 1 */
        const int3 step = make_int3(grid_ray.sign_dir);
        /* Indicates how far we must move (in units of t) to equal the width of a voxel */
        const float3 delta = fabs(grid_ray.r_dir);
        /* Determine t at which the ray crosses the first voxel boundary */
        float3 tmax = ((float3(cell) - entry) + fmaxf(step, float3(0))) * grid_ray.r_dir;

        f32 t = 0.0f;
        u32 axis = 0;
        for (hit.steps = 0; hit.steps < MAX_STEPS; ++hit.steps) {
            /* Fetch the active cell */
            const Brick512* brick = get_brick(cell);
            if (brick->voxcnt > 0) {
                const f32 brick_entry_t = hit.depth + t * rbpu;
                const f32 dist =
                    traverse_brick(brick, cell, grid_ray, brick_entry_t, rbpu, axis, hit);
                if (dist != BIG_F32) {
                    hit.depth = brick_entry_t + dist;
                    if (hit.steps == 0) return hit;

                    hit.normal = 0, hit.normal[axis] = 1;
                    hit.normal = -hit.normal * step;
                    hit.normal = normalize(TransformVector(hit.normal, bb.model));
                    return hit;
                }
            } else if (ray.medium_id) {
                hit.depth = hit.depth + t * rbpu;
                if (hit.steps == 0) return hit;

                hit.normal = 0, hit.normal[axis] = 1;
                hit.normal = -hit.normal * step;
                hit.normal = normalize(TransformVector(hit.normal, bb.model));
                return hit;
            }

            /* Amanatides & Woo */
            /* <http://www.cse.yorku.ca/~amana/research/grid.pdf> */
            if (tmax.x < tmax.y) {
                if (tmax.x < tmax.z) {
                    cell.x += step.x;
                    if (cell.x < 0 || cell.x >= brickmap.size.x) break;
                    axis = 0, t = tmax.x;
                    tmax.x += delta.x;
                } else {
                    cell.z += step.z;
                    if (cell.z < 0 || cell.z >= brickmap.size.z) break;
                    axis = 2, t = tmax.z;
                    tmax.z += delta.z;
                }
            } else {
                if (tmax.y < tmax.z) {
                    cell.y += step.y;
                    if (cell.y < 0 || cell.y >= brickmap.size.y) break;
                    axis = 1, t = tmax.y;
                    tmax.y += delta.y;
                } else {
                    cell.z += step.z;
                    if (cell.z < 0 || cell.z >= brickmap.size.z) break;
                    axis = 2, t = tmax.z;
                    tmax.z += delta.z;
                }
            }
        }

        if (ray.medium_id) {
            hit.material = 0x00; /* Air */
            if (tmax.x < tmax.y) {
                if (tmax.x < tmax.z) {
                    axis = 0;
                } else {
                    axis = 2;
                }
            } else {
                if (tmax.y < tmax.z) {
                    axis = 1;
                } else {
                    axis = 2;
                }
            }
            hit.depth = bb.exit_t(ray);
            hit.normal = 0, hit.normal[axis] = 1;
            hit.normal = -hit.normal * step;
            hit.normal = normalize(TransformVector(hit.normal, bb.model));
            return hit;
        }

        /* No hit occured! */
        hit.depth = BIG_F32;
    } else if (ray.medium_id) {
        hit.material = 0x00; /* Air */
        hit.depth = 0;
    }

    return hit;
}

f32 OVoxelVolume::traverse_brick(const Brick512* brick, const int3& pos, const Ray& ray,
                                 const f32 entry_t, const f32 rbpu, u32& axis, HitInfo& hit) const {
    /* Brick minimum position */
    const float3 bmin = float3(pos) * rbpu;
    /* Brick entry position */
    const float3 entry = ((ray.origin + ray.dir * (entry_t)) - bmin) * vpu;
    /* Clamp the entry point inside the volume grid */
    int3 cell = clamp(floori(entry), 0, 8 - 1);

    /* Which direction each axis will step in -1 or 1 */
    const float3 step = ray.sign_dir;
    /* Indicates how far we must move (in units of t) to equal the width of a voxel */
    const float3 delta = fabs(ray.r_dir);
    /* Determine t at which the ray crosses the first voxel boundary */
    float3 tmax = ((float3(cell) - entry) + fmaxf(step, float3(0))) * ray.r_dir;

    bool exited = false;

    f32 t = 0.0f;
    for (hit.steps; hit.steps < MAX_STEPS; ++hit.steps) {
        /* Fetch the active cell */
        const u8 voxel = get_voxel(brick, cell);
        // TODO: all this logic is not great...
#if USE_BITPACKING
        if (ray.medium_id) {
            const int3 hitc = (pos << 3) + cell;
            const u32 i = (hitc.z * grid_size.y * grid_size.x) + (hitc.y * grid_size.x) + hitc.x;
            const u8 voxel_id = voxels[i];
            if (voxel_id != ray.medium_id) {
                hit.albedo = RGB8_to_RGBF32(palette[voxel_id]);
                hit.material = voxel_id;
                return t / vpu;
            }
        } else if (voxel) {
            const int3 hitc = (pos << 3) + cell;
            const u32 i = (hitc.z * grid_size.y * grid_size.x) + (hitc.y * grid_size.x) + hitc.x;
            const u8 voxel_id = voxels[i];
            if (ray.shadow_ray) {
                // TODO: this is not great...
                if (voxel_id > 16) {
                    hit.albedo = RGB8_to_RGBF32(palette[voxel_id]);
                    hit.material = voxel_id;
                    return t / vpu;
                } else {
                    /* Hacky stohastic absorption */
                    if (RandomFloat() > 0.85f) {
                        hit.albedo = RGB8_to_RGBF32(palette[voxel_id]);
                        hit.material = voxel_id;
                        return t / vpu;
                    }
                }
            } else if (exited || voxel_id != ray.ignore_medium) {
                hit.albedo = RGB8_to_RGBF32(palette[voxel_id]);
                hit.material = voxel_id;
                return t / vpu;
            }
        } else if (ray.ignore_medium) {
            exited = true;
        }
#else
        if (ray.medium_id) {
            // const int3 hitc = (pos << 3) + cell;
            // const u32 i = (hitc.z * grid_size.y * grid_size.x) + (hitc.y * grid_size.x) + hitc.x;
            // const u8 voxel_id = voxels[i];
            if (voxel != ray.medium_id) {
                // const int3 hitc = (pos << 3) + cell;
                // const u32 i =
                //     (hitc.z * grid_size.y * grid_size.x) + (hitc.y * grid_size.x) + hitc.x;
                // const u8 voxel_id = voxels[i];
                hit.albedo = RGB8_to_RGBF32(palette[voxel]);
                hit.material = voxel;
                return t / vpu;
            }
        } else if (voxel) {
            // const int3 hitc = (pos << 3) + cell;
            // const u32 i = (hitc.z * grid_size.y * grid_size.x) + (hitc.y * grid_size.x) + hitc.x;
            // const u8 voxel_id = voxels[i];
            if (ray.shadow_ray) {
                // TODO: this is not great...
                if (voxel > 16) {
                    hit.albedo = RGB8_to_RGBF32(palette[voxel]);
                    hit.material = voxel;
                    return t / vpu;
                } else {
                    /* Hacky stohastic absorption */
                    if (RandomFloat() > 0.85f) {
                        hit.albedo = RGB8_to_RGBF32(palette[voxel]);
                        hit.material = voxel;
                        return t / vpu;
                    }
                }
            } else if (exited || voxel != ray.ignore_medium) {
                hit.albedo = RGB8_to_RGBF32(palette[voxel]);
                hit.material = voxel;
                return t / vpu;
            }
        } else if (ray.ignore_medium) {
            exited = true;
        }
#endif

        /* Amanatides & Woo */
        /* <http://www.cse.yorku.ca/~amana/research/grid.pdf> */
        if (tmax.x < tmax.y) {
            if (tmax.x < tmax.z) {
                cell.x += step.x;
                if (cell.x < 0 || cell.x >= 8) break;
                axis = 0, t = tmax.x;
                tmax.x += delta.x;
            } else {
                cell.z += step.z;
                if (cell.z < 0 || cell.z >= 8) break;
                axis = 2, t = tmax.z;
                tmax.z += delta.z;
            }
        } else {
            if (tmax.y < tmax.z) {
                cell.y += step.y;
                if (cell.y < 0 || cell.y >= 8) break;
                axis = 1, t = tmax.y;
                tmax.y += delta.y;
            } else {
                cell.z += step.z;
                if (cell.z < 0 || cell.z >= 8) break;
                axis = 2, t = tmax.z;
                tmax.z += delta.z;
            }
        }
    }

    /* No hit occured! */
    return BIG_F32;
}

/**
 * @brief Set a voxel in the volume.
 *
 * @param pos Voxel position in the grid.
 * @param value Material index to assign.
 */
void OVoxelVolume::set_voxel(const int3& pos, const u8 value) {
    /* Sanity checks */
    assert(pos.x >= 0 && pos.x < grid_size.x && "Voxel out of range!");
    assert(pos.y >= 0 && pos.y < grid_size.y && "Voxel out of range!");
    assert(pos.z >= 0 && pos.z < grid_size.z && "Voxel out of range!");
    const bool removing = value == 0x00;

    // #if USE_BITPACKING
    { /* Raw voxel data */
        const u32 i = (pos.z * grid_size.y * grid_size.x) + (pos.y * grid_size.x) + pos.x;
        voxels[i] = value;
    }
    // #endif

    { /* Brick map */
        const int3 bpos = pos >> 3;
        const u32 w = brickmap.size.x, h = brickmap.size.y;
        const u32 i = (bpos.z * h * w) + (bpos.y * w) + bpos.x;
        Brick512* brick = brickmap.bricks + i;

        /* Do nothing because the voxel is already empty */
        if (removing && brick->voxels == nullptr) return;

        /* Allocate voxel bits if necessary */
        if (not removing && brick->voxels == nullptr) {
#if USE_BITPACKING
            brick->voxels = (u8*)MALLOC64(64);
            if (not brick->voxels) return;
            memset(brick->voxels, 0x00, 64);
#else
            brick->voxels = (u8*)MALLOC64(512);
            if (not brick->voxels) return;
            memset(brick->voxels, 0x00, 512);
#endif
        }

        /* Which bit do we need to update? */
        const int3 inner_pos = pos - (bpos << 3);

        /* Read the voxel bit, and update it */
        const u8 voxel = get_voxel(brick, inner_pos);
        if (removing) {
            if (voxel != 0x00) {
                /* Set the voxel to 0 (empty) */
                brick->voxcnt--;
                set_voxel(brick, inner_pos, 0);
            }
        } else {
            if (voxel == 0x00) {
                /* Set the voxel to 1 (solid) */
                brick->voxcnt++;
                set_voxel(brick, inner_pos, value);
            }
        }
    }
}

static f32 axis_delta(const f32 axis_dir) { return fabs(1.0f / axis_dir); }

static f32 entry(const f32 ro, const f32 rd, const f32 min, const f32 max) {
    const bool sign = (rd < 0);
    const f32 bmin = sign ? max : min;
    return (bmin - ro) / rd;  // TODO: maybe try getting rid of division here?
}

static f32 safe_entry(const f32 ro, const f32 rd, const f32 min, const f32 max) {
    const bool sign = (rd < 0);
    const f32 bmin = sign ? max : min;
    const f32 tmin = (bmin - ro) / rd;  // TODO: maybe try getting rid of division here?
    return fmaxf(0, tmin);
}

PacketHit8x8 OVoxelVolume::intersect(const RayPacket8x8& packet, const bool debug) const {
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
    const f32 k_min = 0, k_max = bb.size[k] - upv;  // TODO: make sure these are correct!

    /* Compute entry time of the corner rays along major axis K */
    const f32 entry_t = safe_entry(origin[k], ray_tl[k], k_min, k_max);

    /* Entry time along major axis K (floor + 0.5 is to align k_t to the grid) */
    f32 k_t = floorf((origin[k] + ray_tl[k] * entry_t) * vpu + 0.001f);
    const f32 k_o = (k_t - origin[k] * vpu) * k_sign;

    if (debug) {
        db::draw_line(origin, origin + packet.rays[0].dir * 5.5f, 0xFFFF0000);
        db::draw_line(origin, origin + packet.rays[63].dir * 5.5f, 0xFFFF0000);
        db::draw_line(origin, origin + packet.rays[7].dir * 5.5f, 0xFFFF0000);
        db::draw_line(origin, origin + packet.rays[56].dir * 5.5f, 0xFFFF0000);
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
        // if (u_max < 0.0f || v_max < 0.0f) continue;
        // if (u_min >= grid_size[u] || v_min >= grid_size[v]) continue;

        // for (u32 r = 0; r < 8 * 8; r++) {
        //     hit.hits[r].steps += 1;
        // }

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
        if (u_max < 0.0f || v_max < 0.0f) continue;
        if (u_min >= grid_size[u] || v_min >= grid_size[v]) continue;

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
                    // if (debug) {
                    //     float3 min_p = i, max_p = i + 1;

                    //    /* Draw floating point grid slice */
                    //    db::draw_aabb(min_p * upv, max_p * upv, 0xFFFF00FF);
                    //}

                    /* Find which rays intersect this voxel */
                    const float3 grid_o = origin * vpu;
                    for (u32 r = 0; r < 8 * 8; r++) {
                        // TODO: maybe don't do this transform every time?
                        const float3 rd = TransformVector(packet.rays[r].dir, bb.imodel);
                        const float3 ird = 1.0f / rd;

                        /* Entry point */
                        const f32 entry_k = k_t - fminf(k_sign, 0.0f);
                        const f32 entry_t = entry(grid_o[k], rd[k], entry_k, entry_k) + 0.001f;
                        float3 entry_p = grid_o + rd * entry_t;
                        int3 entry_c = floori(entry_p);

                        // if (debug) {
                        //     float3 min_p = floorf(entry_p), max_p = min_p + 1;

                        //    /* Draw floating point grid slice */
                        //    db::draw_aabb(min_p * upv, max_p * upv, 0xFFFF0000);
                        //}

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
                                hit.hits[r].albedo =
                                    RGB8_to_RGBF32(palette[voxels[i.z * grid_size.y * grid_size.x +
                                                                  i.y * grid_size.x + i.x]]);
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
