#include "vv.h"

#define OGT_VOX_IMPLEMENTATION
#include <ogt/ogt_vox.h>

constexpr u32 MAX_STEPS = 256;

inline u32 merge_u8_u32(ogt_vox_rgba c) { return (c.r << 16) | (c.g << 8) | c.b; }

OVoxelVolume::OVoxelVolume(const float3& pos, const char* vox_path, const f32 vpu) : vpu(vpu) {
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
    const ogt_vox_model* model = scene->models[0];

    /* Initialize the volume */
    grid_size = int3(model->size_y, model->size_z, model->size_x);
    brickmap = Brickmap(grid_size);
    voxels = new u8[grid_size.x * grid_size.y * grid_size.z]{};
    bb = OBB(pos, float3(grid_size) / vpu);
    pivot = bb.size * 0.5f; /* center pivot */
    set_rotation(0, 0);

    /* Format the grid */
    for (u32 z = 0; z < model->size_z; z++) {
        for (u32 y = 0; y < model->size_y; y++) {
            for (u32 x = 0; x < model->size_x; x++) {
                const u32 mi = (z * model->size_y * model->size_x) + (y * model->size_x) + x;
                // const u32 i = (x * grid_size.y * grid_size.x) + (z * grid_size.x) + y;

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
      voxels(new u8[grid_size.x * grid_size.y * grid_size.z]{}),
      grid_size(grid_size),
      vpu(vpu) {
    /* Fill the grid with some noise */
    for (u32 z = 0; z < grid_size.z; z++) {
        for (u32 y = 0; y < grid_size.y; y++) {
            for (u32 x = 0; x < grid_size.x; x++) {
                const u32 i = (z * grid_size.y * grid_size.x) + (y * grid_size.x) + x;
                const f32 noise =
                    noise3D((f32)x / grid_size.x, (f32)y / grid_size.y, (f32)z / grid_size.z);
                if (noise > 0.09f) {
                    set_voxel(int3(x, y, z), 0xFF);
                } else {
                    set_voxel(int3(x, y, z), 0x00);
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

void OVoxelVolume::set_rotation(const float3& axis, const f32 angle) {
    bb.set_rotation_pivot(pivot, axis, angle);
}

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
                    hit.normal = TransformVector(hit.normal, bb.model);
                    return hit;
                }
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

        /* No hit occured! */
        hit.depth = BIG_F32;
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

    f32 t = 0.0f;
    for (hit.steps; hit.steps < MAX_STEPS; ++hit.steps) {
        /* Fetch the active cell */
        const u8 voxel = get_voxel(brick, cell);
        if (voxel) {
            const int3 hit_cell = (pos << 3) + cell;
            const u32 i =
                (hit_cell.z * grid_size.y * grid_size.x) + (hit_cell.y * grid_size.x) + hit_cell.x;
            hit.albedo = RGB8_to_RGBF32(palette[voxels[i]]);
            return t / vpu;
        }

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

    { /* Raw voxel data */
        const u32 i = (pos.z * grid_size.y * grid_size.x) + (pos.y * grid_size.x) + pos.x;
        voxels[i] = value;
    }

    { /* Brick map */
        const int3 bpos = pos >> 3;
        const u32 w = brickmap.size.x, h = brickmap.size.y;
        const u32 i = (bpos.z * h * w) + (bpos.y * w) + bpos.x;
        Brick512* brick = brickmap.bricks + i;

        /* Do nothing because the voxel is already empty */
        if (removing && brick->voxels == nullptr) return;

        /* Allocate voxel bits if necessary */
        if (not removing && brick->voxels == nullptr) {
            brick->voxels = (u8*)MALLOC64(64);
            if (not brick->voxels) return;
            memset(brick->voxels, 0x00, 64);
        }

        /* Which bit do we need to update? */
        const int3 inner_pos = pos - (bpos << 3);
        const u32 ii = (inner_pos.z * 8 * 8) + (inner_pos.y * 8) + inner_pos.x;
        const u32 byte = ii >> 3;
        const u8 bitmask = 0b1 << (ii - (byte << 3));

        /* Read the voxel bit, and update it */
        const u8 voxel = brick->voxels[byte] & bitmask;
        if (removing) {
            if (voxel != 0x00) {
                /* Set the voxel to 0 (empty) */
                brick->voxcnt--;
                brick->voxels[byte] ^= bitmask;
            }
        } else {
            if (voxel == 0x00) {
                /* Set the voxel to 1 (solid) */
                brick->voxcnt++;
                brick->voxels[byte] |= bitmask;
            }
        }
    }
}
