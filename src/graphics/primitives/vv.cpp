#include "precomp.h"
#include "vv.h"

// OVoxelVolume::OVoxelVolume(const float3& pos, const float3& scale, const int3& grid_size)
//     : bb(OBB(pos, scale)),
//       grid(new u8[grid_size.x * grid_size.y * grid_size.z]{}),
//       grid_size(grid_size) {
//     /* Fill the grid with some random voxels */
//     for (u32 z = 0; z < grid_size.z; z++) {
//         for (u32 y = 0; y < grid_size.y; y++) {
//             for (u32 x = 0; x < grid_size.x; x++) {
//                 u32 i = (z * grid_size.y * grid_size.x) + (y * grid_size.x) + x;
//                 if (RandomFloat() > 0.2f) {
//                     grid[i] = 0xFF;
//                 } else {
//                     grid[i] = 0x00;
//                 }
//             }
//         }
//     }
// }

AABB OVoxelVolume::get_aabb() const { return bb.get_aabb(); }
float3 OVoxelVolume::center() const { return bb.center(); }

void OVoxelVolume::set_rotation(const float3& axis, const f32 angle) {
    bb.set_rotation(axis, angle);
}

HitInfo OVoxelVolume::intersect(const Ray& ray) const {
    HitInfo hit = bb.intersect(ray);

    if (hit.depth != BIG_F32) {
        const Ray grid_ray = bb.world_to_local(ray);

        /* Voxels per unit */
        const float3 vpu = grid_size / bb.size;
        const float3 entry_pos = (grid_ray.origin + grid_ray.dir * hit.depth) * vpu;

        // TODO: do voxel tracing...
    }

    return hit;
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
        const u32 i = (bpos.z * brickmap.h * brickmap.w) + (bpos.y * brickmap.w) + bpos.x;
        Brickmap::Brick512* brick = brickmap.bricks + i;

        /* Do nothing because the voxel is already empty */
        if (removing && brick->voxels == nullptr) return;

        /* Allocate voxel bits if necessary */
        if (not removing && brick->voxels == nullptr) {
            brick->voxels = (u8*)MALLOC64(64);
            memset(brick->voxels, 0x00, 64);
        }

        const int3 inner_pos = pos - (bpos << 3);
        const u32 ii = (inner_pos.z * 8 * 8) + (inner_pos.y * 8) + inner_pos.x;
        const u32 byte = ii >> 3;
        const u8 bitmask = 0b1 << (ii - (byte << 3));

        // TODO: set the voxel and update voxcnt!
        const u8 voxel = brick->voxels[byte] & bitmask;
        if (removing) {
            if (voxel) {
            
            }
        }
        brick->voxels[byte] |= bitmask;
    }
}
