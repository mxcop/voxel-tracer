#include "precomp.h"
#include "colliders.h"

float3 SphereCollider::furthest_point(const Transform& t, const float3& dir) const {
    return (t.position + center) + normalize(dir) * radius;
}

float3 BoxCollider::furthest_point(const Transform& t, const float3& dir) const {
    float3 maxp;
    f32 maxd = -BIG_F32;

    for (u32 i = 0; i < 8; i++) {
        const f32 x = (i & 1) ? max.x : min.x;
        const f32 y = (i & 2) ? max.y : min.y;
        const f32 z = (i & 4) ? max.z : min.z;
        const float3 corner = t.position + float3(x, y, z);

        const f32 d = dot(corner, dir);
        if (d > maxd) {
            maxd = d;
            maxp = corner;
        }
    }

    return maxp;
}

VoxelCollider::VoxelCollider(const int3& grid_size, const f32 vpu, const u8* voxel_data)
    : Collider(VOXELS), size(grid_size), extend(float3(grid_size) / vpu) {
    /* 
     * 0 = Empty or surrounded. (never generates a contact point, and never needs to be checked)
     * 1 = Face voxel.          (never generates a contact point)
     * 2 = Edge voxel.          (can generate contact point, with "edge" or "corner")
     * 3 = Corner voxel.        (can generate contact point, with "face", "edge", or "corner")
     */

    // TODO: convert "voxel_data" to a collision friendly format.

}
