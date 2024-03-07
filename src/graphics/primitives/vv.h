#pragma once

/**
 * @brief A grid of voxels with rotation.
 */
class OVoxelVolume : public Traceable {
    /* Bounding box */
    OBB bb;

    struct Brickmap {
        /* 8x8x8 voxels, 512 bits, 64 bytes, 1 bit each. (solid | empty) */
        struct Brick512 {
            u8* voxels; /* Can be null! */
            u16 voxcnt = 0;
        };

        Brick512* bricks;
        /* Brick map size (x:w, y:h, z:d) */
        u16 w, h, d;
    } brickmap;

    /* Voxel material indices, 1 byte each. */
    u8* voxels;
    int3 grid_size;

   public:
    OVoxelVolume() = default;
    OVoxelVolume(const float3& pos, const float3& scale, const int3& grid_size);
    ~OVoxelVolume() {
        delete[] voxels;
        delete[] brickmap.bricks;
    }

    /* Copy */
    OVoxelVolume(const OVoxelVolume&) = delete;
    OVoxelVolume& operator=(const OVoxelVolume&) = delete;
    /* Move */
    OVoxelVolume(OVoxelVolume&&) = default;
    OVoxelVolume& operator=(OVoxelVolume&&) = default;

    /* Trace-able functions */
    AABB get_aabb() const override;
    float3 center() const override;
    HitInfo intersect(const Ray& ray) const override;

    void set_rotation(const float3& axis, const f32 angle);

    /**
     * @brief Set a voxel in the volume.
     * 
     * @param pos Voxel position in the grid.
     * @param value Material index to assign.
     */
    void set_voxel(const int3& pos, const u8 value);
};
