#pragma once

/**
 * @brief A grid of voxels with rotation.
 */
class OVoxelVolume : public Traceable {
    /* Bounding box */
    OBB bb;

    /* 8x8x8 voxels, 512 bits, 64 bytes, 1 bit each. (solid | empty) */
    struct Brick512 {
        u8* voxels = nullptr; /* Can be null! */
        u16 voxcnt = 0;
    };

    struct Brickmap {
        Brick512* bricks = nullptr;
        /* Brick map size in bricks. */
        int3 size = 0;

        Brickmap() = default;
        explicit Brickmap(const int3& grid_size) {
            size = grid_size >> 3;
            bricks = new Brick512[size.z * size.y * size.x]{};
        }
    } brickmap;

    /* Voxel material indices, 1 byte each. */
    u8* voxels = nullptr;
    int3 grid_size = 0;

    /* Voxel material palette, 8 bit rgba. */
    u32* palette = nullptr;

    /* Voxels per unit. */
    f32 vpu;

    /**
     * @brief Get a brick from the brickmap by position.
     */
    __forceinline Brick512* get_brick(const int3& brick) const {
        const int3& bs = brickmap.size;
        const u32 i = (brick.z * bs.x * bs.y) + (brick.y * bs.x) + brick.x;
        return brickmap.bricks + i;
    };

    /**
     * @brief Get a voxel from a brick by position.
     */
    __forceinline u8 get_voxel(const Brick512* brick, const int3& voxel) const {
        const u32 i = (voxel.z * 8 * 8) + (voxel.y * 8) + voxel.x;
        const u32 byte = i >> 3;
        const u8 bitmask = 0b1 << (i - (byte << 3));
        return brick->voxels[byte] & bitmask;
    };

    /**
     * @brief Traverse a ray through a brick from the brickmap.
     */
    f32 traverse_brick(const Brick512* brick, const int3& pos, const Ray& ray, const f32 entry_t,
                       const f32 rbpu, u32& axis, HitInfo& hit) const;

   public:
    OVoxelVolume() = default;
    /* Create voxel volume from a .vox file */
    OVoxelVolume(const float3& pos, const char* vox_path, const f32 vpu = 16.0f);
    /* Create voxel volume of certain size and fill it with noise */
    OVoxelVolume(const float3& pos, const int3& grid_size, const f32 vpu = 16.0f);
    ~OVoxelVolume() {
        delete[] voxels;
        delete[] palette;
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
