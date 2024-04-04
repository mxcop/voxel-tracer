#pragma once

/**
 * @brief A grid of voxels with rotation.
 */
class CoherentVoxelVolume {
    /* Bounding box */
    OBB bb;
    float3 pivot;

    /* Voxel material indices, 1 byte each. */
    u8* voxels = nullptr;
    int3 grid_size = 0;

    /* Voxel material palette, 8 bit rgba. */
    u32* palette = nullptr;

    /* Voxels per unit. */
    f32 vpu;

    __forceinline i32 getsign(const f32 f) const { return (i32)(((u32&)f) >> 31) * 2 - 1; }

    void set_voxel(const u32 x, const u32 y, const u32 z, const u8 mat);

   public:
    CoherentVoxelVolume() = default;
    /* Create voxel volume of certain size and fill it with noise */
    CoherentVoxelVolume(const float3& pos, const int3& grid_size, const f32 vpu = 20.0f);
    ~CoherentVoxelVolume() {
        delete[] voxels;
        delete[] palette;
    }

    /* Copy */
    CoherentVoxelVolume(const CoherentVoxelVolume&) = delete;
    CoherentVoxelVolume& operator=(const CoherentVoxelVolume&) = delete;
    /* Move */
    CoherentVoxelVolume(CoherentVoxelVolume&&) = default;
    CoherentVoxelVolume& operator=(CoherentVoxelVolume&&) = default;

    /* Trace-able functions */
    PacketHit8x8 intersect(const RayPacket8x8& packet, const bool debug = false) const;

    inline void set_pivot(const float3 pivot) { bb.pivot = pivot, this->pivot = pivot; };
    inline void set_rotation(const quat& rot) { bb.set_rotation_pivot(pivot, rot); };
    inline void set_position(const float3& pos) { bb.set_position(pos); };
};
