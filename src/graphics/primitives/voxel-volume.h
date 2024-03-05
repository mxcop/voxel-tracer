#pragma once

#include "graphics/rays/hit.h"
#include "morton.h"

/* 16 voxels per unit of space */
constexpr f32 DEFAULT_VOXEL_SCALE = 16.0f;

/* Use morton ordering for voxels (slightly faster) */
#define USE_MORTON 1
/* Use the AVX2 gather instructions to load voxels (worse on my AMD CPU) */
#define USE_AVX2_GATHER 1

#define USE_BRICKMAP 1
constexpr u32 BRICK_LEVELS = 4;
constexpr f32 BRICK_LEVEL_MUL[BRICK_LEVELS] = {2.0f, 2.0f, 2.0f, 2.0f};
constexpr f32 BRICK_LEVEL_REC[BRICK_LEVELS] = {0.5f, 0.5f, 0.5f, 0.5f};
constexpr u32 MAX_LEVEL_SIZE = 16;  // 1u << BRICK_LEVELS;
// constexpr u32 BRICK_SIZE = 8;
// constexpr f32 R_BRICK_SIZE = 1.0f / BRICK_SIZE;

struct VoxelVolume {
    f32 scale = DEFAULT_VOXEL_SCALE;
    /* Bounding box */
    union {
        f128 bmin;
        struct {
            float3 bbmin;
            f32 _;
        };
    };
    union {
        f128 bmax;
        struct {
            float3 bbmax;
            f32 _;
        };
    };
    union {
        f128 vsize;
        struct {
            float3 voxel_size;
            f32 _;
        };
    };
    // TODO: Change from vector[] to u8**!
    u8* voxels[BRICK_LEVELS + 1];

    VoxelVolume() = default;
    ~VoxelVolume() {
        for (u32 i = 0; i < BRICK_LEVELS + 1; i++) {
            delete[] voxels[i];
        }
        // delete[] voxels;
    }

    VoxelVolume(const VoxelVolume&) = delete;
    VoxelVolume(VoxelVolume&&) = default;
    VoxelVolume& operator=(const VoxelVolume&) = delete;
    VoxelVolume& operator=(VoxelVolume&&) = default;

    /**
     * @param pos Position of the voxel volume in world space.
     * @param size Size of the voxel volume in voxels.
     * @param scale The world space scale of 1 voxel. (default 16 voxels per unit)
     */
    VoxelVolume(float3 pos, int3 size, f32 scale = DEFAULT_VOXEL_SCALE)
        : scale(scale), bbmin(pos), bbmax(pos + float3(size) / scale), voxel_size(size) {
#if USE_MORTON
        /* Morton only works for perfect cubes */
        u32 bytes = pow(max(max(size.x, size.y), size.z), 3);
#else
        u32 bytes = size.x * size.y * size.z;
#endif
        voxels[0] = new u8[bytes];
        for (u32 i = 0, scale = 1; i < BRICK_LEVELS; i++) {
            scale *= BRICK_LEVEL_MUL[i];
            const u32 size = (scale * scale * scale);
            voxels[i + 1] = new u8[bytes / size];
        }

        /* Generate the highest level of detail */
        f32 rx = 1.0f / 128.0f;
#pragma omp parallel for schedule(dynamic)
        for (i32 z = 0; z < size.z; ++z) {
            const f32 fz = (f32)z / 128.0f;
            for (u32 y = 0; y < size.y; ++y) {
                const f32 fy = (f32)y / 128.0f;
                f32 fx = 0;
                for (u32 x = 0; x < size.x; ++x, fx += rx) {
                    const f32 noise = noise3D(fx, fy, fz);
                    const f32 noise2 = noise3D(fx * 8.0f, fy * 8.0f, fz * 8.0f) + .5f;
#if USE_MORTON
                    const u64 i = morton_encode(x, y, z);
#else
                    const u64 i = (z * size.x * size.y) + (y * size.x) + x;
#endif
                    voxels[0][i] = noise > 0.09f ? (noise2 * 0xFF) : 0x00;
                }
            }
        }

        int3 bsize = size;
        for (u32 b = 0; b < BRICK_LEVELS; b++) {
            const f32 bmul = BRICK_LEVEL_MUL[b];
            /* See if this brick is solid or not */
            for (u32 z = 0; z < bsize.z; z += bmul) {
                for (u32 y = 0; y < bsize.y; y += bmul) {
                    for (u32 x = 0; x < bsize.x; x += bmul) {
                        voxels[b + 1][morton_encode(x * BRICK_LEVEL_REC[b], y * BRICK_LEVEL_REC[b],
                                                    z * BRICK_LEVEL_REC[b])] = 0;
#if USE_MORTON
                        for (u32 bz = 0; bz < bmul; bz++) {
                            for (u32 by = 0; by < bmul; by++) {
                                for (u32 bx = 0; bx < bmul; bx++) {
                                    voxels[b + 1][morton_encode(x * BRICK_LEVEL_REC[b],
                                                                y * BRICK_LEVEL_REC[b],
                                                                z * BRICK_LEVEL_REC[b])] |=
                                        voxels[b][morton_encode(x + bx, y + by, z + bz)];
                                }
                            }
                        }
#else
                        // TODO: this stuff here...
                        const u64 i = (z * size.x * size.y) + (y * size.x) + x;
#endif
                    }
                }
            }

            bsize = floori(float3(bsize) * BRICK_LEVEL_REC[b]); /* 2, 4, 8, 16 */
        }

#if 0
        u32 brick_indices[BRICK_LEVELS] = {};
        for (u32 i = 0; i < BRICK_LEVELS + 1; i++) {
            const u32 scale = 1u << i;
            const u32 size = (scale * scale * scale);
            voxels[i] = std::vector<u8>();
            voxels[i].resize(bytes / size);
            
            if (i > 0) {
                brick_indices[i - 1] = size;
            }
        }

        /* Populate the voxel data */
        u8 bricks[BRICK_LEVELS] = {};
        for (size_t i = 0; i < bytes; ++i) {
            u32 x, y, z;
            morton_decode(i, x, y, z);
            const f32 fx = (f32)x / 128.0f;
            const f32 fy = (f32)y / 128.0f;
            const f32 fz = (f32)z / 128.0f;
            const f32 noise = noise3D(fx, fy, fz);

            /* Set the voxel */
            if (noise > 0.09f) {
                voxels[0][i] = 0xFF;
                for (u32 j = 0; j < BRICK_LEVELS; j++) {
                    bricks[j] = 0xFF;
                }
            }

            /* Update the brick levels */
            for (u32 j = 0; j < BRICK_LEVELS; j++) {
                if (i == brick_indices[j]) {
                    const u32 scale = 1u << (j + 1);
                    const u32 size = (scale * scale * scale);
                    brick_indices[j] += size;

                    if (bricks[j] == 0xFF) voxels[j + 1][i / size] = 0xFF;
                    bricks[j] = 0x00; /* reset */
                }
            }
        }
#endif

        // #if USE_BRICKMAP
        //        voxels[1] = std::vector<u8>();
        //        voxels[1].resize(bytes / (BRICK_SIZE * BRICK_SIZE * BRICK_SIZE));
        //
        //        float3 bsize = size / BRICK_SIZE;
        //
        //        f32 rx = 1.0f / 128.0f;
        //        for (u32 z = 0; z < size.z; z += BRICK_SIZE) {
        //            for (u32 y = 0; y < size.y; y += BRICK_SIZE) {
        //                for (u32 x = 0; x < size.x; x += BRICK_SIZE) {
        //                    u8 brick = 0x00;
        //                    for (u32 bz = 0; bz < BRICK_SIZE; bz++) {
        //                        for (u32 by = 0; by < BRICK_SIZE; by++) {
        //                            for (u32 bx = 0; bx < BRICK_SIZE; bx++) {
        //                                const f32 fx = (f32)(x + bx) / 128.0f;
        //                                const f32 fy = (f32)(y + by) / 128.0f;
        //                                const f32 fz = (f32)(z + bz) / 128.0f;
        //                                const f32 noise = noise3D(fx, fy, fz);
        //
        // #if USE_MORTON
        //                                const u64 i = morton_encode((x + bx), (y + by), (z + bz));
        // #else
        //                                const u64 i =
        //                                    ((z + bz) * size.x * size.y) + ((y + by) * size.x) +
        //                                    (x + bx);
        // #endif
        //                                if (noise > 0.09f) {
        //                                    brick = 0xFF; /* brick contains at least 1 voxel */
        //
        //                                    const f32 color =
        //                                        noise3D(fx * 8.0f, fy * 8.0f, fz * 8.0f) + .5f;
        //                                    voxels[0][i] = color * 0xFF;
        //                                } else {
        //                                    voxels[0][i] = 0x00;
        //                                }
        //                            }
        //                        }
        //                    }
        // #if USE_MORTON
        //                    const u64 i = morton_encode(x / BRICK_SIZE, y / BRICK_SIZE, z /
        //                    BRICK_SIZE);
        // #else
        //                    const u64 i = ((z / BRICK_SIZE) * bsize.x * bsize.y) +
        //                                  ((y / BRICK_SIZE) * bsize.x) + (x / BRICK_SIZE);
        // #endif
        //                    voxels[1][i] = brick;
        //                }
        //            }
        //        }
        // #else
        //        f32 rx = 1.0f / 128.0f;
        //        for (u32 z = 0; z < size.z; ++z) {
        //            const f32 fz = (f32)z / 128.0f;
        //            for (u32 y = 0; y < size.y; ++y) {
        //                const f32 fy = (f32)y / 128.0f;
        //                f32 fx = 0;
        //                for (u32 x = 0; x < size.x; ++x, fx += rx) {
        //                    const f32 noise = noise3D(fx, fy, fz);
        //                    const f32 noise2 = noise3D(fx * 8.0f, fy * 8.0f, fz * 8.0f) + .5f;
        // #if USE_MORTON
        //                    const u64 i = morton_encode(x, y, z);
        // #else
        //                    const u64 i = (z * size.x * size.y) + (y * size.x) + x;
        // #endif
        //                    voxels[i] = noise > 0.09f ? (noise2 * 0xFF) : 0x00;
        //                }
        //            }
        //        }
        // #endif
    }

    /**
     * @brief Intersect the voxel volume with a ray.
     */
    HitInfo intersect(const Ray& ray) const;
    bool is_occluded(const Ray& ray, u32* steps = nullptr) const;

    /**
     * @brief Intersect the voxel volume with a packet of 4 rays.
     */
    PacketHitInfo intersect(const RayPacket& packet) const;

   private:
    /* Fast ray to aabb check, using SSE and FMA. */
    inline f32 ray_vs_aabb(const Ray& ray) const {
        /* Idea to use fmsub to save 1 instruction came from
         * <http://www.joshbarczak.com/blog/?p=787> */
        const f128 ord = _mm_mul_ps(ray.o, ray.rd);
        const f128 t1 = _mm_fmsub_ps(bmin, ray.rd, ord);
        const f128 t2 = _mm_fmsub_ps(bmax, ray.rd, ord);

        /* Find the near and far intersection point */
        const f128 vmax4 = _mm_max_ps(t1, t2), vmin4 = _mm_min_ps(t1, t2);
        const f32 tmax = _min(_min(vmax4.m128_f32[0], vmax4.m128_f32[1]), vmax4.m128_f32[2]);
        /* The last max(0) here handles being inside of the AABB */
        const f32 tmin =
            _max(_max(_max(vmin4.m128_f32[0], vmin4.m128_f32[1]), vmin4.m128_f32[2]), 0.0f);
        return tmin <= tmax ? tmin : BIG_F32;
    }

    inline void ray_vs_aabb(const Ray& ray, f32& tmin_out, f32& tmax_out) const {
        /* Idea to use fmsub to save 1 instruction came from
         * <http://www.joshbarczak.com/blog/?p=787> */
        const f128 ord = _mm_mul_ps(ray.o, ray.rd);
        const f128 t1 = _mm_fmsub_ps(bmin, ray.rd, ord);
        const f128 t2 = _mm_fmsub_ps(bmax, ray.rd, ord);

        /* Find the near and far intersection point */
        const f128 vmax4 = _mm_max_ps(t1, t2), vmin4 = _mm_min_ps(t1, t2);
        const f32 tmax = _min(_min(vmax4.m128_f32[0], vmax4.m128_f32[1]), vmax4.m128_f32[2]);
        /* The last max(0) here handles being inside of the AABB */
        const f32 tmin =
            _max(_max(_max(vmin4.m128_f32[0], vmin4.m128_f32[1]), vmin4.m128_f32[2]), 0.0f);
        tmin_out = tmin, tmax_out = tmax;
    }

    /* Fast ray packet to aabb check, using SSE */
    void ray4_vs_aabb(const RayPacket& packet, PacketHitInfo& out) const {
        f128 tmin = _mm_setzero_ps();
        f128 tmax = _mm_set_ps1(BIG_F32);

        for (u32 a = 0; a < 3; ++a) {
#if 0
            const f128 aabb_min = _mm_set_ps1(packet.signs[a] > 0.0f ? bbmin[a] : bbmax[a]);
            const f128 aabb_max = _mm_set_ps1(packet.signs[a] > 0.0f ? bbmax[a] : bbmin[a]);

            const f128 dmin = _mm_mul_ps(_mm_sub_ps(aabb_min, packet.ro[a]), packet.ird[a]);
            const f128 dmax = _mm_mul_ps(_mm_sub_ps(aabb_max, packet.ro[a]), packet.ird[a]);

            tmin = _mm_max_ps(tmin, dmin);
            tmax = _mm_min_ps(tmax, dmax);
#else

            // TODO: optimize function check more...
            const f128 t1 =
                _mm_mul_ps(_mm_sub_ps(_mm_set_ps1(bbmin[a]), packet.ro[a]), packet.ird[a]);
            const f128 t2 =
                _mm_mul_ps(_mm_sub_ps(_mm_set_ps1(bbmax[a]), packet.ro[a]), packet.ird[a]);

            tmin = _mm_max_ps(tmin, _mm_min_ps(t1, t2));
            tmax = _mm_min_ps(tmax, _mm_max_ps(t1, t2));
// tmin = _mm_min_ps(_mm_max_ps(t1, tmin), _mm_max_ps(t2, tmin));
// tmax = _mm_max_ps(_mm_min_ps(t1, tmax), _mm_min_ps(t2, tmax));
#endif

            // float t1 = (box->min[d] - ray->origin[d]) * ray->dir_inv[d];
            // float t2 = (box->max[d] - ray->origin[d]) * ray->dir_inv[d];

            // tmin = max(tmin, min(t1, t2));
            // tmax = min(tmax, max(t1, t2));
        }
        // tmin = _mm_max_ps(_mm_setzero_ps(), tmin);

        /* Use a mask to remove non-intersections (tmin > tmax) */
        const f128 mask = _mm_cmple_ps(tmin, tmax);
        out.depth = _mm_blendv_ps(BIG_F128, tmin, mask);
        out.exit_t = tmax;
    }

    /* Fetch a voxel from this volumes voxel data. */
    __forceinline u8 fetch_voxel(const int3& idx, u32 lod = 0) const {
#if USE_MORTON
        const u64 i = morton_encode(idx.x, idx.y, idx.z);
#else
        const float3 size = voxel_size;
        size_t i = ((size_t)idx.z * size.x * size.y) + ((size_t)idx.y * size.x) + idx.x;
#endif
        return voxels[lod][i];
    }

    __forceinline const u8* fetch_voxel(const int3& idx, const u8* lod) const {
#if USE_MORTON
        const u32 i = morton_encode(idx.x, idx.y, idx.z);
#else
        const float3 size = voxel_size;
        size_t i = ((size_t)idx.z * size.x * size.y) + ((size_t)idx.y * size.x) + idx.x;
#endif
        return lod + i;
    }

    /* Fetch a brick from this volumes brick data. */
    //    __forceinline u8 fetch_brick(const int3& idx) const {
    // #if USE_MORTON
    //        const u64 i = morton_encode(idx.x, idx.y, idx.z);
    // #else
    //        const float3 size = voxel_size / BRICK_SIZE;
    //        size_t i = ((size_t)idx.z * size.x * size.y) + ((size_t)idx.y * size.x) + idx.x;
    // #endif
    //        return bricks[0][i];
    //    }

    /* Fetch 4 voxels using the indices inside an SSE register. */
    __forceinline i128 fetch_voxels(const i128 indices) const {
#if USE_AVX2_GATHER
        return gather_voxels(indices);
#else
        return _mm_set_epi32(voxels[indices.m128i_u32[3]], voxels[indices.m128i_u32[2]],
                             voxels[indices.m128i_u32[1]], voxels[indices.m128i_u32[0]]);
#endif
    }

    /* Gather 4 voxels from this volumes voxel data using SIMD. */
    __forceinline i128 gather_voxels(const i128 indices) const {
        const i128 u8mask = _mm_set1_epi32(0x000000FF);
        return _mm_and_epi32(_mm_i32gather_epi32((int*)voxels[0], indices, sizeof(u8)),
                             u8mask);
    }
};
