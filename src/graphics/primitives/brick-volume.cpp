#include "brick-volume.h"

constexpr u32 MAX_STEPS = 256;

BrickVolume::BrickVolume(const float3 pos, const int3 size, const f32 vpu)
    : bbmin(pos),
      bbmax(pos + float3(size) / vpu),
      vsize(size),
      bsize(size >> 3),
      vpu(vpu),
      bpu(vpu * 0.125f) {
    /* Volume size MUST BE a multiple of 8 */
    assert(size.x % 8 == 0 && size.y % 8 == 0 && size.z % 8 == 0);

    /* Setup the voxels and bricks */
    // voxels.resize((u64)vsize.x * vsize.y * vsize.z);
    const u64 bricks = (u64)bsize.x * bsize.y * bsize.z;
    brickmap = new Brick512[bricks];

#if 0
    /* TEMP: Random bricks */
    u32 seed = 18942663u;
    u32 seed2 = 74928367u;
    for (u32 b = 0; b < bricks; ++b) {
        Brick512* brick = brickmap + b;

        if (RandomFloat(seed) < 0.2f) {
            /* Allocate voxel packets */
            brick->packets = new u8[64];
            memset(brick->packets, 0x00, 64);

            /* Randomly assign voxels */
            for (u16 i = 0; i < 64; i++) {
                brick->packets[i] = 0u;
                for (u8 b = 0; b < 8; b++) {
                    if (RandomFloat(seed2) < 0.3f) {
                        brick->packets[i] |= (0b1 << b);
                        brick->popcnt++;
                    }
                }
            }
        }
    }
#else
    const float3 offset = make_float3(RandomFloat(), RandomFloat(), RandomFloat());
#pragma omp parallel for schedule(dynamic)
    for (i32 bz = 0; bz < bsize.z; bz++) {
        for (u32 by = 0; by < bsize.y; by++) {
            for (u32 bx = 0; bx < bsize.x; bx++) {
                Brick512* brick = &brickmap[(bz * bsize.y * bsize.x) + (by * bsize.x) + bx];
                for (u32 z = 0; z < 8; z++) {
                    for (u32 y = 0; y < 8; y++) {
                        for (u32 x = 0; x < 8; x++) {
                            const f32 noise = noise3D(offset.x + (bx * 8 + x) * 0.01f,
                                                      offset.y + (by * 8 + y) * 0.01f,
                                                      offset.z + (bz * 8 + z) * 0.01f);

                            if (noise > 0.03f) {
                                if (not brick->packets) {
                                    brick->packets = (u8*)MALLOC64(64);
                                    brick->voxels = new u32[8 * 8 * 8];
                                    memset(brick->packets, 0x00, 64);
                                    memset(brick->voxels, 0x00, sizeof(u32) * 8 * 8 * 8);
                                }

                                /* Noise value above */
                                const f32 noise_up = noise3D(offset.x + (bx * 8 + x) * 0.01f,
                                                             offset.y + (by * 8 + y + 1) * 0.01f,
                                                             offset.z + (bz * 8 + z) * 0.01f);

                                const u32 i = (z * 8 * 8) + (y * 8) + x;
                                // const u32 i = morton_encode(x, y, z);
                                if (noise_up < 0.03f) {
                                    /* Grass */
                                    const f32 noise_c =
                                        pow((noise3D(offset.x + (bx * 8 + x) * 0.2f,
                                                     offset.y + (by * 8 + y) * 0.2f,
                                                     offset.z + (bz * 8 + z) * 0.2f) +
                                             0.3f),
                                            2) *
                                        2.0f;
                                    const float4 color = float4(0.55f, 1.0f, 0.1f, 0.0f) * noise_c;
                                    brick->voxels[i] = RGBF32_to_RGB8(&color);
                                } else {
                                    /* Dirt */
                                    const f32 noise_c =
                                        pow((noise3D(offset.x + (bx * 8 + x) * 0.2f,
                                                     offset.y + (by * 8 + y) * 0.2f,
                                                     offset.z + (bz * 8 + z) * 0.2f) +
                                             0.3f),
                                            2) *
                                        2.0f;
                                    const float4 color =
                                        float4(1.0f, 0.7f, 0.65f, 0.0f) * (noise_c * 0.95f + 0.05f);
                                    brick->voxels[i] = RGBF32_to_RGB8(&color);
                                }

                                const u32 byte = i >> 3;
                                const u8 bitmask = 0b1 << (i - (byte << 3));
                                brick->packets[byte] |= bitmask;
                                brick->popcnt++;
                            }
                        }
                    }
                }
            }
        }
    }
#endif
}

BrickVolume::~BrickVolume() { delete[] brickmap; }

/**
 * @brief Intersect the volume bounding box. (tmin > tmax, means no intersection)
 */
void BrickVolume::intersect_bb(const Ray& ray, f32& tmin_out, f32& tmax_out) const {
#if 1
    /* Idea to use fmsub to save 1 instruction came from
     * <http://www.joshbarczak.com/blog/?p=787> */
    const f128 ord = _mm_mul_ps(ray.o, ray.rd);
    const f128 t1 = _mm_fmsub_ps(bmin, ray.rd, ord);
    const f128 t2 = _mm_fmsub_ps(bmax, ray.rd, ord);

    /* Find the near and far intersection point */
    const f128 vmax4 = _mm_max_ps(t1, t2), vmin4 = _mm_min_ps(t1, t2);

    /* Set the 4th element to 0 */
    const f128 tmin4 = (f128&)_mm_slli_si128((i128&)vmin4, 4);

    /* Get the horizontal minimum and maximum "t" */
    tmax_out = _mm_hmin3_ps(vmax4), tmin_out = _mm_hmax_ps(tmin4);
#else
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
#endif
}

/**
 * @brief Intersect the volume with a ray.
 */
HitInfo BrickVolume::intersect(const Ray& ray) const {
    HitInfo hit;

    /* Exit if the ray misses the volume */
    f32 tmin, tmax;
    intersect_bb(ray, tmin, tmax);
    if (tmin > tmax - 0.01f) {
        hit.depth = BIG_F32;
        return hit;
    }
    hit.depth = tmin;

    /* Volume entry position */
    const float3 entry = ((ray.origin + ray.dir * (hit.depth)) - bbmin) * bpu;
    /* Clamp the entry point inside the volume grid */
    int3 cell = clamp(floori(entry), int3(0), bsize - 1);

    /* Which direction each axis will step in -1 or 1 */
    const int3 step = make_int3(ray.sign_dir);
    /* Indicates how far we must move (in units of t) to equal the width of a voxel */
    const float3 delta = fabs(ray.r_dir);
    /* Determine t at which the ray crosses the first voxel boundary */
    float3 side = ((float3(cell) - entry) + fmaxf(step, float3(0))) * ray.r_dir;

    f32 t = 0.0f;
    const f32 rbpu = 1.0f / bpu;
    for (hit.steps = 0; hit.steps < MAX_STEPS; ++hit.steps) {
        /* Fetch the active cell */
        const Brick512* brick = get_brick(cell);
        if (brick->popcnt > 0) {
            const f32 brick_entry_t = hit.depth + t * rbpu;
            const f32 dist = traverse_brick(brick, cell, ray, brick_entry_t, hit);
            if (dist != BIG_F32) {
                if (hit.steps == 0) {
                    const float3 p = (ray.origin + ray.dir * tmin);
                    if (fabs(p.x - bbmin.x) < 0.01f) {
                        hit.normal = float3(-1, 0, 0);
                    } else if (fabs(p.x - bbmax.x) < 0.01f) {
                        hit.normal = float3(1, 0, 0);
                    } else if (fabs(p.y - bbmin.y) < 0.01f) {
                        hit.normal = float3(0, -1, 0);
                    } else if (fabs(p.y - bbmax.y) < 0.01f) {
                        hit.normal = float3(0, 1, 0);
                    } else if (fabs(p.z - bbmin.z) < 0.01f) {
                        hit.normal = float3(0, 0, -1);
                    } else if (fabs(p.z - bbmax.z) < 0.01f) {
                        hit.normal = float3(0, 0, 1);
                    }
                } else {
                    hit.normal = -hit.normal * step;
                }
                hit.depth = brick_entry_t + dist;

                return hit;
            }
        }

        /* Amanatides & Woo */
        /* <http://www.cse.yorku.ca/~amana/research/grid.pdf> */
        if (side.x < side.y) {
            if (side.x < side.z) {
                cell.x += step.x;
                if (cell.x < 0 || cell.x >= bsize.x) break;
                hit.normal = float3(1, 0, 0);
                t = side.x;
                side.x += delta.x;
            } else {
                cell.z += step.z;
                if (cell.z < 0 || cell.z >= bsize.z) break;
                hit.normal = float3(0, 0, 1);
                t = side.z;
                side.z += delta.z;
            }
        } else {
            if (side.y < side.z) {
                cell.y += step.y;
                if (cell.y < 0 || cell.y >= bsize.y) break;
                hit.normal = float3(0, 1, 0);
                t = side.y;
                side.y += delta.y;
            } else {
                cell.z += step.z;
                if (cell.z < 0 || cell.z >= bsize.z) break;
                hit.normal = float3(0, 0, 1);
                t = side.z;
                side.z += delta.z;
            }
        }
    }

    /* No hit occured! */
    hit.depth = BIG_F32;
    return hit;
}

bool BrickVolume::is_occluded(const Ray& ray, u32* steps) const {
    /* Exit if the ray misses the volume */
    f32 tmin, tmax;
    intersect_bb(ray, tmin, tmax);
    tmax = min(tmax, 1.0f);
    if (tmin > tmax - 0.001f) {
        return false;
    }

    /* Volume entry position */
    const float3 entry = ((ray.origin + ray.dir * tmin) - bbmin) * bpu;
    /* Clamp the entry point inside the volume grid */
    int3 cell = clamp(floori(entry), 0, bsize - 1);

    /* Which direction each axis will step in -1 or 1 */
    const float3 step = ray.sign_dir;
    /* Indicates how far we must move (in units of t) to equal the width of a voxel */
    const float3 delta = fabs(ray.r_dir);
    /* Determine t at which the ray crosses the first voxel boundary */
    float3 side = ((float3(cell) - entry) + fmaxf(step, float3(0))) * ray.r_dir;

    f32 t = 0.0f;
    const f32 rbpu = 1.0f / bpu;
#ifdef DEV
    u32 j = 0;
    u32& i = steps ? *steps : j;
#else
    u32 i = 0;
#endif
    for (; i < MAX_STEPS; ++i) {
        /* Fetch the active cell */
        const Brick512* brick = get_brick(cell);
        if (brick->popcnt > 0) {
            const f32 brick_entry_t = tmin + t * rbpu;
            if (traverse_brick(brick, cell, ray, brick_entry_t, tmax)) {
                return true;
            }
        }

        /* Amanatides & Woo */
        /* <http://www.cse.yorku.ca/~amana/research/grid.pdf> */
        if (side.x < side.y) {
            if (side.x < side.z) {
                cell.x += step.x;
                if (cell.x < 0 || cell.x >= bsize.x) break;
                t = side.x;
                side.x += delta.x;
            } else {
                cell.z += step.z;
                if (cell.z < 0 || cell.z >= bsize.z) break;
                t = side.z;
                side.z += delta.z;
            }
        } else {
            if (side.y < side.z) {
                cell.y += step.y;
                if (cell.y < 0 || cell.y >= bsize.y) break;
                t = side.y;
                side.y += delta.y;
            } else {
                cell.z += step.z;
                if (cell.z < 0 || cell.z >= bsize.z) break;
                t = side.z;
                side.z += delta.z;
            }
        }

        if (tmin + t * rbpu >= tmax) break;
    }

    /* No hit occured! */
    return false;
}

bool BrickVolume::traverse_brick(const Brick512* brick, const int3& pos, const Ray& ray,
                                 const f32 entry_t, const f32 tmax) const {
    /* Brick minimum position */
    const float3 bmin = bbmin + float3(pos) / bpu;
    /* Brick entry position */
    const float3 entry = ((ray.origin + ray.dir * entry_t) - bmin) * vpu;
    /* Clamp the entry point inside the volume grid */
    int3 cell = clamp(floori(entry), 0, 8 - 1);

    /* Which direction each axis will step in -1 or 1 */
    const float3 step = ray.sign_dir;
    /* Indicates how far we must move (in units of t) to equal the width of a voxel */
    const float3 delta = fabs(ray.r_dir);
    /* Determine t at which the ray crosses the first voxel boundary */
    float3 side = ((float3(cell) - entry) + fmaxf(step, float3(0))) * ray.r_dir;

    f32 t = 0.0f;
    const f32 rvpu = 1.0f / vpu;
    for (u32 i = 0; i < MAX_STEPS; ++i) {
        /* Fetch the active cell */
        const u8 voxel = get_voxel(brick, floori(cell));
        if (voxel) {
            return true;
        }

        /* Amanatides & Woo */
        /* <http://www.cse.yorku.ca/~amana/research/grid.pdf> */
        if (side.x < side.y) {
            if (side.x < side.z) {
                cell.x += step.x;
                if (cell.x < 0 || cell.x >= 8) break;
                t = side.x;
                side.x += delta.x;
            } else {
                cell.z += step.z;
                if (cell.z < 0 || cell.z >= 8) break;
                t = side.z;
                side.z += delta.z;
            }
        } else {
            if (side.y < side.z) {
                cell.y += step.y;
                if (cell.y < 0 || cell.y >= 8) break;
                t = side.y;
                side.y += delta.y;
            } else {
                cell.z += step.z;
                if (cell.z < 0 || cell.z >= 8) break;
                t = side.z;
                side.z += delta.z;
            }
        }

        if (entry_t + t * rvpu >= tmax) break;
    }

    /* No hit occured! */
    return false;
}

f32 BrickVolume::traverse_brick(const Brick512* brick, const int3& pos, const Ray& ray,
                                const f32 entry_t, HitInfo& hit) const {
    /* Brick minimum position */
    const float3 bmin = bbmin + float3(pos) / bpu;
    /* Brick entry position */
    const float3 entry = ((ray.origin + ray.dir * (entry_t)) - bmin) * vpu;
    /* Clamp the entry point inside the volume grid */
    int3 cell = clamp(floori(entry), 0, 8 - 1);

    /* Which direction each axis will step in -1 or 1 */
    const float3 step = ray.sign_dir;
    /* Indicates how far we must move (in units of t) to equal the width of a voxel */
    const float3 delta = fabs(ray.r_dir);
    /* Determine t at which the ray crosses the first voxel boundary */
    float3 side = ((float3(cell) - entry) + fmaxf(step, float3(0))) * ray.r_dir;

    f32 t = 0.0f;
    for (hit.steps; hit.steps < MAX_STEPS; ++hit.steps) {
        /* Fetch the active cell */
        const u8 voxel = get_voxel(brick, cell);
        if (voxel) {
            hit.albedo = RGB8_to_RGBF32(brick->voxels[(cell.z * 8 * 8) + (cell.y * 8) + cell.x]);
            return t / vpu;
        }

        /* Amanatides & Woo */
        /* <http://www.cse.yorku.ca/~amana/research/grid.pdf> */
        if (side.x < side.y) {
            if (side.x < side.z) {
                cell.x += step.x;
                if (cell.x < 0 || cell.x >= 8) break;
                hit.normal = float3(1, 0, 0);
                t = side.x;
                side.x += delta.x;
            } else {
                cell.z += step.z;
                if (cell.z < 0 || cell.z >= 8) break;
                hit.normal = float3(0, 0, 1);
                t = side.z;
                side.z += delta.z;
            }
        } else {
            if (side.y < side.z) {
                cell.y += step.y;
                if (cell.y < 0 || cell.y >= 8) break;
                hit.normal = float3(0, 1, 0);
                t = side.y;
                side.y += delta.y;
            } else {
                cell.z += step.z;
                if (cell.z < 0 || cell.z >= 8) break;
                hit.normal = float3(0, 0, 1);
                t = side.z;
                side.z += delta.z;
            }
        }
    }

    /* No hit occured! */
    return BIG_F32;
}

void BrickVolume::place_voxel(const Ray& ray) {
    const HitInfo hit = intersect(ray);
    if (hit.depth == BIG_F32) return;

    /* Get the cell at the intersection point */
    const int3 hit_cell =
        floori((ray.origin + ray.dir * hit.depth + hit.normal * 0.0001f - bbmin) * vpu);
    const int3 hit_brick = hit_cell >> 3;

    /* Exit if the brick is outside the volume */
    if (hit_brick.x < 0 || hit_brick.y < 0 || hit_brick.z < 0) return;
    if (hit_brick.x >= bsize.x || hit_brick.y >= bsize.y || hit_brick.z >= bsize.z) return;

    Brick512* brick = get_brick(hit_brick);
    if (brick->packets == nullptr) {
        brick->packets = (u8*)MALLOC64(64);
        brick->voxels = new u32[8 * 8 * 8];
        memset(brick->packets, 0x00, 64);
        memset(brick->voxels, 0x00, sizeof(u32) * 8 * 8 * 8);
    }

    /* Set the voxel to 1 */
    const int3 i = hit_cell - (hit_brick << 3);
    const u32 idx = (i.z * 8 * 8) + (i.y * 8) + i.x;
    const u32 byte = idx >> 3;
    const u8 bitmask = 0b1 << (idx - (byte << 3));
    brick->packets[byte] |= bitmask;
    brick->popcnt++;

    /* Noise color */
    const f32 noise_c = (noise3D((hit_cell.x) * 0.2f, (hit_cell.y) * 0.2f, (hit_cell.z) * 0.2f) + 0.3f) * 2.0f;

    /* Insert voxel data */
    const float4 color = float4(0.25f, 0.5f, 1.0f, 0.1f) * noise_c * noise_c;
    brick->voxels[idx] = RGBF32_to_RGB8(&color);
}
