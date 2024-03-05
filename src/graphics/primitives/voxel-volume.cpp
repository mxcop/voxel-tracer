#include "voxel-volume.h"

constexpr u32 MAX_STEPS = 256;

#if USE_BRICKMAP

HitInfo VoxelVolume::intersect(const Ray& ray) const {
    HitInfo hit = {};

    /* Early return if bounding box was not hit */
    f32 tmin, tmax;
    ray_vs_aabb(ray, tmin, tmax);
    if (tmin > tmax - 0.01f) {
        hit.depth = BIG_F32;
        return hit;
    }

    /* Calculate the ray start, end, and extend */
    const float3 p0 = ((ray.origin + ray.dir * (tmin + 0.001f)) - bbmin) * scale;
    const float3 p1 = ((ray.origin + ray.dir * (tmax - 0.001f)) - bbmin) * scale;
    const float3 extend = p1 - p0;
    const float3 inv_extend = 1.0f / extend;
    const float3 volume = (bbmax - bbmin) * scale;

    /* Setup for DDA traversal */
    f32 lod = MAX_LEVEL_SIZE, rlod = 1.0f / MAX_LEVEL_SIZE;
    u32 level = BRICK_LEVELS;
    u8* lod_ptr = voxels[level];
    float3 pos = clamp(floorf(p0 * rlod) * lod, float3(0), volume - 1.0f);
    float3 step = ray.sign_dir * lod;
    float3 delta = inv_extend * step;
    float3 side = (pos + fmaxf(step, 0.0f) - p0) * inv_extend;

    /* Start traversing */
    f32 t = 0.0f;
    u32 j = 0;
    for (hit.steps; hit.steps < MAX_STEPS; ++hit.steps) {
        const u8* voxel = fetch_voxel(floori(pos * rlod), lod_ptr);

        if (*voxel) {
            /* Stop if we hit something in L1 */
            if (level == 0) {
                hit.depth = tmin + (tmax - tmin) * t;
                hit.albedo = float3((f32)*voxel / 256.0f);
                /* If we hit on the first step, we have to compute the normal */
                if (hit.steps <= BRICK_LEVELS) {
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
                return hit;
            }
            level--, j = 0;

            /* Move down one level */
            lod *= BRICK_LEVEL_REC[level], rlod *= BRICK_LEVEL_MUL[level];
            step *= BRICK_LEVEL_REC[level], delta *= BRICK_LEVEL_REC[level];
            lod_ptr = voxels[level];

            pos =
                clamp(floorf((p0 + (t - 0.0001f) * extend) * rlod) * lod, float3(0), volume - 1.0f);
            side = (pos + fmaxf(step, float3(0.0f)) - p0) * inv_extend;
            continue;
        }

        /* (try) move up one level after 12 steps */
        //if (j++ >= 12 && level < BRICK_LEVELS) {
        //    level++, j = 0;

        //    /* Move up one level */
        //    lod *= BRICK_LEVEL_MUL[level - 1], rlod *= BRICK_LEVEL_REC[level - 1];
        //    step *= BRICK_LEVEL_MUL[level - 1], delta *= BRICK_LEVEL_MUL[level - 1];
        //    lod_ptr = voxels[level];

        //    pos =
        //        clamp(floorf((p0 + (t + 0.001f) * extend) * rlod) * lod, float3(0), volume - 1.0f);
        //    side = (pos + fmaxf(step, float3(0.0f)) - p0) * inv_extend;
        //    continue;
        //}

        /* Move to the next voxel */
        if (side.x < side.y) {
            if (side.x < side.z) {
                pos.x += step.x;
                if (pos.x < 0 || pos.x >= volume.x) break;
                hit.normal = float3(1, 0, 0);
                t = side.x;
                side.x += delta.x;
            } else {
                pos.z += step.z;
                if (pos.z < 0 || pos.z >= volume.z) break;
                hit.normal = float3(0, 0, 1);
                t = side.z;
                side.z += delta.z;
            }
        } else {
            if (side.y < side.z) {
                pos.y += step.y;
                if (pos.y < 0 || pos.y >= volume.y) break;
                hit.normal = float3(0, 1, 0);
                t = side.y;
                side.y += delta.y;
            } else {
                pos.z += step.z;
                if (pos.z < 0 || pos.z >= volume.z) break;
                hit.normal = float3(0, 0, 1);
                t = side.z;
                side.z += delta.z;
            }
        }
    }

    hit.depth = BIG_F32;
    hit.normal = float3(0);
    return hit;
}

#else

HitInfo VoxelVolume::intersect(const Ray& ray) const {
    HitInfo hit = {};

    /* Early return if bounding box was not hit */
    hit.depth = ray_vs_aabb(ray);
    if (hit.depth == BIG_F32) return hit;

    /* Calculate the size of the voxel volume in voxels */
    int3 size = floori(voxel_size);

    /* Calculate the starting voxel index */
    float3 ray_pos = ((ray.origin + ray.dir * (hit.depth + 0.001f)) - bbmin) * scale;
    int3 idx = floori(ray_pos);

    /* Clamp the starting voxel index inside the voxel volume */
    idx = clamp(idx, int3(0, 0, 0), size - 1);

    /* Which direction each axis will step in -1 or 1 */
    float3 step = ray.sign_dir;
    /* Indicates how far we must move (in units of t) to equal the width of a voxel */
    float3 delta = fabs(ray.r_dir);
    /* Determine t at which the ray crosses the first voxel boundary */
    float3 tmax = ((float3(idx) - ray_pos) + fmaxf(step, float3(0.0f))) * ray.r_dir;

    u32 i = 0;
    f32 t = 0.0f;
    float3 normal = {};
    for (; i < MAX_STEPS; ++i) {
        /* Fetch the current voxel */
        u8 voxel = fetch_voxel(idx);
        if (voxel != 0x00) {
            hit.depth += t / scale;
            hit.normal = -normal * step;
            hit.albedo = float3((f32)voxel / 256.0f);
            hit.steps = i;
            return hit;
        }

        /* Amanatides & Woo */
        /* <http://www.cse.yorku.ca/~amana/research/grid.pdf> */
        if (tmax.x < tmax.y) {
            if (tmax.x < tmax.z) {
                ray_pos.x += step.x;
                idx.x += step.x;
                if (idx.x < 0 || idx.x >= size.x) break;
                normal = float3(1, 0, 0);
                t = tmax.x;
                tmax.x += delta.x;
            } else {
                ray_pos.z += step.z;
                idx.z += step.z;
                if (idx.z < 0 || idx.z >= size.z) break;
                normal = float3(0, 0, 1);
                t = tmax.z;
                tmax.z += delta.z;
            }
        } else {
            if (tmax.y < tmax.z) {
                ray_pos.y += step.y;
                idx.y += step.y;
                if (idx.y < 0 || idx.y >= size.y) break;
                normal = float3(0, 1, 0);
                t = tmax.y;
                tmax.y += delta.y;
            } else {
                ray_pos.z += step.z;
                idx.z += step.z;
                if (idx.z < 0 || idx.z >= size.z) break;
                normal = float3(0, 0, 1);
                t = tmax.z;
                tmax.z += delta.z;
            }
        }
    }

    /* No hit occured! */
    hit.depth = BIG_F32;
    hit.steps = i;
    return hit;
}

#endif

#if USE_BRICKMAP

bool VoxelVolume::is_occluded(const Ray& ray, u32* steps) const {
    /* Early return if bounding box was not hit */
    f32 tmin, tmax;
    ray_vs_aabb(ray, tmin, tmax);
    // tmax = min(tmax, 1.0f);
    if (tmin > tmax - 0.01f || tmin > 1.0f) {
        return false;
    }
    tmax = min(tmax, 1.0f);

    /* Calculate the ray start, end, and extend */
    const float3 p0 = ((ray.origin + ray.dir * (tmin + 0.0001f)) - bbmin) * scale;
    const float3 p1 = ((ray.origin + ray.dir * (tmax - 0.0001f)) - bbmin) * scale;
    const float3 extend = p1 - p0;
    const float3 inv_extend = 1.0f / extend;
    const float3 volume = (bbmax - bbmin) * scale;

    /* Exit if the end-point is outside the volume */
    //if (p1.x < 0 || p1.y < 0 || p1.z < 0 || p1.x > volume.x - 1 || p1.y > volume.y - 1 ||
    //    p1.z > volume.z - 1) {
    //    return false;
    //}

    /* Setup for DDA traversal */
    f32 lod = MAX_LEVEL_SIZE, rlod = 1.0f / MAX_LEVEL_SIZE;
    u32 level = BRICK_LEVELS;
    u8* lod_ptr = voxels[level];
    float3 pos = clamp(floorf(p0 * rlod) * lod, float3(0), volume - 1.0f);
    float3 step = ray.sign_dir * lod;
    float3 delta = inv_extend * step;
    float3 side = (pos + fmaxf(step, float3(0.0f)) - p0) * inv_extend;

    /* Start traversing */
    float3 normal = {};
    f32 t = 0.0f;
#ifdef DEV
    u32 j = 0;
    u32& i = steps ? *steps : j;
#else
    u32 i = 0;
#endif
    for (i, j = 0; i < MAX_STEPS; ++i) {
        const u8* voxel = fetch_voxel(floori(pos * rlod), lod_ptr);

        if (*voxel) {
            /* Stop if we hit something in L1 */
            if (level == 0) {
                return true;
            }
            level--, j = 0;

            /* Move down one level */
            lod *= BRICK_LEVEL_REC[level], rlod *= BRICK_LEVEL_MUL[level];
            step *= BRICK_LEVEL_REC[level], delta *= BRICK_LEVEL_REC[level];
            lod_ptr = voxels[level];

            pos =
                clamp(floorf((p0 + (t - 0.0001f) * extend) * rlod) * lod, float3(0), volume - 1.0f);
            side = (pos + fmaxf(step, float3(0.0f)) - p0) * inv_extend;
            continue;
        }

        /* (try) move up one level after 12 steps */
        //if (j++ >= 12 && level < BRICK_LEVELS) {
        //    level++, j = 0;

        //    /* Move up one level */
        //    lod *= BRICK_LEVEL_MUL[level - 1], rlod *= BRICK_LEVEL_REC[level - 1];
        //    step *= BRICK_LEVEL_MUL[level - 1], delta *= BRICK_LEVEL_MUL[level - 1];
        //    lod_ptr = voxels[level];

        //    pos =
        //        clamp(floorf((p0 + (t + 0.001f) * extend) * rlod) * lod, float3(0), volume - 1.0f);
        //    side = (pos + fmaxf(step, float3(0.0f)) - p0) * inv_extend;
        //    continue;
        //}

        /* Move to the next voxel */
        if (side.x < side.y) {
            if (side.x < side.z) {
                pos.x += step.x;
                if (pos.x < 0 || pos.x >= volume.x) break;
                normal = float3(1, 0, 0);
                t = side.x;
                side.x += delta.x;
            } else {
                pos.z += step.z;
                if (pos.z < 0 || pos.z >= volume.z) break;
                normal = float3(0, 0, 1);
                t = side.z;
                side.z += delta.z;
            }
        } else {
            if (side.y < side.z) {
                pos.y += step.y;
                if (pos.y < 0 || pos.y >= volume.y) break;
                normal = float3(0, 1, 0);
                t = side.y;
                side.y += delta.y;
            } else {
                pos.z += step.z;
                if (pos.z < 0 || pos.z >= volume.z) break;
                normal = float3(0, 0, 1);
                t = side.z;
                side.z += delta.z;
            }
        }
        if (t > 1.0f) break;
    }

    return false;
}

#else

bool VoxelVolume::is_occluded(const Ray& ray) const {
    HitInfo hit = {};

    /* Early return if bounding box was not hit */
    hit.depth = ray_vs_aabb(ray);
    if (hit.depth == BIG_F32) return false;

    /* Calculate the size of the voxel volume in voxels */
    const int3 size = floori(voxel_size);

    /* Calculate the starting voxel index */
    const float3 ray_pos = ((ray.origin + ray.dir * (hit.depth + 0.001f)) - bbmin) * scale;
    int3 idx = floori(ray_pos);

    /* Clamp the starting voxel index inside the voxel volume */
    idx = clamp(idx, int3(0, 0, 0), size - 1);

    /* Which direction each axis will step in -1 or 1 */
    const float3 step = ray.sign_dir;
    /* Indicates how far we must move (in units of t) to equal the width of a voxel */
    const float3 delta = fabs(ray.r_dir);
    /* Determine t at which the ray crosses the first voxel boundary */
    float3 tmax = ((float3(idx) - ray_pos) + fmaxf(step, float3(0.0f))) * ray.r_dir;

    for (u32 i = 0; i < MAX_STEPS; ++i) {
        /* Fetch the current voxel */
        u8 voxel = fetch_voxel(idx);
        if (voxel != 0x00) {
            return true;
        }

        /* Amanatides & Woo */
        /* <http://www.cse.yorku.ca/~amana/research/grid.pdf> */
        if (tmax.x < tmax.y) {
            if (tmax.x < tmax.z) {
                idx.x += step.x;
                if (idx.x < 0 || idx.x >= size.x) break;
                tmax.x += delta.x;
                // if (tmax.x >= 1.0f) return false;
            } else {
                idx.z += step.z;
                if (idx.z < 0 || idx.z >= size.z) break;
                tmax.z += delta.z;
                // if (tmax.z >= 1.0f) return false;
            }
        } else {
            if (tmax.y < tmax.z) {
                idx.y += step.y;
                if (idx.y < 0 || idx.y >= size.y) break;
                tmax.y += delta.y;
                // if (tmax.y >= 1.0f) return false;
            } else {
                idx.z += step.z;
                if (idx.z < 0 || idx.z >= size.z) break;
                tmax.z += delta.z;
                // if (tmax.z >= 1.0f) return false;
            }
        }
    }
    return false;
}

#endif

static __forceinline i128 normal_encode(const i128 x, const i128 y, const i128 z, const u32 sx,
                                        const u32 sy) {
    const i128 msx = _mm_set1_epi32((i32)sx);
    const i128 msy = _mm_set1_epi32((i32)sy);
    const i128 fc = _mm_mul_epi32(_mm_mul_epi32(z, msx), msy);
    const i128 sc = _mm_add_epi32(_mm_mul_epi32(y, msx), x);
    return _mm_add_epi32(fc, sc);
    // return (z * sx * sy) + (y * sx) + x;
}

/* Amanatides & Woo */
/* <http://www.cse.yorku.ca/~amana/research/grid.pdf> */
PacketHitInfo VoxelVolume::intersect(const RayPacket& packet) const {
    PacketHitInfo hit = {};

    /* Early return if no rays hit the bounding box */
    ray4_vs_aabb(packet, hit);
    const f128 miss_mask = _mm_cmpeq_ps(hit.depth, BIG_F128);
    if (_mm_movemask_ps(miss_mask) == 0b1111) return hit;

    /* Voxels per unit of space, and it's reciprocal */
    const f128 ms = _mm_set_ps1(scale), is = _mm_set_ps1(1.0f / scale);

    /* Determine the ray entry point into the volume */
    const f128 entry_t = _mm_add_ps(hit.depth, _mm_set_ps1(0.0001f));
    const f128 entry_x = _mm_fmadd_ps(packet.rd_x, entry_t, packet.ro_x);
    const f128 entry_y = _mm_fmadd_ps(packet.rd_y, entry_t, packet.ro_y);
    const f128 entry_z = _mm_fmadd_ps(packet.rd_z, entry_t, packet.ro_z);

    /* Ray position (floating) */
    const f128 rp_x = _mm_mul_ps(_mm_sub_ps(entry_x, _mm_set_ps1(bbmin.x)), ms);
    const f128 rp_y = _mm_mul_ps(_mm_sub_ps(entry_y, _mm_set_ps1(bbmin.y)), ms);
    const f128 rp_z = _mm_mul_ps(_mm_sub_ps(entry_z, _mm_set_ps1(bbmin.z)), ms);

    /* Ray step, for moving the ray in voxel space (-1 or 1) */
    const i128 step_x = _mm_set1_epi32(packet.signs.x);
    const i128 step_y = _mm_set1_epi32(packet.signs.y);
    const i128 step_z = _mm_set1_epi32(packet.signs.z);

    /* Distance (in units of t) equal to the size of a voxel */
    const f128 delta_x = _mm_abs_ps(packet.ird_x);
    const f128 delta_y = _mm_abs_ps(packet.ird_y);
    const f128 delta_z = _mm_abs_ps(packet.ird_z);

    /* Ray index, the 3D index of the current voxel */
    const i128 lwb = _mm_set1_epi32(0);
    i128q ri;
    ri.x = _mm_floorclamp_ps(rp_x, lwb, _mm_set1_epi32(voxel_size.x - 1));
    ri.y = _mm_floorclamp_ps(rp_y, lwb, _mm_set1_epi32(voxel_size.y - 1));
    ri.z = _mm_floorclamp_ps(rp_z, lwb, _mm_set1_epi32(voxel_size.z - 1));
    // ri.x = _mm_cvtps_epi32(_mm_floor_ps(rp_x));
    // ri.y = _mm_cvtps_epi32(_mm_floor_ps(rp_y));
    // ri.z = _mm_cvtps_epi32(_mm_floor_ps(rp_z));

    /* Next position of the ray (in voxel space) */
    const f128 next_x = _mm_sub_ps(_mm_cvtepi32_ps(ri.x), rp_x);
    const f128 next_y = _mm_sub_ps(_mm_cvtepi32_ps(ri.y), rp_y);
    const f128 next_z = _mm_sub_ps(_mm_cvtepi32_ps(ri.z), rp_z);

    /* Get the starting offset for tmax, if this isn't done a visual glitch will occur */
    const f128 offset_x = _mm_max_ps(_mm_cvtepi32_ps(step_x), _mm_setzero_ps());
    const f128 offset_y = _mm_max_ps(_mm_cvtepi32_ps(step_y), _mm_setzero_ps());
    const f128 offset_z = _mm_max_ps(_mm_cvtepi32_ps(step_z), _mm_setzero_ps());

    /* "t" at which the ray crosses the next voxel boundary */
    f128 tmax_x = _mm_mul_ps(_mm_add_ps(next_x, offset_x), packet.ird_x);
    f128 tmax_y = _mm_mul_ps(_mm_add_ps(next_y, offset_y), packet.ird_y);
    f128 tmax_z = _mm_mul_ps(_mm_add_ps(next_z, offset_z), packet.ird_z);

    /* Mask that tracks status of the rays (done = 0xFFFFFFFF) */
    i128 status_mask = (i128&)miss_mask;
    i128 exit_mask = (i128&)miss_mask;
    constexpr u32 RAYS = 4;

    /* Traverse */
    u32 s = 0;
    f128 t = _mm_setzero_ps();
    for (; s < MAX_STEPS; ++s) {
        /* Get the voxel indices */
#if USE_MORTON
        const i128 indices = morton_encode(ri.x, ri.y, ri.z);
#else
        const i128 indices = normal_encode(ri.x, ri.y, ri.z, voxel_size.x, voxel_size.y);
#endif

        /* Mark any rays outside the volume as done */
        exit_mask = _mm_or_epi32(exit_mask, _mm_cmplt_epi32(ri.x, _mm_set1_epi32(0)));
        exit_mask = _mm_or_epi32(exit_mask, _mm_cmplt_epi32(ri.y, _mm_set1_epi32(0)));
        exit_mask = _mm_or_epi32(exit_mask, _mm_cmplt_epi32(ri.z, _mm_set1_epi32(0)));
        exit_mask =
            _mm_or_epi32(exit_mask, _mm_cmpgt_epi32(ri.x, _mm_set1_epi32(voxel_size.x - 1)));
        exit_mask =
            _mm_or_epi32(exit_mask, _mm_cmpgt_epi32(ri.y, _mm_set1_epi32(voxel_size.y - 1)));
        exit_mask =
            _mm_or_epi32(exit_mask, _mm_cmpgt_epi32(ri.z, _mm_set1_epi32(voxel_size.z - 1)));
        status_mask = _mm_or_epi32(status_mask, exit_mask);

        /* Early exit if all rays are done */
        if (_mm_movemask_ps((f128&)status_mask) == 0b1111) {
            hit.depth = _mm_blendv_ps(_mm_fmadd_ps(t, is, hit.depth), BIG_F128, (f128&)exit_mask);
            break;
        }

        /* Mark rays that hit a voxel as done */
        const i128 safe_indices = _mm_andnot_epi32(status_mask, indices);
        const i128 vox = fetch_voxels(safe_indices);
        status_mask = _mm_or_epi32(status_mask, _mm_cmpgt_epi32(vox, _mm_setzero_si128()));

        /* Find the smallest axis for each ray */
        const f128 min_yz = _mm_min_ps(tmax_y, tmax_z);
        const f128 min_zx = _mm_min_ps(tmax_z, tmax_x);
        const f128 min_xy = _mm_min_ps(tmax_x, tmax_y);
        /* And make a mask */
        const f128 cmp_x = _mm_cmple_ps(tmax_x, min_yz);
        const f128 cmp_y = _mm_cmple_ps(tmax_y, min_zx);
        const f128 cmp_z = _mm_cmple_ps(tmax_z, min_xy);

        /* Record the current "t" of each ray */
        const f128 tmin = _mm_min_ps(_mm_min_ps(tmax_x, tmax_y), tmax_z);
        t = _mm_blendv_ps(tmin, t, (f128&)status_mask);

        /* Step to the next voxel index */
        ri.x = _mm_add_epi32(ri.x, _mm_and_epi32((i128&)cmp_x, step_x));
        ri.y = _mm_add_epi32(ri.y, _mm_and_epi32((i128&)cmp_y, step_y));
        ri.z = _mm_add_epi32(ri.z, _mm_and_epi32((i128&)cmp_z, step_z));

        /* Update the ray tmax */
        tmax_x = _mm_add_ps(tmax_x, _mm_and_ps(cmp_x, delta_x));
        tmax_y = _mm_add_ps(tmax_y, _mm_and_ps(cmp_y, delta_y));
        tmax_z = _mm_add_ps(tmax_z, _mm_and_ps(cmp_z, delta_z));
    }

    /* DEBUG: Save the step count */
    hit.steps = s;
    return hit;
}
