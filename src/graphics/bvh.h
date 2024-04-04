#pragma once

#include "rays/ray.h"
#include "primitives/basic/aabb.h"
#include "primitives/basic/obb.h"

/*
 * The BVH implementation is inspired by a series of articles from Jacco Bikker.
 * Source: https://jacco.ompf2.com/2022/04/13/how-to-build-a-bvh-part-1-basics/
 */

class Bvh {
   public:
    /* NOTE: a cache line is usually 64 bytes */
    /* This node is 32 bytes, so 2 perfectly fit into a cache line */
    struct alignas(32) Node {
        union {
            f128 aabb_min4;
            struct {
                float3 aabb_min;
                u32 left_first;
            };
        };
        union {
            f128 aabb_max4;
            struct {
                float3 aabb_max;
                u32 prim_count;
            };
        };

        bool is_leaf() const { return prim_count; }
    };

   private:
    Node* nodes = nullptr;
    /* Skip the second node, for better child node cache alignment */
    u16 root_idx = 0, nodes_used = 2;
    u16 size = 2;

    // TODO: use a vector of indices instead.
    // So we don't need to modify the original vector.
    Traceable** prims = nullptr;

    void subdivide(Bvh::Node& node, int lvl);

   public:
    Bvh() = default;
    Bvh(u32 size, Traceable** new_prims);
    ~Bvh() {
        if (nodes) delete[] nodes;
    };

    /* Copy */
    Bvh(const Bvh&) = delete;
    Bvh& operator=(const Bvh&) = delete;
    /* Move */
    Bvh(Bvh&&) = default;
    Bvh& operator=(Bvh&&) = default;

    void build(const u32 size, Traceable** prims);

    /**
     * @brief Evaluate the surface area heuristic of a node along an axis with a certain split
     * position.
     */
    f32 evaluate_sah(const Node& node, i32 axis, f32 pos) const;
    f32 find_best_split_plane(const Node& node, i32& axis, f32& pos) const;

    /**
     * @brief Intersect the BVH and return information about a potential hit.
     */
    HitInfo intersect(const Ray& ray) const;
    PacketHit8x8 coherent_intersect(const RayPacket8x8& packet, bool debug = false) const;

    PacketHitInfo intersect(const RayPacket128& packet) const;

    /**
     * @brief Check if a ray hits anything in its path.
     */
    bool is_occluded(const Ray& ray, const f32 tmax = BIG_F32, u32* steps = nullptr) const;
};
