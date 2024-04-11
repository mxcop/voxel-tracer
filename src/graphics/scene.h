#pragma once

#include "primitives/basic/capsule.h"

class Scene {
    /* Bounding volume hierarchy */
    Bvh* bvh = nullptr;
    SkyDome skydome;

   public:
    /* Lighting */
    vector<SphereLight> lights;
    Capsule laser_segments[8];

    /* Models */
    OVoxelVolume* enemies[4];
    // TODO: remove polymorphism
    Traceable** shapes;
    u32 shapes_len = 0;

    /* Direction pointing towards the sun */
    const float3 sun_dir = {-0.619501f, 0.465931f, -0.631765f};
    float3 sun_light = {0.95f, 0.93f, 0.875f};

    Scene();
    ~Scene();

#ifdef PROFILING
    void set_shapes(Traceable** new_shapes, const u32 len);
#endif

    /**
     * @brief Called each frame.
     */
    void tick(const f32 dt);

    /**
     * @brief Intersect the scene with a ray.
     * @return Information about a potential hit.
     */
    HitInfo intersect(const Ray& ray) const;

    /**
     * @brief Intersect the scene with a ray.
     * @return Information about a potential hit.
     */
    PacketHit8x8 coherent_intersect(const RayPacket8x8& packet, bool debug = false) const;

    /**
     * @return True if the ray hit something before it reached tmax.
     */
    bool is_occluded(const Ray& ray, const f32 tmax = BIG_F32) const;

    /**
     * @brief Get the color of the sky for a ray.
     */
    float3 sample_sky(const Ray& ray) const;
};
