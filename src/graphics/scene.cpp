#include "scene.h"

constexpr u32 MODELS_NUM = 5;

Scene::Scene() {
    /* Material testing shape */
    shapes_len = MODELS_NUM + 8;
    shapes = new Traceable* [shapes_len] {};

    /* Render preview */
    auto cube = new OVoxelVolume(0, "assets/vox/testing/glass-box.vox", 0);
    shapes[0] = cube;

    /* Enemy models */
    for (u32 i = 0; i < 4; i++) {
        enemies[i] = new OVoxelVolume({(f32)i, 2, 0}, "assets/vox/enemy-drone.vox", 0);
        shapes[1 + i] = enemies[i];
    }

    /* Laser segment capsules */
    for (u32 i = 0; i < 8; i++) {
        laser_segments[i] = Capsule(1000, 1000, 0.025f);
        shapes[MODELS_NUM + i] = &laser_segments[i];
    }

    /* Initialize the BVH */
    bvh = new Bvh(shapes_len, shapes);

    /* Load the HDR skydome */
    skydome = SkyDome("assets/kiara_1_dawn_4k.hdr");
}

Scene::~Scene() {
    delete bvh;
    for (u32 i = 0; i < shapes_len; i++) {
        delete shapes[i];
    }
}

void Scene::tick(const f32) {
    /* Reconstruct the BVH in case something moved */
    bvh->build(shapes_len, shapes);
}

/**
 * @brief Intersect the scene with a ray.
 * @return Information about a potential hit.
 */
HitInfo Scene::intersect(const Ray& ray) const {
    HitInfo hit = bvh->intersect(ray);
    /* If the ray misses the scene, use the skydome color */
    if (hit.depth == BIG_F32) hit.albedo = skydome.sample_dir(ray.dir);
    return hit;
}

PacketHit8x8 Scene::coherent_intersect(const RayPacket8x8& packet, bool debug) const {
    PacketHit8x8 hits = bvh->coherent_intersect(packet, debug);
    for (u32 r = 0; r < 8 * 8; r++) {
        /* If the ray misses the scene, use the skydome color */
        if (hits.hits[r].depth == BIG_F32)
            hits.hits[r].albedo = skydome.sample_dir(packet.rays[r].dir);
    }
    return hits;
}

/**
 * @return True if the ray hit something before it reached tmax.
 */
bool Scene::is_occluded(const Ray& ray, const f32 tmax) const {
    Ray o_ray = ray;
    o_ray.shadow_ray = true;
    return bvh->is_occluded(o_ray, tmax);
}

/**
 * @brief Get the color of the sky for a ray.
 */
float3 Scene::sample_sky(const Ray& ray) const { return skydome.sample_dir(ray.dir); }

#ifdef PROFILING
void Scene::set_shapes(Traceable** new_shapes, const u32 len) {
    shapes_len = len;
    // TODO: delete old shapes...
    shapes = new Traceable* [shapes_len] {};
    for (u32 i = 0; i < len; i++) {
        shapes[i] = new_shapes[i];
    }
    bvh->build(shapes_len, shapes);
}
#endif
