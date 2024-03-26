#include "scene.h"

Scene::Scene() {
    /* Material testing shape */
    shapes_len = 6;
    shapes = new Traceable* [shapes_len] {};

    const char* robot_m = "assets/vox/small-robot-walker.vox";
    /* Head & body */
    auto head = new OVoxelVolume({VOXEL * 4, VOXEL * 11, VOXEL * 4}, robot_m, 0);
    head->set_pivot(VOXEL * 4);
    head->set_rotation(quat::from_axis_angle({0, 1, 0}, 0));
    shapes[0] = head;
    shapes[1] = new OVoxelVolume({0, 0, 0}, robot_m, 1);

    /* Legs */
    shapes[2] = new OVoxelVolume({VOXEL * 5, VOXEL * -4, VOXEL * -5}, robot_m, 3);
    shapes[3] = new OVoxelVolume({VOXEL * 5, VOXEL * -4, VOXEL * 5}, robot_m, 4);
    shapes[4] = new OVoxelVolume({VOXEL * -5, VOXEL * -4, VOXEL * 5}, robot_m, 5);
    shapes[5] = new OVoxelVolume({VOXEL * -5, VOXEL * -4, VOXEL * -5}, robot_m, 6);

    // shapes[3] = new OVoxelVolume({2, 0, 0}, "assets/vox/small-robot-walker.vox", 2);

    // shapes[1] = new AABB({-50, -2, -50}, {50, -1, 50}, {1, 1, 1});

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

void Scene::tick(const f32 dt) {
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
