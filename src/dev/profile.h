#pragma once

#ifdef PROFILING

typedef OVoxelVolume vv;

/**
 * @brief Initialize the profiling scene.
 */
void profile_init(Renderer& r) {
    /* Crates */
    //r.shapes[0] = new vv({-3.0f, 2.5f + VOXEL * 6, -0.5f}, "assets/vox/crate-10.vox");
    //r.shapes[1] = new vv({-3.0f, 2.5f - VOXEL * 3, -0.5f}, "assets/vox/crate-16h.vox");
    //r.shapes[2] = new vv({1.0f + VOXEL * 17, 2.5f, -0.5f}, "assets/vox/crate-16.vox");

    /* Try load the camera settings */
    FILE* f = fopen("camera_profiling.bin", "rb");
    if (f) {
        fread(&r.camera, 1, sizeof(Camera), f);
        fclose(f);
    }

    /* Crates */
    Traceable** shapes = new Traceable*[8 * 8 * 8]{};
    u32 i = 0;
    const char* models[2] = {"assets/vox/crate-16.vox", "assets/vox/crate-10.vox"};
    for (u32 z = 0; z < 8; z++) {
        const char* model = models[z % 2];
        for (u32 y = 0; y < 8; y++) {
            for (u32 x = 0; x < 8; x++, i++) {
                shapes[i] = new vv({VOXEL * 32.0f * x, VOXEL * 32.0f * y, VOXEL * 32.0f * z}, model);
            }
        }
    }

    r.bvh = new Bvh(8 * 8 * 8, shapes);
}

/**
 * @brief Update the profiling scene.
 */
void profile_tick(Renderer& r, const f32 dt) {}

#endif
