#pragma once

class Player {
    Camera& camera;
    OVoxelVolume* model;

    f32 pitch = 0, yaw = 0;

    float3 velocity = 0;

   public:
    Player() = delete;
    Player(Camera& camera, OVoxelVolume* model);

    /* Returns the depth delta! */
    f32 tick(const f32 dt, const int2 mouse_delta);
};
