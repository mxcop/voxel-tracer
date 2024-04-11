#pragma once

class Player {
    Camera& camera;

    f32 pitch = 0, yaw = 0;

    float3 velocity = 0;

   public:
    Player() = delete;
    Player(Camera& camera);

    /* Returns the depth delta! */
    f32 tick(const f32 dt, const int2 mouse_delta);
};
