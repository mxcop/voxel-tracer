#pragma once

class Enemy {
    OVoxelVolume* model;

    f32 pitch = 0, yaw = 0;

    float3 pos;
    float3 velocity = 0;

   public:
    Enemy() = delete;
    Enemy(OVoxelVolume* model);

    void tick(const f32 dt, const float3 player, Enemy** enemies, const u32 enemies_len);
};
