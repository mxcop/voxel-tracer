#pragma once

class Enemy {
    OVoxelVolume* model;

    float3 pos;

   public:
    Enemy() = delete;
    Enemy(OVoxelVolume* model);

    void tick(const f32 dt, const float3 player);
};
