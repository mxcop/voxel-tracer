#pragma once

class Enemy {
    OVoxelVolume* model;

    f32 pitch = 0, yaw = 0;

    float3 pos;
    float3 velocity = 0;

    i32 health = 32;

   public:
    Enemy() = delete;
    Enemy(OVoxelVolume* model);

    void respawn() { pos = RandomFloat3() * 32.0f - 16.0f; };

    bool tick(const f32 dt, const float3& player, Enemy** enemies, const u32 enemies_len);

    bool process_hit(const Ray& laser);
};
