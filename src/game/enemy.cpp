#include "enemy.h"

constexpr f32 ENEMY_SPEED = 10.0f;

constexpr f32 PLAYER_WEIGHT = 2.0f;
constexpr f32 ENEMY_WEIGHT = 2.0f;

Enemy::Enemy(OVoxelVolume* model) : model(model) { pos = RandomFloat3() * 32.0f - 16.0f; }

bool Enemy::tick(const f32 dt, const float3& player, Enemy** enemies, const u32 enemies_len) {
    /* Target the player */
    float3 target_dir = normalize(player - pos) * PLAYER_WEIGHT;

    /* Avoid other enemies */
    for (u32 i = 0; i < enemies_len; i++) {
        const float3 extend = pos - enemies[i]->pos;
        const f32 dist = length(extend);
        if (dist == 0) continue; /* Ignore self */

        const f32 factor = fmaxf((5.0f - dist) / 5.0f, 0.0f) * ENEMY_WEIGHT;
        target_dir += factor * (extend / dist);
    }
    target_dir = normalize(target_dir); /* Final target */

    /* Accelerate towards the player */
    velocity += target_dir * dt * ENEMY_SPEED;
    velocity *= powf(0.3f, dt);
    pos += velocity * dt;

    /* Look towards the player */
    const float3 look_dir = normalize(velocity);
    yaw = atan2(look_dir.x, look_dir.z);
    quat new_rot = quat::from_axis_angle({0, 1, 0}, yaw);

    /* Update the model */
    model->set_position(pos);
    model->set_rotation(new_rot);

    if (length(player - pos) < 1.0f) {
        return true;
    }
    return false;
}

bool Enemy::process_hit(const Ray& laser) {
    const HitInfo hit = model->intersect(laser);
    if (hit.no_hit()) return false;

    /* Laser hit */
    const float3 hit_point = laser.intersection(hit) - hit.normal * 0.001f;
    const int3 voxel_point = make_int3(model->to_grid(hit_point));
    model->set_voxel(voxel_point, 0x00);
    health--;

    /* Died */
    if (health <= 0) {
        pos = RandomFloat3() * 32.0f - 16.0f;
        velocity = 0;
        health = 32;

        model->reload_model("assets/vox/enemy-drone.vox");
        return true;
    }
    return false;
}
