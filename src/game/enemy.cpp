#include "enemy.h"

constexpr f32 ENEMY_SPEED = 10.0f;

constexpr f32 PLAYER_WEIGHT = 2.0f;
constexpr f32 ENEMY_WEIGHT = 2.0f;

Enemy::Enemy(OVoxelVolume* model) : model(model) { pos = RandomFloat3() * 32.0f - 16.0f; }

void Enemy::tick(const f32 dt, const float3 player, Enemy** enemies, const u32 enemies_len) { 
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
    pitch = asin(-look_dir.y);
    yaw = atan2(look_dir.x, look_dir.z);
	quat new_rot = quat::from_axis_angle({0, 1, 0}, yaw);
    new_rot = new_rot * quat::from_axis_angle({1, 0, 0}, pitch);

    /* Update the model */
    model->set_position(pos);
	model->set_rotation(new_rot);
}
