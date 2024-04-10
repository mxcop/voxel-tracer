#include "enemy.h"

Enemy::Enemy(OVoxelVolume* model) : model(model) { pos = RandomFloat3() * 5.0f; }

void Enemy::tick(const f32 dt, const float3 player) { 
	model->set_position(pos);
}
