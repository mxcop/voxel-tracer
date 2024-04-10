#include "player.h"

#define DRONE_CONTROLS 0

Player::Player(Camera& camera, OVoxelVolume* model) : camera(camera), model(model) {}

f32 Player::tick(const f32 dt, const int2 mouse_delta) { 
    Camera& c = camera;

    /* Input handling */
    if (IsKeyDown(GLFW_KEY_UP)) pitch += dt;
    if (IsKeyDown(GLFW_KEY_DOWN)) pitch -= dt;
    if (IsKeyDown(GLFW_KEY_LEFT)) yaw -= dt;
    if (IsKeyDown(GLFW_KEY_RIGHT)) yaw += dt;

    const f32 mouse_u = mouse_delta.x * 0.05f;
    const f32 mouse_v = mouse_delta.y * 0.05f;
    yaw += mouse_u * dt;
    pitch -= mouse_v * dt;
    if (pitch < -1.5f) pitch = -1.5f;
    if (pitch > 0.4f) pitch = 0.4f;

    quat new_rot = quat::from_axis_angle({0, 1, 0}, yaw);
    new_rot = new_rot * quat::from_axis_angle({1, 0, 0}, pitch);
    new_rot = new_rot * quat::from_axis_angle({0, 0, 1}, roll);

    const float3 up = new_rot.rotate_vec({0, 1, 0});
    const float3 ahead = new_rot.rotate_vec({0, 0, -1});
    const float3 side = new_rot.rotate_vec({1, 0, 0});

    constexpr f32 MOVE_SPEED = 20.0f;
    if (IsKeyDown(GLFW_KEY_W)) velocity += MOVE_SPEED * ahead * dt;
    if (IsKeyDown(GLFW_KEY_A)) velocity += MOVE_SPEED * side * dt;
    if (IsKeyDown(GLFW_KEY_S)) velocity -= MOVE_SPEED * ahead * dt;
    if (IsKeyDown(GLFW_KEY_D)) velocity -= MOVE_SPEED * side * dt;

    constexpr f32 VMOVE_SPEED = 35.0f;
    if (IsKeyDown(GLFW_KEY_SPACE)) velocity += VMOVE_SPEED * up * dt;
    if (IsKeyDown(GLFW_KEY_LEFT_SHIFT)) velocity -= VMOVE_SPEED * up * dt;

    /* Move camera */
    velocity *= powf(0.3f, dt);
    const float3 prev_pos = camera.pos;
    camera.pos += velocity * dt;

    const f32 depth_delta = dot(ahead, camera.pos) - dot(ahead, prev_pos);

    /* Rotate camera */
    camera.set_rot(new_rot);

    /* Update model transform */
    model->set_position(camera.pos + float3(-0.025f, 0.2f, 0));
    model->set_rotation(quat::from_axis_angle({0, 1, 0}, yaw));

	return depth_delta;
}
