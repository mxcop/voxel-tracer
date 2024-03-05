#include "camera.h"

bool Camera::update(const f32 t) {
    if (not WindowHasFocus()) return false;
    f32 speed = 0.5f * t;

    /* Determine the camera directions */
    float3 ahead = normalize(target - pos);
    float3 right = normalize(cross(UP, ahead));
    float3 up = normalize(cross(ahead, right));

    bool changed = false;
    speed = speed * (IsKeyDown(GLFW_KEY_LEFT_CONTROL) ? 4.0f : 1.0f);
    /* Apply any user inputs */
    if (IsKeyDown(GLFW_KEY_A)) pos -= speed * 4 * right, changed = true;
    if (IsKeyDown(GLFW_KEY_D)) pos += speed * 4 * right, changed = true;
    if (IsKeyDown(GLFW_KEY_W)) pos += speed * 4 * ahead, changed = true;
    if (IsKeyDown(GLFW_KEY_S)) pos -= speed * 4 * ahead, changed = true;
    if (IsKeyDown(GLFW_KEY_SPACE)) pos += speed * 4 * up, changed = true;
    if (IsKeyDown(GLFW_KEY_LEFT_SHIFT)) pos -= speed * 4 * up, changed = true;
    speed = speed / (IsKeyDown(GLFW_KEY_LEFT_CONTROL) ? 4.0f : 1.0f);
    target = pos + ahead;
    if (IsKeyDown(GLFW_KEY_UP)) target += speed * up, changed = true;
    if (IsKeyDown(GLFW_KEY_DOWN)) target -= speed * up, changed = true;
    if (IsKeyDown(GLFW_KEY_LEFT)) target -= speed * right, changed = true;
    if (IsKeyDown(GLFW_KEY_RIGHT)) target += speed * right, changed = true;
    if (not changed) return false;

    /* Update the camera state */
    ahead = normalize(target - pos);
    up = normalize(cross(ahead, right));
    right = normalize(cross(up, ahead));
    tl = pos + 2.0f * ahead - ASPECT_RATIO * right + up;
    tr = pos + 2.0f * ahead + ASPECT_RATIO * right + up;
    bl = pos + 2.0f * ahead - ASPECT_RATIO * right - up;

    return true;
}
