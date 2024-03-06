#include "camera.h"

bool Camera::update(const f32 t) {
    if (not WindowHasFocus()) return false;
    f32 speed = 3.0f * t;

    /* Determine the camera directions */
    float3 ahead = normalize(target - pos);
    float3 right = normalize(cross(UP, ahead));
    float3 up = normalize(cross(ahead, right));

    bool changed = false;
    /* Apply any user inputs */
    if (IsKeyDown(GLFW_KEY_UP)) target += speed * 0.33f * up, changed = true;
    if (IsKeyDown(GLFW_KEY_DOWN)) target -= speed * 0.33f * up, changed = true;
    if (IsKeyDown(GLFW_KEY_LEFT)) target -= speed * 0.33f * right, changed = true;
    if (IsKeyDown(GLFW_KEY_RIGHT)) target += speed * 0.33f * right, changed = true;
    ahead = normalize(target - pos);
    right = normalize(cross(UP, ahead));
    up = normalize(cross(ahead, right));
    speed = speed * (IsKeyDown(GLFW_KEY_LEFT_CONTROL) ? 4.0f : 1.0f);
    if (IsKeyDown(GLFW_KEY_A)) pos -= speed * right, changed = true;
    if (IsKeyDown(GLFW_KEY_D)) pos += speed * right, changed = true;
    if (GetAsyncKeyState('W')) pos += speed * ahead, changed = true;
    if (IsKeyDown(GLFW_KEY_S)) pos -= speed * ahead, changed = true;
    if (IsKeyDown(GLFW_KEY_SPACE)) pos += speed * up, changed = true;
    if (IsKeyDown(GLFW_KEY_LEFT_SHIFT)) pos -= speed * up, changed = true;
    speed = speed / (IsKeyDown(GLFW_KEY_LEFT_CONTROL) ? 4.0f : 1.0f);
    target = pos + ahead;

    /* Update the camera state */
    ahead = normalize(target - pos);
    up = normalize(cross(ahead, right));
    right = normalize(cross(up, ahead));
    tl = pos + 2.0f * ahead - ASPECT_RATIO * right + up;
    tr = pos + 2.0f * ahead + ASPECT_RATIO * right + up;
    bl = pos + 2.0f * ahead - ASPECT_RATIO * right - up;

    return changed;
}
