#include "camera.h"

f32 Camera::update(const f32 t) {
    if (not WindowHasFocus()) return false;
    f32 speed = 1.5f * t;

    prev_pyramid = Pyramid(pos, get_dir(0, 0), get_dir(1, 0), get_dir(0, 1));

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
    bool moved_forward = false;
    bool forward = true;
    if (IsKeyDown(GLFW_KEY_A)) pos -= speed * right, changed = true;
    if (IsKeyDown(GLFW_KEY_D)) pos += speed * right, changed = true;
    if (GetAsyncKeyState('W'))
        pos += speed * ahead, changed = true, moved_forward = true;
    else if (IsKeyDown(GLFW_KEY_S))
        pos -= speed * ahead, changed = true, moved_forward = true, forward = false;
    if (IsKeyDown(GLFW_KEY_SPACE)) pos += speed * up, changed = true;
    if (IsKeyDown(GLFW_KEY_LEFT_SHIFT)) pos -= speed * up, changed = true;
    target = pos + ahead;

    /* Update the camera state */
    ahead = normalize(target - pos);
    up = normalize(cross(ahead, right));
    right = normalize(cross(up, ahead));
    tl = pos + 2.0f * ahead - ASPECT_RATIO * right + up;
    tr = pos + 2.0f * ahead + ASPECT_RATIO * right + up;
    bl = pos + 2.0f * ahead - ASPECT_RATIO * right - up;
    
    if (moved_forward) {
        /* Return the forward delta (for reprojection) */
        return forward ? speed : -speed;
    }
    return 0;
}
