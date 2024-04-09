#include "game.h"

#include "dev/gui.h"

/**
 * ===== Game Initialize =====
 */
void Game::init() {
    if (renderer) delete renderer;
    renderer = new Renderer();
    renderer->init();

    player = new Player(renderer->camera);

    DisableCursor();
}

/**
 * ===== Game Tick / Update =====
 */
void Game::tick(const f32 dt) { 
    renderer->tick(screen, dt);

    renderer->scene.tick(dt);

    /* Update the camera */
#if 0
    renderer->depth_delta = renderer->camera.freecam(dt);

    if (not escaped) {
        const int2 delta = mouse_pos - mouse_old;
        mouse_old = mouse_pos;
        renderer->camera.look(delta, dt);
    }
#else
    if (escaped) return;
    const int2 delta = mouse_pos - mouse_old;
    mouse_old = mouse_pos;
    renderer->depth_delta = player->tick(dt, delta);
#endif
    renderer->camera.tick();
}

/**
 * ===== Game GUI Drawing =====
 */
void Game::gui(const f32 dt) {
    /* Development GUI */
    devgui_stats(dt);
    devgui_control();

#ifdef DEV
    /* Fast mode switch */
    static bool f_down = false;
    if (IsKeyDown(GLFW_KEY_F) && f_down == false) {
        if (dev::display_mode == dev::DM::FINAL) {
            dev::display_mode = dev::DM::ALBEDO;
        } else {
            dev::display_mode = dev::DM::FINAL;
        }
        f_down = true;
    }
    if (!IsKeyDown(GLFW_KEY_F) && f_down == true) {
        f_down = false;
    }

    /* Dev gui switch */
    static bool tilde_down = false;
    if (IsKeyDown(GLFW_KEY_GRAVE_ACCENT) && !tilde_down) {
        dev::hide_devgui = !dev::hide_devgui;
        tilde_down = true;
    }
    if (!IsKeyDown(GLFW_KEY_GRAVE_ACCENT) && tilde_down) {
        tilde_down = false;
    }
#endif

    /* Escape logic */
    static bool escape_down = false;
    if (IsKeyDown(GLFW_KEY_ESCAPE) && !escape_down) {
        if (escaped) {
            running = false;
        } else {
            escaped = true;
            EnableCursor();
        }
        escape_down = true;
    }
    if (!IsKeyDown(GLFW_KEY_ESCAPE) && escape_down) {
        escape_down = false;
    }

    renderer->gui(running, dt);
}

/**
 * ===== Game Shutdown =====
 */
void Game::shutdown() { renderer->shutdown(); }

/**
 * ===== Input Handling =====
 */
void Game::MouseDown(int button) {
#ifdef DEV
    /* Don't listen to mouse if it's over ImGui windows */
    if (ImGui::GetIO().WantCaptureMouse) return;
#endif

    /* Escape logic */
    if (escaped) {
#ifdef DEV
        if (button == 0) {
            const Camera& cam = renderer->camera;
            const Ray ray = cam.get_primary_ray(mouse_pos.x, mouse_pos.y);
            dev::debug_ray = ray;
            dev::debug_ray.debug = true;

            const int2 mp = floori(float2(mouse_pos) / 16.0f) * 16;
            const Ray ray_tl = cam.get_primary_ray(mp.x, mp.y);
            const Ray ray_tr = cam.get_primary_ray(mp.x + 16, mp.y);
            const Ray ray_bl = cam.get_primary_ray(mp.x, mp.y + 16);
            dev::debug_py = Pyramid(cam.pos, normalize(cam.target - cam.pos), ray_tl.dir,
                                    ray_tr.dir, ray_bl.dir);

            dev::debug_packet = cam.get_packet8x8(mouse_pos.x, mouse_pos.y);
        }
#endif

        mouse_old = mouse_pos, escaped = false;
        DisableCursor();
        return;
    } else {
        vector<float3> laser_path = renderer->path(
            Ray(renderer->camera.pos, normalize(renderer->camera.target - renderer->camera.pos)));

        for (u32 i = 0; i < 8; i++) {
            if (i < laser_path.size() - 1) {
                const float3 a = laser_path[i];
                const float3 b = laser_path[i + 1];
                renderer->scene.laser_segments[i].a = a;
                renderer->scene.laser_segments[i].b = b;
            } else {
                renderer->scene.laser_segments[i].a = 1000;
                renderer->scene.laser_segments[i].b = 1000;
            }
        }
    }
}
