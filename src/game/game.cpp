#include "game.h"

#include "gui.h"

#include "dev/gui.h"
#include "dev/debug.h"

/**
 * ===== Game Initialize =====
 */
void Game::init() {
    if (renderer) delete renderer;
    renderer = new Renderer();
    renderer->init();

    player = new Player(renderer->camera);

    for (u32 i = 0; i < 4; i++) {
        enemies[i] = new Enemy(renderer->scene.enemies[i]);
    }

    DisableCursor();
}

/**
 * ===== Game Tick / Update =====
 */
void Game::tick(const f32 dt) {
    if (state == GameState::MENU || state == GameState::GAMEOVER) {
        renderer->camera.pos = {-1.86f, 1.455f, 2.69f};
        renderer->camera.target = {-1.096f, 1.1667f, 2.1135f};
        renderer->camera.tick();

        renderer->scene.tick(dt);
        renderer->tick(screen, dt);

        EnableCursor();
        return;
    }

    /* Update enemies */
    for (u32 i = 0; i < 4; i++) {
        if (enemies[i]->tick(dt, renderer->camera.pos, &enemies[0], 4)) {
            state = GameState::GAMEOVER;
        }
    }

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

    /* Update laser */
    for (u32 i = 0; i < 8; i++) {
        renderer->scene.laser_segments[i].a = 1000;
        renderer->scene.laser_segments[i].b = 1000;
    }
    renderer->scene.tick(dt);

    vector<float3> laser_path =
        renderer->path(Ray(renderer->camera.pos - float3(0, 0.1f, 0),
                           normalize(renderer->camera.target - renderer->camera.pos)));

    for (u32 i = 0; i < 8; i++) {
        if (i < laser_path.size() - 1) {
            const float3 a = laser_path[i];
            const float3 b = laser_path[i + 1];
            renderer->scene.laser_segments[i].a = a;
            renderer->scene.laser_segments[i].b = b;
        }
    }

    const float3 last_point = laser_path[laser_path.size() - 1];
    const float3 second_last_point = laser_path[laser_path.size() - 2];
    const float3 laser_dir = normalize(last_point - second_last_point);
    const Ray laser_ray = Ray(second_last_point, laser_dir);
    
    for (u32 i = 0; i < 4; i++) {
        if (enemies[i]->process_hit(laser_ray)) {
            score += 100;
        }
    }

    renderer->scene.tick(dt);
    renderer->tick(screen, dt);
}

/**
 * ===== Game GUI Drawing =====
 */
void Game::gui(const f32 dt) {
    /* Development GUI */
    devgui_stats(dt);

    switch (state) {
        case GameState::MENU: {
            /* Window composition */
            float menu_width = centered_overlay();
            const ImVec2 BUTTON_SIZE = ImVec2(menu_width * 0.5f, 32.0f);

            /* Draw content */
            if (ImGui::Begin("Menu", nullptr, overlay_flags)) {
                centered_text("Main Menu");
                ImGui::Separator();

                if (aligned_button_h("Play", 0.5f, BUTTON_SIZE)) {
                    state = GameState::GAME;
                    score = 0;
                    for (u32 i = 0; i < 4; i++) {
                        enemies[i]->respawn();
                    }
                    DisableCursor();
                }
                if (aligned_button_h("Quit", 0.5f, BUTTON_SIZE)) {
                    running = false;
                }
            }
            ImGui::End();

            break;
        }
        case GameState::GAME: {
            /* Window composition */
            centered_score();

            /* Draw content */
            if (ImGui::Begin("Score", nullptr, overlay_flags)) {
                centered_text("Score: %i", score);
            }
            ImGui::End();

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

            break;
        }
        case GameState::GAMEOVER: {
            /* Window composition */
            float menu_width = centered_overlay();
            const ImVec2 BUTTON_SIZE = ImVec2(menu_width * 0.5f, 32.0f);

            /* Draw content */
            if (ImGui::Begin("Menu", nullptr, overlay_flags)) {
                centered_text("Game Over");
                ImGui::Separator();

                centered_text("Score: %i", score);

                if (aligned_button_h("Try again", 0.5f, BUTTON_SIZE)) {
                    state = GameState::GAME;
                    score = 0;
                    for (u32 i = 0; i < 4; i++) {
                        enemies[i]->respawn();
                    }
                    DisableCursor();
                }
                if (aligned_button_h("Quit", 0.5f, BUTTON_SIZE)) {
                    running = false;
                }
            }
            ImGui::End();

            break;
        }
        default:
            break;
    }

    /* Development GUI */
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
            const Ray ray = cam.get_primary_ray((f32)mouse_pos.x, (f32)mouse_pos.y);
            dev::debug_ray = ray;
            dev::debug_ray.debug = true;

            const float2 mp = floorf(float2(mouse_pos) / 16.0f) * 16.0f;
            const Ray ray_tl = cam.get_primary_ray(mp.x, mp.y);
            const Ray ray_tr = cam.get_primary_ray(mp.x + 16.0f, mp.y);
            const Ray ray_bl = cam.get_primary_ray(mp.x, mp.y + 16.0f);
            dev::debug_py = Pyramid(cam.pos, normalize(cam.target - cam.pos), ray_tl.dir,
                                    ray_tr.dir, ray_bl.dir);

            dev::debug_packet = cam.get_packet8x8((f32)mouse_pos.x, (f32)mouse_pos.y);
        }
#endif

        if (button == 0) {
            mouse_old = mouse_pos, escaped = false;
            DisableCursor();
            return;
        }
    }
}
