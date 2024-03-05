#include "gui.h"

#ifdef DEV

/* Internal development variables */
namespace inter {
float3 light_color = float3(1);
f32 radius = 0.1f;
f32 power = 16.0f;
}  // namespace inter

/**
 * @brief Update the development statistics GUI.
 */
void devgui_stats(const f32 dt) {
    if (dev::hide_devgui) return;

    /* Window position */
    constexpr float PADDING = 10.0f;
    ImVec2 work_pos = ImGui::GetMainViewport()->WorkPos;
    ImGui::SetNextWindowPos(ImVec2(work_pos.x + PADDING, work_pos.y + PADDING), ImGuiCond_Always,
                            ImVec2(0.0f, 0.0f));

    constexpr ImGuiWindowFlags overlay_flags =
        ImGuiWindowFlags_NoDecoration | ImGuiWindowFlags_AlwaysAutoResize |
        ImGuiWindowFlags_NoSavedSettings | ImGuiWindowFlags_NoFocusOnAppearing |
        ImGuiWindowFlags_NoNav | ImGuiWindowFlags_NoMove;

    /* Display performance */
    ImGui::SetNextWindowBgAlpha(0.35f);
    if (ImGui::Begin("Stats", nullptr, overlay_flags)) {
        f32 frame_time_us = dev::frame_time * 1'000'000.0f;
        f64 win_size = (f64)WIN_WIDTH * WIN_HEIGHT;

        static f32 avg = 10, alpha = 1;
        avg = (1 - alpha) * avg + alpha * dev::frame_time * 1000;
        if (alpha > 0.05f) alpha *= 0.5f;
        f32 fps = 1000.0f / avg, rps = (WIN_WIDTH * WIN_HEIGHT) / avg;

        ImGui::Text("Perf overlay\n");
        ImGui::Separator();
        ImGui::Text("FPS: %.1f", fps);
        ImGui::Separator();
        ImGui::Text("Ray/s: %.2fM", rps / 1000);
        ImGui::Text("Ray time (mean): %.2fns", (frame_time_us / win_size) * 1'000.0);
        ImGui::Text("Ray time (goal): %.2fns", (0.0166666 / win_size) * 1.0e+9);
        ImGui::Separator();
        ImGui::Text("Frame time: %.2fms", dev::frame_time * 1'000.0f);
    }
    ImGui::End();
}

/**
 * @brief Update the development control GUI.
 */
void devgui_control() {
    if (dev::hide_devgui) return;

    if (ImGui::Begin("Control")) {
        /* Change what part of the render process is displayed */
        const char* modes[] = {"Final",         "Albedo",          "Normals", "Depth",
                               "Primary Steps", "Secondary Steps", "Ambient Steps"};
        if (ImGui::Combo("Display Mode", (i32*)&dev::display_mode, modes, IM_ARRAYSIZE(modes))) {
            dev::renderer->reset_accu();
        }

        /* Create sphere lights from the GUI */
        if (ImGui::CollapsingHeader("Sphere Light")) {
            ImGui::ColorEdit3("Color", &inter::light_color.x);
            ImGui::InputFloat("Radius", &inter::radius);
            ImGui::InputFloat("Power", &inter::power);
            if (ImGui::Button("Spawn on camera")) {
                const float3 origin = dev::renderer->camera.pos;
                dev::renderer->area_lights.emplace_back(origin, inter::radius, inter::light_color,
                                                        inter::power);
                dev::renderer->reset_accu();
            }
        }
    }
    ImGui::End();
}

#endif
