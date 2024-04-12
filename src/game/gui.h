#pragma once

#include <imgui.h>

constexpr ImGuiWindowFlags overlay_flags =
    ImGuiWindowFlags_NoDecoration | ImGuiWindowFlags_AlwaysAutoResize |
    ImGuiWindowFlags_NoSavedSettings | ImGuiWindowFlags_NoFocusOnAppearing |
    ImGuiWindowFlags_NoNav | ImGuiWindowFlags_NoMove;

/* Create a horizontally centered ImGui text label */
void centered_text(const char* txt) {
    auto win_w = ImGui::GetWindowSize().x;
    auto txt_w = ImGui::CalcTextSize(txt).x;

    ImGui::SetCursorPosX((win_w - txt_w) * 0.5f);
    ImGui::Text(txt);
}

void centered_text(const char* txt, const u32 num) {
    auto win_w = ImGui::GetWindowSize().x;
    auto txt_w = ImGui::CalcTextSize(txt).x;

    ImGui::SetCursorPosX((win_w - txt_w) * 0.5f);
    ImGui::Text(txt, num);
}

/* Horizontally aligned ImGui button */
/* Source : <https://github.com/ocornut/imgui/discussions/3862> */
bool aligned_button_h(const char* label, float alignment, ImVec2 size) {
    ImGuiStyle& style = ImGui::GetStyle();

    float calc_size = ImGui::CalcTextSize(label).x + style.FramePadding.x * 2.0f;
    if (size.x != 0.0f) calc_size = size.x;
    float avail = ImGui::GetContentRegionAvail().x;

    float off = (avail - calc_size) * alignment;
    if (off > 0.0f) ImGui::SetCursorPosX(ImGui::GetCursorPosX() + off);

    if (size.x != 0.0f) {
        return ImGui::Button(label, size);
    } else {
        return ImGui::Button(label);
    }
}

/* Setup a centered overlay window, returns the overlay width. */
float centered_overlay() {
    const ImVec2 work_pos = ImGui::GetMainViewport()->WorkPos;
    const ImVec2 work_size = ImGui::GetMainViewport()->WorkSize;
    const ImVec2 menu_size = ImVec2(work_size.x * 0.2f, work_size.y * 0.5f);
    ImGui::SetNextWindowPos(ImVec2(work_pos.x + work_size.x * 0.5f - menu_size.x * 0.5f,
                            work_pos.y + work_size.y * 0.5f - menu_size.y * 0.5f));
    ImGui::SetNextWindowSize(menu_size);
    ImGui::SetNextWindowBgAlpha(0.35f);
    return menu_size.x;
}

float centered_score() {
    const ImVec2 work_pos = ImGui::GetMainViewport()->WorkPos;
    const ImVec2 work_size = ImGui::GetMainViewport()->WorkSize;
    const ImVec2 menu_size = ImVec2(128.0f, 32.0f);
    ImGui::SetNextWindowPos(ImVec2(work_pos.x + work_size.x * 0.5f - menu_size.x * 0.5f,
                                   16.0f + menu_size.y));
    ImGui::SetNextWindowSize(menu_size);
    ImGui::SetNextWindowBgAlpha(0.35f);
    return menu_size.x;
}
