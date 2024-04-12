#pragma once

#include "graphics/renderer.h"

#include "player.h"
#include "enemy.h"

enum class GameState {
    MENU,
    GAME,
    GAMEOVER
};

class Game : public TheApp {
    /* User input */
    int2 mouse_old = 0, mouse_pos = 0;
    bool escaped = false;

    Player* player = nullptr;
    Enemy* enemies[4];

    GameState state = GameState::MENU;
    u32 score = 0;

   public:
    /* Ray tracer */
    Renderer* renderer = nullptr;

    /** @brief Initialize the application */
    void init();
    /** @brief Called as often as possible */
    void tick(const f32 dt);
    /** @brief Called after "tick" for ImGUI rendering */
    void gui(const f32 dt);
    /** @brief Called before the application closes */
    void shutdown();

    /* (Tmpl8) User input */
    void MouseUp(int) {}
    void MouseDown(int button);
    void MouseMove(int x, int y) { mouse_pos.x = x, mouse_pos.y = y; }
    void MouseWheel(float) {}
    void KeyUp(int) {}
    void KeyDown(int) {}
};
