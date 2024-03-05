#pragma once

#ifdef DEV

/**
 * @brief Update the development statistics GUI.
 */
extern void devgui_stats(const f32 dt);

/**
 * @brief Update the development control GUI.
 */
extern void devgui_control();

#else

/**
 * @brief DOES NOTHING IF "DEV" IS UNDEFINED.
 */
inline void devgui_stats(const f32 dt){};

/**
 * @brief DOES NOTHING IF "DEV" IS UNDEFINED.
 */
inline void devgui_control(){};

#endif
