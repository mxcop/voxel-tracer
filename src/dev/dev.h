#pragma once

class Renderer;
class Camera;
struct TraceResult;
struct HitInfo;
struct Ray;

namespace Tmpl8 {
class Surface;
}

/* Development namespace (only accessible if "DEV" is defined) */
namespace dev {

#ifdef DEV

/* Dev GUI switch */
extern bool hide_devgui;

extern bool use_projection;

/* Renderer, Camera, and Debug Screen */
extern Renderer* renderer;
extern Camera* main_camera;
extern Tmpl8::Surface* db_screen;

extern Ray debug_ray;

extern f32 frame_time;

enum DisplayMode : int {
	FINAL,
	ALBEDO,
	NORMALS,
	DEPTH,
	PRIMARY_STEPS,
	SECONDARY_STEPS,
	AMBIENT_STEPS
};
typedef DisplayMode DM;
extern DisplayMode display_mode;

/**
 * @brief Display development display modes.
 */
extern bool display_modes(const HitInfo& hit, TraceResult& r);

#else

/**
 * @brief DOES NOTHING IF "DEV" IS UNDEFINED.
 */
inline bool display_modes(const HitInfo& hit) { return false; };

#endif

}  // namespace dev
