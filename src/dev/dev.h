#pragma once

#ifdef DEV

class Renderer;
class Camera;

namespace Tmpl8 {
class Surface;
}

/* Development namespace (only accessible if "DEV" is defined) */
namespace dev {
/* Dev GUI switch */
extern bool hide_devgui;

/* Renderer, Camera, and Debug Screen */
extern Renderer* renderer;
extern Camera* main_camera;
extern Tmpl8::Surface* db_screen;

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
}  // namespace dev

#endif
