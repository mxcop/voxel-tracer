#pragma once

#ifdef DEV

class Renderer;

/* Development namespace (only accessible if "DEV" is defined) */
namespace dev {
extern bool hide_devgui;
extern Renderer* renderer;
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
