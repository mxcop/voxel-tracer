#include "dev.h"

#ifdef DEV

namespace dev {
bool hide_devgui = false;
Renderer* renderer = nullptr;
f32 frame_time = 1;

DisplayMode display_mode = FINAL;
}  // namespace dev

#endif
