#include "dev.h"

#ifdef DEV

namespace dev {
bool hide_devgui = false;
bool use_projection = true;
Renderer* renderer = nullptr;
Camera* main_camera = nullptr;
Tmpl8::Surface* db_screen = nullptr;
f32 frame_time = 1;

Ray debug_ray;
CoherentPacket4x4 debug_packet;

DisplayMode display_mode = FINAL;

/**
 * @brief Display development display modes.
 */
bool display_modes(const HitInfo& hit, TraceResult& r) {
    /* Handle special display modes, for debugging */
    if (dev::display_mode != dev::DM::FINAL) {
        if (dev::display_mode != dev::DM::PRIMARY_STEPS) {
            if (hit.depth >= BIG_F32) {
                r.albedo = float4(0.06f); /* 0xFF101010 */
                r.no_reproject();
                return true;
            }
        }

        switch (dev::display_mode) {
            case dev::DM::ALBEDO:
                r.albedo = hit.albedo;
                r.no_reproject();
                return true;
            case dev::DM::NORMALS:
                r.albedo = float4(fmaxf(hit.normal, 0.0f), 1.0f);
                r.no_reproject();
                return true;
            case dev::DM::DEPTH:
                r.albedo = hit.depth / 16.0f;
                r.no_reproject();
                return true;
            case dev::DM::PRIMARY_STEPS:
                r.albedo = hit.steps / 64.0f;
                r.no_reproject();
                return true;
        }
    }

    return false;
}

}  // namespace dev

#endif
