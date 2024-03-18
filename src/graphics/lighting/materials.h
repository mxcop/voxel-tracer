#pragma once

/**
 * @brief Material row table.
 */
enum class MaterialRow : u32 {
    GLASS = 0,   /* ID : 0..8 */
    MIRROR = 1,  /* ID : 8..16 */
    METAL = 2, /* ID : 16..24 */
    UNUSED_3 = 3, /* ID : 24..32 */
};
