#pragma once

/* 8x8x8 Hilbert curve LUT. */
const u32 HILBERT_512[512] = {
    0x0,   0x7,   0x8,   0xB,   0xD4,  0xD3,  0xCC,  0xCB,  0x3,   0x4,   0x9,   0xA,   0xD7,
    0xD0,  0xCF,  0xC8,  0x3C,  0x3B,  0x36,  0x35,  0xD8,  0xDB,  0xC6,  0xC7,  0x3F,  0x38,
    0x37,  0x34,  0xD9,  0xDA,  0xC1,  0xC0,  0x40,  0x43,  0x7C,  0x7F,  0x80,  0x83,  0xBC,
    0xBF,  0x41,  0x42,  0x7D,  0x7E,  0x81,  0x82,  0xBD,  0xBE,  0x5A,  0x5D,  0x62,  0x65,
    0x9A,  0x9D,  0xA2,  0xA5,  0x59,  0x5E,  0x61,  0x66,  0x99,  0x9E,  0xA1,  0xA6,  0x1,
    0x6,   0xF,   0xC,   0xD5,  0xD2,  0xCD,  0xCA,  0x2,   0x5,   0xE,   0xD,   0xD6,  0xD1,
    0xCE,  0xC9,  0x3D,  0x3A,  0x31,  0x32,  0xDF,  0xDC,  0xC5,  0xC4,  0x3E,  0x39,  0x30,
    0x33,  0xDE,  0xDD,  0xC2,  0xC3,  0x47,  0x44,  0x7B,  0x78,  0x87,  0x84,  0xBB,  0xB8,
    0x46,  0x45,  0x7A,  0x79,  0x86,  0x85,  0xBA,  0xB9,  0x5B,  0x5C,  0x63,  0x64,  0x9B,
    0x9C,  0xA3,  0xA4,  0x58,  0x5F,  0x60,  0x67,  0x98,  0x9F,  0xA0,  0xA7,  0x1A,  0x1B,
    0x10,  0x13,  0xEA,  0xED,  0xF2,  0xF5,  0x1D,  0x1C,  0x11,  0x12,  0xE9,  0xEE,  0xF1,
    0xF6,  0x22,  0x23,  0x2E,  0x2D,  0xE0,  0xE3,  0xFA,  0xFB,  0x25,  0x24,  0x2F,  0x2C,
    0xE1,  0xE2,  0xFD,  0xFC,  0x48,  0x49,  0x76,  0x77,  0x88,  0x89,  0xB6,  0xB7,  0x4F,
    0x4E,  0x71,  0x70,  0x8F,  0x8E,  0xB1,  0xB0,  0x50,  0x51,  0x6E,  0x6F,  0x90,  0x91,
    0xAE,  0xAF,  0x57,  0x56,  0x69,  0x68,  0x97,  0x96,  0xA9,  0xA8,  0x19,  0x18,  0x17,
    0x14,  0xEB,  0xEC,  0xF3,  0xF4,  0x1E,  0x1F,  0x16,  0x15,  0xE8,  0xEF,  0xF0,  0xF7,
    0x21,  0x20,  0x29,  0x2A,  0xE7,  0xE4,  0xF9,  0xF8,  0x26,  0x27,  0x28,  0x2B,  0xE6,
    0xE5,  0xFE,  0xFF,  0x4B,  0x4A,  0x75,  0x74,  0x8B,  0x8A,  0xB5,  0xB4,  0x4C,  0x4D,
    0x72,  0x73,  0x8C,  0x8D,  0xB2,  0xB3,  0x53,  0x52,  0x6D,  0x6C,  0x93,  0x92,  0xAD,
    0xAC,  0x54,  0x55,  0x6A,  0x6B,  0x94,  0x95,  0xAA,  0xAB,  0x1E6, 0x1E7, 0x1E8, 0x1EB,
    0x114, 0x113, 0x10C, 0x10B, 0x1E1, 0x1E0, 0x1E9, 0x1EA, 0x117, 0x110, 0x10F, 0x108, 0x1DE,
    0x1DF, 0x1D6, 0x1D5, 0x118, 0x11B, 0x106, 0x107, 0x1D9, 0x1D8, 0x1D7, 0x1D4, 0x119, 0x11A,
    0x101, 0x100, 0x1B4, 0x1B5, 0x18A, 0x18B, 0x174, 0x175, 0x14A, 0x14B, 0x1B3, 0x1B2, 0x18D,
    0x18C, 0x173, 0x172, 0x14D, 0x14C, 0x1AC, 0x1AD, 0x192, 0x193, 0x16C, 0x16D, 0x152, 0x153,
    0x1AB, 0x1AA, 0x195, 0x194, 0x16B, 0x16A, 0x155, 0x154, 0x1E5, 0x1E4, 0x1EF, 0x1EC, 0x115,
    0x112, 0x10D, 0x10A, 0x1E2, 0x1E3, 0x1EE, 0x1ED, 0x116, 0x111, 0x10E, 0x109, 0x1DD, 0x1DC,
    0x1D1, 0x1D2, 0x11F, 0x11C, 0x105, 0x104, 0x1DA, 0x1DB, 0x1D0, 0x1D3, 0x11E, 0x11D, 0x102,
    0x103, 0x1B7, 0x1B6, 0x189, 0x188, 0x177, 0x176, 0x149, 0x148, 0x1B0, 0x1B1, 0x18E, 0x18F,
    0x170, 0x171, 0x14E, 0x14F, 0x1AF, 0x1AE, 0x191, 0x190, 0x16F, 0x16E, 0x151, 0x150, 0x1A8,
    0x1A9, 0x196, 0x197, 0x168, 0x169, 0x156, 0x157, 0x1FE, 0x1F9, 0x1F0, 0x1F3, 0x12A, 0x12D,
    0x132, 0x135, 0x1FD, 0x1FA, 0x1F1, 0x1F2, 0x129, 0x12E, 0x131, 0x136, 0x1C2, 0x1C5, 0x1CE,
    0x1CD, 0x120, 0x123, 0x13A, 0x13B, 0x1C1, 0x1C6, 0x1CF, 0x1CC, 0x121, 0x122, 0x13D, 0x13C,
    0x1B8, 0x1BB, 0x184, 0x187, 0x178, 0x17B, 0x144, 0x147, 0x1B9, 0x1BA, 0x185, 0x186, 0x179,
    0x17A, 0x145, 0x146, 0x1A4, 0x1A3, 0x19C, 0x19B, 0x164, 0x163, 0x15C, 0x15B, 0x1A7, 0x1A0,
    0x19F, 0x198, 0x167, 0x160, 0x15F, 0x158, 0x1FF, 0x1F8, 0x1F7, 0x1F4, 0x12B, 0x12C, 0x133,
    0x134, 0x1FC, 0x1FB, 0x1F6, 0x1F5, 0x128, 0x12F, 0x130, 0x137, 0x1C3, 0x1C4, 0x1C9, 0x1CA,
    0x127, 0x124, 0x139, 0x138, 0x1C0, 0x1C7, 0x1C8, 0x1CB, 0x126, 0x125, 0x13E, 0x13F, 0x1BF,
    0x1BC, 0x183, 0x180, 0x17F, 0x17C, 0x143, 0x140, 0x1BE, 0x1BD, 0x182, 0x181, 0x17E, 0x17D,
    0x142, 0x141, 0x1A5, 0x1A2, 0x19D, 0x19A, 0x165, 0x162, 0x15D, 0x15A, 0x1A6, 0x1A1, 0x19E,
    0x199, 0x166, 0x161, 0x15E, 0x159};