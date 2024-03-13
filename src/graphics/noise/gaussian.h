#pragma once

/* Credit : <https://github.com/bfraboni/FastGaussianBlur/blob/main/fast_gaussian_blur_template.h> */
inline void horizontal_blur(const float4* in, float4* out, const i32 w, const i32 h,
                                        const i32 r) {
    const f32 iarr = 1.f / (r + r + 1);
    const f32 iwidth = 1.f / w;
#pragma omp parallel for
    for (i32 i = 0; i < h; i++) {
        const i32 begin = i * w;
        const i32 end = begin + w;
        f32 acc[3] = {0};

        // current index, left index, right index
        i32 ti = begin, li = begin - r - 1, ri = begin + r;

        // initial accumulation
        for (i32 j = ti; j < ri; j++)
            for (i32 ch = 0; ch < 3; ++ch) {
                acc[ch] += in[j][ch];
            }

        // 1. left side out and right side in
        for (; li < begin; ri++, ti++, li++)
            for (i32 ch = 0; ch < 3; ++ch) {
                acc[ch] += in[ri][ch];
                const f32 inorm = 1.f / f32(ri + 1 - begin);
                out[ti][ch] = acc[ch] * inorm;
            }

        // 2. left side in and right side in
        for (; ri < end; ri++, ti++, li++)
            for (i32 ch = 0; ch < 3; ++ch) {
                acc[ch] += in[ri][ch] - in[li][ch];
                out[ti][ch] = acc[ch] * iarr;
            }

        // 3. left side in and right side out
        for (; ti < end; ti++, li++)
            for (i32 ch = 0; ch < 3; ++ch) {
                acc[ch] -= in[li][ch];
                const f32 inorm = 1.f / f32(end - li - 1);
                out[ti][ch] = acc[ch] * inorm;
            }
    }
}

/**
 * @brief This function performs a 2D tranposition of an image.
 * Credit : <https://github.com/bfraboni/FastGaussianBlur/blob/main/fast_gaussian_blur_template.h>
 */
inline void flip_block(const float4* in, float4* out, const i32 w, const i32 h) {
    constexpr i32 block = 256 / 4;
#pragma omp parallel for collapse(2)
    for (i32 x = 0; x < w; x += block)
        for (i32 y = 0; y < h; y += block) {
            const float4* p = in + y * w + x;
            float4* q = out + y + x * h;

            const i32 blockx = std::min(w, x + block) - x;
            const i32 blocky = std::min(h, y + block) - y;
            for (i32 xx = 0; xx < blockx; xx++) {
                for (i32 yy = 0; yy < blocky; yy++) {
                    for (i32 k = 0; k < 4; k++) *q = *p;
                    p += w;
                    q++;
                }
                p += -blocky * w + 1;
                q += -blocky * 1 + h;
            }
        }
}

inline f32 sigma_to_box_radius(i32 boxes[], const f32 sigma, const i32 n) {
    // ideal filter width
    float wi = std::sqrt((12 * sigma * sigma / n) + 1);
    int wl = wi;  // no need std::floor
    if (wl % 2 == 0) wl--;
    int wu = wl + 2;

    float mi = (12 * sigma * sigma - n * wl * wl - 4 * n * wl - 3 * n) / (-4 * wl - 4);
    int m = mi + 0.5f;  // avoid std::round by adding 0.5f and cast to integer type

    for (int i = 0; i < n; i++) boxes[i] = ((i < m ? wl : wu) - 1) / 2;

    return std::sqrt((m * wl * wl + (n - m) * wu * wu - n) / 12.f);
}

inline void fast_gaussian_blur(float4*& in, float4*& out, const i32 w, const i32 h, const i32 n, const f32 sigma) {
    // compute box kernel sizes
    i32* boxes = new i32[n];
    sigma_to_box_radius(boxes, sigma, n);

    // perform N horizontal blur passes
    for (u32 i = 0; i < n; ++i) {
        horizontal_blur(in, out, w, h, boxes[i]);
        std::swap(in, out);
    }

    // flip buffer
    flip_block(in, out, w, h);
    std::swap(in, out);

    // perform N horizontal blur passes on flipped image
    for (u32 i = 0; i < n; ++i) {
        horizontal_blur(in, out, h, w, boxes[i]);
        std::swap(in, out);
    }

    // flip buffer
    flip_block(in, out, h, w);
    delete[] boxes;
}
