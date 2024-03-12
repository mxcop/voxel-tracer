#pragma once

class Simplex {
    array<float3, 4> m_points = {};
    i32 size;

   public:
    Simplex() : size(0) {}

    Simplex& operator=(std::initializer_list<float3> list) {
        for (auto v = list.begin(); v != list.end(); v++) {
            m_points[std::distance(list.begin(), v)] = *v;
        }
        size = list.size();

        return *this;
    }

    void push_front(const float3& point) {
        m_points = {point, m_points[0], m_points[1], m_points[2]};
        size = std::min(size + 1, 4);
    }

    float3& operator[](int i) { return m_points[i]; }
    size_t get_size() const { return size; }

    /* Iterator */
    auto begin() const { return m_points.begin(); }
    auto end() const { return m_points.end() - (4 - size); }
};
