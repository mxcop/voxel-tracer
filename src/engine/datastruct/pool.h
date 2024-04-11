/**
 * @file pool.h
 * @author Max (mxcop)
 * @brief Memory pool.
 * @date 2023-07-28
 */
#pragma once
#include <functional>

/**
 * @brief (Fixed size) Memory pool implementation.
 */
template <typename T>
class Pool {
    T* data = nullptr;
    bool* status = nullptr;

    size_t capacity = 0;
    size_t size = 0;

   public:
    Pool() = default;
    Pool(size_t capacity = 16) {
        data = (T*)std::malloc(capacity * sizeof T);
        status = new bool[capacity]{0};
        this->capacity = capacity;
        size = 0;
    }

    ~Pool() {
        delete[] data;
        delete[] status;
    }

    /**
     * @brief Add a new element into the pool.
     *
     * @param element The element to add.
     * @return The index of the added element. (-1 if pool is full)
     */
    i32 add(T element) {
        /* Check if there's an inactive element we can override */
        for (size_t i = 0; i < size; i++) {
            if (status[i]) continue;

            data[i] = element;
            status[i] = true;
            return (i32)i;
        }

        /* Exit if the pool is full */
        if (size >= capacity) {
            return -1;
        }

        /* Insert the new element */
        data[size] = element;
        status[size] = true;

        size++;
        return (i32)size - 1;
    }

    /**
     * @brief Get an element by index.
     *
     * @param idx The index of the element.
     * @return A reference to the element.
     */
    T& get(size_t idx) {
        if (idx < size) return data[idx];
        return data[0]; /* Silly?! */
    }

    /**
     * @brief Deactivate an element in the pool.
     *
     * @param idx The index of the element to deactivate.
     * @return True if the element was deactivated.
     */
    bool remove(size_t idx) {
        if (idx >= size) return false;

        status[idx] = false;
        return true;
    }

    /**
     * @brief Iterate over all active elements in the pool.
     */
    void iterate(std::function<void(T&)> iterator) {
        for (size_t i = 0; i < size; i++) {
            if (status[i] == false) continue;

            iterator(data[i]);
        }
    }

    size_t get_size() { return size; }
    size_t get_capacity() { return capacity; }

    /* Iterator */
    struct Iterator {
        // Iterator tags here...

        // Iterator constructors here...
        Iterator(T* ptr) { m_ptr = ptr; }

        T& operator*() const { return *m_ptr; }
        T* operator->() { return m_ptr; }

        // Prefix increment
        Iterator& operator++() {
            m_ptr++;
            return *this;
        }

        // Postfix increment
        Iterator operator++(int) {
            Iterator tmp = *this;
            ++(*this);
            return tmp;
        }

        friend bool operator==(const Iterator& a, const Iterator& b) { return a.m_ptr == b.m_ptr; };
        friend bool operator!=(const Iterator& a, const Iterator& b) { return a.m_ptr != b.m_ptr; };

       private:
        T* m_ptr;
    };

    Iterator begin() { return Iterator(data); };
    Iterator end() { return Iterator(data + size); };
};
