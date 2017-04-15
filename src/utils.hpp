#ifndef __UTILS_H__
#define __UTILS_H__

#include <vector>
#include <cmath>
#include <cassert>

class Vector
{
private:
    std::vector<double> data;
public:
    Vector(size_t size) {
        data = std::vector<double>(size);
    };

    void set(size_t idx, double val) {
        data[idx] = val;
    }

    double get(size_t idx) const {
        return data[idx];
    }

    size_t size() const {
        return data.size();
    }

    double norm2() const {
        double res = .0f;

        for (double num: data) {
            res += num * num;
        }
        return sqrt(res);
    }

    double dotV(Vector other) const {
        double sum = .0f;

        assert(data.size() == other.size());

        for (size_t i = 0; i < data.size(); i++) {
            sum += data[i] * other.get(i);
        }

        return sum;
    }

    Vector add(Vector other) const {
        assert(data.size() == other.size());

        Vector res(data.size());

        for (size_t i = 0; i < data.size(); i++) {
            res.set(i, data[i] + other.get(i));
        }

        return res;
    }

    // same as add, but inplace
    Vector iadd(Vector other) {
        assert(data.size() == other.size());

        for (size_t i = 0; i < data.size(); i++) {
            data[i] += other.get(i);
        }

        return *this;
    }

    Vector sub(Vector other) const {
        assert(data.size() == other.size());

        Vector res(data.size());

        for (size_t i = 0; i < data.size(); i++) {
            res.set(i, data[i] - other.get(i));
        }

        return res;
    }

    // same as sub, but inplace
    Vector isub(Vector other) {
        assert(data.size() == other.size());

        for (size_t i = 0; i < data.size(); i++) {
            data[i] -= other.get(i);
        }

        return *this;
    }

    Vector mulS(double scaler) const {
        Vector res(data.size());

        for (size_t i = 0; i < data.size(); i++) {
            res.set(i, data[i] * scaler);
        }

        return res;
    }

    // same as mulS, but inplace
    Vector imulS(double scaler) {
        for (size_t i = 0; i < data.size(); i++) {
            data[i] *= scaler;
        }

        return *this;
    }

    ~Vector() {};
};

class Matrix
{
public:
    Matrix() {};
    ~Matrix() {};
};

#endif  // UTILS_H
