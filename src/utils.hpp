#ifndef __UTILS_H__
#define __UTILS_H__

#include <vector>
#include <cmath>

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

    ~Vector() {};
};

class Matrix
{
public:
    Matrix() {};
    ~Matrix() {};
};

#endif  // UTILS_H
