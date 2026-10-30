#pragma once

#include <cstddef>
#include <vector>

namespace advdiff {

// 2D scalar field with one layer of ghost (halo) cells. See docs/model.md, section 4.
//
// Interior cells: i in [0, nx), j in [0, ny).
// Ghost cells:    i = -1 or i = nx, j = -1 or j = ny.
// Row-major storage, x is the fastest index.
class Field {
public:
    Field(int nx, int ny)
        : nx_(nx), ny_(ny), stride_(nx + 2),
          data_(static_cast<std::size_t>(nx + 2) * static_cast<std::size_t>(ny + 2), 0.0) {}

    int nx() const { return nx_; }
    int ny() const { return ny_; }

    // Distance in memory between (i, j) and (i, j + 1).
    int stride() const { return stride_; }

    std::size_t index(int i, int j) const {
        return static_cast<std::size_t>(j + 1) * static_cast<std::size_t>(stride_)
             + static_cast<std::size_t>(i + 1);
    }

    double& operator()(int i, int j) { return data_[index(i, j)]; }
    double operator()(int i, int j) const { return data_[index(i, j)]; }

    // Raw access, for kernels that want to work on pointers.
    double* data() { return data_.data(); }
    const double* data() const { return data_.data(); }

private:
    int nx_;
    int ny_;
    int stride_;
    std::vector<double> data_;
};

}  // namespace advdiff
