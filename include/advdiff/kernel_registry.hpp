#pragma once

#include <advdiff/kernel.hpp>

#include <array>
#include <optional>
#include <string_view>
#include <tuple>

namespace advdiff {

using StepFunction = decltype(&step_incremental);

struct KernelInfo {
    std::string_view name;
    StepFunction step;
};

inline constexpr std::array kernels {
    KernelInfo{.name = "incremental", .step = step_incremental},
    KernelInfo{.name = "omp", .step = step_incremental_omp},
    KernelInfo{.name = "naive", .step = step_naive},
    KernelInfo{.name = "weighted", .step = step_weighted},
};

inline constexpr std::string_view default_kernel = "incremental";

inline constexpr auto kernel_names = std::apply([](const auto&... kernel) {
    return std::array{kernel.name...};
}, kernels);

// Returns the kernel with this name, or nothing if there is none.
std::optional<StepFunction> find_kernel(std::string_view name);

} // namespace advdiff
