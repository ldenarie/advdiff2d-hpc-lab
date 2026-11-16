#include <advdiff/kernel_registry.hpp>

#include <algorithm>

std::optional<advdiff::StepFunction> advdiff::find_kernel(std::string_view name) {
    if (auto it = std::ranges::find(advdiff::kernels, name, &KernelInfo::name); it != advdiff::kernels.cend()) {
        return it->step;
    }
    return std::nullopt;
}
