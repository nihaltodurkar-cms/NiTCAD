// Copies a device into the consistency layer's plain data (Unit 55, TVR-1; DECISIONS 11.2). The benchmark layer
// depends on everything; this is the one place a device::Device meets verify/consistency.
#pragma once

#include "NiTCAD/device/device.hpp"
#include "NiTCAD/verify/consistency/device_data.hpp"

namespace NiTCAD::verify::benchmark {

[[nodiscard]] consistency::DeviceData device_data(const device::Device& device);

}  // namespace NiTCAD::verify::benchmark
