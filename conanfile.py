from conan import ConanFile
from conan.tools.cmake import CMakeDeps, CMakeToolchain, cmake_layout


class RapideDependencies(ConanFile):
    """Resolve native dependencies and generate the CMake build configuration."""

    settings = "os", "arch", "compiler", "build_type"
    options = {"with_coal": [True, False], "with_tests": [True, False]}
    default_options = {
        "with_coal": True,
        "with_tests": True,
        "coal/*:with_octomap": False,
        "coal/*:with_qhull": False,
        "assimp/*:with_3mf_exporter": False,
    }

    def requirements(self):
        self.requires("eigen/3.4.0")
        if self.options.with_coal:
            self.requires("coal/3.0.4")
        if self.options.with_tests:
            self.requires("catch2/3.8.1")

    def layout(self):
        cmake_layout(self)

    def generate(self):
        CMakeDeps(self).generate()
        toolchain = CMakeToolchain(self)
        toolchain.variables["RAPIDE_WITH_COAL"] = bool(self.options.with_coal)
        toolchain.variables["BUILD_TESTING"] = bool(self.options.with_tests)
        toolchain.generate()
