# Development

## Current build dependencies

- **Eigen 3.4.0** — header-only linear algebra
- **Coal 3.0.4** — shared collision library (enabled by default)
- **Catch2 3.8.1** — unit tests
- **CMake 3.23+**, Ninja, and a C++20 compiler

Conan 2 supplies the native dependencies. `conanfile.py` declares versions and
options; `conan.lock` fixes the resolved dependency versions and recipe revisions,
including transitive dependencies. Coal's Python bindings, OctoMap, and Qhull are
not needed by the current adapter and are disabled. Assimp's unused 3MF exporter
is also disabled, avoiding its kuba-zip dependency with a failing source checksum.
Viewer Python dependencies remain separate in `tools/viewer/requirements.txt`.

## Building

Install the build tools in a local environment and detect your compiler once:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r tools/requirements-build.txt
conan profile detect
```

Inspect the detected profile before building. The repository's release/debug
profiles inherit its OS, architecture, and compiler settings and select C++20.
On this Mac, the compiler is Apple Clang 17 targeting ARM64. Linux users should
use their own detected compiler profile rather than copying the Mac settings.

```bash
conan install . -pr:h=conan/profiles/release -pr:b=default --build=missing
cmake --preset conan-release
cmake --build --preset conan-release -j 2
ctest --preset conan-release --output-on-failure
```

Conan automatically reads `conan.lock`, downloads compatible binaries when
available, and builds missing ones. The first install can take several minutes.
Executables are under `build/Release/`, for example
`build/Release/bike_rollout`. No cmeel environment or manually configured
`CMAKE_PREFIX_PATH` is needed. Generated presets, build output, and the local
Python environment are ignored by Git. Use a fresh build directory when switching
from the previous dependency setup.

For Debug, use `conan/profiles/debug` and the `conan-debug` presets. To omit
collision support or tests, add `-o '&:with_coal=False'` or
`-o '&:with_tests=False'` to `conan install`; these options also configure CMake.
Tests are enabled by default and require Catch2 rather than silently skipping it.

Coal remains a shared library. Development builds use the Conan cache at runtime;
keep that cache available. Distributing an executable requires packaging its
shared-library dependencies separately.

To deliberately refresh the dependency lock after editing the manifest:

```bash
conan lock create . -pr:h=conan/profiles/release -pr:b=default --lockfile="" --lockfile-out=conan.lock
```

Review the lockfile change and repeat the build/tests. A lockfile fixes dependency
recipes, not the local compiler or SDK; those remain part of the build environment.

## C++ formatting

Use clang-format with the repository's `.clang-format` configuration. Control
statements require braces, including single-statement loops and conditionals.
The configuration requires clang-format 15 or newer; the current formatting pass
used 23.1.5.

```bash
# Apply formatting (requires ripgrep).
rg --files -0 rapide examples test -g '*.hpp' -g '*.cpp' | xargs -0 clang-format -i

# Check formatting without changing files.
rg --files -0 rapide examples test -g '*.hpp' -g '*.cpp' | xargs -0 clang-format --dry-run --Werror
```
