# Aerodynamics experiments

A set of small MuJoCo demos exploring fluid/aerodynamic forces: wind, wings,
lift/drag, gliding, and a higher-fidelity airfoil-polar model injected via the
`mjcb_passive` callback.

## Models (MJCF)

| File | What it is |
|---|---|
| `wind_demo.xml` | A flat plate in a steady wind (used by `wind_demo.c`). |
| `wind_gui.xml` | A plate tethered to a mast by a ball joint, flutters in wind. |
| `wing.xml` | A single slender wing (used by the sweep programs). |
| `glider.xml` | A pitch-stable glider (wing + tail + fin + nose ballast), tuned `fluidcoef`. |
| `glider_ellipsoid.xml` | Glider using MuJoCo's built-in ellipsoid fluid model. |
| `glider_polar.xml` | Same airframe with **no** built-in fluid; aero supplied by the polar callback. |
| `fly_compare.xml` | Two gliders side by side (ellipsoid vs polar) for the live A/B viewer. |

## Programs

| File | What it does |
|---|---|
| `wind_demo.c` | Reports drag/lift on a plate at several orientations in wind. |
| `wing_sweep.c` | Angle-of-attack and wind-speed sweep for a single wing (ellipsoid model). |
| `glide_test.c` | Simulates the glider trajectory headless. |
| `shape_test.c` | Compares how wing shape (aspect ratio, thickness) affects each model. |
| `airfoil_compare.c` | Built-in ellipsoid vs. a real airfoil polar (with stall): static sweep + glide. |
| `fly_compare.cc` | **Live GLFW viewer**: red glider (ellipsoid) vs blue glider (airfoil polar), HUD. |

## The higher-fidelity aero model

`airfoil_compare.c` and `fly_compare.cc` implement quasi-steady **blade-element
(strip) theory** via `mjcb_passive`: for each lifting surface, compute the local
airflow → angle of attack → look up `Cl(α), Cd(α)` from an analytic airfoil
polar (linear lift → stall at ~14° → flat-plate post-stall, plus induced drag
`Cd0 + Cl^2/(pi*AR*e)`), then apply `L,D = 1/2 rho V^2 S {Cl,Cd}` with
`mj_applyFT`. Surfaces are aerodynamically independent (no downwash/wake).

## Building

These were built against an in-tree MuJoCo build at `<repo>/build` (see the main
docs for `cmake -S . -B build && cmake --build build`). From the repo root:

```bash
MJ=$PWD
# headless C demos:
clang experiments/aero/airfoil_compare.c -o /tmp/airfoil_compare \
  -I $MJ/include -L $MJ/build/lib -lmujoco \
  -Wl,-rpath,$MJ/build/lib -lm

# live viewer (needs the GLFW that the build fetched):
clang++ experiments/aero/fly_compare.cc -o /tmp/fly_compare -std=c++17 \
  -I $MJ/include -I $MJ/build/_deps/glfw3-src/include \
  -L $MJ/build/lib -lmujoco $MJ/build/lib/libglfw3.a \
  -Wl,-rpath,$MJ/build/lib \
  -framework Cocoa -framework IOKit -framework OpenGL -framework CoreVideo
```

> **Note:** the programs currently load their `.xml` models by **absolute
> `/tmp/...` paths**. To run them as-is, copy the XML files to `/tmp` first
> (`cp experiments/aero/*.xml /tmp/`), or edit the paths in the source. macOS
> linker flags shown; on Linux drop the `-framework` flags and link
> `-lGL -lglfw` (or the fetched static lib) instead.
