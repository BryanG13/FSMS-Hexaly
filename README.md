## FSMS-Hexaly: Feeder Service with Madatory Stops with Hexaly

This repository contains a C++17 implementation of a demand‑responsive feeder/shuttle scheduling model solved with the Hexaly optimizer (formerly LocalSolver). It reads passenger and stop data, builds an optimization model, and writes solution artifacts such as bus routes, timetables, and passenger assignments. The code is used in this [academic paper](https://doi.org/10.1002/net.22185).

Outputs are written to `data/output/` (e.g., `xsol_*.txt`, `ysol_*.txt`, `dsol_*.txt`).

## Requirements

- macOS with a modern C++ compiler (Apple Clang recommended)
- CMake >= 3.22
- Hexaly SDK 14.0+ installed
	- Default expected location on macOS: `/opt/hexaly_14_0`
	- Or set `HEXALY_HOME` to your SDK root

## Getting started

### 1) Prepare the Hexaly SDK

Install Hexaly and ensure headers and libraries are available. On macOS the default package usually installs under `/opt/hexaly_14_0`.

If you installed somewhere else, set an environment variable before configuring the build:

```bash
export HEXALY_HOME=/path/to/hexaly_14_0
```

The build system will also accept a CMake hint:

```bash
cmake -DHEXALY_ROOT=/path/to/hexaly_14_0 ..
```

### 2) Configure and build

```bash
cmake -S . -B build
cmake --build build --config Debug
```

Notes:
- On macOS, the build sets rpath to the detected Hexaly library directory to help the executable find the `.dylib` at runtime.
- Ensure your Hexaly license is set up according to Hexaly’s documentation (license file or environment variable).

### 3) Run

Run the executable from the repository root (or from the `build/` directory):

```bash
./build/Hexaly
```

The program expects input files under `data/input/` and writes results to `data/output/`.

## Data layout

Input files (examples provided in `data/input/`):

- `passengers20.txt`
- `arrivals20.txt`
- `departures20.txt`
- `mandatory.txt`
- `optional5.txt`

Output files (written to `data/output/`):

- `xsol_<instance>.txt` — bus routes (stop sequences)
- `ysol_<instance>.txt` — passenger to bus/trip/stop assignment
- `dsol_<instance>.txt` — bus departure timetable

## Project structure

- `src/main.cpp` — optimization model and program entry point
- `data/input/` — sample input data files
- `data/output/` — outputs written by the program
- `CMakeLists.txt` — build configuration, including Hexaly SDK detection
