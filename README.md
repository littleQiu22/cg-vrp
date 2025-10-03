# CG-VRP

**Column Generation for the Vehicle Routing Problem**

`cg-vrp` is a Python library that leverages the **Column Generation** algorithm to solve the **Linear Programming (LP) Relaxation** of Vehicle Routing Problem (VRP). Currently considered VRP variants: VRPTW.

The computationally intensive pricing subproblem, is implemented in C++ and exposed to Python using `nanobind` for maximum performance.

## Installation

You can install `cg-vrp` in two ways: from the Python Package Index (PyPI) or by building it from the source code.

### 1. From PyPI (Recommended)

The easiest way to install the library is directly from PyPI. This will download a pre-compiled binary wheel, so you won't need a C++ compiler on your system.

```bash
pip install cg-vrp
```

### 2. From Source (For Development)

If you want to modify the code, you should build it from the source.

**Step 1: Clone the repository**
```bash
git clone https://github.com/littleQiu22/cg-vrp.git
cd cg-vrp
```

**Step 2: Install the package**
This command will compile the C++ extension and install the Python package.
```bash
pip install .
```

#### Potential Build Issues

Building from source requires a C++ toolchain and CMake. Here are some common issues you might encounter:

*   **CMake Not Found**: The build system, `scikit-build-core`, uses CMake to configure and build the C++ extension. If CMake is not installed or not found in your system's `PATH`, the build will fail.
    *   **Solution**: Install CMake. You can download it from the [official CMake website](https://cmake.org/download/) or install it via a package manager like `brew install cmake` (macOS) or `sudo apt-get install cmake` (Debian/Ubuntu).

*   **Incorrect CMake Generator on Windows**: By default, Python extensions on Windows are built with the MSVC (Microsoft Visual C++) compiler. If you have only installed the MinGW toolchain, `scikit-build-core` might default to an incompatible generator like "NMake Makefiles", causing errors.
    *   **Solution**: If you must use MinGW, you need to tell CMake to use the correct generator by setting the `CMAKE_GENERATOR` environment variable **before** running `pip install`.

    **For PowerShell:**
    ```powershell
    $env:CMAKE_GENERATOR="MinGW Makefiles"
    pip install .
    ```

## Running the Solomon Benchmark Example

The repository includes an example script to solve instances from the well-known Solomon VRPTW benchmark suite.

**Step 1: Get the source code**
If you haven't already, clone the repository:
```bash
git clone https://github.com/littleQiu22/cg-vrp.git
cd cg-vrp
```

**Step 2: Run the main script**
The script is located in the `examples` directory. You can run it as a module from the root directory of the project. The command format is:

```bash
python -m examples.Solomon.main <instance_name>
```

For a concrete example, to solve the `c101` instance, run:
```bash
python -m examples.Solomon.main c101
```

You can find all available instance names in the `examples/Solomon/data/` directory.

### Potential Runtime Issues on Windows with MinGW

If you built the library from source using **MinGW on Windows**, you might encounter a `ModuleNotFoundError` or a DLL loading error when you try to run the script.

*   **The Problem**: The compiled C++ extension (`.pyd` file) will depend on MinGW's runtime DLLs (like `libwinpthread-1.dll`). For security reasons, Python 3.8 and later no longer search the system `Path` environment variable for DLL dependencies when loading extensions. Therefore, even if you have added your MinGW `bin` directory to `Path`, Python will not find the necessary DLLs.

*   **The Solution**: You must tell the library where to find the MinGW DLLs by setting environment variable `MINGW_DLL_PATH` before running the script.

    **For PowerShell:**
    ```powershell
    # Replace the path with the location of your MinGW bin directory
    $env:MINGW_DLL_PATH = "E:\mingw64\bin"
    python -m examples.Solomon.main c101
    ```    

## How It Works: Column Generation

The Vehicle Routing Problem can be formulated and **linearly relaxed** to a massive linear program where each variable (or "column") represents a possible vehicle route. The number of such routes is large, making it impossible to solve directly.

**Column Generation** is an iterative algorithm that elegantly solves this by decomposing the problem:

1.  **Restricted Master Problem (RMP)**: We start by solving a simplified version of the problem that considers only a small subset of all possible routes. This RMP is a small linear program that can be solved efficiently. The solution provides dual prices for each customer.

2.  **Pricing Subproblem**: The core of the algorithm. We use the dual prices from the RMP to "price out" new routes. The goal is to find a route with a **negative reduced cost**. Such a route, if added to the RMP, has the potential to improve the overall solution. This subproblem is equivalent to a Resource Constrained Shortest Path Problem(RCSPP).

3.  **Iteration**: If the pricing subproblem finds one or more routes with negative reduced costs, we add them to the RMP and solve it again. This loop continues until the pricing subproblem can no longer find any routes with negative reduced cost, at which point we have proven that the solution to our RMP is optimal for the original, full problem.

To obtain a final, practical solution where each route is either used or not used, an additional step is required. There are two main approaches:

### 1. The Exact Method: Branch-and-Price

To find the provably optimal integer solution, the Column Generation process is embedded within a **Branch-and-Bound** search tree. This combined algorithm is known as **Branch-and-Price**.

At each node of the tree, Column Generation is used to solve the LP relaxation. If the solution is fractional, a branching decision is made (e.g., forcing an edge to be either included or excluded from all routes), creating two new subproblems (nodes) in the tree. This process continues until an integer solution is found and its optimality is proven.

Implementing a full Branch-and-Price framework is currently beyond the scope of this library.

### 2. Heuristic Methods

A common and practical approach is to use a heuristic to construct an integer solution from the fractional result. Simple methods include:

*   **Feasibility Pump**: These are techniques that try to "pump" the fractional solution towards an integer feasible one.

*   **Solving a final Integer Program**: A very straightforward heuristic is as follows:
    1.  **Collect Final Routes**: Take all the routes (columns) that have a non-zero value in the final LP solution.
    2.  **Solve a Set Covering Problem**: Formulate a new, smaller Integer Program (like a Set Covering Problem) using only this pool of "good" routes. The goal is to select a subset of these routes that covers all customers, at minimum cost.

Currently, `cg-vrp` focuses on solving the root node of the Branch-and-Price tree, providing the essential components (strong bounds and a pool of promising routes) needed for these final integer programming steps.

### Pricing Algorithms

This library implements two algorithms for the pricing subproblem (RCSPP):

*   **Pulsing Algorithm**:
    *   Based on a **Depth-First Search (DFS)** approach that recursively explores paths.
    *   It is a lightweight and memory-efficient algorithm.
    *   Uses **primal bound pruning**: if a partial path's cost, combined with a heuristic estimate to the end, is already worse than the best-known solution, the path is abandoned early.

*   **Labeling Algorithm (A\* variant)**:
    *   Based on an **A\* search** algorithm, which uses a priority queue to explore the most promising paths first.
    *   In addition to **primal bound pruning**, it introduces **dominance rules**: If two partial paths arrive at the same customer, and one path is better or equal in all resources (e.g., lower cost, earlier arrival time, less demand consumed) than the other, the dominated (worse) path can be safely discarded. This reduces the search space.

## Further Reading & References

To learn more about the algorithms used in the pricing subproblem, please refer to the academic papers:

1.  **Pulsing**: [https://dx.doi.org/10.1287/trsc.2014.0582](https://dx.doi.org/10.1287/trsc.2014.0582)

2.  **Labeling**: [https://doi.org/10.1287/trsc.2018.0878](https://doi.org/10.1287/trsc.2018.0878)