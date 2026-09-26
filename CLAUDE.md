# Claude Notes

## Meta
- When the user says "take a note" or similar, update this file (CLAUDE.md)
- Only record higher-order ideas and principles
- Never include specific code details or implementation specifics
- Never include historically contingent information (e.g., "X was renamed to Y")
- When adding or modifying style rules in this file, also update `.cursorrules` to mirror the changes
- Maximize use of whitelisted commands to avoid requiring user approval for each command

## Testing and Development

### Expensive and Integration Tests
Tests that are statistical, slow, or integration-focused should use the `EXPENSIVE_` prefix with environment variable gating:
- Prefix expensive test names with `EXPENSIVE_`
- Use `REQUIRE_EXPENSIVE()` macro at the start of the test (defined in `src/geometry_art/testing/geometric_assertions.hpp`)
- All tests are always compiled
- Skipped by default (normal runs skip them)
- Opt-in to run by setting environment variable: `EXPENSIVE=1`
- No recompilation needed

Pattern:
```cpp
TEST(ComponentTest, EXPENSIVE_StatisticalValidation) {
    REQUIRE_EXPENSIVE();

    // Test with 100k+ samples, multiple runs, etc.
}
```

Running tests:
```bash
# Normal run (skips expensive tests automatically)
ctest

# Run expensive tests too
EXPENSIVE=1 ctest
```

What makes a test expensive:
- Large sample counts (10k+ samples)
- Statistical validation requiring multiple runs
- Integration tests involving multiple subsystems
- Performance benchmarks

### Pre-Commit Checklist
Before making a commit, always verify:
1. The code builds successfully: `cmake --build .` (or `make`)
2. All tests pass, including expensive tests: `EXPENSIVE=1 ctest`

These checks ensure that expensive tests are not accidentally broken and remain executable.

### Build System Behavior
- CMakeLists.txt automatically picks up all `*_test.cpp` files using `GLOB_RECURSE`
- New test files are discovered at cmake configuration time, not build time
- After adding a new test file, run `cmake ..` from the build directory to reconfigure
- The build system will then automatically compile and register the new test

### File Organization for Related Types
When a concept has multiple concrete implementations, organize them in a directory named after the concept:
- Directory name matches the concept (e.g., `interval_sampler/`, `circular_interval_sampler/`)
- Concept file and its implementations live together in the directory
- The directory is a sibling to related types (e.g., `math/interval_sampler/` is sibling to `math/interval.hpp`)

When a class has internal helper classes that exist purely to support its implementation:
- Put the main class and its helpers in a directory named after the main class
- The helpers are implementation details, not part of the public interface
- Tests for the helpers can be separate files within the directory

### Dependency Injection via Concepts
**ALL injected dependencies must use concepts** - not just random/slow ones. Any dependency injected via template parameter should be constrained by a concept:
- Use concepts to constrain ALL injected template parameters
- Template parameters must use the concept as the constraint, never `typename`
- Provide concrete default template arguments for production use
- This enables dependency injection while maintaining type safety and clear interfaces
- Benefits: testability (inject test doubles), type safety (compiler enforces interface), documentation (concept shows requirements)

Pattern:
```cpp
template<
    ScalarField ScalarFieldType,
    SpherePointGenerator GeneratorType = RandomSpherePointGenerator<>,
    IntervalSampler IntervalSamplerType = UniformIntervalSampler
>
class MyClass {
 public:
    MyClass() = default;

    MyClass(
        ScalarFieldType field,
        GeneratorType generator,
        IntervalSamplerType sampler = IntervalSamplerType()
    );

 private:
    ScalarFieldType _field;
    GeneratorType _generator;
    IntervalSamplerType _sampler;
};
```

Key principles:
- **Concept as constraint**: Use the concept name (e.g., `IntervalSampler`) not `typename` for ALL injected template parameters
- **Concrete default**: Provide a concrete type as the default (e.g., `UniformIntervalSampler`)
- **Constructor injection**: Allow injecting dependencies via constructor with default-constructed fallback
- **Default constructor**: Provide a default constructor that uses default-constructed dependencies
- **Never create internally**: Never instantiate dependencies as uninjected private members - ALL dependencies must be injectable

### Unit Test Design Principles
Unit tests should be minimal, deterministic, and fast:
- Use 1-3 data points to create contrived situations that exhibit specific behaviors
- Use test doubles (mocks, sequences, seeded samplers) instead of real random dependencies
- Each test should verify a single behavior with a minimal scenario
- Avoid loops and large sample counts in unit tests
- Place pure unit tests at the top of test files

What belongs in unit tests:
- Edge cases with minimal data (empty inputs, single elements, boundary values)
- Specific behaviors with contrived scenarios (e.g., "filters outside points" with 1 inside + 1 outside)
- Error conditions and validation
- Deterministic calculations with known inputs/outputs

### Integration Test Design Principles
Integration tests should use real collaborators and realistic scenarios:
- Mark with `EXPENSIVE_` prefix and `REQUIRE_EXPENSIVE()` macro
- Use real random generators (optionally seeded for reproducibility)
- Use realistic sample counts (10k+ samples)
- Test convergence, statistical properties, and system-level behaviors
- Place integration tests at the bottom of test files

What belongs in integration tests:
- Statistical validation requiring many samples
- Convergence tests for numerical algorithms
- System-level behaviors involving multiple components working together
- Performance benchmarks
- End-to-end workflows with realistic data volumes

### Terminal Testing Flags
The `tessellate` executable supports command-line flags for terminal/AI-assisted testing:
- `--points <N>` / `-p`: Number of sites (default: 10)
- `--density-field <type>` / `-f`: 'constant', 'linear', 'quadratic', 'noise' (piecewise polynomial mesh), 'noise-smooth' (C1 spline projection), 'noise-fit' (global polynomial fit) or 'image' (with `--image <path>`, darker is denser) (default: quadratic); flat domains carry 'constant', 'noise', 'image' and, where the height is walled, 'linear'
- `--lloyd-passes <N>`: Density-weighted Lloyd warm-start passes (default: 5)
- `--max-outer-iterations <N>` / `--max-inner-iterations <N>`: Augmented Lagrangian budget
- `--capacity-tolerance <x>`: Relative RMS capacity error at which optimization stops (default: 1e-7)
- `--density-tolerance <x>`: Refine the density's spline (sphere) or sampling grid (flat) until its relative representation error is below this (0 keeps the cell-scale resolution)
- `--snapshot <path>`: Write the tessellation as `<path>.json` and `<path>.svg`; `--snapshot-interval <s>` also rewrites it every `s` seconds while running
- `--geometry <g>` / `-g`: 'sphere' (default), or 'torus' / 'cylinder' / 'plane' for the flat family; the torus wraps both axes, the cylinder walls its rims, the plane walls all four sides; flat runs use the newton inner solver
- `--scale <k>`: Output units per model unit (default: 50); the solver works in model units, the snapshot records the scale, and the sphere's radius is one model unit
- `--width <w>` / `--height <h>`: The flat domain's rectangle in output units (cylinder circumference and height; default: 100 x 50)
- `--contrast <r>`: Densest-to-sparsest density ratio of the 'linear' field, which rises from bottom to top (every geometry but the torus)
- `--help`: Show usage information

Examples:
- Fast terminal test with constant density: `./tessellate -f constant -p 10`
- Larger run with the quadratic field: `./tessellate -p 200`

The `skeletonize` executable thickens a snapshot's edge graph into a printable solid:
- `skeletonize <snapshot.json>`: Reads the sites and domain back from the snapshot and writes the model beside it
- `--format <f>` / `-f`: 'stl' (default, binary), 'obj', 'ply' or 'off'; `--output <path>` / `-o` names the file instead, its extension picking the format
- `--scale <k>` / `-s`: Output units per model unit (the sphere's radius, or one unit of a flat domain; default: the snapshot's, else 50)
- `--bar-width <w>` / `-w`, `--bar-thickness <t>` / `-t`: The bars along and off the surface, in output units (default: 1.5)
- `--resolution <r>` / `-r`: Longest facet edge along a curved surface, in output units (default: 1.0)
- `--window <w>`: Cut a flat geometry's model down to a centred square this wide, laid flat, in output units; the torus shows its tessellation with no boundary effects at all
- `--frame-width <f>`: The frame around the window, outside it, in output units (default: the bar width)

Visual inspection goes through the web studio -- `./studio`, documented in
`viewer/README.md` -- which runs the executables and renders their snapshots
live. A program the studio offers declares its parameters in one place there;
a new flag reaches the form by being declared, not by being wired up twice.

## CGAL Considerations

### Order-Agnostic Operations
Important geometric computations should produce identical results regardless of equivalent but different orderings of connectivity/combinatorics. For example, polygon area and moment calculations should not depend on edge iteration order. Use exact arithmetic where needed to achieve this.

### Checkpoint/Restore Pattern
When implementing checkpoint/restore for optimization (saving best state to return to later), rebuild the VoronoiSphere from scratch rather than using `update_site()` in a loop:

```cpp
// Preferred: rebuild from scratch
void restore_checkpoint(const std::vector<Point3>& checkpoint) {
    auto new_voronoi = std::make_unique<VoronoiSphere>();
    for (const auto& site : checkpoint) {
        new_voronoi->insert(site);
    }
    _voronoi_sphere = std::move(new_voronoi);
}
```

The rebuild approach is preferred for performance - bulk insertion is faster than updating each site individually.

## Code Style Preferences

### Naming Conventions
- Use full, descriptive variable names - never abbreviations
- Prefer `bounding_box` over `bbox`, `point` over `pt`, `index` over `idx`, etc.
- Parameter names in concept definitions should be descriptive (e.g., `const BoundingBox& bounding_box`)

### Comments
- Do not use verbose documentation comments
- Avoid TODOs in code comments
- Avoid inline explanatory comments
- Keep code clean and minimal - let it speak for itself
- This applies even when writing CGAL-style code

### Control Flow
- Avoid the "V antipattern" - large nested conditional blocks that create a V-shape indentation
- Prefer private helper functions with early returns for exceptional cases
- Keep main functions flat and linear when possible

### Function and Method Organization
- Always organize functions/methods with highest abstraction first, followed by lower abstraction
- In implementation files (.cpp), place main() or the highest-level entry point first
- Use forward declarations when necessary to maintain this ordering (never reorganize to avoid forward declarations)
- Within each section, order by abstraction level (high-level orchestration methods before low-level utility methods)

### Class Member Ordering
- Always put `public:` section first in class definitions
- Public interface is what users care about - show it first
- Private implementation details (methods, type aliases, members) come after
- This applies to all class members: methods, type aliases, using declarations, data members

### Formatting Multi-line Expressions
- Continuation operators (like `<<`) go at the end of the line, not the beginning
- Opening elements (parenthesis, bracket, etc.) stay at end of line
- Arguments/content on new lines, indented once (not aligned with previous line)
- Closing elements (parenthesis, bracket, etc.) on their own line for nested structures
- Never align with elements from previous lines
- Separate consecutive multi-line blocks with a blank line before and after
- This applies to all multi-line expressions including chained stream operations (`<<`)

Examples:
```cpp
// Function call
double length = std::sqrt(
    candidate.x() * candidate.x() +
    candidate.y() * candidate.y() +
    candidate.z() * candidate.z()
);

Point3 normalized(
    candidate.x() / length,
    candidate.y() / length,
    candidate.z() / length
);

// Stream operations - simple case (no nested structure)
std::cout << "Usage: " << program_name << " [options]\n"
    << "Options:\n"
    << "  --points <N>  Number of points\n"
    << "  --help        Show help\n"
    << std::endl;

// Stream operations - with nested structure (opening/closing on own lines)
std::cout << "Centroid: (" <<
    info.centroid.x() << ", " << info.centroid.y() << ", " << info.centroid.z() <<
    ")" << std::endl;
```

### File Endings
- Always end files with a single trailing newline