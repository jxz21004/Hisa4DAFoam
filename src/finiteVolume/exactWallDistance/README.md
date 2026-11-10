# Geometric wall distance with shape derivatives

`wallDist { method exactWallDistance; updateInterval 1; }` selects nearest distance to the global triangulated wall. It is opt-in and is compiled into HisaSolver by src/Make/files. Do not also load the standalone plugin when using that integrated build (duplicate runtime registration).

For standalone testing, this directory supports `wmake libso`, with `FOAM_USER_LIBBIN` selecting the output directory. Native and ADR variants must be built against their respective OpenFOAM libraries. The main HISA build follows WM_AD_MODE.

## Geometry and cache

Polygon centre-fan decomposition matches OpenFOAM face::nearestPoint; triangular faces are used directly. All selected wall patches are gathered into replicated geometry and an indexed search tree. Memory scales with global wall faces per rank. There is no per-iteration search on a static primal mesh.

The cache stores integer triangle associations and passive coordinate snapshots only. With identical geometry and query coordinates, a repeated recording reuses associations but always recomputes the projection with current active coordinates. Across actual mesh changes, previous triangles give search upper bounds and a global search verifies the closest surface. This is deliberately stricter than freezing an association for arbitrary deformation; moving geometry still incurs tree reconstruction/search. Topology changes clear the cache.

Typed scalar MPI reduction carries wall coordinate dependencies across ranks. Nested pointField serialization must not replace it: that version passed value checks but failed parallel reverse-mode derivatives.

Wall patch distances retain native SMALL protection; internal geometric distances are not floored. Normals point toward the wall. Search failure or absence of walls is fatal.

## DAFoam coupling

DAInputVolCoord calls mesh.movePoints while the reverse tape is active. wallDist updates yWall, which DASpalartAllmaras references. Set updateInterval=1 so every coordinate update recomputes distances, even when optimization updates occur at the same time index. Do not use meshWaveFrozen for this method: its stored numerical distance loses the shape dependence. Existing forceMeshWaveFrozen only replaces method=meshWave and leaves this method unchanged. No shared SA equation modification is needed.

Nearest-face/edge selection is piecewise differentiable. AD records the current projection branch; derivative checks crossing a change in nearest feature need one-sided interpretation. This is not a claim of globally smooth distance.

## Validation

- L1p5 original geometry: 181152 cells; serial/8-rank maximum difference 1.14e-13 m; 52 near-wall queries agree with brute-force polygons.
- Warped-wall 60-cell test: rigid translation, uniform scaling, and smooth deformation; native and ADR, serial, 2 and 8 ranks, FD steps 1e-4,1e-5,1e-6. The 2-rank partition is wall-normal so one rank has no wall faces, exercising remote surface derivatives.
- Repeated recordings exercise cached associations. Local deformation AD/FD absolute discrepancy about 5e-9 at h=1e-5.

Project test launcher and outputs are under validation/hisa_adflow_comparison/wall_distance_shape and reports/wall_distance_shape_test_build. Full aerodynamic adjoint gradients, topology-changing AD and large moving-mesh performance are not yet certified.
