# Mesh of the 3D model

`test_05_main.m` looks here for the file named in `cfg.mesh_file`. Two formats
are accepted, told apart by the extension.

## `.mat`

Must contain two variables:

| variable   | size              | content                                       |
|------------|-------------------|-----------------------------------------------|
| `nodes`    | `[n_nodes x 3]`   | nodal coordinates, **in micrometres**         |
| `elements` | `[n_elem x n_pn]` | connectivity, 1-based indices into `nodes`     |

`FeStructure` also accepts the alternative names `Nodes`/`coord`/`Coord`/`XYZ`
and `Elements`/`elem`/`Elem`/`conn`, and casts both to `double`: the TDK meshes
store the connectivity as `int64`.

The number of columns of `elements` must match `cfg.element_type`: WED15 → 15,
HEX8 → 8, HEX20 → 20, TET4 → 4, TET10 → 10. **The node ordering inside an
element must be the one YaFEc expects**: if the frequencies come out badly
wrong while the mesh plots correctly, this is the first thing to suspect.

## Units

The TDK mesh is **not in metres**. `Accelerometer_3D.m` works in the consistent
system µm / MPa / kg, and `test_05_main.m` is configured to match:

| quantity     | unit      | value in the config |
|--------------|-----------|---------------------|
| length       | µm        | node coordinates    |
| Young        | MPa       | `cfg.E = 168e3`     |
| density      | kg/µm³    | `cfg.rho = 2.33e-15`|
| gravity      | µm/s²     | `cfg.g_value = 9.81e6` |

In that system the assembled `K` comes out in N/m and `M` in kg, so
**frequencies are in Hz and times in seconds**, while displacements are in µm
and forces in µN. This was verified by building the same structure twice, once
in SI and once in this system: the frequencies agree to 1e-13.

## `.inp`

Abaqus file, read through `abqmesh` plus the node-set parser of `FeStructure`,
which also expands the `*Nset, ..., generate` form.

## Choosing the contact and anchor nodes

When the mesh carries no node sets, which is the case for a `.mat`, they are
declared geometrically in `cfg.node_sets`, and `describe_node_sets()` prints
the count and the extent of each one. That printout is the check that caught a
node-set bug on the 2D model at the first attempt. Two things to watch:

- **anchors**: only the nodes actually bonded to the substrate. On the dummy
  model these were two pads of 3 nodes each, not a surface: a single box
  containing both caught 846 nodes instead of 6. Disjoint regions need the
  union of several boxes (a struct array), not one wide box.
- **contact nodes**: the nodes of the surface that can hit the stopper. Worth
  checking with `plot_node_sets()` that the set is a surface and not a
  thickness, and that the nodes are coplanar.
