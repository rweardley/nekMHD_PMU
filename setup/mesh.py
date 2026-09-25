# Script created for Coreform Cubit 2025.8

import cubit
import sys
import numpy as np

sys.path.append("../nekrs_mhd_examples/python")
import boundary_layer_high_order_convert as bl

cubit.cmd("reset")

### User Settings ###

# Set file names
in_stepname = "'PMU_SMALL_MOD_MHD_origin.stp'"
out_meshname = "pmu"

nondimensionalise = True
fluid_only = True

polynomial_order = 9
bl_growth_rate = 1.15
MHD_first_layer_multiplier = 0.5
# MHD_first_layer_multiplier = 0.3
num_qps_in_first_layer = 1

# Set reference length for nondimensionalisation
# Applied after mesh generation
ref_length = 0.231

Ha = 18468.36

surf_external = "9 54 61 62 68"
surf_inlet = 132
surf_outlet = 131
surf_first_wall = 56


### Mesh calculations ###


kolmogorov_length = 1.83097376e-4  # m
kolmogorov_multiplier = 6
y1 = 4.5e-4


### Estimate mesh parameters as for a low-order mesh

# MHD requirements:
MHD_HL = ref_length / Ha  # Hartmann layer
MHD_SL = ref_length / Ha**0.5  # Side layer
MHD_FLT = MHD_first_layer_multiplier * MHD_HL  # First layer thickness
MHD_BL = 5 * MHD_SL  # BL total thickness
MHD_n = np.ceil(np.log(1 - (MHD_BL * (1 - bl_growth_rate)) / MHD_FLT) / np.log(bl_growth_rate))  # num layers
MHD_NLT = MHD_FLT * bl_growth_rate ** (MHD_n - 1)  # Nth layer thickness
MHD_core_sz = MHD_NLT  # core size
print(f"MHD: Hartmann Layer Thickness (HL): {MHD_HL:.6e} m")
print(f"MHD: Side Layer Thickness (SL): {MHD_SL:.6e} m")
print(f"MHD: First Layer Thickness (FLT): {MHD_FLT:.6e} m")
print(f"MHD: Boundary Layer Total Thickness (BL): {MHD_BL:.6e} m")
print(f"MHD: Growth Rate (gr): {bl_growth_rate:.2f}")
print(f"MHD: Number of Layers (n): {int(MHD_n)}")
print(f"MHD: Nth Layer Thickness (NLT): {MHD_NLT:.6e} m")
print(f"MHD: Suggested Core Sizing: {MHD_core_sz:.6e} m")

# Turbulence requirements:
print(f"HD: Kolmogorov length (η): {kolmogorov_length:.3e} m")
HD_res = kolmogorov_length * kolmogorov_multiplier
print(f"HD: Turbulence resolution ({kolmogorov_multiplier}η): {HD_res:.6e} m")
print(f"HD: y1: {y1:.3e} m")

# Compare 1st layer requirements:
FLT_options = {"MHD_FLT": MHD_FLT, "HD_y1": y1, f"{kolmogorov_multiplier}η": HD_res}
options = ", ".join(FLT_options.keys())
print(f"Comparing FLT options (minimum): {options}")
FLT_choice = min(FLT_options, key=FLT_options.get)
FLT = FLT_options[FLT_choice]
print(f"Selecting FLT = {FLT_choice} = {FLT} m")

# Compare core resolution requirements:
core_sz_options = {"MHD_core_sz": MHD_core_sz, f"{kolmogorov_multiplier}η": HD_res}
options = ", ".join(core_sz_options.keys())
print(f"Comparing core size options (minimum): {options}")
core_sz_choice = min(core_sz_options, key=core_sz_options.get)
core_sz = core_sz_options[core_sz_choice]
print(f"Selecting core_sz = {core_sz_choice} = {core_sz} m")

# calculate settings for geometric boundary layer
# estimate n from current growth rate, FLT and final layer thickness (xy_core)
num_layers = (np.log(core_sz / FLT)) / np.log(bl_growth_rate) + 1
print(f"BL: a_1st={FLT:.3e} m, a_final={core_sz:.3e} m, r={bl_growth_rate:.3f}, n={num_layers:.3f}")
# estimate BL thickness
BL = bl.geometric_series_total(FLT, bl_growth_rate, num_layers)
print(f"BL: Total thickness L = {BL:.3e} m = {BL/ref_length:.3e} * ref_length")

# Compare BL thickness requirements:
BL_L_options = {"MHD_BL": MHD_BL, "BL_L": BL}
options = ", ".join(BL_L_options.keys())
print(f"Comparing BL total thickness options (maximum): {options}")
BL_L_choice = max(BL_L_options, key=BL_L_options.get)
BL_L = BL_L_options[BL_L_choice]
print(f"Selecting BL total thickness = {BL_L_choice} = {BL_L} m")

# Calculate growth rate required for given FLT and xy_core
print(f"Recalculating for new BL total thickness")
bl_growth_rate = (FLT - BL_L) / (core_sz - BL_L)
print(f"Recalculated growth rate: {bl_growth_rate}")
num_layers = bl.geometric_series_nterms(FLT, bl_growth_rate, BL_L, initial_guess=num_layers)[0]
print(f"Recalculated number of layers: {num_layers}")
num_layers = int(np.floor(num_layers))
print(f"Round number of layers down to nearest int: {num_layers}")
print(f"BL: a_1st={FLT:.3e} m, r={bl_growth_rate:.3f}, n={num_layers}, total={BL_L:.3e} m")
print(f"Target final layer thickness: {core_sz:.3e}")
final_layer_thickness = FLT * bl_growth_rate ** (num_layers - 1)
print(f"Calculated final layer thickness: {final_layer_thickness:.3e}")

lo_delta_core = core_sz

lo_fluid_first_row = FLT
lo_fluid_growth_factor = bl_growth_rate
lo_fluid_num_layers = num_layers

# Print low-order mesh parameters
print("\n--- Low-Order Mesh Parameters ---")
print(f"Core mesh size (lo_delta_core): {lo_delta_core:.6e} m")
print(f"Fluid first row thickness (lo_fluid_first_row): {lo_fluid_first_row:.6e} m")
print(f"Fluid growth factor (lo_fluid_growth_factor): {lo_fluid_growth_factor:.2f}")
print(f"Fluid number of layers (lo_fluid_num_layers): {int(lo_fluid_num_layers)}")

# Convert mesh parameters to the coarser settings required
# for a comparable high-order mesh

fluid_hi_solution = bl.boundary_layer_low_to_high_order(lo_fluid_first_row, lo_fluid_growth_factor, lo_fluid_num_layers, polynomial_order, num_qps_in_first_layer)
# Print low-order mesh parameters
print(f"\n--- Converting Mesh Parameters to High-Order (N={polynomial_order})---")
fluid_first_row, fluid_growth_factor, fluid_num_layers = fluid_hi_solution[0]
fluid_lo_total_thickness, fluid_hi_total_thickness = fluid_hi_solution[1]

delta_core = lo_delta_core * polynomial_order

print(f"Core mesh size (delta_core): {delta_core:.6e} m")
print(f"Fluid first row thickness (fluid_first_row): {fluid_first_row:.6e} m")
print(f"Fluid growth factor (fluid_growth_factor): {fluid_growth_factor:.2f}")
print(f"Fluid number of layers (fluid_num_layers): {int(fluid_num_layers)}")

# Set up mesh specs based on the above calculations

bl_y1_size = fluid_first_row  # wedges don't get split in wall normal direction
bl_fluid_growth = fluid_growth_factor
bl_fluid_layers = 3
bl_solid_growth = fluid_growth_factor
bl_solid_layers = 2
wall_size = fluid_first_row * 70
exterior_size = wall_size * 7
surf_layers = 2  # not quite sure what this one does
bulk_size = delta_core * 2  # edges of tets get split in half
bulk_gradient = 1.1
bulk_min_num_layers_3d = 3
bulk_min_num_layers_2d = 2
bulk_min_num_layers_1d = 1
trimesher_surf_gradation = 1.3
trimesher_vol_gradation = 1.3
tetmesher_growth_factor = 1.3

tetmesh_optimise = True

### Load and Prepare Geometry ###

cubit.cmd(f"import step {in_stepname} noheal")

cubit.cmd("imprint volume all")
cubit.cmd("merge volume all")

cubit.cmd('create group "vol_solid"')
cubit.cmd("vol_solid add volume 1")

cubit.cmd('create group "vol_fluid"')
cubit.cmd("vol_fluid add volume 2")

cubit.cmd('create group "boundary_surfs"')
cubit.cmd("boundary_surfs add surface all")
cubit.cmd(f"boundary_surfs remove surface {surf_inlet}")
cubit.cmd(f"boundary_surfs remove surface {surf_outlet}")
cubit.cmd(f"boundary_surfs remove surface {surf_external}")
cubit.cmd(f"boundary_surfs remove surface {surf_first_wall}")

cubit.cmd(f'group "inlet" add surf {surf_inlet}')
cubit.cmd(f'group "outlet" add surf {surf_outlet}')

cubit.cmd('create group "fs_interface"')
cubit.cmd("fs_interface add surface all")
cubit.cmd(f"fs_interface remove surface {surf_inlet}")
cubit.cmd(f"fs_interface remove surface {surf_outlet}")
cubit.cmd(f"fs_interface remove surface {surf_external}")
cubit.cmd(f"fs_interface remove surface {surf_first_wall}")

cubit.cmd(f'group "exterior" add surf {surf_external}')
cubit.cmd(f'group "first_wall" add surf {surf_first_wall}')

# add surfaces to sidesets

cubit.cmd("sideset 1 add surface in inlet")
cubit.cmd('sideset 1 name "inlet"')
cubit.cmd("sideset 2 add surface in outlet")
cubit.cmd('sideset 2 name "outlet"')
cubit.cmd("sideset 3 add surface in fs_interface")
cubit.cmd('sideset 3 name "fs_interface"')
cubit.cmd("sideset 4 add surface in exterior")
cubit.cmd('sideset 4 name "exterior"')
cubit.cmd("sideset 5 add surface in first_wall")
cubit.cmd('sideset 5 name "first_wall"')

# create element blocks
cubit.cmd("block 1 add volume in vol_fluid")
cubit.cmd("block 2 add volume in vol_solid")

# Force nodes to follow curved surfaces (walls)
cubit.cmd("set node constraint on")

### Set mesh sizing

# set surface mesh size

cubit.cmd(f"surface in boundary_surfs size {wall_size}")

# add boundary layers

cubit.cmd("create boundary_layer 1")
cubit.cmd(f"modify boundary_layer 1 uniform height {bl_y1_size} growth {bl_fluid_growth} layers {bl_fluid_layers}")
cubit.cmd(f"modify boundary_layer 1 add surface in boundary_surfs volume in vol_fluid")
cubit.cmd("modify boundary_layer 1 continuity on")

cubit.cmd("create boundary_layer 2")
cubit.cmd(f"modify boundary_layer 2 uniform height {bl_y1_size} growth {bl_solid_growth} layers {bl_solid_layers}")
cubit.cmd(f"modify boundary_layer 2 add surface in boundary_surfs volume in vol_solid")
cubit.cmd("modify boundary_layer 2 continuity on")

# Mesh the boundary surfaces

cubit.cmd("surf all scheme trimesh")

# # these don't seem to do anything?
# cubit.cmd("set trimesher coarse off")
# cubit.cmd(f"set trimesher surface gradation {trimesher_surf_gradation}")
# cubit.cmd(f"set trimesher volume gradation {trimesher_vol_gradation}")

cubit.cmd("mesh surface in boundary_surfs")

# test: coarsely mesh exterior surfaces
cubit.cmd(f"surface {surf_external} {surf_first_wall} size {exterior_size}")
cubit.cmd(f"mesh surface {surf_external} {surf_first_wall}")

# tetmesher settings

if tetmesh_optimise:
    cubit.cmd("set tetmesher optimize overconstrained tetrahedra on")
    cubit.cmd("set tetmesher optimize overconstrained edges on")
    cubit.cmd("set tetmesher optimize sliver on")
    cubit.cmd("tetmesher optimize level 6")

### Mesh the volumes ###

# fluid setup

cubit.cmd(f"volume in vol_fluid sizing function type skeleton min_size auto max_size {bulk_size} max_gradient {bulk_gradient} min_num_layers_3d {bulk_min_num_layers_3d} min_num_layers_2d {bulk_min_num_layers_2d} min_num_layers_1d {bulk_min_num_layers_1d}")
cubit.cmd(f"volume in vol_fluid sizing function type skeleton add size_source surface in boundary_surfs size {wall_size} num_layers {surf_layers}")

cubit.cmd("volume in vol_fluid scheme tetmesh")
cubit.cmd(f"volume in vol_fluid tetmesh growth_factor {tetmesher_growth_factor}")

# solid setup

cubit.cmd(f"volume in vol_solid sizing function type skeleton min_size auto max_size {bulk_size} max_gradient {bulk_gradient} min_num_layers_3d {bulk_min_num_layers_3d} min_num_layers_2d {bulk_min_num_layers_2d} min_num_layers_1d {bulk_min_num_layers_1d}")
cubit.cmd(f"volume in vol_solid sizing function type skeleton add size_source surface in boundary_surfs size {wall_size} num_layers {surf_layers}")

cubit.cmd("volume in vol_solid scheme tetmesh")
cubit.cmd(f"volume in vol_solid tetmesh growth_factor {tetmesher_growth_factor}")

# mesh solid then fluid

cubit.cmd("mesh volume in vol_solid")
cubit.cmd("mesh volume in vol_fluid")

if fluid_only: cubit.cmd("del volume in vol_solid") # mesh it but delete it to ensure fluid mesh is consistent

### Nondimensionalise geometry and mesh length scale ###

if nondimensionalise:
    cubit.cmd(f"volume all scale {1.0/ref_length}")

### Save mesh (exodus) ###

cubit.cmd("set exodus netcdf4 off")
cubit.cmd("set large exodus file on")
cubit.cmd(f'export mesh "{out_meshname}_fluid.exo" block 1 overwrite')
if not fluid_only:
    cubit.cmd(f'export mesh "{out_meshname}_solid.exo" block 2 overwrite')

### Calculate number of hex elements this will create in .re2 mesh ###

num_tets = len(cubit.parse_cubit_list("tet", "all"))
num_wedges = len(cubit.parse_cubit_list("wedge", "all"))

print(f"Total number of tets:  \t{num_tets:.3e}")
print(f"Total number of wedges:\t{num_wedges:.3e}")

num_re2_hexes = 3 * num_wedges + 4 * num_tets
num_QPs = num_re2_hexes * polynomial_order**3

print(f"Expected num .re2 hexes:\t{num_re2_hexes:.3e}")
print(f"Expected NekRS QPs:     \t{num_QPs:.3e}")
print(f"30M DOFs/GPU requires {np.ceil(num_QPs/(30e6))} GPUs")
print(f"8 GPUs per node requires {np.ceil(num_QPs/(8*30e6))} nodes")
