import sys
import numpy as np

sys.path.append(r"../nekrs_mhd_examples/python")
import nekrs_mhd_calc as calc


def kolmogorov_estimate(Dh, Re):
    eta = ((0.07) ** (1 / 4)) * Dh * (Re ** (-3 / 4))
    return eta


def y1_estimate(y1p, Dh, Re, Cf):
    y1 = (y1p * Dh) / (Re * np.sqrt(0.5 * Cf))
    return y1


def Cf_prandtl(Re):
    Cf = 0.074 * (Re ** (-1 / 5))
    return Cf


def Cf_schlichting(Re):
    Cf = (2 * np.log10(Re) - 0.65) ** (-2.3)
    return Cf


def Cf_schultz_grunov(Re):
    Cf = 0.370 * (np.log10(Re)) ** (-2.584)
    return Cf


PbLi_density = 9.814e3  # [kg/m3]
PbLi_dyn_viscosity = 1.981e-3  # [kg/m.s]
PbLi_elec_conductivity = 8.769e5  # [S/m]
SS316_elec_conductivity = 1.057e6  # [S/m]
permeability = 4e-7 * np.pi

a_BZ = 0.231  # breeder zone half-width in B direction; use for Ha and Re
d = 1.2  # approximate flow path length between inlet and outlet [m]

# Inlet flow settings

a_in = 0.0750 / 2  # inlet half-width (in B-field direction)
b_in = 0.0615 / 2  # inlet half-height (in B-field direction)
Mdot = np.array([0.56, 1, 2, 5.6])  # mass flow rate [kg/s]
Q = Mdot / PbLi_density  # volumetric flow rate [m^3/h]
A = (2 * a_in) * (2 * b_in)  # inlet area
U = Q / A  # m/s
Bfield_mag = 3.8  # T
Bfield_orientation = (0, 0, 1)

print("Inlet:")
print(f"Mass flow rate = {Mdot} kg/s")
print(f"Volumetric flow rate = {Q} m3/s")
print(f"Average velocity = {U} m/s\n")

case = calc.NekRescaleMHD(
    U,
    Bfield_mag,
    Bfield_orientation,
    a_BZ,
    PbLi_density,
    PbLi_dyn_viscosity,
    permeability,
    PbLi_elec_conductivity,
    ref_conductivity_solid=SS316_elec_conductivity,
    ref_length_is_mhd=True,
    ref_viscosity_is_dynamic=True,
    inductionless_NJxB=True,
    # missing some arguments here, need to set it up for a case with finite wall conductivity; this script is based on smalllab with insulating walls
)
case.print_scaling_summary()
print(f"Flow-through time (dimensional): {d/U} s")
print(f"Flow-through time (nondimensional): {(d/U)/case.scale_time} [ND]\n")
# estimate scales to be resolved
print(f"Hartmann layer thickness: {a_BZ/case.Ha:.3e}")
eta = kolmogorov_estimate(a_BZ * 2, case.Re)
print(f"Kolmogorov scale: {eta} m")
y1_prandtl = y1_estimate(1, a_BZ * 2, case.Re, Cf_prandtl(case.Re))
y1_schlichting = y1_estimate(1, a_BZ * 2, case.Re, Cf_schlichting(case.Re))
y1_sg = y1_estimate(1, a_BZ * 2, case.Re, Cf_schultz_grunov(case.Re))
print(f"Turbulence y1 (Prandtl, y1+ = 1): {y1_prandtl} m")
print(f"Turbulence y1 (Schlichting, y1+ = 1): {y1_schlichting} m")
print(f"Turbulence y1 (Schultz-Grunov, y1+ = 1): {y1_sg} m\n")

case.par_inputs()
