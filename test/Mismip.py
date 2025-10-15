import numpy as np
from SetIceShelfBC import SetIceShelfBC
from mismipbasalforcings import mismipbasalforcings

# -----------------------------------------------------------------------------
# MISMIP parameterization block — fully commented
# -----------------------------------------------------------------------------
# Assumes an ISSM model object `md` already exists with a mesh (md.mesh.x/y)
# and all standard fields allocated. This file fills in geometry, rheology,
# basal friction, boundary conditions, and basic forcings for a canonical
# MISMIP-style setup.
# -----------------------------------------------------------------------------

# ------------------------------
# Geometry: bed, surface, base, thickness
# ------------------------------
print('      creating thickness')
# Longitudinal bed shape (polynomial in nondimensional x): controls fjord-like slope
bx = (
    -150
    - 728.8 * (md.mesh.x / 300000) ** 2
    + 343.91 * (md.mesh.x / 300000) ** 4
    - 50.57 * (md.mesh.x / 300000) ** 6
)
# Transverse valley walls via two opposing logistic ramps in y (buttress-like)
by = (
    500.0 / (1 + np.exp(-2.0 / 4000.0 * (md.mesh.y - 80000.0 / 2.0 - 24000.0)))
    + 500.0 / (1 + np.exp( 2.0 / 4000.0 * (md.mesh.y - 80000.0 / 2.0 + 24000.0)))
)
# Same transverse feature evaluated at y=0 (centerline reference), used for surface
by0 = (
    500.0 / (1 + np.exp(-2.0 / 4000.0 * (0.0 - 80000.0 / 2.0 - 24000.0)))
    + 500.0 / (1 + np.exp( 2.0 / 4000.0 * (0.0 - 80000.0 / 2.0 + 24000.0)))
)

# Bed elevation (clipped to -720 m) — ensures no unrealistically deep troughs
md.geometry.bed = np.maximum(bx + by, -720)

# Prescribed surface elevation (centerline + 100 m buffer, floor at 10 m)
md.geometry.surface = np.maximum(bx + by0 + 100, 10)

# Ice base cannot go below -90 m (simple flotation limiter for initial state)
md.geometry.base = np.maximum(md.geometry.bed, -90)

# Thickness from surface-base; positive by construction
md.geometry.thickness = md.geometry.surface - md.geometry.base

# ------------------------------
# Basal friction (Weertman/Glen friction parameters)
# ------------------------------
print('      creating drag')
# Scalar drag coefficient C (SI-consistent). Using sqrt here mirrors classic ISSM
# MISMIP parameters (C^2 ~ 3.160e6 Pa m^-1/3 s^1/3 when n=3). Adjust as needed.
md.friction.coefficient = np.sqrt(3.160e6) * np.ones(md.mesh.numberofvertices)

# p, q exponents define sliding law: tau_b = C^2 * N^q * |u|^(p-1) * u
md.friction.p = 3 * np.ones(md.mesh.numberofelements)  # Glen exponent
md.friction.q = 0 * np.ones(md.mesh.numberofelements)  # no effective-pressure dependence

# ------------------------------
# Flow law (ice rheology)
# ------------------------------
print('      creating flow law parameter')
# Convert Arrhenius softness A to B = A^{-1/n}. Here A = 6.338e-25 Pa^-3 s^-1 (n=3).
# B has units of Pa s^(1/3) and is vertex-based in ISSM.
md.materials.rheology_B = (1 / (6.338e-25) ** (1 / 3)) * np.ones(md.mesh.numberofvertices)

# Glen exponent n and law selector
md.materials.rheology_n = 3 * np.ones(md.mesh.numberofelements)
md.materials.rheology_law = 'None'  # use the simple isothermal Glen law

# ------------------------------
# Boundary conditions — ice shelf front & masks
# ------------------------------
print('      boundary conditions for diagnostic model')
# Apply shelf front boundary from a contour (Neumann conditions at calving front)
md = SetIceShelfBC(md, './Front.exp')

# Initialize level-set masks: -1 = inside ice/ocean (no front), +1 = outside
md.mask.ice_levelset[:] = -1
md.mask.ocean_levelset[:] = -1

# Define a vertical calving/front line near x = 640 km (level-set = 0)
pos = np.where((md.mesh.x < 640000.1) & (md.mesh.x > 639999.9))[0]
md.mask.ice_levelset[pos] = 0

# Velocity Dirichlet BCs: NaN means "free" (no SPC), set along domain edges
md.stressbalance.spcvx[:] = np.nan
md.stressbalance.spcvy[:] = np.nan

# No-slip in y at north/south boundaries (y≈0 and y≈80 km)
pos = np.where(
    ((md.mesh.y < 80000.1) & (md.mesh.y > 79999.9)) |
    ((md.mesh.y <     0.1) & (md.mesh.y >    -0.1))
)[0]
md.stressbalance.spcvy[pos] = 0  # clamp along those top/bottom lines

# No-flow at the left wall (x≈0): vx=vy=0 to anchor solution
pos2 = np.where((md.mesh.x < 0.1) & (md.mesh.x > -0.1))[0]
md.stressbalance.spcvx[pos2] = 0
md.stressbalance.spcvy[pos2] = 0

# ------------------------------
# Forcings (melt, SMB, geothermal, etc.)
# ------------------------------
print('      forcing conditions')
# Basal melt & buttressing parameterization container
# NOTE: In Python you *assign* the instance to md.basalforcings. The comment in
# the original snippet ("← remove md") likely referred to MATLAB (where one would
# write `md.basalforcings = mismipbasalforcings();`). Here this line is correct.
md.basalforcings = mismipbasalforcings()

# Basal melt scaling (dimensionless factor applied inside the parameterization)
md.basalforcings.meltrate_factor = 0
# Thickness below which the shelf "thins" rapidly (m)
md.basalforcings.threshold_thickness = 75
# Depth above which melt is applied (m, negative = below sea level)
md.basalforcings.upperdepth_melt = -100

# Surface mass balance (m ice eq. / yr)
md.smb.mass_balance = 0.3 * np.ones(md.mesh.numberofvertices)

# Geothermal heat flux (W m^-2)
md.basalforcings.geothermalflux = 0.5 * np.ones(md.mesh.numberofvertices)

# Grounded-ice basal melt rate (m/yr; often set to 0)
md.basalforcings.groundedice_melting_rate = 0.0 * np.ones(md.mesh.numberofvertices)

# ------------------------------
# Physical constants & model switches
# ------------------------------
md.thermal.spctemperature = np.nan * np.ones(md.mesh.numberofvertices)  # no fixed T BCs
md.groundingline.migration = 'SubelementMigration'  # accurate GL tracking

# Densities (kg/m^3) and gravity (m/s^2); yts = seconds per year
md.materials.rho_ice = 918
md.materials.rho_water = 1028
md.constants.g = 9.8
md.constants.yts = 31556926

# Transient model toggles
md.transient.isthermal = 0        # isothermal
md.transient.isgroundingline = 1  # simulate GL motion
md.stressbalance.isnewton = 0     # disable Newton linearization for now

# ------------------------------
# Initial conditions (simple uniform state)
# ------------------------------
# Start with a small, non-zero speed to avoid exact-stagnation numerical quirks
md.initialization.vx = np.ones(md.mesh.numberofvertices)
md.initialization.vy = np.ones(md.mesh.numberofvertices)
md.initialization.vz = np.ones(md.mesh.numberofvertices)
md.initialization.vel = np.sqrt(2) * np.ones(md.mesh.numberofvertices)

# Hydrostatic overburden as a first-guess pressure (Pa)
md.initialization.pressure = md.constants.g * md.materials.rho_ice * md.geometry.thickness

# Uniform isothermal field (K)
md.initialization.temperature = 273 * np.ones(md.mesh.numberofvertices)
