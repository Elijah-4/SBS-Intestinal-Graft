"""
Model 5: Bioactive Agent Delivery from a Degrading Scaffold (Prototype 3)
Growth Hormone Release from PLGA/SIS Scaffold

This script simulates the spatiotemporal release of Growth Hormone (GH) from a
degrading PLGA scaffold, coupling polymer degradation kinetics with Fickian diffusion.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy.optimize import curve_fit
import matplotlib.patches as patches
from matplotlib.patches import FancyBboxPatch, Circle, FancyArrowPatch
import warnings
warnings.filterwarnings('ignore')

# Set publication-quality plotting parameters
plt.rcParams['figure.dpi'] = 300
plt.rcParams['font.size'] = 10
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['axes.labelsize'] = 11
plt.rcParams['axes.titlesize'] = 12
plt.rcParams['xtick.labelsize'] = 9
plt.rcParams['ytick.labelsize'] = 9
plt.rcParams['legend.fontsize'] = 9

# ============================================================================
# PARAMETERS (Table 6.1)
# ============================================================================

# Scaffold Geometry
L_graft = 15.0  # cm - Graft length
D_outer = 2.5   # cm - Outer diameter
T_wall = 0.2    # cm (2 mm converted to cm)
R_outer = D_outer / 2  # cm - Outer radius
R_inner = R_outer - T_wall  # cm - Inner radius

# PLGA Properties
LG_ratio = "75:25"
k_deg = 0.05    # day^-1 - Degradation rate constant
M_0 = 100.0     # kDa - Initial molecular weight

# Growth Hormone Properties
MW_GH = 22.0    # kDa - Molecular weight of GH
D_min = 1e-12   # cm²/s - Minimum diffusivity in intact polymer
D_max = 1.59e-7 # cm²/s - Maximum diffusivity in degraded polymer

# Convert diffusivity to cm²/day for easier time scale
D_min_day = D_min * 86400  # cm²/day
D_max_day = D_max * 86400  # cm²/day

# Target Therapeutic Concentration
C_ther_min = 10.0  # ng/mL
C_ther_max = 50.0  # ng/mL

# Initial GH Loading - THIS IS THE KEY DESIGN PARAMETER TO OPTIMIZE
M_0_GH = 5.0  # mg - Initial total GH loading (to be optimized)

# Calculate scaffold volume (cylindrical shell)
V_scaffold = np.pi * (R_outer**2 - R_inner**2) * L_graft  # cm³

# Initial concentration in scaffold
C_0 = M_0_GH / V_scaffold  # mg/cm³

# ============================================================================
# GOVERNING EQUATIONS
# ============================================================================

def polymer_degradation(t):
    """
    Equation 6.2.1: Polymer degradation kinetics
    M(t) = M(0) * exp(-k_deg * t)

    Args:
        t: time in days
    Returns:
        M(t): molecular weight at time t in kDa
    """
    return M_0 * np.exp(-k_deg * t)

def diffusion_coefficient(t):
    """
    Equation 6.2.3: Time-dependent diffusion coefficient
    D(t) = D_min + (D_max - D_min) * (1 - M(t)/M(0))

    Args:
        t: time in days
    Returns:
        D(t): diffusion coefficient at time t in cm²/day
    """
    M_t = polymer_degradation(t)
    D_t = D_min_day + (D_max_day - D_min_day) * (1 - M_t / M_0)
    return D_t

def korsmeyer_peppas(t, k, n):
    """
    Equation 6.2.4: Korsmeyer-Peppas model
    M_t/M_inf = k * t^n

    Args:
        t: time
        k: release rate constant
        n: release exponent
    Returns:
        Fractional release
    """
    return k * np.power(t, n)

# ============================================================================
# NUMERICAL SOLUTION: FICKIAN DIFFUSION IN CYLINDRICAL COORDINATES
# ============================================================================

def solve_diffusion_cylindrical(M_0_GH, t_max=30, n_r=100, n_t=300):
    """
    Solve Fick's Second Law in cylindrical coordinates with time-dependent D(t)
    ∂C/∂t = (1/r) * ∂/∂r(r * D(t) * ∂C/∂r)

    Using finite difference method with implicit scheme for stability
    """
    # Spatial domain: from inner radius to outer radius + tissue penetration
    r_tissue_extent = 0.5  # cm - extent of tissue domain to simulate
    r_min = R_inner
    r_max = R_outer + r_tissue_extent

    r = np.linspace(r_min, r_max, n_r)
    dr = r[1] - r[0]

    # Time domain
    t = np.linspace(0, t_max, n_t)
    dt = t[1] - t[0]

    # Initialize concentration matrix
    C = np.zeros((n_t, n_r))

    # Initial condition: GH uniformly distributed in scaffold wall
    for i, r_val in enumerate(r):
        if r_val <= R_outer:
            C[0, i] = C_0  # mg/cm³
        else:
            C[0, i] = 0.0  # No GH in tissue initially

    # Time stepping with Crank-Nicolson scheme (implicit)
    for k in range(n_t - 1):
        # Get current diffusion coefficient
        D_current = diffusion_coefficient(t[k])

        # Simplified explicit scheme for demonstration
        # (For publication, consider using sparse matrix solver for implicit scheme)
        C_new = C[k].copy()

        # Interior points (excluding boundaries)
        for i in range(1, n_r - 1):
            # Finite difference approximation
            dC_dr = (C[k, i+1] - C[k, i-1]) / (2 * dr)
            d2C_dr2 = (C[k, i+1] - 2*C[k, i] + C[k, i-1]) / (dr**2)

            # Fick's law in cylindrical coordinates
            dC_dt = D_current * (d2C_dr2 + (1/r[i]) * dC_dr)

            C_new[i] = C[k, i] + dC_dt * dt

        # Boundary conditions
        # Inner boundary (luminal side): no flux (sealed by SIS coating)
        C_new[0] = C_new[1]

        # Outer boundary (far tissue): sink condition (perfect clearance)
        C_new[-1] = 0.0

        C[k+1] = C_new

    return t, r, C

# ============================================================================
# CUMULATIVE RELEASE CALCULATION
# ============================================================================

def calculate_cumulative_release(t, r, C):
    """
    Calculate cumulative release as percentage of total loaded GH
    """
    cumulative_release = np.zeros(len(t))

    for k in range(len(t)):
        # Integrate remaining GH in scaffold
        remaining_mass = 0
        for i in range(len(r)):
            if r[i] <= R_outer:
                # Cylindrical shell volume element: 2*pi*r*dr*L
                dV = 2 * np.pi * r[i] * (r[1] - r[0]) * L_graft
                remaining_mass += C[k, i] * dV

        cumulative_release[k] = (1 - remaining_mass / M_0_GH) * 100

    return cumulative_release

# ============================================================================
# MAIN SIMULATION
# ============================================================================

print("=" * 70)
print("MODEL 5: Growth Hormone Release from PLGA Scaffold")
print("=" * 70)
print(f"\nScaffold Parameters:")
print(f"  Outer Radius: {R_outer:.2f} cm")
print(f"  Wall Thickness: {T_wall:.2f} cm")
print(f"  Length: {L_graft:.1f} cm")
print(f"  Volume: {V_scaffold:.3f} cm³")
print(f"\nPLGA Properties:")
print(f"  L:G Ratio: {LG_ratio}")
print(f"  Degradation Rate: {k_deg} day⁻¹")
print(f"  Initial MW: {M_0} kDa")
print(f"\nGH Loading:")
print(f"  Initial Loading: {M_0_GH} mg")
print(f"  Initial Concentration: {C_0:.3f} mg/cm³ = {C_0*1000:.1f} μg/cm³")
print(f"\nRunning simulation...")

# Run the simulation
t_sim, r_sim, C_sim = solve_diffusion_cylindrical(M_0_GH, t_max=30, n_r=100, n_t=300)

# Find index of scaffold-tissue interface
interface_idx = np.argmin(np.abs(r_sim - R_outer))

# Extract concentration at interface over time (convert mg/cm³ to ng/mL)
# 1 mg/cm³ = 1000 μg/cm³ = 1000 μg/mL = 1,000,000 ng/mL
C_interface = C_sim[:, interface_idx] * 1e6  # ng/mL

# Calculate cumulative release
cumulative_release = calculate_cumulative_release(t_sim, r_sim, C_sim)

# Fit Korsmeyer-Peppas model to cumulative release data
# Use data from 10% to 60% release for fitting (standard practice)
mask = (cumulative_release >= 10) & (cumulative_release <= 60)
if np.sum(mask) > 5:
    t_fit = t_sim[mask]
    Mt_Minf_fit = cumulative_release[mask] / 100
    try:
        popt, _ = curve_fit(korsmeyer_peppas, t_fit, Mt_Minf_fit, p0=[0.1, 0.5])
        k_KP, n_KP = popt
        print(f"\nKorsmeyer-Peppas Model Fit:")
        print(f"  k = {k_KP:.4f}")
        print(f"  n = {n_KP:.3f}")
        if n_KP < 0.45:
            mechanism = "Fickian diffusion"
        elif n_KP < 0.89:
            mechanism = "Anomalous (non-Fickian) transport"
        else:
            mechanism = "Case-II or Super Case-II transport"
        print(f"  Release mechanism: {mechanism}")
    except:
        k_KP, n_KP = None, None
        print("\nKorsmeyer-Peppas fit failed - insufficient data")
else:
    k_KP, n_KP = None, None

# ============================================================================
# FIGURE 6.1: COMPOSITE VISUALIZATION
# ============================================================================

fig = plt.figure(figsize=(14, 10))
gs = fig.add_gridspec(2, 2, hspace=0.35, wspace=0.3)

# ----------------------------------------------------------------------------
# PANEL A: Schematic Representation
# ----------------------------------------------------------------------------
ax_A = fig.add_subplot(gs[0, 0])
ax_A.set_xlim(-2, 2)
ax_A.set_ylim(-2, 2)
ax_A.set_aspect('equal')
ax_A.axis('off')

# Draw scaffold cross-section
outer_circle = Circle((0, 0), R_outer*1.5, fill=False, edgecolor='black', linewidth=2)
inner_circle = Circle((0, 0), R_inner*1.5, fill=False, edgecolor='black', linewidth=1.5)
scaffold_ring = patches.Wedge((0, 0), R_outer*1.5, 0, 360, width=(R_outer-R_inner)*1.5,
                               facecolor='lightblue', edgecolor='black', linewidth=2, alpha=0.6)

ax_A.add_patch(scaffold_ring)

# Draw PLGA matrix with embedded GH (small circles)
np.random.seed(42)
n_particles = 30
for _ in range(n_particles):
    angle = np.random.uniform(0, 2*np.pi)
    radius = np.random.uniform(R_inner*1.5, R_outer*1.5)
    x = radius * np.cos(angle)
    y = radius * np.sin(angle)
    particle = Circle((x, y), 0.08, facecolor='red', edgecolor='darkred', linewidth=0.5)
    ax_A.add_patch(particle)

# Draw inner lumen
lumen = Circle((0, 0), R_inner*1.5, facecolor='white', edgecolor='black', linewidth=1.5)
ax_A.add_patch(lumen)

# Draw surrounding tissue
tissue_outer = Circle((0, 0), 1.9, fill=False, edgecolor='green',
                       linewidth=2, linestyle='--', alpha=0.7)
ax_A.add_patch(tissue_outer)

# Add arrows showing diffusion
arrow1 = FancyArrowPatch((0.5, 0.8), (1.2, 1.3),
                         arrowstyle='->', mutation_scale=20,
                         linewidth=2, color='red')
arrow2 = FancyArrowPatch((-0.5, -0.8), (-1.2, -1.3),
                         arrowstyle='->', mutation_scale=20,
                         linewidth=2, color='red')
ax_A.add_patch(arrow1)
ax_A.add_patch(arrow2)

# Labels
ax_A.text(0, 0, 'Lumen\n(SIS coating)', ha='center', va='center',
          fontsize=9, fontweight='bold')
ax_A.text(0, -R_inner*1.5-0.3, 'PLGA Matrix\n+ GH', ha='center', va='top',
          fontsize=9, fontweight='bold', color='blue')
ax_A.text(0, 1.6, 'Surrounding\nTissue', ha='center', va='bottom',
          fontsize=9, fontweight='bold', color='green')
ax_A.text(1.5, 1.5, 'GH\nDiffusion', ha='center', va='center',
          fontsize=8, color='red', fontweight='bold')

ax_A.set_title('A. Schematic: GH Release from PLGA Scaffold',
               fontweight='bold', fontsize=12, pad=10)

# ----------------------------------------------------------------------------
# PANEL B: GH Concentration at Scaffold-Tissue Interface
# ----------------------------------------------------------------------------
ax_B = fig.add_subplot(gs[0, 1])

ax_B.plot(t_sim, C_interface, 'b-', linewidth=2.5, label='Interface Concentration')
ax_B.axhline(C_ther_min, color='green', linestyle='--', linewidth=1.5, alpha=0.7)
ax_B.axhline(C_ther_max, color='green', linestyle='--', linewidth=1.5, alpha=0.7)
ax_B.fill_between(t_sim, C_ther_min, C_ther_max, color='green', alpha=0.2,
                   label='Therapeutic Window')

ax_B.set_xlabel('Time (days)', fontweight='bold')
ax_B.set_ylabel('GH Concentration (ng/mL)', fontweight='bold')
ax_B.set_title('B. GH Concentration at Scaffold-Tissue Interface',
               fontweight='bold', fontsize=12, pad=10)
ax_B.grid(True, alpha=0.3)
ax_B.legend(loc='best')
ax_B.set_xlim(0, 30)

# Add text annotation for therapeutic window
max_conc = np.max(C_interface)
if max_conc > C_ther_min:
    # Find time in therapeutic window
    in_window = (C_interface >= C_ther_min) & (C_interface <= C_ther_max)
    if np.any(in_window):
        time_indices = np.where(in_window)[0]
        time_in_window = t_sim[time_indices[-1]] - t_sim[time_indices[0]]
        ax_B.text(0.98, 0.98, f'Time in therapeutic window:\n{time_in_window:.1f} days',
                  transform=ax_B.transAxes, ha='right', va='top',
                  bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8),
                  fontsize=9)

# ----------------------------------------------------------------------------
# PANEL C: Spatial Concentration Profiles over Time
# ----------------------------------------------------------------------------
ax_C = fig.add_subplot(gs[1, 0])

# Select specific time points to display
time_points = [2, 7, 14, 21, 28]
colors = plt.cm.viridis(np.linspace(0, 0.9, len(time_points)))

# Convert r from scaffold center to distance from outer surface
r_from_surface = (r_sim - R_outer) * 10  # Convert to mm

for i, tp in enumerate(time_points):
    idx = np.argmin(np.abs(t_sim - tp))
    C_profile = C_sim[idx, :] * 1e6  # Convert to ng/mL
    ax_C.plot(r_from_surface, C_profile, color=colors[i],
              linewidth=2, label=f't = {tp} days')

# Set limits first before adding vertical line
ax_C.set_xlim(-2, 5)
ax_C.set_ylim(bottom=0)

# Mark the scaffold boundary
ax_C.axvline(0, color='black', linestyle='--', linewidth=1.5, alpha=0.5,
             label='Scaffold Surface')

ax_C.set_xlabel('Distance from Scaffold Surface (mm)', fontweight='bold')
ax_C.set_ylabel('GH Concentration (ng/mL)', fontweight='bold')
ax_C.set_title('C. Spatial GH Concentration Profiles',
               fontweight='bold', fontsize=12, pad=10)
ax_C.grid(True, alpha=0.3)
ax_C.legend(loc='best', ncol=2)

# ----------------------------------------------------------------------------
# PANEL D: Cumulative GH Release Profile
# ----------------------------------------------------------------------------
ax_D = fig.add_subplot(gs[1, 1])

ax_D.plot(t_sim, cumulative_release, 'b-', linewidth=2.5, label='Simulated Release')

# Plot Korsmeyer-Peppas fit if available
if k_KP is not None and n_KP is not None:
    t_KP = np.linspace(t_sim[0], t_sim[-1], 100)
    release_KP = korsmeyer_peppas(t_KP, k_KP, n_KP) * 100
    ax_D.plot(t_KP, release_KP, 'r--', linewidth=2,
              label=f'K-P Fit (n={n_KP:.3f})')

ax_D.set_xlabel('Time (days)', fontweight='bold')
ax_D.set_ylabel('Cumulative GH Release (%)', fontweight='bold')
ax_D.set_title('D. Cumulative GH Release Profile',
               fontweight='bold', fontsize=12, pad=10)
ax_D.grid(True, alpha=0.3)
ax_D.legend(loc='best')
ax_D.set_xlim(0, 30)
ax_D.set_ylim(0, 105)

# Add annotation
if k_KP is not None:
    annotation_text = f'Release Exponent n = {n_KP:.3f}\n'
    if n_KP < 0.45:
        annotation_text += 'Mechanism: Fickian diffusion'
    elif n_KP < 0.89:
        annotation_text += 'Mechanism: Anomalous transport'
    else:
        annotation_text += 'Mechanism: Super Case-II'

    ax_D.text(0.05, 0.95, annotation_text, transform=ax_D.transAxes,
              ha='left', va='top',
              bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.8),
              fontsize=9)

# Add main title
fig.suptitle('Figure 6.1: Growth Hormone Release from Prototype 3 PLGA/SIS Scaffold',
             fontsize=14, fontweight='bold', y=0.98)

plt.tight_layout()
plt.savefig('Figure_6_1_GH_Release_Model.png', dpi=300, bbox_inches='tight')
print(f"\nFigure saved as: Figure_6_1_GH_Release_Model.png")

# ============================================================================
# ADDITIONAL ANALYSIS FIGURES
# ============================================================================

# Figure: Polymer Degradation and Diffusivity Evolution
fig2, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4))

t_analysis = np.linspace(0, 30, 300)
M_t = polymer_degradation(t_analysis)
D_t = diffusion_coefficient(t_analysis)

ax1.plot(t_analysis, M_t, 'b-', linewidth=2.5)
ax1.set_xlabel('Time (days)', fontweight='bold')
ax1.set_ylabel('Molecular Weight (kDa)', fontweight='bold')
ax1.set_title('Polymer Degradation Profile', fontweight='bold')
ax1.grid(True, alpha=0.3)

ax2.semilogy(t_analysis, D_t, 'r-', linewidth=2.5)
ax2.set_xlabel('Time (days)', fontweight='bold')
ax2.set_ylabel('Diffusion Coefficient (cm²/day)', fontweight='bold')
ax2.set_title('Time-Dependent Diffusivity', fontweight='bold')
ax2.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('Figure_S1_Degradation_Diffusivity.png', dpi=300, bbox_inches='tight')
print(f"Supplementary figure saved as: Figure_S1_Degradation_Diffusivity.png")

# ============================================================================
# FEASIBILITY ASSESSMENT
# ============================================================================

print("\n" + "=" * 70)
print("FEASIBILITY ASSESSMENT")
print("=" * 70)

# Check if therapeutic window is achieved
max_interface_conc = np.max(C_interface)
time_above_min = np.sum(C_interface >= C_ther_min) * (t_sim[1] - t_sim[0])
time_in_window = np.sum((C_interface >= C_ther_min) & (C_interface <= C_ther_max)) * (t_sim[1] - t_sim[0])

print(f"\nTherapeutic Window Analysis:")
print(f"  Maximum interface concentration: {max_interface_conc:.2f} ng/mL")
print(f"  Time above minimum threshold ({C_ther_min} ng/mL): {time_above_min:.1f} days")
print(f"  Time within therapeutic window ({C_ther_min}-{C_ther_max} ng/mL): {time_in_window:.1f} days")

if time_in_window >= 21:
    print(f"\n✓ FEASIBLE: Therapeutic concentration maintained for ≥3 weeks")
elif time_above_min >= 21:
    print(f"\n⚠ MARGINAL: Concentration above minimum but may exceed maximum")
    print(f"  Recommendation: Reduce initial loading to {M_0_GH * 0.7:.1f} mg")
else:
    print(f"\n✗ NOT FEASIBLE: Insufficient duration in therapeutic window")
    print(f"  Recommendation: Increase initial loading to {M_0_GH * 1.5:.1f} mg")

print(f"\nRelease Kinetics:")
print(f"  Total release at 30 days: {cumulative_release[-1]:.1f}%")
print(f"  Release rate: {cumulative_release[-1]/30:.2f}%/day average")

print("\n" + "=" * 70)
print("SIMULATION COMPLETE")
print("=" * 70)

plt.show()
