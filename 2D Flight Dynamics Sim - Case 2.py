# imports
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
import os
from datetime import datetime

# Create timestamped output folder for this run
run_timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
output_dir = os.path.join("simulation_outputs", f"run_{run_timestamp}")
os.makedirs(output_dir, exist_ok=True)
print(f"Saving plots to: {output_dir}")

# Gen Variables
g = 9.81
rho = 1.225
H = 8500

# F-16 Inspired parameters( Double check )(not after burner at military thrust)
m_v = 8500      # kg dry mass (actual F-16 is ~8,500 kg)
m_f = 3300      # kg fuel (internal fuel load)
mm = 10      # kg/s fuel flow at military thrust
ve = 3000       # m/s effective exhaust velocity
A = 27.87       # m² wing reference area
Aw = 27.87      # same reference area for lift
CL = 0.5      # cruise lift coefficient
Cd = 0.02       # clean configuration drag coefficient

# Initial conditions
# case 1
# starts in the air so lift will counteract gravity property
y0 = 100
x0 =0
# helps to prevent sinking caused by gravity at start. Sacrifice since this is a particle ,and we are focusing on
# the relationship between thrust vectoring and non-thrust vectoring
# have at min velocity to produce thrust modeling the takeoff point from a runway
vy0= 10
vx0 = 120
F0 = mm * ve
m_f0 = m_f
state0 = [y0, x0,vy0,vx0, m_f0, np.radians(45)]
t_span = (0, 2500)
theta_b = np.radians(45)
theta_v = 0
target = 2500
max_angle = np.radians(60)
#identifying phase
phase = [1]
# to record what happening at different phases
phase_log = []

# Separate phase tracking for case2b
phase_2b = [1]
phase_log_2b = []
min_angle = np.radians(45)

def save_and_show(filename, dpi=150):
    filepath = os.path.join(output_dir, filename)
    plt.savefig(filepath, dpi=dpi, bbox_inches='tight')
    print(f"Saved: {filename}")
    plt.show()
    plt.close()  # prevents memory buildup across many plots

#Non-Thrust Vectoring Case
def case2a(t, state):
    y, x, vy, vx, m_fuel, theta_b = state

    # ---------- Mass & Thrust ----------
    if m_fuel > 0:
        Ft = F0
        mt = m_v + m_fuel
        dm_fuel_dt = -mm
    else:
        Ft = 0
        mt = m_v
        dm_fuel_dt = 0

    v_mag = np.hypot(vx, vy)

    # ---------- Atmosphere ----------
    rho_a = rho * np.exp(-y / H)
    q = max(0.5 * rho_a * v_mag ** 2, 1.0)

    # ---------- Phase Transitions (with logging) ----------
    if phase[0] == 1 and y >= target / 2:
        phase[0] = 2
        if 2 > (phase_log[-1][0] if phase_log else 0):
            phase_log.append((2, t, y, v_mag))
            print(f"Phase 2 at y={y:.1f}m, v={v_mag:.1f}m/s, t={t:.1f}s")

    if phase[0] == 2 and y >= target:
        phase[0] = 3
        if 3 > (phase_log[-1][0] if phase_log else 0):
            phase_log.append((3, t, y, v_mag))
            print(f"Phase 3 at y={y:.1f}m, v={v_mag:.1f}m/s, t={t:.1f}s")

    if phase[0] == 3 and m_fuel <= mm:
        phase[0] = 4
        if 4 > (phase_log[-1][0] if phase_log else 0):
            phase_log.append((4, t, y, v_mag))
            print(f"Phase 4 at y={y:.1f}m, v={v_mag:.1f}m/s, t={t:.1f}s")

    # ---------- Flight Path Angle ----------
    if v_mag < 20.0 or vx <= 0:
        theta_v = theta_b
    else:
        theta_v = np.arctan2(vy, vx)

    D = q * A * Cd
    L = q * Aw * CL
    #---------- Phase Logic -------
    if phase[0] == 1:
        theta_b_target = min_angle + (max_angle - min_angle) * (y / (target / 2))
    elif phase[0] == 2:
        progress = (y - target / 2) / (target / 2)
        theta_b_target = max(max_angle * (1 - progress), 0)
    elif phase[0] == 3:
        theta_b_target = np.radians(5)
    else:
        theta_b_target = theta_b  # freeze during glide

    # ---------- Rate Limiter ----------
    max_rate = np.radians(10)
    dtheta_b_dt = np.clip(theta_b_target - theta_b, -max_rate, max_rate)


    # ---------- Equations of Motion ----------
    dydt = vy
    dxdt = vx
    dvydt = (Ft * np.sin(theta_b) - mt * g - D * np.sin(theta_v) + L * np.cos(theta_v)) / mt
    dvxdt = (Ft * np.cos(theta_b) - D * np.cos(theta_v) - L * np.sin(theta_v)) / mt

    return [dydt, dxdt, dvydt, dvxdt, dm_fuel_dt, dtheta_b_dt]


# Thrust Vectoring Case
def case2b(t, state):
    y, x, vy, vx, m_fuel, theta_b = state

    # ---------- Mass & Thrust ----------
    if m_fuel > 0:
        Ft = F0
        mt = m_v + m_fuel
        dm_fuel_dt = -mm
    else:
        Ft = 0
        mt = m_v
        dm_fuel_dt = 0

    v_mag = np.hypot(vx, vy)

    # ---------- Atmosphere ----------
    rho_a = rho * np.exp(-y / H)
    q = max(0.5 * rho_a * v_mag ** 2, 1.0)

    # ---------- Phase Transitions (with logging) ----------
    if phase_2b[0] == 1 and y >= target / 2:
        phase_2b[0] = 2
        if 2 > (phase_log_2b[-1][0] if phase_log_2b else 0):
            phase_log_2b.append((2, t, y, v_mag))
            print(f"2B Phase 2 at y={y:.1f}m, v={v_mag:.1f}m/s, t={t:.1f}s")

    if phase_2b[0] == 2 and y >= target:
        phase_2b[0] = 3
        if 3 > (phase_log_2b[-1][0] if phase_log_2b else 0):
            phase_log_2b.append((3, t, y, v_mag))
            print(f"2B Phase 3 at y={y:.1f}m, v={v_mag:.1f}m/s, t={t:.1f}s")

    if phase_2b[0] == 3 and m_fuel <= mm:
        phase_2b[0] = 4
        if 4 > (phase_log_2b[-1][0] if phase_log_2b else 0):
            phase_log_2b.append((4, t, y, v_mag))
            print(f"2B Phase 4 at y={y:.1f}m, v={v_mag:.1f}m/s, t={t:.1f}s")

    # ---------- Flight Path Angle ----------
    if v_mag < 20.0 or vx <= 0:
        theta_v = theta_b
    else:
        theta_v = np.arctan2(vy, vx)

    D = q * A * Cd
    L = q * Aw * CL

    # ---------- Thrust Vectoring — only difference from case2a ----------
    # theta_b continuously follows theta_v rather than a prescribed schedule
    # No phase schedule needed — thrust always aligns with velocity vector
    theta_b_target = theta_v

    # ---------- Rate Limiter ----------
    max_rate = np.radians(10)
    dtheta_b_dt = np.clip(theta_b_target - theta_b, -max_rate, max_rate)

    # ---------- Equations of Motion ----------
    dydt = vy
    dxdt = vx
    dvydt = (Ft * np.sin(theta_b) - mt * g - D * np.sin(theta_v) + L * np.cos(theta_v)) / mt
    dvxdt = (Ft * np.cos(theta_b) - D * np.cos(theta_v) - L * np.sin(theta_v)) / mt

    return [dydt, dxdt, dvydt, dvxdt, dm_fuel_dt, dtheta_b_dt]




#*Events*
# Check for when the ground is hit after fuel is totally consumed
def hit_ground(t, state):
    y, x, vy, vx, m_f,theta_b  = state
    return y - 0.1
    # event occurs when y = 0


hit_ground.terminal = True  # stop integration
hit_ground.direction = -1  # only trigger when y is decreasing

phase[0] = 1
phase_log = []
solution = solve_ivp(
    case2a,
    t_span,
    state0,
    events=[hit_ground],
    max_step=1.0,
    method='RK45'
)

# Run case2b
phase_2b[0] = 1
phase_log_2b = []
solution_2b = solve_ivp(
    case2b,
    t_span,
    state0,
    events=[hit_ground],
    max_step=1.0,
    method='RK45'
)

# Diagnostics
print(f"Termination status: {solution.status}")
print(f"Termination message: {solution.message}")
print(f"Final time: {solution.t[-1]:.2f}s")
print(f"Final altitude: {solution.y[0][-1]:.2f}m")
print(f"Final vx: {solution.y[3][-1]:.2f}m/s")

# Extract states
t = solution.t
y = solution.y[0]
x = solution.y[1]
vy = solution.y[2]
vx = solution.y[3]
m_fuel_sol = solution.y[4]
theta_b_sol = solution.y[5]  # directly from state — no reconstruction needed

# Extract case2b states
t_2b = solution_2b.t
y_2b = solution_2b.y[0]
x_2b = solution_2b.y[1]
vy_2b = solution_2b.y[2]
vx_2b = solution_2b.y[3]
m_fuel_2b = solution_2b.y[4]
theta_b_2b = solution_2b.y[5]

v_mag_2b = np.hypot(vx_2b, vy_2b)

phase_times_2b = [entry[1] for entry in phase_log_2b]
phase_alts_2b  = [entry[2] for entry in phase_log_2b]
phase_nums_2b  = [entry[0] for entry in phase_log_2b]

theta_v_2b = np.arctan2(vy_2b, vx_2b)
theta_b_deg_2b = np.degrees(theta_b_2b)
theta_v_deg_2b = np.degrees(theta_v_2b)

# Solving v_mag over time
v_mag_sol = np.hypot(vx, vy)

# Extract phase transition points for markers
phase_times = [entry[1] for entry in phase_log]
phase_alts  = [entry[2] for entry in phase_log]
phase_nums  = [entry[0] for entry in phase_log]
phase_vmags = [entry[3] for entry in phase_log]

# Theta_v recomputed from solution velocity components
# Theta_b now comes directly from state — reconstruction loop removed
theta_v_history = np.arctan2(vy, vx)

# Convert to degrees for readability
theta_b_deg = np.degrees(theta_b_sol)
theta_v_deg = np.degrees(theta_v_history)

# ============================================================
# COMPARISON PLOTS — Case 2A vs Case 2B
# ============================================================

# Plot 1 - Altitude vs Time (Comparison)
plt.figure(figsize=(10, 6))
plt.plot(t, y, label='Case 2A (Non-Thrust Vectoring)', color='blue')
plt.plot(t_2b, y_2b, label='Case 2B (Thrust Vectoring)', color='red')
for i, (pt, py, pn) in enumerate(zip(phase_times, phase_alts, phase_nums)):
    plt.axvline(x=pt, color='blue', linestyle='--', alpha=0.3)
for i, (pt, py, pn) in enumerate(zip(phase_times_2b, phase_alts_2b, phase_nums_2b)):
    plt.axvline(x=pt, color='red', linestyle='--', alpha=0.3)
plt.xlabel("Time (s)")
plt.ylabel("Height (m)")
plt.title("Altitude vs Time — Case 2A vs Case 2B")
plt.grid()
plt.legend()
save_and_show("01_altitude_comparison.png")

# Plot 2 - Trajectory (Comparison)
plt.figure(figsize=(10, 6))
plt.plot(x, y, label='Case 2A (Non-Thrust Vectoring)', color='blue')
plt.plot(x_2b, y_2b, label='Case 2B (Thrust Vectoring)', color='red')
for i, (pt, py, pn) in enumerate(zip(phase_times, phase_alts, phase_nums)):
    x_phase = x[np.argmin(np.abs(t - pt))]
    plt.scatter(x_phase, py, color='blue', s=50, zorder=5)
for i, (pt, py, pn) in enumerate(zip(phase_times_2b, phase_alts_2b, phase_nums_2b)):
    x_phase = x_2b[np.argmin(np.abs(t_2b - pt))]
    plt.scatter(x_phase, py, color='red', s=50, zorder=5)
plt.xlabel("Downrange Distance (m)")
plt.ylabel("Altitude (m)")
plt.title("Trajectory — Case 2A vs Case 2B")
plt.grid()
plt.legend()
save_and_show("02_trajectory_comparison.png")

# Plot 3 - Velocity Components (Comparison)
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 10))

# Horizontal velocity
ax1.plot(t, vx, label='Case 2A vx', color='green')
ax1.plot(t_2b, vx_2b, label='Case 2B vx', color='darkgreen', linestyle='--')
for i, (pt, py, pn) in enumerate(zip(phase_times, phase_alts, phase_nums)):
    ax1.axvline(x=pt, color='green', linestyle='--', alpha=0.3)
for i, (pt, py, pn) in enumerate(zip(phase_times_2b, phase_alts_2b, phase_nums_2b)):
    ax1.axvline(x=pt, color='darkgreen', linestyle='--', alpha=0.3)
ax1.set_xlabel("Time (s)")
ax1.set_ylabel("Horizontal Velocity (m/s)")
ax1.set_title("Horizontal Velocity vs Time — Case 2A vs Case 2B")
ax1.grid()
ax1.legend()

# Vertical velocity
ax2.plot(t, vy, label='Case 2A vy', color='orange')
ax2.plot(t_2b, vy_2b, label='Case 2B vy', color='darkorange', linestyle='--')
for i, (pt, py, pn) in enumerate(zip(phase_times, phase_alts, phase_nums)):
    ax2.axvline(x=pt, color='orange', linestyle='--', alpha=0.3)
for i, (pt, py, pn) in enumerate(zip(phase_times_2b, phase_alts_2b, phase_nums_2b)):
    ax2.axvline(x=pt, color='darkorange', linestyle='--', alpha=0.3)
ax2.set_xlabel("Time (s)")
ax2.set_ylabel("Vertical Velocity (m/s)")
ax2.set_title("Vertical Velocity vs Time — Case 2A vs Case 2B")
ax2.grid()
ax2.legend()

plt.tight_layout()
save_and_show("03_velocity_components_comparison.png")

# Plot 4 - Total Speed (Comparison)
plt.figure(figsize=(10, 6))
plt.plot(t, v_mag_sol, label='Case 2A (Non-Thrust Vectoring)', color='blue')
plt.plot(t_2b, v_mag_2b, label='Case 2B (Thrust Vectoring)', color='red')
plt.axhline(y=75, color='black', linestyle='--',
            label='Min cruise speed (~75 m/s)', alpha=0.7)
for i, (pt, py, pn) in enumerate(zip(phase_times, phase_alts, phase_nums)):
    plt.axvline(x=pt, color='blue', linestyle='--', alpha=0.3)
for i, (pt, py, pn) in enumerate(zip(phase_times_2b, phase_alts_2b, phase_nums_2b)):
    plt.axvline(x=pt, color='red', linestyle='--', alpha=0.3)
plt.xlabel("Time (s)")
plt.ylabel("Speed (m/s)")
plt.title("Total Speed vs Time — Case 2A vs Case 2B")
plt.grid()
plt.legend()
save_and_show("04_total_speed_comparison.png")

# Plot 5 - Theta Comparison (Subplots for each case)
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 10))

# Case 2A
ax1.plot(t, theta_b_deg, label='θ_b (body angle)', color='blue')
ax1.plot(t, theta_v_deg, label='θ_v (velocity angle)', color='orange')
for i, (pt, py, pn) in enumerate(zip(phase_times, phase_alts, phase_nums)):
    ax1.axvline(x=pt, color='red', linestyle='--', alpha=0.5)
    ax1.annotate(f'Phase {pn}',
                xy=(pt, 0),
                xytext=(pt + 5, 10),
                fontsize=8,
                color='red')
ax1.axhline(y=0, color='black', linestyle='-', alpha=0.3)
ax1.set_xlabel("Time (s)")
ax1.set_ylabel("Angle (degrees)")
ax1.set_title("θ_b vs θ_v — Case 2A (Non-Thrust Vectoring)")
ax1.grid()
ax1.legend()

# Case 2B
ax2.plot(t_2b, theta_b_deg_2b, label='θ_b (body angle)', color='blue')
ax2.plot(t_2b, theta_v_deg_2b, label='θ_v (velocity angle)', color='orange')
for i, (pt, py, pn) in enumerate(zip(phase_times_2b, phase_alts_2b, phase_nums_2b)):
    ax2.axvline(x=pt, color='red', linestyle='--', alpha=0.5)
    ax2.annotate(f'Phase {pn}',
                xy=(pt, 0),
                xytext=(pt + 5, 10),
                fontsize=8,
                color='red')
ax2.axhline(y=0, color='black', linestyle='-', alpha=0.3)
ax2.set_xlabel("Time (s)")
ax2.set_ylabel("Angle (degrees)")
ax2.set_title("θ_b vs θ_v — Case 2B (Thrust Vectoring)")
ax2.grid()
ax2.legend()

plt.tight_layout()
save_and_show("05_theta_comparison_subplots.png")

# Plot 6 - Linear Momentum Components (Comparison)
# Compute mass for Case 2B
mass_sol_2b = m_v + m_fuel_2b

# Momentum components for Case 2B

p_x_2b = mass_sol_2b * vx_2b
p_y_2b = mass_sol_2b * vy_2b

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 10))

# Horizontal momentum
ax1.plot(t, p_x, label='Case 2A p_x', color='green')
ax1.plot(t_2b, p_x_2b, label='Case 2B p_x', color='darkgreen', linestyle='--')
for i, (pt, py_alt, pn) in enumerate(zip(phase_times, phase_alts, phase_nums)):
    ax1.axvline(x=pt, color='green', linestyle='--', alpha=0.3)
for i, (pt, py_alt, pn) in enumerate(zip(phase_times_2b, phase_alts_2b, phase_nums_2b)):
    ax1.axvline(x=pt, color='darkgreen', linestyle='--', alpha=0.3)
ax1.axhline(y=0, color='black', linestyle='-', alpha=0.3)
ax1.set_xlabel("Time (s)")
ax1.set_ylabel("Horizontal Momentum (kg·m/s)")
ax1.set_title("Horizontal Momentum vs Time — Case 2A vs Case 2B")
ax1.grid()
ax1.legend()

# Vertical momentum
ax2.plot(t, p_y, label='Case 2A p_y', color='orange')
ax2.plot(t_2b, p_y_2b, label='Case 2B p_y', color='darkorange', linestyle='--')
for i, (pt, py_alt, pn) in enumerate(zip(phase_times, phase_alts, phase_nums)):
    ax2.axvline(x=pt, color='orange', linestyle='--', alpha=0.3)
for i, (pt, py_alt, pn) in enumerate(zip(phase_times_2b, phase_alts_2b, phase_nums_2b)):
    ax2.axvline(x=pt, color='darkorange', linestyle='--', alpha=0.3)
ax2.axhline(y=0, color='black', linestyle='-', alpha=0.3)
ax2.set_xlabel("Time (s)")
ax2.set_ylabel("Vertical Momentum (kg·m/s)")
ax2.set_title("Vertical Momentum vs Time — Case 2A vs Case 2B")
ax2.grid()
ax2.legend()

plt.tight_layout()
save_and_show("06_momentum_comparison.png")


# ============================================================
# Run Summary with Case Comparisons
# ============================================================
summary_path = os.path.join(output_dir, "run_summary.txt")
with open(summary_path, 'w') as f:
    f.write(f"F-16 SIMULATION COMPARISON — Case 2A vs Case 2B\n")
    f.write(f"Timestamp: {run_timestamp}\n")
    f.write(f"{'=' * 60}\n\n")

    # Common Parameters
    f.write(f"COMMON PARAMETERS:\n")
    f.write(f"{'─' * 40}\n")
    f.write(f"  Dry mass (m_v)          = {m_v:>8.0f} kg\n")
    f.write(f"  Initial fuel (m_f)      = {m_f:>8.0f} kg\n")
    f.write(f"  Fuel flow rate (mm)     = {mm:>8.0f} kg/s\n")
    f.write(f"  Exhaust velocity (ve)   = {ve:>8.0f} m/s\n")
    f.write(f"  Initial thrust (F0)     = {F0:>8.0f} N\n")
    f.write(f"  Target altitude         = {target:>8.0f} m\n")
    f.write(f"  Initial vx              = {vx0:>8.0f} m/s\n")
    f.write(f"  Initial vy              = {vy0:>8.0f} m/s\n")
    f.write(f"  Initial altitude        = {y0:>8.0f} m\n")
    f.write(f"  Lift coefficient (CL)   = {CL:>8.2f}\n")
    f.write(f"  Drag coefficient (Cd)   = {Cd:>8.2f}\n")
    f.write(f"  Initial theta_b         = {np.degrees(state0[5]):>8.1f} deg\n")
    f.write(f"  Theoretical burn time   = {m_f / mm:>8.1f} s\n\n")

    # Case 2A Results
    f.write(f"CASE 2A — NON-THRUST VECTORING\n")
    f.write(f"{'─' * 40}\n")
    f.write(f"Phase Transitions:\n")
    for entry in phase_log:
        pn, pt, py, pv = entry
        f.write(f"  Phase {pn}: t={pt:8.1f}s  y={py:8.1f}m  v={pv:8.1f}m/s\n")

    # Calculate overshoot for Case 2A
    max_alt_2a = max(y)
    overshoot_2a = ((max_alt_2a - target) / target) * 100

    f.write(f"\nKey Results:\n")
    f.write(f"  Max altitude             = {max_alt_2a:>10.1f} m\n")
    f.write(f"  Target overshoot         = {overshoot_2a:>10.1f}%\n")
    f.write(f"  Max speed                = {max(v_mag_sol):>10.1f} m/s\n")
    f.write(f"  Min speed                = {min(v_mag_sol):>10.1f} m/s\n")
    f.write(f"  Max horizontal range     = {max(x):>10.1f} m\n")
    f.write(f"  Flight time              = {t[-1]:>10.1f} s\n")
    f.write(f"  Final altitude           = {y[-1]:>10.1f} m\n")
    f.write(f"  Final speed              = {v_mag_sol[-1]:>10.1f} m/s\n")
    f.write(f"  Final theta_b            = {theta_b_deg[-1]:>10.1f} deg\n")
    f.write(f"  Final theta_v            = {theta_v_deg[-1]:>10.1f} deg\n")
    f.write(f"  Final fuel remaining     = {m_fuel_sol[-1]:>10.1f} kg\n")

    # Altitude at specific time milestones for Case 2A
    f.write(f"\n  Altitude Milestones:\n")
    time_milestones = [50, 100, 200, 300]
    for tm in time_milestones:
        if tm <= t[-1]:
            idx = np.argmin(np.abs(t - tm))
            f.write(f"    t={tm:4.0f}s:  y={y[idx]:8.1f}m  v={v_mag_sol[idx]:7.1f}m/s\n")

    f.write(f"\n{'=' * 60}\n\n")

    # Case 2B Results
    f.write(f"CASE 2B — THRUST VECTORING\n")
    f.write(f"{'─' * 40}\n")
    f.write(f"Phase Transitions:\n")
    for entry in phase_log_2b:
        pn, pt, py, pv = entry
        f.write(f"  Phase {pn}: t={pt:8.1f}s  y={py:8.1f}m  v={pv:8.1f}m/s\n")

    # Calculate overshoot for Case 2B
    max_alt_2b = max(y_2b)
    overshoot_2b = ((max_alt_2b - target) / target) * 100

    f.write(f"\nKey Results:\n")
    f.write(f"  Max altitude             = {max_alt_2b:>10.1f} m\n")
    f.write(f"  Target overshoot         = {overshoot_2b:>10.1f}%\n")
    f.write(f"  Max speed                = {max(v_mag_2b):>10.1f} m/s\n")
    f.write(f"  Min speed                = {min(v_mag_2b):>10.1f} m/s\n")
    f.write(f"  Max horizontal range     = {max(x_2b):>10.1f} m\n")
    f.write(f"  Flight time              = {t_2b[-1]:>10.1f} s\n")
    f.write(f"  Final altitude           = {y_2b[-1]:>10.1f} m\n")
    f.write(f"  Final speed              = {v_mag_2b[-1]:>10.1f} m/s\n")
    f.write(f"  Final theta_b            = {theta_b_deg_2b[-1]:>10.1f} deg\n")
    f.write(f"  Final theta_v            = {theta_v_deg_2b[-1]:>10.1f} deg\n")
    f.write(f"  Final fuel remaining     = {m_fuel_2b[-1]:>10.1f} kg\n")

    # Altitude at specific time milestones for Case 2B
    f.write(f"\n  Altitude Milestones:\n")
    for tm in time_milestones:
        if tm <= t_2b[-1]:
            idx = np.argmin(np.abs(t_2b - tm))
            f.write(f"    t={tm:4.0f}s:  y={y_2b[idx]:8.1f}m  v={v_mag_2b[idx]:7.1f}m/s\n")

    f.write(f"\n{'=' * 60}\n\n")

    # Direct Comparisons
    f.write(f"DIRECT COMPARISONS — Case 2A vs Case 2B\n")
    f.write(f"{'─' * 40}\n")

    # Altitude comparison
    f.write(f"\n  Altitude Performance:\n")
    f.write(f"    {'Metric':<30} {'Case 2A':>12} {'Case 2B':>12} {'Difference':>12}\n")
    f.write(f"    {'─' * 30} {'─' * 12} {'─' * 12} {'─' * 12}\n")
    f.write(f"    {'Max altitude (m)':<30} {max_alt_2a:>12.1f} {max_alt_2b:>12.1f} {max_alt_2b - max_alt_2a:>+12.1f}\n")
    f.write(
        f"    {'Target overshoot (%)':<30} {overshoot_2a:>11.1f}% {overshoot_2b:>11.1f}% {overshoot_2b - overshoot_2a:>+11.1f}%\n")
    f.write(f"    {'Final altitude (m)':<30} {y[-1]:>12.1f} {y_2b[-1]:>12.1f} {y_2b[-1] - y[-1]:>+12.1f}\n")

    # Speed comparison
    f.write(f"\n  Speed Performance:\n")
    f.write(f"    {'Metric':<30} {'Case 2A':>12} {'Case 2B':>12} {'Difference':>12}\n")
    f.write(f"    {'─' * 30} {'─' * 12} {'─' * 12} {'─' * 12}\n")
    f.write(
        f"    {'Max speed (m/s)':<30} {max(v_mag_sol):>12.1f} {max(v_mag_2b):>12.1f} {max(v_mag_2b) - max(v_mag_sol):>+12.1f}\n")
    f.write(
        f"    {'Min speed (m/s)':<30} {min(v_mag_sol):>12.1f} {min(v_mag_2b):>12.1f} {min(v_mag_2b) - min(v_mag_sol):>+12.1f}\n")
    f.write(
        f"    {'Final speed (m/s)':<30} {v_mag_sol[-1]:>12.1f} {v_mag_2b[-1]:>12.1f} {v_mag_2b[-1] - v_mag_sol[-1]:>+12.1f}\n")

    # Range and endurance comparison
    f.write(f"\n  Range and Endurance:\n")
    f.write(f"    {'Metric':<30} {'Case 2A':>12} {'Case 2B':>12} {'Difference':>12}\n")
    f.write(f"    {'─' * 30} {'─' * 12} {'─' * 12} {'─' * 12}\n")
    f.write(f"    {'Max range (m)':<30} {max(x):>12.1f} {max(x_2b):>12.1f} {max(x_2b) - max(x):>+12.1f}\n")
    f.write(f"    {'Flight time (s)':<30} {t[-1]:>12.1f} {t_2b[-1]:>12.1f} {t_2b[-1] - t[-1]:>+12.1f}\n")

    # Angle comparison
    f.write(f"\n  Final Angle Comparison:\n")
    f.write(f"    {'Metric':<30} {'Case 2A':>12} {'Case 2B':>12} {'Difference':>12}\n")
    f.write(f"    {'─' * 30} {'─' * 12} {'─' * 12} {'─' * 12}\n")
    f.write(
        f"    {'Final theta_b (deg)':<30} {theta_b_deg[-1]:>12.1f} {theta_b_deg_2b[-1]:>12.1f} {theta_b_deg_2b[-1] - theta_b_deg[-1]:>+12.1f}\n")
    f.write(
        f"    {'Final theta_v (deg)':<30} {theta_v_deg[-1]:>12.1f} {theta_v_deg_2b[-1]:>12.1f} {theta_v_deg_2b[-1] - theta_v_deg[-1]:>+12.1f}\n")

    # Time to reach target altitude
    f.write(f"\n  Time to Reach Target ({target} m):\n")
    # Find when each case first reaches target altitude
    idx_target_2a = np.argmax(y >= target) if np.any(y >= target) else None
    idx_target_2b = np.argmax(y_2b >= target) if np.any(y_2b >= target) else None

    if idx_target_2a is not None:
        time_to_target_2a = t[idx_target_2a]
        f.write(f"    Case 2A: {time_to_target_2a:>8.1f}s at y={y[idx_target_2a]:.1f}m\n")
    else:
        f.write(f"    Case 2A: Did not reach target altitude\n")
        time_to_target_2a = None

    if idx_target_2b is not None:
        time_to_target_2b = t_2b[idx_target_2b]
        f.write(f"    Case 2B: {time_to_target_2b:>8.1f}s at y={y_2b[idx_target_2b]:.1f}m\n")
    else:
        f.write(f"    Case 2B: Did not reach target altitude\n")
        time_to_target_2b = None

    if time_to_target_2a is not None and time_to_target_2b is not None:
        f.write(f"    Difference: {time_to_target_2b - time_to_target_2a:>+8.1f}s\n")

    # Efficiency comparison (altitude gained per kg of fuel)
    f.write(f"\n  Fuel Efficiency (Altitude Gain):\n")
    fuel_used_2a = m_f - m_fuel_sol[-1]
    fuel_used_2b = m_f - m_fuel_2b[-1]
    alt_gain_2a = y[-1] - y0
    alt_gain_2b = y_2b[-1] - y0

    if fuel_used_2a > 0:
        efficiency_2a = alt_gain_2a / fuel_used_2a
        f.write(f"    Case 2A: {efficiency_2a:>6.2f} m/kg fuel\n")
    else:
        f.write(f"    Case 2A: N/A (no fuel used)\n")

    if fuel_used_2b > 0:
        efficiency_2b = alt_gain_2b / fuel_used_2b
        f.write(f"    Case 2B: {efficiency_2b:>6.2f} m/kg fuel\n")
    else:
        f.write(f"    Case 2B: N/A (no fuel used)\n")

    if fuel_used_2a > 0 and fuel_used_2b > 0:
        f.write(f"    Improvement: {((efficiency_2b - efficiency_2a) / efficiency_2a * 100):>+6.1f}%\n")

    f.write(f"\n{'=' * 60}\n")
    f.write(f"\nSUMMARY:\n")
    f.write(
        f"  The thrust vectoring case (2B) {'improved' if max_alt_2b > max_alt_2a else 'reduced'} maximum altitude by {abs(max_alt_2b - max_alt_2a):.1f}m\n")
    f.write(f"  compared to the non-thrust vectoring case (2A).\n")
    f.write(f"  Target altitude overshoot: Case 2A = {overshoot_2a:.1f}%, Case 2B = {overshoot_2b:.1f}%\n")

print(f"Run summary saved to: {summary_path}")
