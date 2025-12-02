# Example

This example generates a trajectory and plots all kinematic derivatives.

```@example trajectory
using PruntTrajectories
using Plots

# Define kinematic limits
accel_max = 50.0    # Maximum acceleration
jerk_max = 100.0    # Maximum jerk
snap_max = 200.0    # Maximum snap
crackle_max = 500.0 # Maximum crackle

# Trajectory parameters
start_vel = 10.0    # Starting velocity
max_vel = 80.0      # Maximum allowed velocity
end_vel = 5.0       # Ending velocity
distance = 500.0    # Total distance to travel

# Generate the optimal profile
profile = optimal_full_profile(
    start_vel, max_vel, end_vel, distance,
    accel_max, jerk_max, snap_max, crackle_max
)

println("Total time: ", total_time(profile), " s")
println("Acceleration phase: ", total_time(profile.accel), " s")
println("Coast phase: ", profile.coast, " s")
println("Deceleration phase: ", total_time(profile.decel), " s")
```

## Plotting All Derivatives

```@example trajectory
# Generate time samples
t_total = total_time(profile)
n_points = 1000
times = range(0, t_total * 0.9999, length=n_points)

# Compute all derivatives at each time point
distances = [distance_at_time(profile, t, crackle_max, start_vel) for t in times]
velocities = [velocity_at_time(profile, t, crackle_max, start_vel) for t in times]
accelerations = [acceleration_at_time(profile, t, crackle_max) for t in times]
jerks = [jerk_at_time(profile, t, crackle_max) for t in times]
snaps = [snap_at_time(profile, t, crackle_max) for t in times]
crackles = [crackle_at_time(profile, t, crackle_max) for t in times]
nothing # hide
```

```@setup trajectory
# Create subplots
p1 = plot(times, distances, label="Distance", xlabel="Time (s)", ylabel="Distance", linewidth=2, color=:blue)
p2 = plot(times, velocities, label="Velocity", xlabel="Time (s)", ylabel="Velocity", linewidth=2, color=:green)
hline!(p2, [max_vel], label="Max velocity", linestyle=:dash, color=:red, alpha=0.5)

p3 = plot(times, accelerations, label="Acceleration", xlabel="Time (s)", ylabel="Acceleration", linewidth=2, color=:orange)
hline!(p3, [accel_max, -accel_max], label="Limits", linestyle=:dash, color=:red, alpha=0.5)

p4 = plot(times, jerks, label="Jerk", xlabel="Time (s)", ylabel="Jerk", linewidth=2, color=:purple)
hline!(p4, [jerk_max, -jerk_max], label="Limits", linestyle=:dash, color=:red, alpha=0.5)

p5 = plot(times, snaps, label="Snap", xlabel="Time (s)", ylabel="Snap", linewidth=2, color=:brown)
hline!(p5, [snap_max, -snap_max], label="Limits", linestyle=:dash, color=:red, alpha=0.5)

p6 = plot(times, crackles, label="Crackle", xlabel="Time (s)", ylabel="Crackle", linewidth=2, color=:teal)
hline!(p6, [crackle_max, -crackle_max], label="Limits", linestyle=:dash, color=:red, alpha=0.5)

# Combine into a single figure
plot(p1, p2, p3, p4, p5, p6, layout=(3, 2), size=(900, 800), legend=:topright)
savefig("trajectory_derivatives.svg")
```

![Trajectory Derivatives](trajectory_derivatives.svg)

## Querying Profile at Specific Times

```@example trajectory
# Query at the midpoint
t_mid = t_total / 2

println("At t = ", t_mid, " s:")
println("  Distance: ", distance_at_time(profile, t_mid, crackle_max, start_vel))
println("  Velocity: ", velocity_at_time(profile, t_mid, crackle_max, start_vel))
println("  Acceleration: ", acceleration_at_time(profile, t_mid, crackle_max))
println("  Jerk: ", jerk_at_time(profile, t_mid, crackle_max))
println("  Snap: ", snap_at_time(profile, t_mid, crackle_max))
println("  Crackle: ", crackle_at_time(profile, t_mid, crackle_max))
```

## Using `distance_at_time_with_accel_flag`

The [`distance_at_time_with_accel_flag`](@ref) function returns both the distance and a boolean indicating whether the query time is past the acceleration phase:

```@example trajectory
# During acceleration
dist1, past_accel1 = distance_at_time_with_accel_flag(profile, 0.1, crackle_max, start_vel)
println("At t=0.1s: distance=", round(dist1, digits=3), ", past_accel=", past_accel1)

# During coast/decel
t_late = total_time(profile.accel) + 0.1
dist2, past_accel2 = distance_at_time_with_accel_flag(profile, t_late, crackle_max, start_vel)
println("At t=", round(t_late, digits=3), "s: distance=", round(dist2, digits=3), ", past_accel=", past_accel2)
```
