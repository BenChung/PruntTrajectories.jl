# PruntTrajectories.jl

A Julia package for generating time-optimal motion profiles with bounded crackle (5th derivative of position). This is useful for smooth motion control in robotics, CNC machines, and 3D printers.
The algorithm is a direct translation from the one used by [Prunt](https://prunt3d.com/).

## Overview

PruntTrajectories generates feedrate profiles that minimize travel time while respecting kinematic constraints on:

- **Crackle** (5th derivative of position)
- **Snap** (4th derivative of position)
- **Jerk** (3rd derivative of position)
- **Acceleration** (2nd derivative of position)
- **Velocity** (1st derivative of position)

The package uses a 15-phase profile structure for acceleration/deceleration segments and a 31-phase structure for complete trajectories (acceleration + coast + deceleration).

## Installation

```julia
using Pkg
Pkg.add("PruntTrajectories")
```

## Quick Start

```julia
using PruntTrajectories

# Define kinematic limits
accel_max = 50.0    # Maximum acceleration
jerk_max = 100.0    # Maximum jerk
snap_max = 200.0    # Maximum snap
crackle_max = 500.0 # Maximum crackle

# Generate a trajectory
start_vel = 0.0     # Starting velocity
max_vel = 100.0     # Maximum allowed velocity
end_vel = 0.0       # Ending velocity
distance = 1000.0   # Total distance to travel

profile = optimal_full_profile(
    start_vel, max_vel, end_vel, distance,
    accel_max, jerk_max, snap_max, crackle_max
)

# Query the trajectory at any time
t = 0.5
vel = velocity_at_time(profile, t, crackle_max, start_vel)
dist = distance_at_time(profile, t, crackle_max, start_vel)
```

See the [Example](@ref) page for a complete example with plotting.

## Types

```@docs
FeedrateProfileTimes
FeedrateProfile
```
