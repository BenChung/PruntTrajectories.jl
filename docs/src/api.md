# API Reference

## Profile Generation

```@docs
optimal_full_profile
optimal_profile_for_distance
optimal_profile_for_delta_v
```

## Time Queries

```@docs
total_time
```

## Kinematic Queries

These functions return the value of each kinematic quantity at a given time.

```@docs
crackle_at_time
snap_at_time
jerk_at_time
acceleration_at_time
velocity_at_time
distance_at_time
distance_at_time_with_accel_flag
```

## Fast Queries

Optimized functions for querying at the end of a profile.

```@docs
fast_velocity_at_max_time
fast_distance_at_max_time
```
