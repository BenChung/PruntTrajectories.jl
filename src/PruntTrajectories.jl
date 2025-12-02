module PruntTrajectories

export FeedrateProfileTimes, FeedrateProfile
export Fast_Distance_At_Max_Time, Fast_Velocity_At_Max_Time
export Total_Time
export Crackle_At_Time, Snap_At_Time, Jerk_At_Time, Acceleration_At_Time, Velocity_At_Time, Distance_At_Time
export Distance_At_Time_With_Accel_Flag
export Optimal_Profile_For_Distance, Optimal_Profile_For_Delta_V, Optimal_Full_Profile

"""
    FeedrateProfileTimes

Represents the timings for segments in a 15-phase motion profile. Note that some times are
used for multiple segments.
"""
struct FeedrateProfileTimes
    t1::Float64
    t2::Float64
    t3::Float64
    t4::Float64
end

FeedrateProfileTimes() = FeedrateProfileTimes(0.0, 0.0, 0.0, 0.0)

"""
    FeedrateProfile

Represents the timings for segments in a 31-phase motion profile.
"""
struct FeedrateProfile
    accel::FeedrateProfileTimes
    coast::Float64
    decel::FeedrateProfileTimes
end

"""
    Fast_Distance_At_Max_Time(profile::FeedrateProfileTimes, max_crackle, start_vel)

Calculates the total distance covered during an acceleration or deceleration phase defined
by `profile` with a given starting velocity and maximum crackle.

For an acceleration phase `max_crackle` should be positive and for a deceleration phase
`max_crackle` should be negative.

This function is an optimised version of [`Distance_At_Time`](@ref) where `T` is equal to
[`Total_Time`](@ref)`(profile)`.
"""
function Fast_Distance_At_Max_Time(profile::FeedrateProfileTimes, max_crackle, start_vel)
    T1 = profile.t1
    T2 = profile.t2
    T3 = profile.t3
    T4 = profile.t4
    Cm = max_crackle
    Vs = start_vel

    return (Vs + Cm * T1 * (T1 + T2) * (2.0 * T1 + T2 + T3) * (4.0 * T1 + 2.0 * T2 + T3 + T4) / 2.0) *
           (8.0 * T1 + 4.0 * T2 + 2.0 * T3 + T4)
end

"""
    Fast_Velocity_At_Max_Time(profile::FeedrateProfileTimes, max_crackle, start_vel)

Calculates the final velocity after an acceleration or deceleration phase defined by
`profile` with a given starting velocity and maximum crackle.

For an acceleration phase `max_crackle` should be positive and for a deceleration phase
`max_crackle` should be negative.

This function is an optimised version of [`Velocity_At_Time`](@ref) where `T` is equal to
[`Total_Time`](@ref)`(profile)`.
"""
function Fast_Velocity_At_Max_Time(profile::FeedrateProfileTimes, max_crackle, start_vel)
    T1 = profile.t1
    T2 = profile.t2
    T3 = profile.t3
    T4 = profile.t4
    Cm = max_crackle
    Vs = start_vel

    return Vs + Cm * T1 * (T1 + T2) * (2.0 * T1 + T2 + T3) * (4.0 * T1 + 2.0 * T2 + T3 + T4)
end

"""
    Total_Time(times::FeedrateProfileTimes)

Calculates the total duration of a single acceleration or deceleration phase. This is not
equivalent to the sum of components as some components are used multiple times.
"""
function Total_Time(times::FeedrateProfileTimes)
    return 8.0 * times.t1 + 4.0 * times.t2 + 2.0 * times.t3 + times.t4
end

"""
    Total_Time(profile::FeedrateProfile)

Calculates the total duration of a complete feedrate profile. This is not equivalent to
the sum of components as some components are used multiple times.
"""
function Total_Time(profile::FeedrateProfile)
    return Total_Time(profile.accel) + profile.coast + Total_Time(profile.decel)
end

"""
    Crackle_At_Time(profile::FeedrateProfileTimes, T, max_crackle)

Returns the crackle at a specific time `T` within a single acceleration or deceleration
phase. The crackle will be either `+max_crackle`, `-max_crackle`, or zero. For an
acceleration phase `max_crackle` should be positive and for a deceleration phase
`max_crackle` should be negative.

The return value may be negative.
"""
function Crackle_At_Time(profile::FeedrateProfileTimes, T, max_crackle)
    @assert T <= Total_Time(profile)
    T1, T2, T3, T4 = profile.t1, profile.t2, profile.t3, profile.t4
    Cm = max_crackle

    if T < T1
        return Cm
    elseif T < T1 + T2
        return 0.0
    elseif T < 2.0 * T1 + T2
        return -Cm
    elseif T < 2.0 * T1 + T2 + T3
        return 0.0
    elseif T < 3.0 * T1 + T2 + T3
        return -Cm
    elseif T < 3.0 * T1 + 2.0 * T2 + T3
        return 0.0
    elseif T < 4.0 * T1 + 2.0 * T2 + T3
        return Cm
    elseif T < 4.0 * T1 + 2.0 * T2 + T3 + T4
        return 0.0
    elseif T < 5.0 * T1 + 2.0 * T2 + T3 + T4
        return -Cm
    elseif T < 5.0 * T1 + 3.0 * T2 + T3 + T4
        return 0.0
    elseif T < 6.0 * T1 + 3.0 * T2 + T3 + T4
        return Cm
    elseif T < 6.0 * T1 + 3.0 * T2 + 2.0 * T3 + T4
        return 0.0
    elseif T < 7.0 * T1 + 3.0 * T2 + 2.0 * T3 + T4
        return Cm
    elseif T < 7.0 * T1 + 4.0 * T2 + 2.0 * T3 + T4
        return 0.0
    else
        return -Cm
    end
end

"""
    Snap_At_Time(profile::FeedrateProfileTimes, T, max_crackle)

Returns the snap at a specific time `T` within a single acceleration or deceleration
phase. For an acceleration phase `max_crackle` should be positive and for a deceleration
phase `max_crackle` should be negative.

The return value may be negative.
"""
function Snap_At_Time(profile::FeedrateProfileTimes, T, max_crackle)
    T1, T2, T3, T4 = profile.t1, profile.t2, profile.t3, profile.t4
    Cm = max_crackle
    function Snap_At_Stage(DT, Stage)
        if Stage == 1
            return Cm * DT
        elseif Stage == 2
            return Snap_At_Stage(T1, 1)
        elseif Stage == 3
            return Snap_At_Stage(T2, 2) - Cm * DT
        elseif Stage == 4
            return Snap_At_Stage(T1, 3)
        elseif Stage == 5
            return Snap_At_Stage(T3, 4) - Cm * DT
        elseif Stage == 6
            return Snap_At_Stage(T1, 5)
        elseif Stage == 7
            return Snap_At_Stage(T2, 6) + Cm * DT
        elseif Stage == 8
            return Snap_At_Stage(T1, 7)
        elseif Stage == 9
            return Snap_At_Stage(T4, 8) - Cm * DT
        elseif Stage == 10
            return Snap_At_Stage(T1, 9)
        elseif Stage == 11
            return Snap_At_Stage(T2, 10) + Cm * DT
        elseif Stage == 12
            return Snap_At_Stage(T1, 11)
        elseif Stage == 13
            return Snap_At_Stage(T3, 12) + Cm * DT
        elseif Stage == 14
            return Snap_At_Stage(T1, 13)
        elseif Stage == 15
            return Snap_At_Stage(T2, 14) - Cm * DT
        end
    end
    @assert T <= Total_Time(profile)

    if T < T1
        return Snap_At_Stage(T, 1)
    elseif T < T1 + T2
        return Snap_At_Stage(T - (T1), 2)
    elseif T < 2.0 * T1 + T2
        return Snap_At_Stage(T - (T1 + T2), 3)
    elseif T < 2.0 * T1 + T2 + T3
        return Snap_At_Stage(T - (2.0 * T1 + T2), 4)
    elseif T < 3.0 * T1 + T2 + T3
        return Snap_At_Stage(T - (2.0 * T1 + T2 + T3), 5)
    elseif T < 3.0 * T1 + 2.0 * T2 + T3
        return Snap_At_Stage(T - (3.0 * T1 + T2 + T3), 6)
    elseif T < 4.0 * T1 + 2.0 * T2 + T3
        return Snap_At_Stage(T - (3.0 * T1 + 2.0 * T2 + T3), 7)
    elseif T < 4.0 * T1 + 2.0 * T2 + T3 + T4
        return Snap_At_Stage(T - (4.0 * T1 + 2.0 * T2 + T3), 8)
    elseif T < 5.0 * T1 + 2.0 * T2 + T3 + T4
        return Snap_At_Stage(T - (4.0 * T1 + 2.0 * T2 + T3 + T4), 9)
    elseif T < 5.0 * T1 + 3.0 * T2 + T3 + T4
        return Snap_At_Stage(T - (5.0 * T1 + 2.0 * T2 + T3 + T4), 10)
    elseif T < 6.0 * T1 + 3.0 * T2 + T3 + T4
        return Snap_At_Stage(T - (5.0 * T1 + 3.0 * T2 + T3 + T4), 11)
    elseif T < 6.0 * T1 + 3.0 * T2 + 2.0 * T3 + T4
        return Snap_At_Stage(T - (6.0 * T1 + 3.0 * T2 + T3 + T4), 12)
    elseif T < 7.0 * T1 + 3.0 * T2 + 2.0 * T3 + T4
        return Snap_At_Stage(T - (6.0 * T1 + 3.0 * T2 + 2.0 * T3 + T4), 13)
    elseif T < 7.0 * T1 + 4.0 * T2 + 2.0 * T3 + T4
        return Snap_At_Stage(T - (7.0 * T1 + 3.0 * T2 + 2.0 * T3 + T4), 14)
    else
        return Snap_At_Stage(T - (7.0 * T1 + 4.0 * T2 + 2.0 * T3 + T4), 15)
    end
end

"""
    Jerk_At_Time(profile::FeedrateProfileTimes, T, max_crackle)

Returns the jerk at a specific time `T` within a single acceleration or deceleration
phase. For an acceleration phase `max_crackle` should be positive and for a deceleration
phase `max_crackle` should be negative.

The return value may be negative.
"""
function Jerk_At_Time(profile::FeedrateProfileTimes, T, max_crackle)
    T1, T2, T3, T4 = profile.t1, profile.t2, profile.t3, profile.t4
    Cm = max_crackle
    
    function Jerk_At_Stage(DT, Stage)
        if Stage == 1
            return Cm * DT^2 / 2.0;
        elseif Stage == 2
            return Jerk_At_Stage(T1, 1) + Cm * DT * T1;
        elseif Stage == 3
            return Jerk_At_Stage(T2, 2) + Cm * DT * (-DT + 2.0 * T1) / 2.0;
        elseif Stage == 4
            return Jerk_At_Stage(T1, 3);
        elseif Stage == 5
            return Jerk_At_Stage(T3, 4) - Cm * DT^2 / 2.0;
        elseif Stage == 6
            return Jerk_At_Stage(T1, 5) - Cm * DT * T1;
        elseif Stage == 7
            return Jerk_At_Stage(T2, 6) + Cm * DT * (DT - 2.0 * T1) / 2.0;
        elseif Stage == 8
            return Jerk_At_Stage(T1, 7);
        elseif Stage == 9
            return Jerk_At_Stage(T4, 8) - Cm * DT^2 / 2.0;
        elseif Stage == 10
            return Jerk_At_Stage(T1, 9) - Cm * DT * T1;
        elseif Stage == 11
            return Jerk_At_Stage(T2, 10) + Cm * DT * (DT - 2.0 * T1) / 2.0;
        elseif Stage == 12
            return Jerk_At_Stage(T1, 11);
        elseif Stage == 13
            return Jerk_At_Stage(T3, 12) + Cm * DT^2 / 2.0;
        elseif Stage == 14
            return Jerk_At_Stage(T1, 13) + Cm * DT * T1;
        elseif Stage == 15
            return Jerk_At_Stage(T2, 14) + Cm * DT * (-DT + 2.0 * T1) / 2.0;
        end
    end
    @assert T <= Total_Time(profile)

    if T < T1
        return Jerk_At_Stage(T, 1)
    elseif T < T1 + T2
        return Jerk_At_Stage(T - (T1), 2)
    elseif T < 2.0 * T1 + T2
        return Jerk_At_Stage(T - (T1 + T2), 3)
    elseif T < 2.0 * T1 + T2 + T3
        return Jerk_At_Stage(T - (2.0 * T1 + T2), 4)
    elseif T < 3.0 * T1 + T2 + T3
        return Jerk_At_Stage(T - (2.0 * T1 + T2 + T3), 5)
    elseif T < 3.0 * T1 + 2.0 * T2 + T3
        return Jerk_At_Stage(T - (3.0 * T1 + T2 + T3), 6)
    elseif T < 4.0 * T1 + 2.0 * T2 + T3
        return Jerk_At_Stage(T - (3.0 * T1 + 2.0 * T2 + T3), 7)
    elseif T < 4.0 * T1 + 2.0 * T2 + T3 + T4
        return Jerk_At_Stage(T - (4.0 * T1 + 2.0 * T2 + T3), 8)
    elseif T < 5.0 * T1 + 2.0 * T2 + T3 + T4
        return Jerk_At_Stage(T - (4.0 * T1 + 2.0 * T2 + T3 + T4), 9)
    elseif T < 5.0 * T1 + 3.0 * T2 + T3 + T4
        return Jerk_At_Stage(T - (5.0 * T1 + 2.0 * T2 + T3 + T4), 10)
    elseif T < 6.0 * T1 + 3.0 * T2 + T3 + T4
        return Jerk_At_Stage(T - (5.0 * T1 + 3.0 * T2 + T3 + T4), 11)
    elseif T < 6.0 * T1 + 3.0 * T2 + 2.0 * T3 + T4
        return Jerk_At_Stage(T - (6.0 * T1 + 3.0 * T2 + T3 + T4), 12)
    elseif T < 7.0 * T1 + 3.0 * T2 + 2.0 * T3 + T4
        return Jerk_At_Stage(T - (6.0 * T1 + 3.0 * T2 + 2.0 * T3 + T4), 13)
    elseif T < 7.0 * T1 + 4.0 * T2 + 2.0 * T3 + T4
        return Jerk_At_Stage(T - (7.0 * T1 + 3.0 * T2 + 2.0 * T3 + T4), 14)
    else
         return Jerk_At_Stage(T - (7.0 * T1 + 4.0 * T2 + 2.0 * T3 + T4), 15)
    end
end

"""
    Acceleration_At_Time(profile::FeedrateProfileTimes, T, max_crackle)

Returns the acceleration at a specific time `T` within a single acceleration or
deceleration phase. For an acceleration phase `max_crackle` should be positive and for
a deceleration phase `max_crackle` should be negative.

The return value may be negative.
"""
function Acceleration_At_Time(profile::FeedrateProfileTimes, T, max_crackle)
    T1, T2, T3, T4 = profile.t1, profile.t2, profile.t3, profile.t4
    Cm = max_crackle
    
    function Acceleration_At_Stage(DT, Stage)
        if Stage == 1
            return Cm * DT^3 / 6.0;
        elseif Stage == 2
            return Acceleration_At_Stage(T1, 1) + Cm * DT * T1 * (DT + T1) / 2.0;
        elseif Stage == 3
               return Acceleration_At_Stage(T2, 2) + Cm * DT * (-DT^2 + 3.0 * DT * T1 + 3.0 * T1 * (T1 + 2.0 * T2)) / 6.0;
        elseif Stage == 4
            return Acceleration_At_Stage(T1, 3) + Cm * DT * T1 * (T1 + T2);
        elseif Stage == 5
            return Acceleration_At_Stage(T3, 4) + Cm * DT * (-DT^2 + 6.0 * T1 * (T1 + T2)) / 6.0;
        elseif Stage == 6
            return Acceleration_At_Stage(T1, 5) + Cm * DT * T1 * (-DT + T1 + 2.0 * T2) / 2.0;
        elseif Stage == 7
            return Acceleration_At_Stage(T2, 6) + Cm * DT * (DT^2 - 3.0 * DT * T1 + 3.0 * T1^2) / 6.0;
        elseif Stage == 8
            return Acceleration_At_Stage(T1, 7);
        elseif Stage == 9
            return Acceleration_At_Stage(T4, 8) - Cm * DT^3 / 6.0;
        elseif Stage == 10
            return Acceleration_At_Stage(T1, 9) + Cm * DT * T1 * (-DT - T1) / 2.0;
        elseif Stage == 11
            return Acceleration_At_Stage(T2, 10) + Cm * DT * (DT^2 - 3.0 * DT * T1 - 3.0 * T1 * (T1 + 2.0 * T2)) / 6.0;
        elseif Stage == 12
            return Acceleration_At_Stage(T1, 11) - Cm * DT * T1 * (T1 + T2);
        elseif Stage == 13
            return Acceleration_At_Stage(T3, 12) + Cm * DT * (DT^2 - 6.0 * T1 * (T1 + T2)) / 6.0;
        elseif Stage == 14
            return Acceleration_At_Stage(T1, 13) + Cm * DT * T1 * (DT - T1 - 2.0 * T2) / 2.0;
        elseif Stage == 15
            return Acceleration_At_Stage(T2, 14) + Cm * DT * (-DT^2 + 3.0 * DT * T1 - 3.0 * T1^2) / 6.0;
        end
    end

    @assert T <= Total_Time(profile)

    if T < T1
        return Acceleration_At_Stage(T, 1);
    elseif T < T1 + T2
        return Acceleration_At_Stage(T - (T1), 2);
    elseif T < 2.0 * T1 + T2
        return Acceleration_At_Stage(T - (T1 + T2), 3);
    elseif T < 2.0 * T1 + T2 + T3
        return Acceleration_At_Stage(T - (2.0 * T1 + T2), 4);
    elseif T < 3.0 * T1 + T2 + T3
        return Acceleration_At_Stage(T - (2.0 * T1 + T2 + T3), 5);
    elseif T < 3.0 * T1 + 2.0 * T2 + T3
        return Acceleration_At_Stage(T - (3.0 * T1 + T2 + T3), 6);
    elseif T < 4.0 * T1 + 2.0 * T2 + T3
        return Acceleration_At_Stage(T - (3.0 * T1 + 2.0 * T2 + T3), 7);
    elseif T < 4.0 * T1 + 2.0 * T2 + T3 + T4
        return Acceleration_At_Stage(T - (4.0 * T1 + 2.0 * T2 + T3), 8);
    elseif T < 5.0 * T1 + 2.0 * T2 + T3 + T4
        return Acceleration_At_Stage(T - (4.0 * T1 + 2.0 * T2 + T3 + T4), 9);
    elseif T < 5.0 * T1 + 3.0 * T2 + T3 + T4
        return Acceleration_At_Stage(T - (5.0 * T1 + 2.0 * T2 + T3 + T4), 10);
    elseif T < 6.0 * T1 + 3.0 * T2 + T3 + T4
        return Acceleration_At_Stage(T - (5.0 * T1 + 3.0 * T2 + T3 + T4), 11);
    elseif T < 6.0 * T1 + 3.0 * T2 + 2.0 * T3 + T4
        return Acceleration_At_Stage(T - (6.0 * T1 + 3.0 * T2 + T3 + T4), 12);
    elseif T < 7.0 * T1 + 3.0 * T2 + 2.0 * T3 + T4
        return Acceleration_At_Stage(T - (6.0 * T1 + 3.0 * T2 + 2.0 * T3 + T4), 13);
    elseif T < 7.0 * T1 + 4.0 * T2 + 2.0 * T3 + T4
        return Acceleration_At_Stage(T - (7.0 * T1 + 3.0 * T2 + 2.0 * T3 + T4), 14);
    else
        return Acceleration_At_Stage(T - (7.0 * T1 + 4.0 * T2 + 2.0 * T3 + T4), 15);
    end
end

"""
    Velocity_At_Time(profile::FeedrateProfileTimes, T, max_crackle, start_vel)

Returns the velocity at a specific time `T` within a single acceleration or deceleration
phase. For an acceleration phase `max_crackle` should be positive and for a deceleration
phase `max_crackle` should be negative.

The return value may be negative.
"""
function Velocity_At_Time(profile::FeedrateProfileTimes, T, max_crackle, start_vel)
    T1, T2, T3, T4 = profile.t1, profile.t2, profile.t3, profile.t4
    Cm = max_crackle
    
    function Velocity_At_Stage(DT, Stage)
        if Stage == 1
            return start_vel + Cm * DT^4 / 24.0;
        elseif Stage == 2 
            return Velocity_At_Stage(T1, 1) + Cm * DT * T1 * (2.0 * DT^2 + 3.0 * DT * T1 + 2.0 * T1^2) / 12.0;
        elseif Stage == 3
            return (Velocity_At_Stage(T2, 2) + Cm * DT * (-DT^3
                      + 4.0 * DT^2 * T1
                      + 6.0 * DT * T1 * (T1 + 2.0 * T2)
                      + 4.0 * T1 * (T1^2 + 3.0 * T1 * T2 + 3.0 * T2^2)) / 24.0)
        elseif Stage == 4
            return (Velocity_At_Stage(T1, 3)
                 + Cm * DT * T1 * (DT * (T1 + T2) + 2.0 * T1^2 + 3.0 * T1 * T2 + T2^2) / 2.0)
        elseif Stage == 5
            return (Velocity_At_Stage(T3, 4)
                 + Cm * DT * (-DT^3
                      + 12.0 * DT * T1 * (T1 + T2)
                      + 12.0 * T1 * (2.0 * T1^2 + 3.0 * T1 * T2 + 2.0 * T1 * T3 + T2^2 + 2.0 * T2 * T3)) / 24.0)
        elseif Stage == 6
            return (Velocity_At_Stage(T1, 5)
                 + Cm * DT * T1 * (-2.0 * DT^2
                      + 3.0 * DT * (T1 + 2.0 * T2)
                      + 22.0 * T1^2
                      + 30.0 * T1 * T2
                      + 12.0 * T1 * T3
                      + 6.0 * T2^2
                      + 12.0 * T2 * T3) / 12.0)
        elseif Stage == 7
            return (Velocity_At_Stage(T2, 6)
                 + Cm * DT * (DT^3
                      - 4.0 * DT^2 * T1
                      + 6.0 * DT * T1^2
                      + 4.0 * T1 * (11.0 * T1^2 + 18.0 * T1 * T2 + 6.0 * T1 * T3 + 6.0 * T2^2 + 6.0 * T2 * T3)) / 24.0)
        elseif Stage == 8
            return (Velocity_At_Stage(T1, 7) + Cm * DT * T1 * (2.0 * T1^2 + 3.0 * T1 * T2 + T1 * T3 + T2^2 + T2 * T3))
        elseif Stage == 9
            return (Velocity_At_Stage(T4, 8)
                 + Cm * DT * (-DT^3 + 24.0 * T1 * (2.0 * T1^2 + 3.0 * T1 * T2 + T1 * T3 + T2^2 + T2 * T3)) / 24.0)
        elseif Stage == 10
            return (Velocity_At_Stage(T1, 9)
                 + Cm * DT * T1 * (-2.0 * DT^2
                      - 3.0 * DT * T1
                      + 22.0 * T1^2
                      + 36.0 * T1 * T2
                      + 12.0 * T1 * T3
                      + 12.0 * T2^2
                      + 12.0 * T2 * T3) / 12.0)
        elseif Stage == 11
            return (Velocity_At_Stage(T2, 10)
                 + Cm * DT * (DT^3
                      - 4.0 * DT^2 * T1
                      - 6.0 * DT * T1 * (T1 + 2.0 * T2)
                      + 4.0 * T1 * (11.0 * T1^2 + 15.0 * T1 * T2 + 6.0 * T1 * T3 + 3.0 * T2^2 + 6.0 * T2 * T3)) / 24.0)
        elseif Stage == 12
            return (Velocity_At_Stage(T1, 11)
                 + Cm * DT * T1 * (-DT * (T1 + T2) + 2.0 * T1^2 + 3.0 * T1 * T2 + 2.0 * T1 * T3 + T2^2 + 2.0 * T2 * T3) / 2.0)
        elseif Stage == 13
            return (Velocity_At_Stage(T3, 12)
                 + Cm * DT * (DT^3 - 12.0 * DT * T1 * (T1 + T2) + 12.0 * T1 * (2.0 * T1^2 + 3.0 * T1 * T2 + T2^2)) / 24.0)
        elseif Stage == 14
            return (Velocity_At_Stage(T1, 13)
                 + Cm * DT * T1 * (2.0 * DT^2 - 3.0 * DT * (T1 + 2.0 * T2) + 2.0 * T1^2 + 6.0 * T1 * T2 + 6.0 * T2^2) / 12.0)
        elseif Stage == 15
            return (Velocity_At_Stage(T2, 14)
                 + Cm * DT * (-DT^3 + 4.0 * DT^2 * T1 - 6.0 * DT * T1^2 + 4.0 * T1^3) / 24.0)
        end
    end
        
    @assert T <= Total_Time(profile)

    if T < T1
        return Velocity_At_Stage(T, 1);
    elseif T < T1 + T2
        return Velocity_At_Stage(T - (T1), 2);
    elseif T < 2.0 * T1 + T2
        return Velocity_At_Stage(T - (T1 + T2), 3);
    elseif T < 2.0 * T1 + T2 + T3
        return Velocity_At_Stage(T - (2.0 * T1 + T2), 4);
    elseif T < 3.0 * T1 + T2 + T3
        return Velocity_At_Stage(T - (2.0 * T1 + T2 + T3), 5);
    elseif T < 3.0 * T1 + 2.0 * T2 + T3
        return Velocity_At_Stage(T - (3.0 * T1 + T2 + T3), 6);
    elseif T < 4.0 * T1 + 2.0 * T2 + T3
        return Velocity_At_Stage(T - (3.0 * T1 + 2.0 * T2 + T3), 7);
    elseif T < 4.0 * T1 + 2.0 * T2 + T3 + T4
        return Velocity_At_Stage(T - (4.0 * T1 + 2.0 * T2 + T3), 8);
    elseif T < 5.0 * T1 + 2.0 * T2 + T3 + T4
        return Velocity_At_Stage(T - (4.0 * T1 + 2.0 * T2 + T3 + T4), 9);
    elseif T < 5.0 * T1 + 3.0 * T2 + T3 + T4
        return Velocity_At_Stage(T - (5.0 * T1 + 2.0 * T2 + T3 + T4), 10);
    elseif T < 6.0 * T1 + 3.0 * T2 + T3 + T4
        return Velocity_At_Stage(T - (5.0 * T1 + 3.0 * T2 + T3 + T4), 11);
    elseif T < 6.0 * T1 + 3.0 * T2 + 2.0 * T3 + T4
        return Velocity_At_Stage(T - (6.0 * T1 + 3.0 * T2 + T3 + T4), 12);
    elseif T < 7.0 * T1 + 3.0 * T2 + 2.0 * T3 + T4
        return Velocity_At_Stage(T - (6.0 * T1 + 3.0 * T2 + 2.0 * T3 + T4), 13);
    elseif T < 7.0 * T1 + 4.0 * T2 + 2.0 * T3 + T4
        return Velocity_At_Stage(T - (7.0 * T1 + 3.0 * T2 + 2.0 * T3 + T4), 14);
    else
        return Velocity_At_Stage(T - (7.0 * T1 + 4.0 * T2 + 2.0 * T3 + T4), 15);
    end
end

"""
    Distance_At_Time(profile::FeedrateProfileTimes, T, max_crackle, start_vel)

Returns the distance from the start point at a specific time `T` within a single
acceleration or deceleration phase. For an acceleration phase `max_crackle` should be
positive and for a deceleration phase `max_crackle` should be negative.

The return value may be negative.
"""
function Distance_At_Time(profile::FeedrateProfileTimes, T, max_crackle, start_vel)
    T1, T2, T3, T4 = profile.t1, profile.t2, profile.t3, profile.t4
    Cm = max_crackle
    
    function Distance_At_Stage(DT, Stage)
        if Stage == 1
            return start_vel * T + Cm * DT^5 / 120.0
        elseif Stage == 2
            return (Distance_At_Stage(T1, 1)
                 + Cm * DT * T1 * (DT^3 + 2.0 * DT^2 * T1 + 2.0 * DT * T1^2 + T1^3) / 24.0)
        elseif Stage == 3
            return (Distance_At_Stage(T2, 2)
                 + Cm * DT * (-DT^4
                      + 5.0 * DT^3 * T1
                      + 10.0 * DT^2 * T1 * (T1 + 2.0 * T2)
                      + 10.0 * DT * T1 * (T1^2 + 3.0 * T1 * T2 + 3.0 * T2^2)
                      + 5.0 * T1 * (T1^3 + 4.0 * T1^2 * T2 + 6.0 * T1 * T2^2 + 4.0 * T2^3)) / 120.0)
        elseif Stage == 4
            return (Distance_At_Stage(T1, 3)
                 + Cm * DT * T1 * (2.0 * DT^2 * (T1 + T2)
                      + 3.0 * DT * (2.0 * T1^2 + 3.0 * T1 * T2 + T2^2)
                      + 7.0 * T1^3
                      + 14.0 * T1^2 * T2
                      + 9.0 * T1 * T2^2
                      + 2.0 * T2^3) / 12.0)
        elseif Stage == 5
            return (Distance_At_Stage(T3, 4)
                 + Cm * DT * (-DT^4
                      + 20.0 * DT^2 * T1 * (T1 + T2)
                      + 30.0 * DT * T1 * (2.0 * T1^2 + 3.0 * T1 * T2 + 2.0 * T1 * T3 + T2^2 + 2.0 * T2 * T3)
                      + 10.0 * T1 * (7.0 * T1^3
                           + 14.0 * T1^2 * T2
                           + 12.0 * T1^2 * T3
                           + 9.0 * T1 * T2^2
                           + 18.0 * T1 * T2 * T3
                           + 6.0 * T1 * T3^2
                           + 2.0 * T2^3
                           + 6.0 * T2^2 * T3
                           + 6.0 * T2 * T3^2)) / 120.0)
        elseif Stage == 6
            return (Distance_At_Stage(T1, 5)
                 + Cm * DT * T1 * (-DT^3
                      + 2.0 * DT^2 * (T1 + 2.0 * T2)
                      + 2.0 * DT * (11.0 * T1^2 + 15.0 * T1 * T2 + 6.0 * T1 * T3 + 3.0 * T2^2 + 6.0 * T2 * T3)
                      + 49.0 * T1^3
                      + 76.0 * T1^2 * T2
                      + 48.0 * T1^2 * T3
                      + 30.0 * T1 * T2^2
                      + 60.0 * T1 * T2 * T3
                      + 12.0 * T1 * T3^2
                      + 4.0 * T2^3
                      + 12.0 * T2^2 * T3
                      + 12.0 * T2 * T3^2) / 24.0)
        elseif Stage == 7
            return (Distance_At_Stage(T2, 6)
                 + Cm * DT * (DT^4
                      - 5.0 * DT^3 * T1
                      + 10.0 * DT^2 * T1^2
                      + 10.0 * DT * T1 * (11.0 * T1^2 + 18.0 * T1 * T2 + 6.0 * T1 * T3 + 6.0 * T2^2 + 6.0 * T2 * T3)
                      + 5.0 * T1 * (49.0 * T1^3
                           + 120.0 * T1^2 * T2
                           + 48.0 * T1^2 * T3
                           + 96.0 * T1 * T2^2
                           + 84.0 * T1 * T2 * T3
                           + 12.0 * T1 * T3^2
                           + 24.0 * T2^3
                           + 36.0 * T2^2 * T3
                           + 12.0 * T2 * T3^2)) / 120.0)
        elseif Stage == 8
            return (Distance_At_Stage(T1, 7)
                 + Cm * DT * T1 * (DT * (2.0 * T1^2 + 3.0 * T1 * T2 + T1 * T3 + T2^2 + T2 * T3)
                      + 8.0 * T1^3
                      + 16.0 * T1^2 * T2
                      + 6.0 * T1^2 * T3
                      + 10.0 * T1 * T2^2
                      + 9.0 * T1 * T2 * T3
                      + T1 * T3^2
                      + 2.0 * T2^3
                      + 3.0 * T2^2 * T3
                      + T2 * T3^2) / 2.0)
        elseif Stage == 9
            return (Distance_At_Stage(T4, 8)
                 + Cm * DT * (-DT^4
                      + 60.0 * DT * T1 * (2.0 * T1^2 + 3.0 * T1 * T2 + T1 * T3 + T2^2 + T2 * T3)
                      + 60.0 * T1 * (8.0 * T1^3
                           + 16.0 * T1^2 * T2
                           + 6.0 * T1^2 * T3
                           + 4.0 * T1^2 * T4
                           + 10.0 * T1 * T2^2
                           + 9.0 * T1 * T2 * T3
                           + 6.0 * T1 * T2 * T4
                           + T1 * T3^2
                           + 2.0 * T1 * T3 * T4
                           + 2.0 * T2^3
                           + 3.0 * T2^2 * T3
                           + 2.0 * T2^2 * T4
                           + T2 * T3^2
                           + 2.0 * T2 * T3 * T4)) / 120.0)
        elseif Stage == 10
            return (Distance_At_Stage(T1, 9)
                 + Cm * DT * T1 * (-DT^3
                      - 2.0 * DT^2 * T1
                      + 2.0 * DT * (11.0 * T1^2 + 18.0 * T1 * T2 + 6.0 * T1 * T3 + 6.0 * T2^2 + 6.0 * T2 * T3)
                      + 143.0 * T1^3
                      + 264.0 * T1^2 * T2
                      + 96.0 * T1^2 * T3
                      + 48.0 * T1^2 * T4
                      + 144.0 * T1 * T2^2
                      + 132.0 * T1 * T2 * T3
                      + 72.0 * T1 * T2 * T4
                      + 12.0 * T1 * T3^2
                      + 24.0 * T1 * T3 * T4
                      + 24.0 * T2^3
                      + 36.0 * T2^2 * T3
                      + 24.0 * T2^2 * T4
                      + 12.0 * T2 * T3^2
                      + 24.0 * T2 * T3 * T4) / 24.0)
        elseif Stage == 11
            return (Distance_At_Stage(T2, 10)
                 + Cm * DT * (DT^4
                      - 5.0 * DT^3 * T1
                      - 10.0 * DT^2 * T1 * (T1 + 2.0 * T2)
                      + 10.0 * DT * T1 * (11.0 * T1^2 + 15.0 * T1 * T2 + 6.0 * T1 * T3 + 3.0 * T2^2 + 6.0 * T2 * T3)
                      + 5.0 * T1 * (143.0 * T1^3
                           + 308.0 * T1^2 * T2
                           + 96.0 * T1^2 * T3
                           + 48.0 * T1^2 * T4
                           + 210.0 * T1 * T2^2
                           + 156.0 * T1 * T2 * T3
                           + 72.0 * T1 * T2 * T4
                           + 12.0 * T1 * T3^2
                           + 24.0 * T1 * T3 * T4
                           + 44.0 * T2^3
                           + 60.0 * T2^2 * T3
                           + 24.0 * T2^2 * T4
                           + 12.0 * T2 * T3^2
                           + 24.0 * T2 * T3 * T4)) / 120.0)
        elseif Stage == 12
            return (Distance_At_Stage(T1, 11)
                 + Cm * DT * T1 * (-2.0 * DT^2 * (T1 + T2)
                      + 3.0 * DT * (2.0 * T1^2 + 3.0 * T1 * T2 + 2.0 * T1 * T3 + T2^2 + 2.0 * T2 * T3)
                      + 89.0 * T1^3
                      + 178.0 * T1^2 * T2
                      + 60.0 * T1^2 * T3
                      + 24.0 * T1^2 * T4
                      + 111.0 * T1 * T2^2
                      + 90.0 * T1 * T2 * T3
                      + 36.0 * T1 * T2 * T4
                      + 6.0 * T1 * T3^2
                      + 12.0 * T1 * T3 * T4
                      + 22.0 * T2^3
                      + 30.0 * T2^2 * T3
                      + 12.0 * T2^2 * T4
                      + 6.0 * T2 * T3^2
                      + 12.0 * T2 * T3 * T4) / 12.0)
        elseif Stage == 13
            return (Distance_At_Stage(T3, 12)
                 + Cm * DT * (DT^4
                      - 20.0 * DT^2 * T1 * (T1 + T2)
                      + 30.0 * DT * T1 * (2.0 * T1^2 + 3.0 * T1 * T2 + T2^2)
                      + 10.0 * T1 * (89.0 * T1^3
                           + 178.0 * T1^2 * T2
                           + 72.0 * T1^2 * T3
                           + 24.0 * T1^2 * T4
                           + 111.0 * T1 * T2^2
                           + 108.0 * T1 * T2 * T3
                           + 36.0 * T1 * T2 * T4
                           + 12.0 * T1 * T3^2
                           + 12.0 * T1 * T3 * T4
                           + 22.0 * T2^3
                           + 36.0 * T2^2 * T3
                           + 12.0 * T2^2 * T4
                           + 12.0 * T2 * T3^2
                           + 12.0 * T2 * T3 * T4)) / 120.0)
        elseif Stage == 14
            return (Distance_At_Stage(T1, 13)
                 + Cm * DT * T1 * (DT^3
                      - 2.0 * DT^2 * (T1 + 2.0 * T2)
                      + 2.0 * DT * (T1^2 + 3.0 * T1 * T2 + 3.0 * T2^2)
                      + 191.0 * T1^3
                      + 380.0 * T1^2 * T2
                      + 144.0 * T1^2 * T3
                      + 48.0 * T1^2 * T4
                      + 234.0 * T1 * T2^2
                      + 216.0 * T1 * T2 * T3
                      + 72.0 * T1 * T2 * T4
                      + 24.0 * T1 * T3^2
                      + 24.0 * T1 * T3 * T4
                      + 44.0 * T2^3
                      + 72.0 * T2^2 * T3
                      + 24.0 * T2^2 * T4
                      + 24.0 * T2 * T3^2
                      + 24.0 * T2 * T3 * T4) / 24.0)
        elseif Stage == 15
            return (Distance_At_Stage(T2, 14)
                 + Cm * DT * (-DT^4
                      + 5.0 * DT^3 * T1
                      - 10.0 * DT^2 * T1^2
                      + 10.0 * DT * T1^3
                      + 5.0 * T1 * (191.0 * T1^3
                           + 384.0 * T1^2 * T2
                           + 144.0 * T1^2 * T3
                           + 48.0 * T1^2 * T4
                           + 240.0 * T1 * T2^2
                           + 216.0 * T1 * T2 * T3
                           + 72.0 * T1 * T2 * T4
                           + 24.0 * T1 * T3^2
                           + 24.0 * T1 * T3 * T4
                           + 48.0 * T2^3
                           + 72.0 * T2^2 * T3
                           + 24.0 * T2^2 * T4
                           + 24.0 * T2 * T3^2
                           + 24.0 * T2 * T3 * T4)) / 120.0)
        end
    end
    @assert T <= Total_Time(profile)

    if T < T1
        return Distance_At_Stage(T, 1)
    elseif T < T1 + T2
        return Distance_At_Stage(T - (T1), 2)
    elseif T < 2.0 * T1 + T2
        return Distance_At_Stage(T - (T1 + T2), 3)
    elseif T < 2.0 * T1 + T2 + T3
        return Distance_At_Stage(T - (2.0 * T1 + T2), 4)
    elseif T < 3.0 * T1 + T2 + T3
        return Distance_At_Stage(T - (2.0 * T1 + T2 + T3), 5)
    elseif T < 3.0 * T1 + 2.0 * T2 + T3
        return Distance_At_Stage(T - (3.0 * T1 + T2 + T3), 6)
    elseif T < 4.0 * T1 + 2.0 * T2 + T3
        return Distance_At_Stage(T - (3.0 * T1 + 2.0 * T2 + T3), 7)
    elseif T < 4.0 * T1 + 2.0 * T2 + T3 + T4
        return Distance_At_Stage(T - (4.0 * T1 + 2.0 * T2 + T3), 8)
    elseif T < 5.0 * T1 + 2.0 * T2 + T3 + T4
        return Distance_At_Stage(T - (4.0 * T1 + 2.0 * T2 + T3 + T4), 9)
    elseif T < 5.0 * T1 + 3.0 * T2 + T3 + T4
        return Distance_At_Stage(T - (5.0 * T1 + 2.0 * T2 + T3 + T4), 10)
    elseif T < 6.0 * T1 + 3.0 * T2 + T3 + T4
        return Distance_At_Stage(T - (5.0 * T1 + 3.0 * T2 + T3 + T4), 11)
    elseif T < 6.0 * T1 + 3.0 * T2 + 2.0 * T3 + T4
        return Distance_At_Stage(T - (6.0 * T1 + 3.0 * T2 + T3 + T4), 12)
    elseif T < 7.0 * T1 + 3.0 * T2 + 2.0 * T3 + T4
        return Distance_At_Stage(T - (6.0 * T1 + 3.0 * T2 + 2.0 * T3 + T4), 13)
    elseif T < 7.0 * T1 + 4.0 * T2 + 2.0 * T3 + T4
        return Distance_At_Stage(T - (7.0 * T1 + 3.0 * T2 + 2.0 * T3 + T4), 14)
    else
        return Distance_At_Stage(T - (7.0 * T1 + 4.0 * T2 + 2.0 * T3 + T4), 15)
    end
end

"""
    Crackle_At_Time(profile::FeedrateProfile, T, max_crackle)

Returns the crackle at a specific time `T` within a feedrate profile. The crackle will
be either `+max_crackle`, `-max_crackle`, or zero.

The return value may be negative.
"""
function Crackle_At_Time(profile::FeedrateProfile, T, max_crackle)
    @assert T <= Total_Time(profile)

    if T <= Total_Time(profile.accel)
        return Crackle_At_Time(profile.accel, T, max_crackle)
    elseif T < Total_Time(profile.accel) + profile.coast
        return 0.0
    else
        return Crackle_At_Time(profile.decel, T - (Total_Time(profile.accel) + profile.coast), -max_crackle)
    end
end

"""
    Snap_At_Time(profile::FeedrateProfile, T, max_crackle)

Returns the snap at a specific time `T` within a feedrate profile. The return value may
be negative.
"""
function Snap_At_Time(profile::FeedrateProfile, T, max_crackle)
    @assert T <= Total_Time(profile)

    if T <= Total_Time(profile.accel)
        return Snap_At_Time(profile.accel, T, max_crackle)
    elseif T < Total_Time(profile.accel) + profile.coast
        return 0.0
    else
        return Snap_At_Time(profile.decel, T - (Total_Time(profile.accel) + profile.coast), -max_crackle)
    end
end

"""
    Jerk_At_Time(profile::FeedrateProfile, T, max_crackle)

Returns the jerk at a specific time `T` within a feedrate profile. The return value may
be negative.
"""
function Jerk_At_Time(profile::FeedrateProfile, T, max_crackle)
    @assert T <= Total_Time(profile)

    if T <= Total_Time(profile.accel)
        return Jerk_At_Time(profile.accel, T, max_crackle)
    elseif T < Total_Time(profile.accel) + profile.coast
        return 0.0
    else
        return Jerk_At_Time(profile.decel, T - (Total_Time(profile.accel) + profile.coast), -max_crackle)
    end
end

"""
    Acceleration_At_Time(profile::FeedrateProfile, T, max_crackle)

Returns the acceleration at a specific time `T` within a feedrate profile. The return
value may be negative.
"""
function Acceleration_At_Time(profile::FeedrateProfile, T, max_crackle)
    @assert T <= Total_Time(profile)

    if T <= Total_Time(profile.accel)
        return Acceleration_At_Time(profile.accel, T, max_crackle)
    elseif T < Total_Time(profile.accel) + profile.coast
        return 0.0
    else
        return Acceleration_At_Time(profile.decel, T - (Total_Time(profile.accel) + profile.coast), -max_crackle)
    end
end

"""
    Velocity_At_Time(profile::FeedrateProfile, T, max_crackle, start_vel)

Returns the velocity at a specific time `T` within a feedrate profile. The return value
may be negative.
"""
function Velocity_At_Time(profile::FeedrateProfile, T, max_crackle, start_vel)
    mid_vel = Velocity_At_Time(profile.accel, Total_Time(profile.accel), max_crackle, start_vel)

    @assert T <= Total_Time(profile)

    if T <= Total_Time(profile.accel)
        return Velocity_At_Time(profile.accel, T, max_crackle, start_vel)
    elseif T < Total_Time(profile.accel) + profile.coast
        return mid_vel
    else
        return Velocity_At_Time(profile.decel, T - (Total_Time(profile.accel) + profile.coast), -max_crackle, mid_vel)
    end
end

"""
    Distance_At_Time(profile::FeedrateProfile, T, max_crackle, start_vel)

Returns the distance from the start point at a specific time `T` within a feedrate
profile. The return value may be negative.
"""
function Distance_At_Time(profile::FeedrateProfile, T, max_crackle, start_vel)
    mid_vel = Velocity_At_Time(profile.accel, Total_Time(profile.accel), max_crackle, start_vel)
    accel_dist = Distance_At_Time(profile.accel, Total_Time(profile.accel), max_crackle, start_vel)
    mid_dist = mid_vel * profile.coast

    @assert T <= Total_Time(profile)

    if T <= Total_Time(profile.accel)
        return Distance_At_Time(profile.accel, T, max_crackle, start_vel)
    elseif T < Total_Time(profile.accel) + profile.coast
        return accel_dist + mid_vel * (T - Total_Time(profile.accel))
    else
        return accel_dist + mid_dist +
               Distance_At_Time(profile.decel, T - (Total_Time(profile.accel) + profile.coast), -max_crackle, mid_vel)
    end
end

"""
    Distance_At_Time_With_Accel_Flag(profile::FeedrateProfile, T, max_crackle, start_vel)

Returns the distance from the start point at a specific time `T` within a feedrate
profile, along with a boolean indicating whether `T` is past the acceleration part.

Returns a tuple `(distance, is_past_accel_part)` where `is_past_accel_part` is `true` if
`T` is in the coasting or deceleration part, otherwise `false`.

The distance value may be negative.
"""
function Distance_At_Time_With_Accel_Flag(profile::FeedrateProfile, T, max_crackle, start_vel)
    mid_vel = Velocity_At_Time(profile.accel, Total_Time(profile.accel), max_crackle, start_vel)
    accel_dist = Distance_At_Time(profile.accel, Total_Time(profile.accel), max_crackle, start_vel)
    mid_dist = mid_vel * profile.coast

    @assert T <= Total_Time(profile)

    if T <= Total_Time(profile.accel)
        return (Distance_At_Time(profile.accel, T, max_crackle, start_vel), false)
    elseif T < Total_Time(profile.accel) + profile.coast
        return (accel_dist + mid_vel * (T - Total_Time(profile.accel)), true)
    else
        return (accel_dist + mid_dist +
                Distance_At_Time(profile.decel, T - (Total_Time(profile.accel) + profile.coast), -max_crackle, mid_vel),
                true)
    end
end

"""
    Optimal_Profile_For_Distance(start_vel, distance, acceleration_max, jerk_max, snap_max, crackle_max)

Compute the acceleration part of a feedrate profile that has the lowest total time to
travel the given distance without violating any of the given constraints. Note that there
is no velocity limit here.
"""
function Optimal_Profile_For_Distance(start_vel, distance, acceleration_max, jerk_max, snap_max, crackle_max)
    D = distance
    Vs = start_vel
    Am = acceleration_max
    Jm = jerk_max
    Sm = snap_max
    Cm = crackle_max

    function solve_distance_at_time(profile::FeedrateProfileTimes, variable::Int)
        lower = 0.0
        upper = 86400.0  # 24 hours should be more than enough

        while true
            mid = lower + (upper - lower) / 2.0
            mid = reinterpret(Float64, reinterpret(UInt64, lower) + (reinterpret(UInt64, upper) - reinterpret(UInt64, lower)) ÷ 2)

            if lower == mid || upper == mid
                break
            end

            result_arr = [profile.t1, profile.t2, profile.t3, profile.t4]
            result_arr[variable] = mid
            test_profile = FeedrateProfileTimes(result_arr[1], result_arr[2], result_arr[3], result_arr[4])

            if Fast_Distance_At_Max_Time(test_profile, Cm, Vs) <= D
                lower = mid
            else
                upper = mid
            end
        end

        result_arr = [profile.t1, profile.t2, profile.t3, profile.t4]
        result_arr[variable] = lower
        return FeedrateProfileTimes(result_arr[1], result_arr[2], result_arr[3], result_arr[4])
    end

    # Build the cases array based on kinematic constraints
    if Sm^2 < Jm * Cm
        if Am >= Jm * (Jm / Sm + Sm / Cm)
            cases = (
                # Case 1: Reachable: None
                FeedrateProfileTimes(0.0, 0.0, 0.0, 0.0),
                # Case 2: Reachable: Sm
                FeedrateProfileTimes(Sm / Cm, 0.0, 0.0, 0.0),
                # Case 3: Reachable: Sm, Jm
                FeedrateProfileTimes(Sm / Cm, Jm / Sm - Sm / Cm, 0.0, 0.0),
                # Case 4: Reachable: Sm, Jm, Am
                FeedrateProfileTimes(Sm / Cm, Jm / Sm - Sm / Cm, Am / Jm - Jm / Sm - Sm / Cm, 0.0),
            )
        elseif Am >= 2.0 * Sm^3 / Cm^2
            cases = (
                # Case 1: Reachable: None
                FeedrateProfileTimes(0.0, 0.0, 0.0, 0.0),
                # Case 2: Reachable: Sm
                FeedrateProfileTimes(Sm / Cm, 0.0, 0.0, 0.0),
                # Case 3: Impossible case
                FeedrateProfileTimes(Sm / Cm, (0.25 * Sm^2 / Cm^2 + Am / Sm)^(1/2) - 1.5 * Sm / Cm, 0.0, 0.0),
                # Case 4: Reachable: Sm, Am
                FeedrateProfileTimes(Sm / Cm, (0.25 * Sm^2 / Cm^2 + Am / Sm)^(1/2) - 1.5 * Sm / Cm, 0.0, 0.0),
            )
        else
            cases = (
                # Case 1: Reachable: None
                FeedrateProfileTimes(0.0, 0.0, 0.0, 0.0),
                # Case 2: Impossible case
                FeedrateProfileTimes((0.5 * Am / Cm)^(1/3), 0.0, 0.0, 0.0),
                # Case 3: Impossible case
                FeedrateProfileTimes((0.5 * Am / Cm)^(1/3), 0.0, 0.0, 0.0),
                # Case 4: Reachable: Am
                FeedrateProfileTimes((0.5 * Am / Cm)^(1/3), 0.0, 0.0, 0.0),
            )
        end
    else
        if Am > 2.0 * Jm * (Jm / Cm)^(1/2)
            cases = (
                # Case 1: Reachable: None
                FeedrateProfileTimes(0.0, 0.0, 0.0, 0.0),
                # Case 2: Impossible case
                FeedrateProfileTimes((Jm / Cm)^(1/2), 0.0, 0.0, 0.0),
                # Case 3: Reachable: Jm
                FeedrateProfileTimes((Jm / Cm)^(1/2), 0.0, 0.0, 0.0),
                # Case 4: Reachable: Jm, Am
                FeedrateProfileTimes((Jm / Cm)^(1/2), 0.0, Am / Jm - 2.0 * (Jm / Cm)^(1/2), 0.0),
            )
        else
            cases = (
                # Case 1: Reachable: None
                FeedrateProfileTimes(0.0, 0.0, 0.0, 0.0),
                # Case 2: Impossible case
                FeedrateProfileTimes((Am / (2.0 * Cm))^(1/3), 0.0, 0.0, 0.0),
                # Case 3: Impossible case
                FeedrateProfileTimes((Am / (2.0 * Cm))^(1/3), 0.0, 0.0, 0.0),
                # Case 4: Reachable: Am
                FeedrateProfileTimes((Am / (2.0 * Cm))^(1/3), 0.0, 0.0, 0.0),
            )
        end
    end

    # Find the appropriate case and solve
    for i in 4:-1:1
        if i == 1 || D > Fast_Distance_At_Max_Time(cases[i], Cm, Vs)
            return solve_distance_at_time(cases[i], i)
        end
    end

    # Unreachable
    error("Unreachable state in Optimal_Profile_For_Distance")
end

"""
    Optimal_Profile_For_Delta_V(delta_v, acceleration_max, jerk_max, snap_max, crackle_max)

Compute the acceleration part of a feedrate profile that achieves the given change in
velocity in the lowest time without violating any of the given constraints. Note that
there is no distance limit here.
"""
function Optimal_Profile_For_Delta_V(delta_v, acceleration_max, jerk_max, snap_max, crackle_max)
    Vd = abs(delta_v)
    Am = acceleration_max
    Jm = jerk_max
    Sm = snap_max
    Cm = crackle_max

    function solve_velocity_at_time(profile::FeedrateProfileTimes, variable::Int, target)
        lower = 0.0
        upper = 86400.0  # 24 hours should be more than enough

        while true
            mid = reinterpret(Float64, reinterpret(UInt64, lower) + (reinterpret(UInt64, upper) - reinterpret(UInt64, lower)) ÷ 2)

            if lower == mid || upper == mid
                break
            end

            result_arr = [profile.t1, profile.t2, profile.t3, profile.t4]
            result_arr[variable] = mid
            test_profile = FeedrateProfileTimes(result_arr[1], result_arr[2], result_arr[3], result_arr[4])

            if Fast_Velocity_At_Max_Time(test_profile, Cm, 0.0) <= target
                lower = mid
            else
                upper = mid
            end
        end

        result_arr = [profile.t1, profile.t2, profile.t3, profile.t4]
        result_arr[variable] = lower
        return FeedrateProfileTimes(result_arr[1], result_arr[2], result_arr[3], result_arr[4])
    end

    if Sm^2 < Jm * Cm
        if Am >= Jm * (Jm / Sm + Sm / Cm)
            if Vd > Am * (Am / Jm + Jm / Sm + Sm / Cm)
                # Reachable: Sm, Jm, Am
                return FeedrateProfileTimes(Sm / Cm, Jm / Sm - Sm / Cm, Am / Jm - Jm / Sm - Sm / Cm,
                                           Vd / Am - Am / Jm - Jm / Sm - Sm / Cm)
            elseif Vd > 2.0 * Jm * (Jm / Sm + Sm / Cm)^2
                # Reachable: Sm, Jm
                return FeedrateProfileTimes(Sm / Cm, Jm / Sm - Sm / Cm,
                                           0.5 * ((Jm / Sm + Sm / Cm)^2 + 4.0 * Vd / Jm)^(1/2) - 1.5 * (Jm / Sm + Sm / Cm),
                                           0.0)
            elseif Vd > 8.0 * Sm^4 / Cm^3
                # Reachable: Sm
                return solve_velocity_at_time(FeedrateProfileTimes(Sm / Cm, 0.0, 0.0, 0.0), 2, Vd)
            else
                # Reachable: None
                return FeedrateProfileTimes((0.125 * Vd / Cm)^(1/4), 0.0, 0.0, 0.0)
            end
        elseif Am >= 2.0 * Sm^3 / Cm^2
            if Vd > Am * (2.0 * (0.25 * Sm^2 / Cm^2 + Am / Sm)^(1/2) + Sm / Cm)
                # Reachable: Sm, Am
                return FeedrateProfileTimes(Sm / Cm,
                                           (0.25 * Sm^2 / Cm^2 + Am / Sm)^(1/2) - 1.5 * Sm / Cm,
                                           0.0,
                                           Vd / Am - Sm / Cm - 2.0 * (0.25 * Sm^2 / Cm^2 + Am / Sm)^(1/2))
            elseif Vd > 8.0 * Sm^4 / Cm^3
                # Reachable: Sm
                return solve_velocity_at_time(FeedrateProfileTimes(Sm / Cm, 0.0, 0.0, 0.0), 2, Vd)
            else
                # Reachable: None
                return FeedrateProfileTimes((0.125 * Vd / Cm)^(1/4), 0.0, 0.0, 0.0)
            end
        else
            if Vd > 8.0 * Cm * (0.5 * Am / Cm)^(4/3)
                # Reachable: Am
                return FeedrateProfileTimes((0.5 * Am / Cm)^(1/3), 0.0, 0.0, Vd / Am - 4.0 * (0.5 * Am / Cm)^(1/3))
            else
                # Reachable: None
                return FeedrateProfileTimes((0.125 * Vd / Cm)^(1/4), 0.0, 0.0, 0.0)
            end
        end
    else
        if Am > 2.0 * Jm * (Jm / Cm)^(1/2)
            if Vd > Am * (Am / Jm + 2.0 * (Jm / Cm)^(1/2))
                # Reachable: Jm, Am
                return FeedrateProfileTimes((Jm / Cm)^(1/2), 0.0,
                                           Am / Jm - 2.0 * (Jm / Cm)^(1/2),
                                           Vd / Am - Am / Jm - 2.0 * (Jm / Cm)^(1/2))
            elseif Vd > 8.0 * Jm^2 / Cm
                # Reachable: Jm
                return FeedrateProfileTimes((Jm / Cm)^(1/2), 0.0,
                                           (Jm / Cm + Vd / Jm)^(1/2) - 3.0 * (Jm / Cm)^(1/2),
                                           0.0)
            else
                # Reachable: None
                return FeedrateProfileTimes((0.125 * Vd / Cm)^(1/4), 0.0, 0.0, 0.0)
            end
        else
            if Vd > 8.0 * Cm * (0.5 * Am / Cm)^(4/3)
                # Reachable: Am
                return FeedrateProfileTimes((0.5 * Am / Cm)^(1/3), 0.0, 0.0, Vd / Am - 4.0 * (0.5 * Am / Cm)^(1/3))
            else
                # Reachable: None
                return FeedrateProfileTimes((0.125 * Vd / Cm)^(1/4), 0.0, 0.0, 0.0)
            end
        end
    end
end

"""
    Optimal_Full_Profile(start_vel, max_vel, end_vel, distance, acceleration_max, jerk_max, snap_max, crackle_max)

Compute the feedrate profile with the minimal time without violating the given
constraints. Throws an error if there is no legal feedrate profile which can meet the
given constraints, specifically regarding `end_vel` being reachable. Also throws an error
if `start_vel` or `end_vel` are higher than `max_vel`.
"""
function Optimal_Full_Profile(start_vel, max_vel, end_vel, distance, acceleration_max, jerk_max, snap_max, crackle_max)
    if max_vel < start_vel
        error("max_vel cannot be smaller than start_vel")
    end

    if max_vel < end_vel
        error("max_vel cannot be smaller than end_vel")
    end

    if distance == 0.0
        return FeedrateProfile(FeedrateProfileTimes(), 0.0, FeedrateProfileTimes())
    end

    # Check if end_vel is reachable
    check_profile = Optimal_Profile_For_Delta_V(start_vel - end_vel, acceleration_max, jerk_max, snap_max, crackle_max)
    check_crackle = start_vel < end_vel ? crackle_max : -crackle_max
    profile_distance = Fast_Distance_At_Max_Time(check_profile, check_crackle, start_vel)

    if distance < profile_distance
        error("end_vel is not reachable under given constraints")
    end

    accel = Optimal_Profile_For_Delta_V(start_vel - max_vel, acceleration_max, jerk_max, snap_max, crackle_max)
    decel = Optimal_Profile_For_Delta_V(end_vel - max_vel, acceleration_max, jerk_max, snap_max, crackle_max)

    accel_distance = Fast_Distance_At_Max_Time(accel, crackle_max, start_vel)
    decel_distance = Fast_Distance_At_Max_Time(decel, -crackle_max, max_vel)

    if accel_distance + decel_distance <= distance
        coast = (distance - accel_distance - decel_distance) / max_vel
        return FeedrateProfile(accel, coast, decel)
    else
        # Binary search for the optimal mid velocity
        coast = 0.0
        upper = max_vel
        lower = max(start_vel, end_vel)

        while true
            mid = reinterpret(Float64, reinterpret(UInt64, lower) + (reinterpret(UInt64, upper) - reinterpret(UInt64, lower)) ÷ 2)

            if lower == mid || upper == mid
                break
            end

            accel = Optimal_Profile_For_Delta_V(start_vel - mid, acceleration_max, jerk_max, snap_max, crackle_max)
            decel = Optimal_Profile_For_Delta_V(end_vel - mid, acceleration_max, jerk_max, snap_max, crackle_max)

            accel_distance = Fast_Distance_At_Max_Time(accel, crackle_max, start_vel)
            decel_distance = Fast_Distance_At_Max_Time(decel, crackle_max, end_vel)

            if accel_distance + decel_distance <= distance
                lower = mid
            else
                upper = mid
            end
        end

        return FeedrateProfile(accel, coast, decel)
    end
end

end
