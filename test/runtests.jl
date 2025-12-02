using PruntTrajectories
using Test
using Random

@testset "PruntTrajectories.jl" begin
    @testset "FeedrateProfileTimes" begin
        p = FeedrateProfileTimes(0.1, 0.2, 0.3, 0.4)
        @test p.t1 == 0.1
        @test p.t2 == 0.2
        @test p.t3 == 0.3
        @test p.t4 == 0.4

        p_default = FeedrateProfileTimes()
        @test p_default.t1 == 0.0
        @test p_default.t2 == 0.0
        @test p_default.t3 == 0.0
        @test p_default.t4 == 0.0
    end

    @testset "total_time" begin
        p = FeedrateProfileTimes(0.1, 0.2, 0.3, 0.4)
        @test total_time(p) == 8.0 * 0.1 + 4.0 * 0.2 + 2.0 * 0.3 + 0.4
        @test total_time(p) == 2.6

        profile = FeedrateProfile(p, 1.0, p)
        @test total_time(profile) == 2.6 + 1.0 + 2.6
    end

    @testset "fast_velocity_at_max_time" begin
        p = FeedrateProfileTimes(0.1, 0.2, 0.3, 0.4)
        v = fast_velocity_at_max_time(p, 1000.0, 10.0)
        @test v == 41.5
    end

    @testset "fast_distance_at_max_time" begin
        p = FeedrateProfileTimes(0.1, 0.2, 0.3, 0.4)
        d = fast_distance_at_max_time(p, 1000.0, 10.0)
        @test d == 66.95
    end

    @testset "crackle_at_time (FeedrateProfileTimes)" begin
        p = FeedrateProfileTimes(0.1, 0.2, 0.3, 0.4)
        Cm = 1000.0

        # Stage 1: T < T1 (0.1) => Cm
        @test crackle_at_time(p, 0.05, Cm) == Cm
        # Stage 2: T1 <= T < T1+T2 (0.3) => 0.0
        @test crackle_at_time(p, 0.15, Cm) == 0.0
        @test crackle_at_time(p, 0.25, Cm) == 0.0
        # Stage 3: T1+T2 <= T < 2*T1+T2 (0.5) => -Cm
        @test crackle_at_time(p, 0.35, Cm) == -Cm
    end

    @testset "snap_at_time (FeedrateProfileTimes)" begin
        p = FeedrateProfileTimes(0.1, 0.2, 0.3, 0.4)
        Cm = 1000.0

        snap = snap_at_time(p, 0.05, Cm)
        @test snap == Cm * 0.05
    end

    @testset "jerk_at_time (FeedrateProfileTimes)" begin
        p = FeedrateProfileTimes(0.1, 0.2, 0.3, 0.4)
        Cm = 1000.0

        jerk = jerk_at_time(p, 0.05, Cm)
        @test jerk ≈ Cm * 0.05^2 / 2.0
    end

    @testset "acceleration_at_time (FeedrateProfileTimes)" begin
        p = FeedrateProfileTimes(0.1, 0.2, 0.3, 0.4)
        Cm = 1000.0

        accel = acceleration_at_time(p, 0.05, Cm)
        @test accel ≈ Cm * 0.05^3 / 6.0
    end

    @testset "velocity_at_time (FeedrateProfileTimes)" begin
        p = FeedrateProfileTimes(0.1, 0.2, 0.3, 0.4)
        Cm = 1000.0
        start_vel = 10.0

        vel = velocity_at_time(p, 0.05, Cm, start_vel)
        @test vel ≈ start_vel + Cm * 0.05^4 / 24.0
    end

    @testset "distance_at_time (FeedrateProfileTimes)" begin
        p = FeedrateProfileTimes(0.1, 0.2, 0.3, 0.4)
        Cm = 1000.0
        start_vel = 10.0

        dist = distance_at_time(p, 0.05, Cm, start_vel)
        @test dist ≈ start_vel * 0.05 + Cm * 0.05^5 / 120.0
    end

    @testset "FeedrateProfile At_Time functions" begin
        times = FeedrateProfileTimes(0.1, 0.2, 0.3, 0.4)
        profile = FeedrateProfile(times, 1.0, times)
        Cm = 1000.0
        start_vel = 10.0

        @test crackle_at_time(profile, 0.05, Cm) == Cm
        @test snap_at_time(profile, 0.05, Cm) == Cm * 0.05
        @test jerk_at_time(profile, 0.05, Cm) ≈ Cm * 0.05^2 / 2.0
        @test acceleration_at_time(profile, 0.05, Cm) ≈ Cm * 0.05^3 / 6.0

        @test crackle_at_time(profile, 3.0, Cm) == 0.0
        @test snap_at_time(profile, 3.0, Cm) == 0.0
        @test jerk_at_time(profile, 3.0, Cm) == 0.0
        @test acceleration_at_time(profile, 3.0, Cm) == 0.0
    end

    @testset "optimal_profile_for_delta_v" begin
        profile = optimal_profile_for_delta_v(100.0, 50.0, 100.0, 200.0, 500.0)
        @test profile isa FeedrateProfileTimes
        @test profile.t1 > 0.0
        @test profile.t4 > 0.0

        achieved_delta_v = fast_velocity_at_max_time(profile, 500.0, 0.0)
        @test achieved_delta_v ≈ 100.0 atol=1e-10
    end

    @testset "optimal_profile_for_distance" begin
        profile = optimal_profile_for_distance(0.0, 1000.0, 50.0, 100.0, 200.0, 500.0)
        @test profile isa FeedrateProfileTimes

        achieved_distance = fast_distance_at_max_time(profile, 500.0, 0.0)
        @test achieved_distance ≈ 1000.0 atol=1e-10
    end

    @testset "optimal_full_profile" begin
        full_profile = optimal_full_profile(0.0, 100.0, 0.0, 1000.0, 50.0, 100.0, 200.0, 500.0)
        @test full_profile isa FeedrateProfile
        @test full_profile.coast >= 0.0

        @test_throws ErrorException optimal_full_profile(150.0, 100.0, 0.0, 1000.0, 50.0, 100.0, 200.0, 500.0)
        @test_throws ErrorException optimal_full_profile(0.0, 100.0, 150.0, 1000.0, 50.0, 100.0, 200.0, 500.0)

        zero_profile = optimal_full_profile(0.0, 100.0, 0.0, 0.0, 50.0, 100.0, 200.0, 500.0)
        @test total_time(zero_profile) == 0.0
    end

    @testset "distance_at_time_with_accel_flag" begin
        times = FeedrateProfileTimes(0.1, 0.2, 0.3, 0.4)
        profile = FeedrateProfile(times, 1.0, times)
        Cm = 1000.0
        start_vel = 10.0

        dist1, flag1 = distance_at_time_with_accel_flag(profile, 0.5, Cm, start_vel)
        @test flag1 == false

        dist2, flag2 = distance_at_time_with_accel_flag(profile, 3.0, Cm, start_vel)
        @test flag2 == true
    end

    @testset "Consistency between Fast_* and *_At_Time" begin
        p = FeedrateProfileTimes(0.1, 0.2, 0.3, 0.4)
        Cm = 1000.0
        start_vel = 10.0
        t_max = total_time(p)

        fast_vel = fast_velocity_at_max_time(p, Cm, start_vel)
        slow_vel = velocity_at_time(p, t_max, Cm, start_vel)
        @test fast_vel ≈ slow_vel

        fast_dist = fast_distance_at_max_time(p, Cm, start_vel)
        slow_dist = distance_at_time(p, t_max, Cm, start_vel)
        @test fast_dist ≈ slow_dist
    end

    @testset "Derivative integration consistency" begin
        # Trapezoidal integration from 0 to t_end
        function integrate_trapezoidal(f, t_end, n)
            h = t_end / n
            result = 0.5 * (f(0.0) + f(t_end))
            for i in 1:(n-1)
                result += f(i * h)
            end
            return result * h
        end

        # Test derivative integration consistency for a given profile
        function test_derivative_integration(p::FeedrateProfileTimes, Cm, start_vel; n_points=10000)
            t_max = total_time(p)
            if t_max == 0
                return true
            end

            # Test points as fractions of t_max (avoid exact boundaries and zero-crossings)
            test_fracs = [0.02, 0.04, 0.15, 0.25]

            # Snap is integral of Crackle
            for frac in test_fracs
                t_test = frac * t_max
                n_test = max(100, round(Int, n_points * frac))
                integrated = snap_at_time(p, 0.0, Cm) +
                             integrate_trapezoidal(t -> crackle_at_time(p, t, Cm), t_test, n_test)
                analytical = snap_at_time(p, t_test, Cm)
                @test isapprox(integrated, analytical, rtol=0.02, atol=Cm * t_max * 0.01)
            end

            # Jerk is integral of Snap
            for frac in test_fracs
                t_test = frac * t_max
                n_test = max(100, round(Int, n_points * frac))
                integrated = jerk_at_time(p, 0.0, Cm) +
                             integrate_trapezoidal(t -> snap_at_time(p, t, Cm), t_test, n_test)
                analytical = jerk_at_time(p, t_test, Cm)
                @test isapprox(integrated, analytical, rtol=0.02, atol=Cm * t_max^2 * 0.01)
            end

            # Acceleration is integral of Jerk
            for frac in test_fracs
                t_test = frac * t_max
                n_test = max(100, round(Int, n_points * frac))
                integrated = acceleration_at_time(p, 0.0, Cm) +
                             integrate_trapezoidal(t -> jerk_at_time(p, t, Cm), t_test, n_test)
                analytical = acceleration_at_time(p, t_test, Cm)
                @test isapprox(integrated, analytical, rtol=0.02, atol=Cm * t_max^3 * 0.01)
            end

            # Velocity is integral of Acceleration
            for frac in [0.25, 0.5, 0.75, 1.0]
                t_test = frac * t_max
                n_test = max(100, round(Int, n_points * frac))
                integrated = velocity_at_time(p, 0.0, Cm, start_vel) +
                             integrate_trapezoidal(t -> acceleration_at_time(p, t, Cm), t_test, n_test)
                analytical = velocity_at_time(p, t_test, Cm, start_vel)
                @test isapprox(integrated, analytical, rtol=0.01, atol=max(1.0, abs(start_vel) * 0.01))
            end

            # Distance is integral of Velocity
            for frac in [0.25, 0.5, 0.75, 1.0]
                t_test = frac * t_max
                n_test = max(100, round(Int, n_points * frac))
                integrated = distance_at_time(p, 0.0, Cm, start_vel) +
                             integrate_trapezoidal(t -> velocity_at_time(p, t, Cm, start_vel), t_test, n_test)
                analytical = distance_at_time(p, t_test, Cm, start_vel)
                @test isapprox(integrated, analytical, rtol=0.01, atol=max(1.0, abs(start_vel) * t_max * 0.01))
            end
        end

        # Test full profile derivative integration
        function test_full_profile_integration(profile::FeedrateProfile, Cm, start_vel; n_points=10000)
            t_max = total_time(profile)
            if t_max == 0
                return true
            end

            for frac in [0.25, 0.5, 0.75, 0.99]  # Use 0.99 instead of 1.0 to avoid floating point boundary issues
                t_test = frac * t_max
                n_test = max(100, round(Int, n_points * frac))
                integrated = distance_at_time(profile, 0.0, Cm, start_vel) +
                             integrate_trapezoidal(t -> velocity_at_time(profile, t, Cm, start_vel), t_test, n_test)
                analytical = distance_at_time(profile, t_test, Cm, start_vel)
                @test isapprox(integrated, analytical, rtol=0.01, atol=max(1.0, abs(start_vel) * t_max * 0.01))
            end
        end

        # Test with fixed profile
        @testset "Fixed profile" begin
            p = FeedrateProfileTimes(0.1, 0.2, 0.3, 0.4)
            test_derivative_integration(p, 1000.0, 10.0)
        end

        # Test with 100 random profiles
        @testset "Random profiles (n=100)" begin
            rng = Random.MersenneTwister(42)  # Fixed seed for reproducibility

            for i in 1:100
                # Random profile times (each between 0.01 and 1.0)
                t1 = 0.01 + 0.99 * rand(rng)
                t2 = 0.01 + 0.99 * rand(rng)
                t3 = 0.01 + 0.99 * rand(rng)
                t4 = 0.01 + 0.99 * rand(rng)
                p = FeedrateProfileTimes(t1, t2, t3, t4)

                # Random parameters
                Cm = 100.0 + 900.0 * rand(rng)  # 100 to 1000
                start_vel = 50.0 * rand(rng)     # 0 to 50

                @testset "Profile $i" begin
                    test_derivative_integration(p, Cm, start_vel)
                end
            end
        end

        # Test that velocity and distance constraints are met
        function test_velocity_and_distance(profile::FeedrateProfile, Cm, start_vel, max_vel, end_vel, distance; n_samples=1000)
            t_max = total_time(profile)
            if t_max == 0
                @test distance == 0.0
                return true
            end

            # Check start velocity
            v_start = velocity_at_time(profile, 0.0, Cm, start_vel)
            @test isapprox(v_start, start_vel, rtol=1e-10, atol=1e-10)

            # Check end velocity
            v_end = velocity_at_time(profile, t_max * 0.9999, Cm, start_vel)
            @test isapprox(v_end, end_vel, rtol=0.01, atol=0.1)

            # Check total distance
            d_total = distance_at_time(profile, t_max * 0.9999, Cm, start_vel)
            @test isapprox(d_total, distance, rtol=0.01, atol=0.1)

            # Check velocity never exceeds max_vel (sample throughout trajectory)
            for i in 0:n_samples
                t = (i / n_samples) * t_max * 0.9999
                v = velocity_at_time(profile, t, Cm, start_vel)
                @test v <= max_vel + max_vel * 0.01 + 1e-9
                @test v >= 0.0 - 1e-9  # Velocity should be non-negative
            end
        end

        # Test that kinematic limits are respected throughout a full profile
        function test_kinematic_limits(profile::FeedrateProfile, Cm, Am, Jm, Sm; n_samples=1000)
            t_max = total_time(profile)
            if t_max == 0
                return true
            end

            # Sample the trajectory and check limits
            for i in 0:n_samples
                t = (i / n_samples) * t_max * 0.9999  # Avoid exact endpoint

                crackle = abs(crackle_at_time(profile, t, Cm))
                snap = abs(snap_at_time(profile, t, Cm))
                jerk = abs(jerk_at_time(profile, t, Cm))
                accel = abs(acceleration_at_time(profile, t, Cm))

                # Allow small tolerance for numerical precision
                tol = 1e-9
                @test crackle <= Cm + tol
                @test snap <= Sm + Sm * 0.01 + tol
                @test jerk <= Jm + Jm * 0.01 + tol
                @test accel <= Am + Am * 0.01 + tol
            end
        end

        # Test full profiles with random parameters
        @testset "Random full profiles (n=100)" begin
            rng = Random.MersenneTwister(123)  # Different seed

            for i in 1:100
                # Random kinematic limits
                Am = 10.0 + 90.0 * rand(rng)   # 10 to 100
                Jm = 50.0 + 150.0 * rand(rng)  # 50 to 200
                Sm = 100.0 + 200.0 * rand(rng) # 100 to 300
                Cm = 200.0 + 500.0 * rand(rng) # 200 to 700

                # Random velocities and distance
                max_vel = 50.0 + 100.0 * rand(rng)  # 50 to 150
                start_vel = max_vel * rand(rng)      # 0 to max_vel
                end_vel = max_vel * rand(rng)        # 0 to max_vel
                distance = 100.0 + 900.0 * rand(rng) # 100 to 1000

                # Try to create profile (may fail for some parameter combinations)
                try
                    profile = optimal_full_profile(start_vel, max_vel, end_vel, distance, Am, Jm, Sm, Cm)
                    @testset "Full profile $i" begin
                        test_full_profile_integration(profile, Cm, start_vel)
                        test_kinematic_limits(profile, Cm, Am, Jm, Sm)
                        test_velocity_and_distance(profile, Cm, start_vel, max_vel, end_vel, distance)
                    end
                catch e
                    # Skip invalid parameter combinations
                    if !(e isa ErrorException)
                        rethrow(e)
                    end
                end
            end
        end
    end
end
