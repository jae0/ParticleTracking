using Test
using ParticleTracking

# Reduction of stored Oceananigans arrays to interior cell centres. The archive layout that
# exposed the original fault: interior 345x245x13, halo (7, 7, 5), u on Bounded x-faces.
@testset "Interior-centre reduction of haloed, staggered archives" begin
    ic = ParticleTracking._interior_centers
    N, H = (5, 4, 3), (2, 2, 1)

    # Centred field with halo: interior values are recovered exactly, halo is discarded.
    T_int = reshape(collect(1.0:prod(N)), N)
    T_raw = fill(NaN, N .+ 2 .* H)
    T_raw[3:7, 3:6, 2:4] = T_int
    @test ic(T_raw, N, H) == T_int

    # x-face field with halo: (Nx + 1 + 2Hx) faces, averaged onto centres.
    faces = collect(0.0:N[1])                                # face i at x = i - 1
    u_raw = fill(NaN, N[1] + 1 + 2H[1], N[2] + 2H[2], N[3] + 2H[3])
    for i in 1:(N[1] + 1)
        u_raw[H[1] + i, 3:6, 2:4] .= faces[i]
    end
    u_c = ic(u_raw, N, H)
    @test size(u_c) == N
    @test u_c[:, 1, 1] ≈ faces[1:N[1]] .+ 0.5

    # Halo-free layouts used by older archives remain supported.
    @test ic(T_int, N, H) == T_int
    @test size(ic(zeros(N[1] + 1, N[2], N[3]), N, H)) == N

    # 2-D (free surface) path.
    η_raw = fill(NaN, N[1] + 2H[1], N[2] + 2H[2])
    η_raw[3:7, 3:6] .= 1.0
    @test ic(η_raw, (N[1], N[2]), (H[1], H[2])) == ones(N[1], N[2])

    # An unrecognised length is an error, not a truncation.
    @test_throws ErrorException ic(zeros(N[1] + 3, N[2], N[3]), N, H)
    @test_throws ErrorException ic(zeros(N[1], N[2]), N, H)
end
