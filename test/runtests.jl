using Delaunay
using LinearAlgebra
using SparseArrays
using Test



const Dmax = 5



# TODO: Move this into the library?
"""
Convert a vector-of-vectors into a matrix
"""
function flatten(A::AbstractVector)::Matrix
    N = length(A)
    @assert N > 0
    D = length(A[1])
    @assert D > 0
    T = typeof(A[1][1])
    B = Array{T}(undef, N, D)
    for n = 1:N
        @assert length(A[n]) == D
        for d = 1:D
            B[n, d] = A[n][d]
        end
    end
    return B
end

"""
Volume of a simplex

See <https://mathworld.wolfram.com/Cayley-MengerDeterminant.html>.
"""
function volume(S::AbstractArray{T,2})::T where {T}
    N, D = size(S)
    @assert N == D + 1
    B = Array{T}(undef, N + 1, N + 1)
    for i = 1:N+1, j = 1:N+1
        if i == N + 1 && j == N + 1
            B[i, j] = 0
        elseif i == N + 1 || j == N + 1
            B[i, j] = 1
        else
            B[i, j] = sum((S[i, :] - S[j, :]) .^ 2)
        end
    end
    return sqrt((-1)^(D + 1) * det(B) / (2^D * factorial(D)^2))
end



"""
Check consistency of result
"""
function check_triangulation(tri::Triangulation)
    nvertices, dim = size(tri.points)

    nsimplices = size(tri.simplices, 1)
    @test size(tri.simplices, 2) == dim + 1
    for i = 1:nsimplices
        si = tri.simplices[i, :]
        for d = 1:dim+1
            @test 1 <= si[d] <= nvertices
        end
    end

    @test size(tri.neighbors) == (nsimplices, dim + 1)
    for i = 1:nsimplices
        si = tri.simplices[i, :]
        for j in tri.neighbors[i, :]
            j == 0 && continue
            @test 1 <= j <= nsimplices
            @test j != i
            sj = tri.simplices[j, :]
            # Neighbouring simplices must differ in exactly 1 vertex
            @test length(setdiff(si, sj)) == 1
            @test length(setdiff(sj, si)) == 1
        end
    end

    @test size(tri.vertex_to_simplex) == (nvertices,)
    for i = 1:nvertices
        j = tri.vertex_to_simplex[i]
        @test 1 <= j <= nsimplices
        sj = tri.simplices[j, :]
        @test i in sj
    end

    vnv = tri.vertex_neighbor_vertices
    @test size(vnv) == (nvertices, nvertices)
    for j = 1:nvertices
        # Find all simplices containing vertex j
        nbs1 = Set{Int}()
        for n = 1:nsimplices
            if j in tri.simplices[n, :]
                for i in tri.simplices[n, :]
                    if i != j
                        push!(nbs1, i)
                    end
                end
            end
        end

        nbs2 = Set(rowvals(vnv)[nzrange(vnv, j)])
        @test nbs1 == nbs2
    end
end



@testset "Orthogonal simplices D=$D" for D = 1:Dmax
    coords = NTuple{D,Float64}[]
    push!(coords, ntuple(d -> 0, D))
    for n = 1:D
        push!(coords, ntuple(d -> d == n ? 1 : 0, D))
    end
    coords = flatten(coords)
    mesh = delaunay(coords)
    mesh::Triangulation
    @test size(mesh.points) == (D + 1, D)
    @test mesh.points == coords
    @test isempty(mesh.coplanar)
    @test size(mesh.simplices) == (1, D + 1)

    check_triangulation(mesh)

    V = 0.0
    for i = 1:size(mesh.simplices, 1)
        s = mesh.simplices[i, :]
        V += volume(mesh.points[s, :])
    end
end

@testset "Hypercubes D=$D" for D = 1:Dmax
    coords = NTuple{D,Float64}[]
    for i = CartesianIndex(ntuple(d -> 0, D)):CartesianIndex(ntuple(d -> 1, D))
        push!(coords, i.I)
    end
    coords = flatten(coords)
    mesh = delaunay(coords)
    mesh::Triangulation
    @test size(mesh.points) == (2^D, D)
    @test mesh.points == coords
    @test isempty(mesh.coplanar)
    # This would be true for the "standard" hypercube triangulation
    # @test size(mesh.simplices) == (factorial(D), D + 1)

    check_triangulation(mesh)

    V = 0.0
    for i = 1:size(mesh.simplices, 1)
        s = mesh.simplices[i, :]
        V += volume(mesh.points[s, :])
    end
    @test V ≈ 1
end

@testset "Random points D=$D" for D = 1:Dmax
    # Embed points into a hypercube
    coords = NTuple{D,Float64}[]
    for i = CartesianIndex(ntuple(d -> 0, D)):CartesianIndex(ntuple(d -> 1, D))
        push!(coords, i.I)
    end
    coords = [flatten(coords); rand(100, D)]
    mesh = delaunay(coords)
    mesh::Triangulation
    @test mesh.points == coords
    @test isempty(mesh.coplanar)

    check_triangulation(mesh)

    V = 0.0
    for i = 1:size(mesh.simplices, 1)
        s = mesh.simplices[i, :]
        V += volume(mesh.points[s, :])
    end
    @test V ≈ 1
end
