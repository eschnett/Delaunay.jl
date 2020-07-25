module Delaunay

using PyCall
using SparseArrays

export Triangulation, delaunay



const spatial = PyNULL()

function __init__()
    copy!(spatial, pyimport_conda("scipy.spatial", "scipy"))
end



"""
   Triangulation

Result of Delaunay triangulation.
"""
struct Triangulation
    "Coordinates of input points (npoints, ndim)"
    points::Array{Float64,2}
    "(One-based) indices of the points forming the simplices (nsimplex, ndim+1)"
    simplices::Array{Int,2}
    "Indices of neighbor simplices for each simplex; 0 at boundary (nsimplex, ndim+1)"
    neighbors::Array{Int,2}
    "[normal, offset] forming the hyperplane equation of the facet on the paraboloid (nsimplex, ndim+2)"
    equations::Array{Float64,2}
    "Affine transform from x to the barycentric coordinates c (nsimplex, ndim+1,ndim)"
    transform::Array{Float64,3}
    "Lookup array, from a vertex, to some simplex which it is a part of (npoints)"
    vertex_to_simplex::Array{Int,1}
    "Vertices of facets forming the convex hull of the point set (nfaces, ndim)"
    convex_hull::Array{Int,2}
    "Indices of coplanar points and the corresponding indices of the nearest facet and the nearest vertex (ncoplanar, 3)"
    coplanar::Array{Int,2}
    "Neighboring vertices of vertices [neighbors, vertex]"
    vertex_neighbor_vertices::SparseMatrixCSC{Nothing,Int}
end



function Base.show(io::IO, tri::Triangulation)
    println(
        io,
        "Triangulation of $(size(tri.points, 1)) vertices in $(size(tri.points, 2)) dimensions",
    )
    println(io, "Coordinates of vertices:")
    println(io, "    ", tri.points)
    println(io, "Indices of simplices:")
    println(io, "    ", tri.simplices)
    println(io, "Lookup array from vertices to some containing simplex:")
    println(io, "    ", tri.vertex_to_simplex)
    println(io, "Indices of facets of convex hull of the vertices:")
    println(io, "    ", tri.convex_hull)
    println(io, "Indices of coplanar points, nearest facet, and nearest vertex:")
    println(io, "    ", isempty(tri.coplanar) ? "[]" : tri.coplanar)
    println(io, "Neighboring vertices of vertices:")
    vnv = tri.vertex_neighbor_vertices
    for j = 1:size(vnv, 2)
        print(io, "    [$j]: [")
        comma = ""
        for i0 in nzrange(vnv, j)
            i = rowvals(vnv)[i0]
            print(io, comma, i)
            comma = ", "
        end
        println(io, "]")
    end
end



"""
    delaunay(vertices::Array{Float64, 2})::Triangulation

Calculate the Delaunay triangulation of the given vertices.

Input format: `size(vertices) = N, D` where `N` is the number of
vertices and `D` the spatial dimension.

Algorithm: Uses SciPy's `Delaunay` function. See
<https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.Delaunay.html>.
"""
function delaunay(vertices::Array{Float64,2})::Triangulation
    nvertices, dim = size(vertices)
    dim == 1 && return delaunay_1d(vertices)
    py = spatial.Delaunay(vertices)
    points = convert(Array{Float64,2}, py."points")
    simplices = inc!(convert(Array{Int,2}, py."simplices"))
    neighbors = inc!(convert(Array{Int,2}, py."neighbors"))
    equations = convert(Array{Float64,2}, py."equations")
    transform = convert(Array{Float64,3}, py."transform")
    vertex_to_simplex = inc!(convert(Array{Int,1}, py."vertex_to_simplex"))
    convex_hull = inc!(convert(Array{Int,2}, py."convex_hull"))
    coplanar = inc!(convert(Array{Int,2}, py."coplanar"))
    indptr = inc!(convert(Array{Int,1}, get(py."vertex_neighbor_vertices", 0)))
    indices = inc!(convert(Array{Int,1}, get(py."vertex_neighbor_vertices", 1)))
    # # Convert to sparse matrix format
    # I = Int[]
    # J = Int[]
    # V = Nothing[]
    # for i = 1:length(indptr)-1
    #     j0min = indptr[i]
    #     j0max = indptr[i+1]
    #     for j in indices[j0min:j0max-1]
    #         push!(I, i)
    #         push!(J, j)
    #         push!(V, nothing)
    #     end
    # end
    # # We transpose I and J here
    # vertex_neighbor_vertices = sparse(J, I, V, nvertices, nvertices)
    vertex_neighbor_vertices = SparseMatrixCSC{Nothing,Int}(
        nvertices,
        nvertices,
        indptr,
        indices,
        fill(nothing, length(indices)),
    )
    return Triangulation(
        points,
        simplices,
        neighbors,
        equations,
        transform,
        vertex_to_simplex,
        convex_hull,
        coplanar,
        vertex_neighbor_vertices,
    )
end

function delaunay_1d(vertices::Array{Float64,2})::Triangulation
    nvertices, dim = size(vertices)
    @assert dim == 1

    points = vertices

    # Sort vertices
    perm = sortperm(reshape(vertices, nvertices))

    nsimplices = nvertices - 1
    simplices = Array{Int}(undef, nsimplices, 2)
    for i = 1:nsimplices
        pi0 = perm[i]
        pi1 = perm[i+1]
        simplices[i, 1] = min(pi0, pi1)
        simplices[i, 2] = max(pi0, pi1)
    end

    neighbors = Array{Int}(undef, nsimplices, dim + 1)
    for i = 1:nsimplices
        neighbors[i, 1] = i == 1 ? 0 : i - 1
        neighbors[i, 2] = i == nsimplices ? 0 : i + 1
    end

    equations = Array{Int}(undef, 0, dim + 2)          # TODO
    transform = Array{Float64}(undef, 0, dim + 1, dim) # TODO

    vertex_to_simplex = Array{Int}(undef, nvertices)
    for i = 1:nvertices
        vertex_to_simplex[perm[i]] = i == nvertices ? nvertices - 1 : i
    end

    convex_hull = Array{Int}(undef, 2, dim)
    convex_hull[1] = perm[1]
    convex_hull[2] = perm[nvertices]

    coplanar = Array{Int}(undef, 0, 3)

    # Convert to sparse matrix format
    I = Int[]
    J = Int[]
    V = Nothing[]
    for i = 1:nvertices
        if i > 1
            push!(I, perm[i])
            push!(J, perm[i-1])
            push!(V, nothing)
        end
        if i < nvertices
            push!(I, perm[i])
            push!(J, perm[i+1])
            push!(V, nothing)
        end
    end
    # We transpose I and J here
    vertex_neighbor_vertices = sparse(J, I, V, nvertices, nvertices)
    return Triangulation(
        points,
        simplices,
        neighbors,
        equations,
        transform,
        vertex_to_simplex,
        convex_hull,
        coplanar,
        vertex_neighbor_vertices,
    )
end



"""
   inc!(A::AbstractArray)

Convert Python zero-based to Julia 1-based array indices
"""
function inc!(A)
    for i in eachindex(A)
        A[i] += 1
    end
    return A
end

end
