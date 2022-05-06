# Delaunay / də lɔˈnɛ /: Find the Delaunay triangulation for a set of points

* [GitHub](https://github.com/eschnett/Delaunay.jl): Source code
  repository
* [![GitHub CI](https://github.com/eschnett/Delaunay.jl/workflows/CI/badge.svg)](https://github.com/eschnett/Delaunay.jl/actions)

This package finds the [Delaunay
triangulation](https://en.wikipedia.org/wiki/Delaunay_triangulation)
for a set of points in arbitrary dimensions. It uses the Python
package
[`scipy.spatial.Delaunay`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.Delaunay.html)
to perform the actual calculation.

This package is inspired by
[QHull.jl](https://github.com/JuliaPolyhedra/QHull.jl), which uses the
same Python library.

## Example in 2D

```Julia
using Delaunay
points = rand(10, 2)
mesh = delaunay(points)

mesh.points                      # the points
mesh.simplices                   # the simplices (triangles in 2d)
mesh.neighbors                   # neighbouring simplices of a simplex
mesh.vertex_to_simplex           # find a simplex for a point
mesh.convex_hull                 # convex hull of the domain
mesh.vertex_neighbor_vertices    # neighbouring vertices of a vertex

using GLMakie # or CairoMakie/WGLMakie/RPRMakie
color = rand(size(mesh.points, 1))
fig, ax, pl = Makie.poly(mesh.points, mesh.simplices, color=color, strokewidth=2, figure=(resolution=(800, 400),))
points = randn(100, 2)
mesh = delaunay(points)
color = rand(size(mesh.points, 1))
poly(fig[1, 2], mesh.points, mesh.simplices, color=color, strokewidth=2)
save("delaunay2d.png", fig) 
```
![delaunay2d](https://user-images.githubusercontent.com/1010467/167169390-4c6b80b5-1370-424c-a495-8413996bdf68.png)


## Example in 3D

```Julia
using Delaunay
points = rand(100, 3)
mesh = delaunay(points)

using GeometryBasics
# Convert to tetrahedra faces
tetras = [GeometryBasics.TetrahedronFace(mesh.simplices[i, :]...) for i in 1:size(mesh.simplices, 1)]
points = Makie.to_vertices(mesh.points) # Use Makie to convert to Vector{Point3f}
m = GeometryBasics.Mesh(points, tetras) # create tetrahedra mesh
# Triangulate it, since Makie's mesh conversion currently doesn't handle tetrahedras itself 
tris = GeometryBasics.triangle_mesh(m)
fig, ax, pl = Makie.mesh(tris, color=rand(length(tris.position)), colormap=(:viridis, 0.5), transparency=true)
# add wireframe plot, which actually supports tetrahedras...
wireframe!(ax, m, color=:white)
save("delaunay3d.png", fig)
```
![delaunay3d](https://user-images.githubusercontent.com/1010467/167169584-447d1c1c-0f9f-4105-8897-cd5ad502e7d2.png)

The test cases contain also examples in higher dimensions.
