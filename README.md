# Delaunay: Find the Delaunay triangulation for a set of points

This package provides at `delaunay` function that finds the [Delaunay
triangulation](https://en.wikipedia.org/wiki/Delaunay_triangulation)
for a set of points in arbitrary dimensions. It uses
[scipy.spatial.Delaunay](https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.Delaunay.html)
to perform the actual calculation. This package resembles
[QHull.jl](https://github.com/JuliaPolyhedra/QHull.jl), which uses the
same Python library.

## Example

```Julia
using Delaunay
points = rand(10, 2)
mesh = delaunay(points)
mesh.points                     # the points
mesh.simplices                  # the simplices (triangles in 2d)
mesh.neighbors                  # neighbouring simplices of a simplex
mesh.vertex_to_simplex          # find a simplex for a point
mesh.convex_hull                # convex hull of the domain
mesh.vertex_neighbor_vertices   # neighbouring vertices of a vertex

using Makie
scene = Makie.mesh(mesh.points, mesh.simplices)
wireframe!(scene[end][1], color = (:black, 0.6), linewidth = 3)
```
