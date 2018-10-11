##################################################
# Helper functions
##################################################

"""
    spline_catmullrom(u, y0, y1, y2, y3)

Smooth interpolation between `y1, y2`, using the Catmull-Rom spline.
The function uses `y0` and `y3` to ensure that the derivative at `y1`
and `y2` is continuous.  `u` is a parameter between `0` and `1`.

"""
spline_catmullrom(u, y0, y1, y2, y3) = y1 + 0.5u*(-y0 + y2 + u*(2y0 - 5y1 + 4y2 - y3 + u*(-y0 + 3y1 - 3y2 + y3)))

_orbitalrank(orbital_symbol) =
    startswith(string(orbital_symbol), "s") ? 1 :
    startswith(string(orbital_symbol), "p") ? 2 :
    startswith(string(orbital_symbol), "d") ? 3 :
    startswith(string(orbital_symbol), "f") ? 4 :
    error("invalid orbital symbol: $orbital_symbol")

##################################################

"""
    hamiltonian(skt::SlaterKosterTable, r::NTuple{3,<:Real}, ϕ1, ϕ2)

Return the hamiltonian matrix element corresponding to hopping from
`ϕ2` to `ϕ1`, where `ϕ1` is located at `r` relative to `ϕ2`, given in
atomic units (i. e. Bohr radii).

`ϕ1`, `ϕ2` are pairs of `element_symbol => orbital_symbol`, where
`element_symbol` is a valid `element_symbol` from `elementsymbols(skt)`,
and `orbital_symbol` are tesseral spherical harmonic orbitals.

Accepted symbols for the orbitals are:

:s, :px, :py, :pz

[!!!Not yet supported: :dx2y2, :dxy, :dxz, :dyz, :dz2, :f...]

"""
function hamiltonian(skt::SlaterKosterTable, r, ϕ1, ϕ2)
    orbital_symbol1, orbital_symbol2 = last(ϕ1), last(ϕ2)
    if _orbitalrank(orbital_symbol2) < _orbitalrank(orbital_symbol1)
        ϕ1, ϕ2 = ϕ2, ϕ1
        r = -r
    end
    table = skt.data[first(ϕ1), first(ϕ2)]
    return hamiltonian(table, r, last(ϕ1), last(ϕ2))
end

function hamiltonian(skt::PrimitiveTable, r, orbital_symbol1::Symbol, orbital_symbol2::Symbol)
    @assert length(r) == 3 "`length(r) = $(length(r)) != 3`: `r` is not a valid 3-dimensional vector: $r"
    min_dist, max_dist = extrema(skt.integral_table[:dist])
    norm_r = norm(r)
    @assert min_dist <= norm_r <= max_dist "`norm(r) = $norm_r` is outside table bounds: [$min_dist, $max_dist]"
    grid_dist = skt.grid_dist
    lo = floor(Int, norm_r/grid_dist)
    hi = norm_r == max_dist ? lo : lo + 1
    # To support automatic differentiation, we cannot use the
    # remainder operator (`%`) here, but subtract instead
    rem = norm_r/grid_dist - lo
    # Interpolate the table for the given distance
    slice = max(1, lo-1):min(skt.n_grid_points, hi+1)
    iprev, ilo = lo > 1 ? (1, 2) : (1, 1)
    ihi, inext = hi < skt.n_grid_points ? (ilo + 1, ilo + 2) : (ilo + 1, ilo + 1)
    table = DataFrame(Dict(
        name => spline_catmullrom(rem, col[iprev], col[ilo], col[ihi], col[inext]) for (name,col) in eachcol(skt.integral_table[slice, :])))
    # Normalize the separation vector
    d = r ./ norm_r
    return hamiltonian(table, d, Val(orbital_symbol1), Val(orbital_symbol2))[]
end

function hamiltonian(df::DataFrame, r, ::Val{:s}, ::Val{:s})
    return df.Hss0
end

for (i, x) in [(1, :(:x)), (2, :(:y)), (3, :(:z))]
    @eval @inline function hamiltonian(df::DataFrame, d, ::Val{:s}, ::Val{Symbol("p", $x)})
        return d[$i]*df.Hsp0
    end

    for (j, y) in [(1, :(:x)), (2, :(:y)), (3, :(:z))]
        if i == j
            @eval @inline function hamiltonian(df::DataFrame, d, ::Val{Symbol("p", $x)}, ::Val{Symbol("p", $x)})
                return d[$i]^2*df.Hpp0 + (1 - d[$i]^2)*df.Hpp1
            end
        else
            @eval @inline function hamiltonian(df::DataFrame, d, ::Val{Symbol("p", $x)}, ::Val{Symbol("p", $y)})
                return d[$i]*d[$j]*(df.Hpp0 - df.Hpp1)
            end
        end
    end
end