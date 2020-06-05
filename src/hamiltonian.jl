"""
    hamiltonian(skt::SlaterKosterTable)

Return a function equivalent to
`(args...) -> hamiltonian(skt, args...)`.

"""
hamiltonian(skt::SlaterKosterTable) = (args...) -> hamiltonian(skt, args...)

"""
    hamiltonian(skt::SlaterKosterTable, r::NTuple{3,<:Real}, ϕ1, ϕ2)

Return the hamiltonian matrix element corresponding to hopping from
`ϕ2` to `ϕ1`, where `ϕ1` is located at `r` relative to `ϕ2`, given in
Angstrom.

`ϕ1`, `ϕ2` are pairs of `element_symbol => orbital_symbol`, where
`element_symbol` is a valid `element_symbol` from `elementsymbols(skt)`,
and `orbital_symbol` are tesseral spherical harmonic orbitals.

Accepted symbols for the orbitals are:

:s, :px, :py, :pz

!!! note
    The following orbitals are not yet supported: :dx2y2, :dxy, :dxz,
    :dyz, :dz2, :f...

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

"""
    hamiltonian(skt::SlaterKosterTable, ϕ)

Return the diagonal hamiltonian matrix element corresponding to `ϕ`.
`ϕ` is a pair of `element_symbol => orbital_symbol`, where
`element_symbol` is a valid `element_symbol` from
`elementsymbols(skt)`, and `orbital_symbol` is a tesseral spherical
harmonic orbital.

Accepted symbols for the orbitals are:

:s, :px, :py, :pz

!!! note
    The following orbitals are not yet supported: :dx2y2, :dxy, :dxz,
    :dyz, :dz2, :f...

"""
hamiltonian(skt::SlaterKosterTable, ϕ) =
    hamiltonian(skt.data[first(ϕ), first(ϕ)], last(ϕ))

function hamiltonian(skt::HomoNuclearTable, orbital_symbol::Symbol)
    orbital = Symbol(first(string(orbital_symbol)))
    df = skt.energy_table
    E = orbital === :s ? df.Es :
        orbital === :p ? df.Ep :
        orbital === :d ? df.Ed :
        orbital === :f ? df.Ef :
        error("invalid orbital symbol: $orbital_symbol")
    return electronvolt*E[]
end

function hamiltonian(skt::PrimitiveTable, r, orbital_symbol1::Symbol, orbital_symbol2::Symbol)
    # Convert distance to atomic units for indexing
    r = r/a0
    @assert length(r) == 3 "`length(r) = $(length(r)) != 3`: `r` is not a valid 3-dimensional vector: $r"
    min_dist, max_dist = extrema(skt.integral_table[!, :dist])
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
        name => spline_catmullrom(rem, col[iprev], col[ilo], col[ihi], col[inext]) for (name,col) in collect(pairs(eachcol(skt.integral_table[slice, :])))
    ))
    # Normalize the separation vector
    d = r ./ norm_r
    return electronvolt*hamiltonian(table, d, Val(orbital_symbol1), Val(orbital_symbol2))[]
end

function hamiltonian(df::DataFrame, d, ::Val{:s}, ::Val{:s})
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
