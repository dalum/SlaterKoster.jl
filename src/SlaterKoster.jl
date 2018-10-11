module SlaterKoster

using DataFrames
using LinearAlgebra

abstract type SlaterKosterTable{T} end

struct HomoNuclearTable{T} <: SlaterKosterTable{T}
    element_name::Symbol
    grid_dist::T
    n_grid_points::Int
    energy_table::DataFrame
    mass_table::DataFrame
    integral_table::DataFrame
end

struct HeteroNuclearTable{T} <: SlaterKosterTable{T}
    first_element_name::Symbol
    second_element_name::Symbol
    grid_dist::T
    n_grid_points::Int
    mass_table::DataFrame
    integral_table::DataFrame
end

"""
    read(filepath::String; [element_names::NTuple{2,String}])

Parse an SKF located at `filepath` and return a Slater-Koster table
derived from its contents.  The element names will be inferred from
the file name if it follows the standard `X-Y.skf` format, where `X`,
and `Y` are element names.  Alternatively, if a tuple of names is
provided to the keyword argument `element_names`, the element names
listed in the tuple will be used instead.

"""
read(filepath; element_names=(nothing, nothing)) = read(Float64, filepath, element_names=element_names)
function read(T::Type, filepath::String; element_names=(nothing, nothing))
    dirname, filename = splitdir(filepath)
    if element_names === (nothing, nothing)
        # Try to derive the element names from the file name by
        # matching file name patterns of known datasets
        m = match(r"^([A-Z][a-z]*)[-]*([A-Z][a-z]*)\.(skf|SKF)$", filename)
        m === nothing && error("could not determine element names from SKF file; please specify them manually using the keyword argument `element_names`")
        first_element_name, second_element_name, _ = m.captures
    else
        first_element_name, second_element_name = element_names
    end
    lines = Base.readlines(filepath)
    if first_element_name == second_element_name
        if startswith(lines[1], "@")
            return HomoNuclearTable(
                Symbol(first_element_name),
                _readextendedformat_homonuclear(T, lines)...)
        else
            return HomoNuclearTable(
                Symbol(first_element_name),
                _readsimpleformat_homonuclear(T, lines)...)
        end
    else
        if startswith(lines[1], "@")
            return HeteroNuclearTable(
                Symbol(first_element_name),
                Symbol(second_element_name),
                _readextendedformat_heteronuclear(T, lines)...)
        else
            return HeteroNuclearTable(
                Symbol(first_element_name),
                Symbol(second_element_name),
                _readsimpleformat_heteronuclear(T, lines)...)
        end
    end
end

function _readsimpleformat_homonuclear(T, lines)
    grid_dist, n_grid_points = parse.((T, Int), split(lines[1]))

    energy_table = DataFrame(
        Ed = T[], Ep = T[], Es = T[],
        SPE = T[],
        Ud = T[], Up = T[], Us = T[],
        fd = T[], fp = T[], fs = T[])
    push!(energy_table, parse.(T, split(lines[2])))

    mass_table = DataFrame(
        mass = T[],
        c2 = T[], c3 = T[], c4 = T[], c5 = T[],
        c6 = T[], c7 = T[], c8 = T[], c9 = T[],
        rcut = T[],
        d1 = T[], d2 = T[], d3 = T[], d4 = T[], d5 = T[],
        d6 = T[], d7 = T[], d8 = T[], d9 = T[], d10 = T[])
    push!(mass_table, parse.(T, split(lines[3])))

    integral_table = DataFrame(
        dist = T[],
        Hdd0 = T[], Hdd1 = T[], Hdd2 = T[], Hpd0 = T[], Hpd1 = T[], Hpp0 = T[],
        Hpp1 = T[], Hsd0 = T[], Hsp0 = T[], Hss0 = T[],
        Sdd0 = T[], Sdd1 = T[], Sdd2 = T[], Spd0 = T[], Spd1 = T[], Spp0 = T[],
        Spp1 = T[], Ssd0 = T[], Ssp0 = T[], Sss0 = T[])

    for (i, line) in enumerate(lines[4:4+n_grid_points-1])
        push!(integral_table, vcat([i*grid_dist], parse.(T, split(line))))
    end
    return grid_dist, n_grid_points, energy_table, mass_table, integral_table
end

function _readsimpleformat_heteronuclear(T, lines)
    grid_dist, n_grid_points = parse.((T, Int), split(lines[1]))

    mass_table = DataFrame(
        mass = T[],
        c2 = T[], c3 = T[], c4 = T[], c5 = T[],
        c6 = T[], c7 = T[], c8 = T[], c9 = T[],
        rcut = T[],
        d1 = T[], d2 = T[], d3 = T[], d4 = T[], d5 = T[],
        d6 = T[], d7 = T[], d8 = T[], d9 = T[], d10 = T[])
    push!(mass_table, parse.(T, split(lines[2])))

    integral_table = DataFrame(
        dist = T[],
        Hdd0 = T[], Hdd1 = T[], Hdd2 = T[], Hpd0 = T[], Hpd1 = T[], Hpp0 = T[],
        Hpp1 = T[], Hsd0 = T[], Hsp0 = T[], Hss0 = T[],
        Sdd0 = T[], Sdd1 = T[], Sdd2 = T[], Spd0 = T[], Spd1 = T[], Spp0 = T[],
        Spp1 = T[], Ssd0 = T[], Ssp0 = T[], Sss0 = T[])

    for (i, line) in enumerate(lines[3:3+n_grid_points-1])
        push!(integral_table, vcat([i*grid_dist], parse.(T, split(line))))
    end
    return grid_dist, n_grid_points, mass_table, integral_table
end

function _readextendedformat_homonuclear(T, lines)
    grid_dist, n_grid_points = parse.((T, Int), split(lines[2]))

    energy_table = DataFrame(
        Ef = T[], Ed = T[], Ep = T[], Es = T[],
        SPE = T[],
        Uf = T[], Ud = T[], Up = T[], Us = T[],
        ff = T[], fd = T[], fp = T[], fs = T[])
    push!(energy_table, parse.(T, split(lines[3])))

    mass_table = DataFrame(
        mass = T[],
        c2 = T[], c3 = T[], c4 = T[], c5 = T[],
        c6 = T[], c7 = T[], c8 = T[], c9 = T[],
        rcut = T[],
        d1 = T[], d2 = T[], d3 = T[], d4 = T[], d5 = T[],
        d6 = T[], d7 = T[], d8 = T[], d9 = T[], d10 = T[])
    push!(mass_table, parse.(T, split(lines[4])))

    integral_table = DataFrame(
        dist = T[],
        Hff0 = T[], Hff1 = T[], Hff2 = T[], Hff3 = T[],
        Hdf0 = T[], Hdf1 = T[], Hdf2 = T[], Hdd0 = T[],
        Hdd1 = T[], Hdd2 = T[], Hpf0 = T[], Hpf1 = T[],
        Hpd0 = T[], Hpd1 = T[], Hpp0 = T[], Hpp1 = T[],
        Hsf0 = T[], Hsd0 = T[], Hsp0 = T[], Hss0 = T[],
        Sdd1 = T[], Sdd2 = T[], Sff0 = T[], Sff1 = T[],
        Sff2 = T[], Sff3 = T[], Sdf0 = T[], Sdf1 = T[],
        Sdf2 = T[], Sdd0 = T[], Spf0 = T[], Spf1 = T[],
        Spd0 = T[], Spd1 = T[], Spp0 = T[], Spp1 = T[],
        Ssf0 = T[], Ssd0 = T[], Ssp0 = T[], Sss0 = T[])

    for (i, line) in enumerate(lines[5:5+n_grid_points-1])
        push!(integral_table, vcat([i*grid_dist], parse.(T, split(line))))
    end
    return grid_dist, n_grid_points, energy_table, mass_table, integral_table
end

function _readextendedformat_heteronuclear(T, lines)
    grid_dist, n_grid_points = parse.((T, Int), split(lines[2]))

    mass_table = DataFrame(
        mass = T[],
        c2 = T[], c3 = T[], c4 = T[], c5 = T[],
        c6 = T[], c7 = T[], c8 = T[], c9 = T[],
        rcut = T[],
        d1 = T[], d2 = T[], d3 = T[], d4 = T[], d5 = T[],
        d6 = T[], d7 = T[], d8 = T[], d9 = T[], d10 = T[])
    push!(mass_table, parse.(T, split(lines[3])))

    integral_table = DataFrame(
        dist = T[],
        Hff0 = T[], Hff1 = T[], Hff2 = T[], Hff3 = T[],
        Hdf0 = T[], Hdf1 = T[], Hdf2 = T[], Hdd0 = T[],
        Hdd1 = T[], Hdd2 = T[], Hpf0 = T[], Hpf1 = T[],
        Hpd0 = T[], Hpd1 = T[], Hpp0 = T[], Hpp1 = T[],
        Hsf0 = T[], Hsd0 = T[], Hsp0 = T[], Hss0 = T[],
        Sdd1 = T[], Sdd2 = T[], Sff0 = T[], Sff1 = T[],
        Sff2 = T[], Sff3 = T[], Sdf0 = T[], Sdf1 = T[],
        Sdf2 = T[], Sdd0 = T[], Spf0 = T[], Spf1 = T[],
        Spd0 = T[], Spd1 = T[], Spp0 = T[], Spp1 = T[],
        Ssf0 = T[], Ssd0 = T[], Ssp0 = T[], Sss0 = T[])

    for (i, line) in enumerate(lines[4:4+n_grid_points-1])
        push!(integral_table, vcat([i*grid_dist], parse.(T, split(line))))
    end
    return grid_dist, n_grid_points, mass_table, integral_table
end

"""
    hamiltonian(skt::SlaterKosterTable, r::NTuple{3,<:Real}, ϕ1::Symbol, ϕ2::Symbol)

Return the matrix element corresponding to `ϕ1(r)'H*ϕ2(0)`, where
`ϕ1`, `ϕ2` are tesseral spherical harmonic orbitals separated by a
vector `r` with entries given in atomic units (i. e. Bohr radii).

Accepted symbols for orbitals are:

:s, :px, :py, :pz

[!!!Not yet supported: :dx2y2, :dxy, :dxz, :dyz, :dz2, :f...]

"""
function hamiltonian(skt::SlaterKosterTable{T}, r::AbstractVector{T}, ϕ1::Symbol, ϕ2::Symbol) where {T<:Real}
    @assert length(r) == 3 "`length(r) = $(length(r)) != 3`: `r` is not a valid 3-dimensional vector: $r"
    min_dist, max_dist = extrema(skt.integral_table[:dist])
    norm_r = norm(r)
    @assert min_dist <= norm_r <= max_dist "`norm(r) = $norm_r` is outside table bounds: [$min_dist, $max_dist]"
    grid_dist = skt.grid_dist
    lo = floor(Int, norm_r/grid_dist)
    hi = lo + 1
    rem = norm_r % grid_dist
    # Interpolate the table for the given distance
    table = DataFrame(Dict(
        name => (1 - rem)*col[1] + rem*col[2] for (name,col) in eachcol(skt.integral_table[lo:hi, :])))
    # Normalize the separation vector
    d = r ./ norm_r
    return hamiltonian(table, d, Val(ϕ1), Val(ϕ2))[]
end

function hamiltonian(df::DataFrame, ::AbstractVector{T}, ::Val{:s}, ::Val{:s}) where {T<:Real}
    return df.Hss0
end

for (i, x) in [(1, :(:x)), (2, :(:y)), (3, :(:z))]
    @eval @inline function hamiltonian(df::DataFrame, d::AbstractVector{T}, ::Val{:s}, ::Val{Symbol("p", $x)}) where {T<:Real}
        return d[$i] * df.Hsp0
    end

    for (j, y) in [(1, :(:x)), (2, :(:y)), (3, :(:z))]
        if i == j
            @eval @inline function hamiltonian(df::DataFrame, d::AbstractVector{T}, ::Val{Symbol("p", $x)}, ::Val{Symbol("p", $x)}) where {T<:Real}
                return d[$i]^2*df.Hpp0 + (1 - d[$i]^2)*df.Hpp1
            end
        else
            @eval @inline function hamiltonian(df::DataFrame, d::AbstractVector{T}, ::Val{Symbol("p", $x)}, ::Val{Symbol("p", $y)}) where {T<:Real}
                return d[$i]*d[$j]*(df.Hpp0 - df.Hpp1)
            end
        end
    end
end

end # module
