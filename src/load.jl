##################################################
# Helper functions
##################################################

_elementtuple(t::HomoNuclearTable) = (t.element_symbol, t.element_symbol)
_elementtuple(t::HeteroNuclearTable) = (t.first_element_symbol, t.second_element_symbol)

##################################################

"""
    loaddir([T], dirpath::String)

Parse a directory at `dirpath` containing SKF files and return a
`SlaterKosterTable` derived from its contents.

See also [`load`](@ref SlaterKoster.load).

"""
loaddir(dirpath) = loaddir(Float64, dirpath)
function loaddir(T::Type, dirpath)
    data = Dict{Tuple{Symbol,Symbol},PrimitiveTable{T}}()
    for file in filter(f -> occursin(r"(skf|SKF)$", f), Base.readdir(dirpath))
        table = load(T, joinpath(dirpath, file))
        data[_elementtuple(table)] = table
    end
    return SlaterKosterTable(data)
end

"""
    load([T], filepath::String; [element_symbols::NTuple{2,String}])

Parse an SKF located at `filepath` and return a primitive table
derived from its contents.  The element symbols will be inferred from
the file name if it follows the standard `X-Y.skf` format, where `X`,
and `Y` are element symbols.  Alternatively, if a tuple of symbols is
provided to the keyword argument `element_symbols`, the element
symbols listed in the tuple will be used instead.  If `T` is supplied,
floating-point data will be parsed using the type `T` instead of the
default `Float64`.

See also [`loaddir`](@ref SlaterKoster.loaddir).

"""
load(filepath; element_symbols=(nothing, nothing)) = load(Float64, filepath, element_symbols=element_symbols)
function load(T::Type, filepath::String; element_symbols=(nothing, nothing))
    dirname, filename = splitdir(filepath)
    if element_symbols === (nothing, nothing)
        # Try to derive the element symbols from the file name by
        # matching file name patterns of known datasets
        m = match(r"^([A-Z][a-z]*)[-]*([A-Z][a-z]*)\.(skf|SKF)$", filename)
        m === nothing && error("could not determine element symbols from SKF file ($filename); please specify them manually using the keyword argument `element_symbols`")
        first_element_symbol, second_element_symbol, _ = m.captures
    else
        first_element_symbol, second_element_symbol = element_symbols
    end
    lines = Base.readlines(filepath)
    if first_element_symbol == second_element_symbol
        if startswith(lines[1], "@")
            return HomoNuclearTable(
                Symbol(first_element_symbol),
                _readextendedformat_homonuclear(T, lines)...)
        else
            return HomoNuclearTable(
                Symbol(first_element_symbol),
                _readsimpleformat_homonuclear(T, lines)...)
        end
    else
        if startswith(lines[1], "@")
            return HeteroNuclearTable(
                Symbol(first_element_symbol),
                Symbol(second_element_symbol),
                _readextendedformat_heteronuclear(T, lines)...)
        else
            return HeteroNuclearTable(
                Symbol(first_element_symbol),
                Symbol(second_element_symbol),
                _readsimpleformat_heteronuclear(T, lines)...)
        end
    end
end

function parseheader(T::Type, line)
    words = filter(word -> !all(isspace, word), split(line, x -> isspace(x) || x === ','))
    return parse(T, words[1]), parse(Int, words[2])
end

function parseline(T::Type, line)
    words = filter(word -> !all(isspace, word), split(line, x -> isspace(x) || x === ','))
    i = 1
    while i <= length(words)
        m = match(r"([0-9]+)\*([0-9\.]+)", words[i])
        if m !== nothing
            deleteat!(words, i)
            for _ in 1:parse(Int, m[1])
                insert!(words, i, m[2])
            end
        end
        i += 1
    end
    return parse.(T, words)
end

function _readsimpleformat_homonuclear(T, lines)
    grid_dist, n_grid_points = parseheader(T, lines[1])

    energy_table = DataFrame(
        Ed = T[], Ep = T[], Es = T[],
        SPE = T[],
        Ud = T[], Up = T[], Us = T[],
        fd = T[], fp = T[], fs = T[])
    push!(energy_table, parseline(T, lines[2]))

    mass_table = DataFrame(
        mass = T[],
        c2 = T[], c3 = T[], c4 = T[], c5 = T[],
        c6 = T[], c7 = T[], c8 = T[], c9 = T[],
        rcut = T[],
        d1 = T[], d2 = T[], d3 = T[], d4 = T[], d5 = T[],
        d6 = T[], d7 = T[], d8 = T[], d9 = T[], d10 = T[])
    push!(mass_table, parseline(T, lines[3]))

    integral_table = DataFrame(
        dist = T[],
        Hdd0 = T[], Hdd1 = T[], Hdd2 = T[], Hpd0 = T[], Hpd1 = T[], Hpp0 = T[],
        Hpp1 = T[], Hsd0 = T[], Hsp0 = T[], Hss0 = T[],
        Sdd0 = T[], Sdd1 = T[], Sdd2 = T[], Spd0 = T[], Spd1 = T[], Spp0 = T[],
        Spp1 = T[], Ssd0 = T[], Ssp0 = T[], Sss0 = T[])

    for (i, line) in enumerate(lines[4:4+n_grid_points-1])
        push!(integral_table, vcat([i*grid_dist], parseline(T, line)))
    end
    return grid_dist, n_grid_points, energy_table, mass_table, integral_table
end

function _readsimpleformat_heteronuclear(T, lines)
    grid_dist, n_grid_points = parseheader(T, lines[1])

    mass_table = DataFrame(
        mass = T[],
        c2 = T[], c3 = T[], c4 = T[], c5 = T[],
        c6 = T[], c7 = T[], c8 = T[], c9 = T[],
        rcut = T[],
        d1 = T[], d2 = T[], d3 = T[], d4 = T[], d5 = T[],
        d6 = T[], d7 = T[], d8 = T[], d9 = T[], d10 = T[])
    push!(mass_table, parseline(T, lines[2]))

    integral_table = DataFrame(
        dist = T[],
        Hdd0 = T[], Hdd1 = T[], Hdd2 = T[], Hpd0 = T[], Hpd1 = T[], Hpp0 = T[],
        Hpp1 = T[], Hsd0 = T[], Hsp0 = T[], Hss0 = T[],
        Sdd0 = T[], Sdd1 = T[], Sdd2 = T[], Spd0 = T[], Spd1 = T[], Spp0 = T[],
        Spp1 = T[], Ssd0 = T[], Ssp0 = T[], Sss0 = T[])

    for (i, line) in enumerate(lines[3:3+n_grid_points-1])
        push!(integral_table, vcat([i*grid_dist], parseline(T, line)))
    end
    return grid_dist, n_grid_points, mass_table, integral_table
end

function _readextendedformat_homonuclear(T, lines)
    grid_dist, n_grid_points = parseheader(T, lines[1])

    energy_table = DataFrame(
        Ef = T[], Ed = T[], Ep = T[], Es = T[],
        SPE = T[],
        Uf = T[], Ud = T[], Up = T[], Us = T[],
        ff = T[], fd = T[], fp = T[], fs = T[])
    push!(energy_table, parseline(T, lines[3]))

    mass_table = DataFrame(
        mass = T[],
        c2 = T[], c3 = T[], c4 = T[], c5 = T[],
        c6 = T[], c7 = T[], c8 = T[], c9 = T[],
        rcut = T[],
        d1 = T[], d2 = T[], d3 = T[], d4 = T[], d5 = T[],
        d6 = T[], d7 = T[], d8 = T[], d9 = T[], d10 = T[])
    push!(mass_table, parseline(T, lines[4]))

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
        push!(integral_table, vcat([i*grid_dist], parseline(T, line)))
    end
    return grid_dist, n_grid_points, energy_table, mass_table, integral_table
end

function _readextendedformat_heteronuclear(T, lines)
    grid_dist, n_grid_points = parseheader(T, lines[1])

    mass_table = DataFrame(
        mass = T[],
        c2 = T[], c3 = T[], c4 = T[], c5 = T[],
        c6 = T[], c7 = T[], c8 = T[], c9 = T[],
        rcut = T[],
        d1 = T[], d2 = T[], d3 = T[], d4 = T[], d5 = T[],
        d6 = T[], d7 = T[], d8 = T[], d9 = T[], d10 = T[])
    push!(mass_table, parseline(T, lines[3]))

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
        push!(integral_table, vcat([i*grid_dist], parseline(T, line)))
    end
    return grid_dist, n_grid_points, mass_table, integral_table
end
