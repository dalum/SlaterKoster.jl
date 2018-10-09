module SlaterKoster

using DataFrames

abstract type SlaterKosterTable end

struct HomoNuclearTable <: SlaterKosterTable
    element_name::String
    energy_table::DataFrame
    mass_table::DataFrame
    integral_table::DataFrame
end

struct HeteroNuclearTable <: SlaterKosterTable
    first_element_name::String
    second_element_name::String
    mass_table::DataFrame
    integral_table::DataFrame
end

"""
    read(filepath; [element_names])

Parse `filepath` and return a Slater-Koster table derived from its
contents.

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
                first_element_name,
                _readextendedformat_homonuclear(T, lines)...)
        else
            return HomoNuclearTable(
                first_element_name,
                _readsimpleformat_homonuclear(T, lines)...)
        end
    else
        if startswith(lines[1], "@")
            return HeteroNuclearTable(
                first_element_name,
                second_element_name,
                _readextendedformat_heteronuclear(T, lines)...)
        else
            return HeteroNuclearTable(
                first_element_name,
                second_element_name,
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
    return energy_table, mass_table, integral_table
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
    return mass_table, integral_table
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
    return energy_table, mass_table, integral_table
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
    return mass_table, integral_table
end

end # module
