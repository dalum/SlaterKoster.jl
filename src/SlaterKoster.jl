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
function read(filepath::String; element_names=(nothing, nothing))
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
                _readextendedformat_homonuclear(lines)...)
        else
            return HomoNuclearTable(
                first_element_name,
                _readsimpleformat_homonuclear(lines)...)
        end
    else
        if startswith(lines[1], "@")
            return HeteroNuclearTable(
                first_element_name,
                second_element_name,
                _readextendedformat_heteronuclear(lines)...)
        else
            return HeteroNuclearTable(
                first_element_name,
                second_element_name,
                _readsimpleformat_heteronuclear(lines)...)
        end
    end
end

function _readsimpleformat_homonuclear(lines)
    grid_dist, n_grid_points = parse.((Float64, Int), split(lines[1]))

    energy_table = DataFrame(
        Ed = Float64[], Ep = Float64[], Es = Float64[],
        SPE = Float64[],
        Ud = Float64[], Up = Float64[], Us = Float64[],
        fd = Float64[], fp = Float64[], fs = Float64[])
    push!(energy_table, parse.(Float64, split(lines[2])))

    mass_table = DataFrame(
        mass = Float64[],
        c2 = Float64[], c3 = Float64[], c4 = Float64[], c5 = Float64[],
        c6 = Float64[], c7 = Float64[], c8 = Float64[], c9 = Float64[],
        rcut = Float64[],
        d1 = Float64[], d2 = Float64[], d3 = Float64[], d4 = Float64[], d5 = Float64[],
        d6 = Float64[], d7 = Float64[], d8 = Float64[], d9 = Float64[], d10 = Float64[])
    push!(mass_table, parse.(Float64, split(lines[3])))

    integral_table = DataFrame(
        dist = Float64[],
        Hdd0 = Float64[], Hdd1 = Float64[], Hdd2 = Float64[], Hpd0 = Float64[], Hpd1 = Float64[], Hpp0 = Float64[],
        Hpp1 = Float64[], Hsd0 = Float64[], Hsp0 = Float64[], Hss0 = Float64[],
        Sdd0 = Float64[], Sdd1 = Float64[], Sdd2 = Float64[], Spd0 = Float64[], Spd1 = Float64[], Spp0 = Float64[],
        Spp1 = Float64[], Ssd0 = Float64[], Ssp0 = Float64[], Sss0 = Float64[])

    for (i, line) in enumerate(lines[4:4+n_grid_points-1])
        push!(integral_table, vcat([i*grid_dist], parse.(Float64, split(line))))
    end
    return energy_table, mass_table, integral_table
end

function _readsimpleformat_heteronuclear(lines)
    grid_dist, n_grid_points = parse.((Float64, Int), split(lines[1]))

    mass_table = DataFrame(
        mass = Float64[],
        c2 = Float64[], c3 = Float64[], c4 = Float64[], c5 = Float64[],
        c6 = Float64[], c7 = Float64[], c8 = Float64[], c9 = Float64[],
        rcut = Float64[],
        d1 = Float64[], d2 = Float64[], d3 = Float64[], d4 = Float64[], d5 = Float64[],
        d6 = Float64[], d7 = Float64[], d8 = Float64[], d9 = Float64[], d10 = Float64[])
    push!(mass_table, parse.(Float64, split(lines[2])))

    integral_table = DataFrame(
        dist = Float64[],
        Hdd0 = Float64[], Hdd1 = Float64[], Hdd2 = Float64[], Hpd0 = Float64[], Hpd1 = Float64[], Hpp0 = Float64[],
        Hpp1 = Float64[], Hsd0 = Float64[], Hsp0 = Float64[], Hss0 = Float64[],
        Sdd0 = Float64[], Sdd1 = Float64[], Sdd2 = Float64[], Spd0 = Float64[], Spd1 = Float64[], Spp0 = Float64[],
        Spp1 = Float64[], Ssd0 = Float64[], Ssp0 = Float64[], Sss0 = Float64[])

    for (i, line) in enumerate(lines[3:3+n_grid_points-1])
        push!(integral_table, vcat([i*grid_dist], parse.(Float64, split(line))))
    end
    return mass_table, integral_table
end

function _readextendedformat_homonuclear(lines)
    grid_dist, n_grid_points = parse.((Float64, Int), split(lines[2]))

    energy_table = DataFrame(
        Ef = Float64[], Ed = Float64[], Ep = Float64[], Es = Float64[],
        SPE = Float64[],
        Uf = Float64[], Ud = Float64[], Up = Float64[], Us = Float64[],
        ff = Float64[], fd = Float64[], fp = Float64[], fs = Float64[])
    push!(energy_table, parse.(Float64, split(lines[3])))

    mass_table = DataFrame(
        mass = Float64[],
        c2 = Float64[], c3 = Float64[], c4 = Float64[], c5 = Float64[],
        c6 = Float64[], c7 = Float64[], c8 = Float64[], c9 = Float64[],
        rcut = Float64[],
        d1 = Float64[], d2 = Float64[], d3 = Float64[], d4 = Float64[], d5 = Float64[],
        d6 = Float64[], d7 = Float64[], d8 = Float64[], d9 = Float64[], d10 = Float64[])
    push!(mass_table, parse.(Float64, split(lines[4])))

    integral_table = DataFrame(
        dist = Float64[],
        Hff0 = Float64[], Hff1 = Float64[], Hff2 = Float64[], Hff3 = Float64[],
        Hdf0 = Float64[], Hdf1 = Float64[], Hdf2 = Float64[], Hdd0 = Float64[],
        Hdd1 = Float64[], Hdd2 = Float64[], Hpf0 = Float64[], Hpf1 = Float64[],
        Hpd0 = Float64[], Hpd1 = Float64[], Hpp0 = Float64[], Hpp1 = Float64[],
        Hsf0 = Float64[], Hsd0 = Float64[], Hsp0 = Float64[], Hss0 = Float64[],
        Sdd1 = Float64[], Sdd2 = Float64[], Sff0 = Float64[], Sff1 = Float64[],
        Sff2 = Float64[], Sff3 = Float64[], Sdf0 = Float64[], Sdf1 = Float64[],
        Sdf2 = Float64[], Sdd0 = Float64[], Spf0 = Float64[], Spf1 = Float64[],
        Spd0 = Float64[], Spd1 = Float64[], Spp0 = Float64[], Spp1 = Float64[],
        Ssf0 = Float64[], Ssd0 = Float64[], Ssp0 = Float64[], Sss0 = Float64[])

    for (i, line) in enumerate(lines[5:5+n_grid_points-1])
        push!(integral_table, vcat([i*grid_dist], parse.(Float64, split(line))))
    end
    return energy_table, mass_table, integral_table
end

function _readextendedformat_heteronuclear(lines)
    grid_dist, n_grid_points = parse.((Float64, Int), split(lines[2]))

    mass_table = DataFrame(
        mass = Float64[],
        c2 = Float64[], c3 = Float64[], c4 = Float64[], c5 = Float64[],
        c6 = Float64[], c7 = Float64[], c8 = Float64[], c9 = Float64[],
        rcut = Float64[],
        d1 = Float64[], d2 = Float64[], d3 = Float64[], d4 = Float64[], d5 = Float64[],
        d6 = Float64[], d7 = Float64[], d8 = Float64[], d9 = Float64[], d10 = Float64[])
    push!(mass_table, parse.(Float64, split(lines[3])))

    integral_table = DataFrame(
        dist = Float64[],
        Hff0 = Float64[], Hff1 = Float64[], Hff2 = Float64[], Hff3 = Float64[],
        Hdf0 = Float64[], Hdf1 = Float64[], Hdf2 = Float64[], Hdd0 = Float64[],
        Hdd1 = Float64[], Hdd2 = Float64[], Hpf0 = Float64[], Hpf1 = Float64[],
        Hpd0 = Float64[], Hpd1 = Float64[], Hpp0 = Float64[], Hpp1 = Float64[],
        Hsf0 = Float64[], Hsd0 = Float64[], Hsp0 = Float64[], Hss0 = Float64[],
        Sdd1 = Float64[], Sdd2 = Float64[], Sff0 = Float64[], Sff1 = Float64[],
        Sff2 = Float64[], Sff3 = Float64[], Sdf0 = Float64[], Sdf1 = Float64[],
        Sdf2 = Float64[], Sdd0 = Float64[], Spf0 = Float64[], Spf1 = Float64[],
        Spd0 = Float64[], Spd1 = Float64[], Spp0 = Float64[], Spp1 = Float64[],
        Ssf0 = Float64[], Ssd0 = Float64[], Ssp0 = Float64[], Sss0 = Float64[])

    for (i, line) in enumerate(lines[4:4+n_grid_points-1])
        push!(integral_table, vcat([i*grid_dist], parse.(Float64, split(line))))
    end
    return mass_table, integral_table
end

end # module
