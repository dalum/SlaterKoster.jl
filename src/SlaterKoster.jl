module SlaterKoster

using DataFrames
using LinearAlgebra

##################################################
# Constants
##################################################

# Atomic mass unit in atomic units (`amu / m_e`)
const amu = 1_822.888_486_192
# Hartree in units of electronvolts (eV)
const electronvolt = 27.211_386_02

##################################################

abstract type PrimitiveTable{T} end

struct HomoNuclearTable{T} <: PrimitiveTable{T}
    element_symbol::Symbol
    grid_dist::T
    n_grid_points::Int
    energy_table::DataFrame
    mass_table::DataFrame
    integral_table::DataFrame
end

struct HeteroNuclearTable{T} <: PrimitiveTable{T}
    first_element_symbol::Symbol
    second_element_symbol::Symbol
    grid_dist::T
    n_grid_points::Int
    mass_table::DataFrame
    integral_table::DataFrame
end

function Base.:(==)(t1::T, t2::T) where {T<:PrimitiveTable}
    fields = fieldnames(T)
    return all(field -> getfield(t1, field) == getfield(t2, field), fields)
end

struct SlaterKosterTable{T}
    data::Dict{Tuple{Symbol,Symbol},PrimitiveTable{T}}
end

Base.getindex(skt::SlaterKosterTable, x) = getindex(skt.data, x)

function Base.show(io::IO, skt::SlaterKosterTable)
    print(io, "SlaterKosterTable with entries for elements: ")
    print(io, join(elementsymbols(skt), ", ", " and "))
    print(io, ":\n")
    print(io, "    ")
    print(io, join(("$(first(key)) => $(last(key))" for key in keys(skt.data) |> collect |> sort!), "\n    "))
end

elementsymbols(skt::HomoNuclearTable) = [skt.element_symbol]
elementsymbols(skt::HeteroNuclearTable) = [skt.first_element_symbol, skt.second_element_symbol]
elementsymbols(skt::SlaterKosterTable) = reduce(vcat, elementsymbols.(values(skt.data))) |> unique! |> sort!

include("load.jl")
include("hamiltonian.jl")

mass(skt::SlaterKosterTable, element_symbol::Symbol) = mass(skt.data[element_symbol, element_symbol])
mass(skt::HomoNuclearTable) = amu*skt.mass_table[:mass][]

end # module
