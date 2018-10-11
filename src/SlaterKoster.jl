module SlaterKoster

using DataFrames
using LinearAlgebra

##################################################
# Constants
##################################################

# Electron mass in atomic mass units
const m_e = 1822.888486192

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

struct SlaterKosterTable{T}
    data::Dict{Tuple{Symbol,Symbol},PrimitiveTable{T}}
end

function Base.show(io::IO, skt::SlaterKosterTable)
    print(io, "SlaterKosterTable with entries for elements: ")
    print(io, join(elementsymbols(skt), ", ", " and "))
    print(io, ":\n")
    print(io, "    ")
    print(io, join(("$(first(key)) => $(last(key))" for key in keys(skt.data) |> collect |> sort!), "\n    "))
end

elementsymbols(t::HomoNuclearTable) = [t.element_symbol]
elementsymbols(t::HeteroNuclearTable) = [t.first_element_symbol, t.second_element_symbol]
elementsymbols(skt::SlaterKosterTable) = reduce(vcat, elementsymbols.(values(skt.data))) |> unique! |> sort!

include("load.jl")
include("hamiltonian.jl")

end # module
