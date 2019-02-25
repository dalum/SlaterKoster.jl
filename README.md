# SlaterKoster.jl

[![Build Status](https://travis-ci.org/dalum/SlaterKoster.jl.svg?branch=master)](https://travis-ci.org/dalum/SlaterKoster.jl)

A Julia package for loading SKF files and calculating Slater-Koster overlaps.

## Usage
```julia
julia> using SlaterKoster

julia> skt = SlaterKoster.loaddir("path/to/skf/files")

julia> # Calculate the hamiltonian matrix element between a px orbital on a
julia> # carbon atom, and an s orbital on a hydrogen atom, separated by one
julia> # Bohr radius along the x axis using the datasets loaded into `skt`:
julia> SlaterKoster.hamiltonian(skt, [1., 0., 0.,], (:C, :px), (:H, :s))
0.5553591362883
```

