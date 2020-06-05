mass(skt::SlaterKosterTable) = (args...) -> mass(skt, args...)
mass(skt::SlaterKosterTable, element_symbol::Symbol) = mass(skt.data[element_symbol, element_symbol])
mass(skt::HomoNuclearTable) = amu*skt.mass_table[!, :mass][]
