module SlaterKosterTests

using Test
using SlaterKoster

const DIR = joinpath(@__DIR__, "datasets")

@testset "Load" begin
    @test SlaterKoster.load(joinpath(DIR, "C-C.skf")) == SlaterKoster.load(Float64, joinpath(DIR, "C-C.skf"))

    for T in (Float32, Float64)
        skt = SlaterKoster.loaddir(T, DIR)
        cc = SlaterKoster.load(T, joinpath(DIR, "C-C.skf"))
        ch = SlaterKoster.load(T, joinpath(DIR, "C-H.skf"), element_symbols=(:C, :H))

        @test skt[(:C, :C)] == cc
        @test skt[(:C, :H)] == ch
    end
end

end
