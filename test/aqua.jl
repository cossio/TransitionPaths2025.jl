import Aqua
import TransitionPaths2024
using Test: @testset

@testset verbose = true "aqua" begin
    Aqua.test_all(TransitionPaths2024; ambiguities = false, stale_deps = false)
end
