import Aqua
import TransitionPaths2025
using Test: @testset

@testset verbose = true "aqua" begin
    Aqua.test_all(TransitionPaths2025; ambiguities = false, stale_deps = false)
end
