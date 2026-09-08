using Test
using AkaghefAlgebra

@testset "AkaghefAlgebra" begin
  include("Core/numerictype_tests.jl")
  include("Core/spacespec_tests.jl")
  include("Core/interfaces_tests.jl")
  include("Polynomial/polalg_tests.jl")
  include("Sparse/sparseex_tests.jl")
  include("Hopf/hopf_tests.jl")
  include("Hopf/qsl2_dual_tests.jl")
  include("Hopf/taft_exact_tests.jl")
  include("Hopf/symbolic_taft_tests.jl")
  include("Cohomology/cochain_tests.jl")
  # include("parity/runtests_parity.jl")  # Phase D-5
end
