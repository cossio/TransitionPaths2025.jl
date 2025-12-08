import Pfam

Pfam.set_pfam_directory(mktempdir())
Pfam.set_pfam_version("35.0")

module aqua_tests include("aqua.jl") end
module eugenio_tests include("eugenio.jl") end
module onehot_tests include("onehot.jl") end
module pfam_tests include("pfam.jl") end
module hamming_tests include("hamming.jl") end
module hu2004_tests include("hu2004.jl") end
module exp2307_tests include("experiments/202307.jl") end
