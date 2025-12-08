### A Pluto.jl notebook ###
# v0.20.21

using Markdown
using InteractiveUtils

# ╔═╡ 2bd83dc9-f133-451e-ab9e-4ca0770969fa
import Pkg, Revise; Pkg.activate(Base.current_project())

# ╔═╡ e720a2e1-af6f-4c8f-a806-405022dafe35
import PlutoUI, Makie, CairoMakie, Logomaker, PythonPlot, ColorSchemes, BioStructures, TransitionPaths2024, FASTX, BioSequences, DataFrames, CSV

# ╔═╡ e3048c41-237e-4817-b119-f08b4d65f4c8
ancestor_1to2_batch4_sequences = BioSequences.LongAA.(
	FASTX.sequence.(FASTX.FASTAReader(open("/Users/jfdcd/projects/2024/Transition_Paths/TransitionPaths2024.jl/data/ancestral/ancestor_1to2_batch4.fa"))));

# ╔═╡ c264b192-3e71-4155-aaad-00efef9b7222
ancestor_1to2_sequences = BioSequences.LongAA.(
	FASTX.sequence.(FASTX.FASTAReader(open("/Users/jfdcd/projects/2024/Transition_Paths/TransitionPaths2024.jl/data/ancestral/ancestor_1to2.fa"))));

# ╔═╡ 840c0178-6247-4655-8ed3-b03d4db9cdf2
ancestor_1to4_sequences = BioSequences.LongAA.(
	FASTX.sequence.(FASTX.FASTAReader(open("/Users/jfdcd/projects/2024/Transition_Paths/TransitionPaths2024.jl/data/ancestral/ancestor_1to4.fa"))));

# ╔═╡ a821dc9e-b172-4f37-9c69-2ec7cc202778
ancestor_1to2_batch4 = TransitionPaths2024.onehot(ancestor_1to2_batch4_sequences);

# ╔═╡ 0e99ba98-bdae-49b5-b9aa-77f2a4efc195
ancestor_1to2 = TransitionPaths2024.onehot(ancestor_1to2_sequences);

# ╔═╡ 449c216f-cc8d-47d3-9868-522989c1cf14
ancestor_1to4 = TransitionPaths2024.onehot(ancestor_1to4_sequences);

# ╔═╡ e32748a4-9db8-4558-a301-5dcd3378b8f1
df_ancestor_1to2_batch4 = DataFrames.DataFrame(;
	seq = string.(ancestor_1to2_batch4_sequences),
	Eglob = TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:global, ancestor_1to2_batch4),
	EI = TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:I, ancestor_1to2_batch4),
	EIII = TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:II, ancestor_1to2_batch4),
	EIV = TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:IV, ancestor_1to2_batch4),
);

# ╔═╡ f525ebd2-0ee1-4d32-a380-deec48160d92
df_ancestor_1to2 = DataFrames.DataFrame(;
	seq = string.(ancestor_1to2_sequences),
	Eglob = TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:global, ancestor_1to2),
	EI = TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:I, ancestor_1to2),
	EIII = TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:II, ancestor_1to2),
	EIV = TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:IV, ancestor_1to2),
);

# ╔═╡ 9c6072c1-e06c-4e79-8333-cd2379f29ecd
df_ancestor_1to4 = DataFrames.DataFrame(;
	seq = string.(ancestor_1to4_sequences),
	Eglob = TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:global, ancestor_1to4),
	EI = TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:I, ancestor_1to4),
	EIII = TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:II, ancestor_1to4),
	EIV = TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:IV, ancestor_1to4),
);

# ╔═╡ 8ef5f5d8-286b-4c75-bf4d-bc2752d576c5
let ancestor_1to2_batch4_csv = tempname()
	CSV.write(ancestor_1to2_batch4_csv, df_ancestor_1to2_batch4)
	println(ancestor_1to2_batch4_csv)
end

# ╔═╡ 19211e68-365f-45d6-ad9c-2ca99e29d727
let ancestor_1to2_csv = tempname()
	CSV.write(ancestor_1to2_csv, df_ancestor_1to2)
	println(ancestor_1to2_csv)
end

# ╔═╡ d149f29b-c01a-4177-9883-e3f647438187
let ancestor_1to4_csv = tempname()
	CSV.write(ancestor_1to4_csv, df_ancestor_1to4)
	println(ancestor_1to4_csv)
end

# ╔═╡ aef24d32-d58a-49fb-90d3-214947678c9b
md"# Sampled paths with RBM scores for P.G."

# ╔═╡ b215dbd6-e594-4ee7-8422-9d919e429bad
dir_path_for_rbm_scores = mktempdir()

# ╔═╡ 9e436d03-2fcb-412b-884f-65a788c35b11
let file_path = joinpath(dir_path_for_rbm_scores, "1to4_batch3_full.csv")
	sampled_path = unique(TransitionPaths2024.sampled_path_1to4_20240703())
	df = DataFrames.DataFrame(sequence = sampled_path)
	for rbm = (:I, :II, :IV, :global)
		scores = TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(rbm,  TransitionPaths2024.onehot(sampled_path))
		df[!, "RBM_score-$rbm"] = scores
	end
	CSV.write(file_path, df)
	println(file_path)
end

# ╔═╡ fc187edf-b1de-474d-871f-7a1e7fb9b223
let dir_path = mktempdir()
	file_path = joinpath(dir_path_for_rbm_scores, "1to2_rep1_batch4_full.csv")
	sampled_path = unique(TransitionPaths2024.sampled_path_1to2rep1_20240703())
	df = DataFrames.DataFrame(sequence = sampled_path)
	for rbm = (:I, :II, :IV, :global)
		scores = TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(rbm,  TransitionPaths2024.onehot(sampled_path))
		df[!, "RBM_score-$rbm"] = scores
	end
	CSV.write(file_path, df)
	println(file_path)
end

# ╔═╡ 6f3d9d22-743c-4c0d-b976-ba37b32dfaf9
let dir_path = mktempdir()
	file_path = joinpath(dir_path_for_rbm_scores, "1to1_batch2_full.csv")
	sampled_path = unique(TransitionPaths2024.sampled_path_1to1_20240703())
	df = DataFrames.DataFrame(sequence = sampled_path)
	for rbm = (:I, :II, :IV, :global)
		scores = TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(rbm,  TransitionPaths2024.onehot(sampled_path))
		df[!, "RBM_score-$rbm"] = scores
	end
	CSV.write(file_path, df)
	println(file_path)
end

# ╔═╡ c30d0876-e337-432f-9a7d-8622832a38b8
file_path

# ╔═╡ Cell order:
# ╠═2bd83dc9-f133-451e-ab9e-4ca0770969fa
# ╠═e720a2e1-af6f-4c8f-a806-405022dafe35
# ╠═e3048c41-237e-4817-b119-f08b4d65f4c8
# ╠═c264b192-3e71-4155-aaad-00efef9b7222
# ╠═840c0178-6247-4655-8ed3-b03d4db9cdf2
# ╠═a821dc9e-b172-4f37-9c69-2ec7cc202778
# ╠═0e99ba98-bdae-49b5-b9aa-77f2a4efc195
# ╠═449c216f-cc8d-47d3-9868-522989c1cf14
# ╠═e32748a4-9db8-4558-a301-5dcd3378b8f1
# ╠═f525ebd2-0ee1-4d32-a380-deec48160d92
# ╠═9c6072c1-e06c-4e79-8333-cd2379f29ecd
# ╠═8ef5f5d8-286b-4c75-bf4d-bc2752d576c5
# ╠═19211e68-365f-45d6-ad9c-2ca99e29d727
# ╠═d149f29b-c01a-4177-9883-e3f647438187
# ╠═aef24d32-d58a-49fb-90d3-214947678c9b
# ╠═b215dbd6-e594-4ee7-8422-9d919e429bad
# ╠═9e436d03-2fcb-412b-884f-65a788c35b11
# ╠═fc187edf-b1de-474d-871f-7a1e7fb9b223
# ╠═6f3d9d22-743c-4c0d-b976-ba37b32dfaf9
# ╠═c30d0876-e337-432f-9a7d-8622832a38b8
