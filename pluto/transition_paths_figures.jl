### A Pluto.jl notebook ###
# v0.20.21

using Markdown
using InteractiveUtils

# ╔═╡ 70eb32c0-f459-403f-a1d8-8db8a3f2ca06
import Pkg, Revise; Pkg.activate(Base.current_project())

# ╔═╡ d301243a-53c7-4e7d-bb5f-056a9c7a77f8
using XLSX: readxlsx

# ╔═╡ 63e09505-9d7c-4db0-8c54-852ee51626be
using Makie: @L_str, @colorant_str, RGBf

# ╔═╡ 1e3060e6-f40e-44a9-937a-6965a2fa2134
using Statistics: cor, mean, cov, var, std

# ╔═╡ d1b1a928-3daf-4c30-8892-b6ced5c2b3a4
using StatsBase: corspearman, countmap

# ╔═╡ 1b8d68ed-ad51-4c3b-9023-e6372868e7d4
using LinearAlgebra: norm

# ╔═╡ 474d459d-083e-4639-9e7f-753bde25c217
using NaNStatistics: nancor

# ╔═╡ ec73dce9-638a-4f34-8e52-29d03c83c9f6
using RestrictedBoltzmannMachines: inputs_h_from_v, sample_from_inputs, free_energy, sample_v_from_v

# ╔═╡ 5f7923a5-47c3-42d3-94ce-89fef8c57871
using BioSequences: LongAA

# ╔═╡ 7809f3a5-e704-43da-af08-4e47f412c919
using ProgressLogging: @progress

# ╔═╡ 0089bfee-3916-11ef-0452-d739816a3d18
md"# Imports"

# ╔═╡ 7ff741cb-ea2d-435a-80fa-9f53391a3117
import PlutoUI, Makie, CairoMakie, Logomaker, PythonPlot, ColorSchemes, BioStructures, TransitionPaths2024

# ╔═╡ 8febf5c6-92bb-4b91-a022-eb6080c8f2a0
PlutoUI.TableOfContents()

# ╔═╡ 4e3e5c24-41df-4fbd-91d9-42ab34f45cbe
md"# Load data"

# ╔═╡ b5677fef-d1cd-46f5-8f74-7b32fdd4b030
thresh_significance = 0.212046695;

# ╔═╡ 6165d60b-bada-4400-b76b-9deb2bd1450e
all_response_data = TransitionPaths2024.response_data_all_20240703();

# ╔═╡ cd1ca226-00f6-4eb4-b10e-bb2a390b236d
Eeff_I = TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:I, TransitionPaths2024.onehot(TransitionPaths2024.Exp_20240703_sequences().sequences));

# ╔═╡ d93cf7b6-05a0-40d1-8f7b-905f4cd72487
Eeff_II = TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:II, TransitionPaths2024.onehot(TransitionPaths2024.Exp_20240703_sequences().sequences));

# ╔═╡ 383e921b-47d6-43ef-86f4-278ec9e0ce29
Eeff_IV = TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:IV, TransitionPaths2024.onehot(TransitionPaths2024.Exp_20240703_sequences().sequences));

# ╔═╡ c3ae9f17-68e0-45d1-a44c-ee4135ea44d7
Eeff_Glob = TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:global, TransitionPaths2024.onehot(TransitionPaths2024.Exp_20240703_sequences().sequences));

# ╔═╡ f07f8761-54e5-41a9-837a-268eb73b5e4e
response_data_sequences = TransitionPaths2024.Exp_20240703_sequences().sequences[indexin(all_response_data.response_data_WW_names, TransitionPaths2024.Exp_20240703_sequences().names)];

# ╔═╡ e3307dc4-638b-4cf8-af1e-b71088e9cedb
msa_local_rbm_likelihoods = stack(TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(rbm, TransitionPaths2024.onehot(TransitionPaths2024.Eugenio_MSA_20230419())) for rbm = (:I, :II, :IV));

# ╔═╡ ae2e577d-dfcb-4bd4-9e5c-41de0ce892a0
#msa_local_rbm_class_prediction = vec(last.(Tuple.(argmax(msa_local_rbm_likelihoods; dims=2))))
msa_local_rbm_class_prediction = TransitionPaths2024.Eugenio_MSA_Predicted_Classes_20240715().predicted_class;

# ╔═╡ d7959c5e-6c90-4995-9c34-b4de1c6765f1
_class_scatter_colors = (:cyan, :darkorange, :green); # :I, :II, :IV

# ╔═╡ ed61d8ed-751f-4a0b-96af-e5168fff09bc
df_20241213 = TransitionPaths2024.load_sequences_file_20241213();

# ╔═╡ 401dddb8-7acd-47d6-ac95-4e2fdce29361
df_20241213.aligned_sequences = TransitionPaths2024.Exp_20240703_sequences().sequences[indexin(df_20241213.ww_id, replace(TransitionPaths2024.Exp_20240703_sequences().names, "WW44" => "WW44/147", "WW147" => "WW44/147", "WW148" => "WW148/58", "WW149" => "WW9/149"))];

# ╔═╡ f93754f9-437c-486e-bf3f-1f4d433a19dc
response_data_sequences_new = TransitionPaths2024.Exp_20240703_sequences().sequences[indexin(df_20241213.ww_id, replace(
	TransitionPaths2024.Exp_20240703_sequences().names, "WW44" => "WW44/147", "WW147" => "WW44/147", "WW148" => "WW148/58", "WW149" => "WW9/149"
))];

# ╔═╡ f5d53e91-81ba-45d1-a92c-99e7813cb04f
unique(TransitionPaths2024.sampled_path_1to2rep2_20240703())

# ╔═╡ 0de36432-13f5-437a-b007-a23b56bed610
"EMGDWQEVWDENTGCYYYWNTQTNEVTWELP" == "EMGDWQEVWDENTGCYYYWNTQT-NEVTWELP"

# ╔═╡ 7a784b46-a25a-4d0e-8053-8286ec8991a8
unique(TransitionPaths2024.sampled_path_1to1_20240703())

# ╔═╡ 7ccab879-dd94-48fa-84ac-e4b0ce916ae1
"LPEGWEIRYTRE-GVRYFVDHNTRTTTFKDP" == "LPEGWEIRYTRE-GVRYFVDHNTRTTTFKDP"

# ╔═╡ 78e9c9a7-9280-4d2b-b495-4a725b5eee1c
md"# Main figures"

# ╔═╡ d5fb2c5a-e2dc-4f7b-9b19-e96787ddff71
md"""
## Fig. 2. I -> I

Path 1 $\rightarrow$ 1.
"""

# ╔═╡ 73169a6f-559e-46d3-af29-0e35bd487fb5
let fig = Makie.Figure()
	sampled_path = collect(unique(TransitionPaths2024.sampled_path_1to1_20240703()))
	path_response_data = TransitionPaths2024.response_data_path_1to1_20240703()

	C1_responses = df_20241213.C1[indexin(path_response_data.ww_names, df_20241213.ww_id)]
	C2_responses = df_20241213.C2[indexin(path_response_data.ww_names, df_20241213.ww_id)]
	C1_responses_err = df_20241213.eC1[indexin(path_response_data.ww_names, df_20241213.ww_id)]
	C2_responses_err = df_20241213.eC2[indexin(path_response_data.ww_names, df_20241213.ww_id)]
	
	probed_seqs = TransitionPaths2024.Exp_20240703_sequences().sequences[indexin(path_response_data.ww_names, TransitionPaths2024.Exp_20240703_sequences().names)]
	probed_seqs_idx::Vector{Int} = indexin(probed_seqs, sampled_path)

	ax = Makie.Axis(fig[1,2][1,1]; width=175, height=75, xgridvisible=false, ygridvisible=false, ylabel=L"\ln P_\mathrm{RBM}(\text{glob.})")
	Makie.scatterlines!(ax, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:global, TransitionPaths2024.onehot(sampled_path)); color=:black)
	Makie.scatter!(ax, probed_seqs_idx, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:global, TransitionPaths2024.onehot(probed_seqs)); color=:red, markersize=10)
	Makie.hidexdecorations!(ax)

	ax_I = Makie.Axis(
		fig[1,2][2,1]; width=175, height=75, xgridvisible=false, ygridvisible=false,
		ylabel=L"\log P_\mathrm{RBM}(\text{I})", #yticklabelcolor=_class_scatter_colors[1], ylabelcolor=_class_scatter_colors[1]
	)
	# ax_II = Makie.Axis(
	# 	fig[1,1][2,1]; width=175, height=75, xgridvisible=false, ygridvisible=false, yaxisposition=:right,
	# 	ylabel=L"\log P_\mathrm{RBM}(\text{II})", yticklabelcolor=_class_scatter_colors[2], ylabelcolor=_class_scatter_colors[2],
	# )
	#Makie.lines!(ax_II, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:II, TransitionPaths2024.onehot(sampled_path)); color=_class_scatter_colors[2], linewidth=0.75)
	Makie.scatterlines!(ax_I, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:I, TransitionPaths2024.onehot(sampled_path)); color=:black)
	Makie.scatter!(ax_I, probed_seqs_idx, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:I, TransitionPaths2024.onehot(probed_seqs)); color=:red, markersize=10)
	Makie.hidexdecorations!(ax)

	ax = Makie.Axis(fig[1,2][3,1]; width=175, height=75, xgridvisible=false, ygridvisible=false, xlabel="step", ylabel="Norm. resp. I", xticks=1:3:20)
	# Makie.scatterlines!(ax, probed_seqs_idx, path_response_data.C1_response; color=:black)
	# Makie.errorbars!(ax, probed_seqs_idx, path_response_data.C1_response, path_response_data.C1_response_err; color=:gray, direction=:y, linewidth=2, whiskerwidth=5)
	# Makie.scatter!(ax, probed_seqs_idx, path_response_data.C1_response; color=:red)

	Makie.scatterlines!(ax, probed_seqs_idx, C1_responses; color=:black)
	Makie.errorbars!(ax, probed_seqs_idx, C1_responses, C1_responses_err; color=:gray, direction=:y, linewidth=2, whiskerwidth=5)
	Makie.scatter!(ax, probed_seqs_idx, C1_responses; color=:red)

	Makie.hspan!(ax, -5, thresh_significance; color=(:red, 0.5))
	Makie.ylims!(ax, -2, 35)

	ax = Makie.Axis(fig[1,1], width=300, height=300, xlabel=L"I_{1}", ylabel=L"I_{2}", xgridvisible=false, ygridvisible=false) # h_38 -> h_1, h_36 -> h_2
	Makie.vlines!(ax, -1; color=:gray, linestyle=:dash, linewidth=1)
	Makie.hlines!(ax, -3; color=:gray, linestyle=:dash, linewidth=1)
	for (class, color) = zip((:I, :II, :IV), _class_scatter_colors)
		Makie.scatter!(
			ax,
			-inputs_h_from_v(TransitionPaths2024.Eugenio_RBM_20230419(:global), TransitionPaths2024.onehot(TransitionPaths2024.Eugenio_MSA_20230419()))[38, msa_local_rbm_class_prediction .== class],
			 inputs_h_from_v(TransitionPaths2024.Eugenio_RBM_20230419(:global), TransitionPaths2024.onehot(TransitionPaths2024.Eugenio_MSA_20230419()))[36, msa_local_rbm_class_prediction .== class];
			color=(color, 0.5), markersize=4
		)
	end

	@show length(sampled_path) length(probed_seqs)

	# Makie.scatterlines!(ax, 
	# 	-inputs_h_from_v(TransitionPaths2024.Eugenio_RBM_20230419(:global), TransitionPaths2024.onehot(sampled_path[[i for i=1:17 if i ∉ probed_seqs_idx]]))[38, :],
	# 	 inputs_h_from_v(TransitionPaths2024.Eugenio_RBM_20230419(:global), TransitionPaths2024.onehot(sampled_path[[i for i=1:17 if i ∉ probed_seqs_idx]]))[36, :];
	# 	markersize=[i ∉ [5, 7, 9, 12, 15] ? 15 : 6 for i = 1:17 if i ∉ probed_seqs_idx], 
	# 	color=[i ∉ [5, 7, 9, 12, 15] ? (:gray, 0.5) : :black for i = 1:17 if i ∉ probed_seqs_idx]
	# )

	Makie.scatterlines!(ax, 
		-inputs_h_from_v(TransitionPaths2024.Eugenio_RBM_20230419(:global), TransitionPaths2024.onehot(sampled_path))[38, :],
		 inputs_h_from_v(TransitionPaths2024.Eugenio_RBM_20230419(:global), TransitionPaths2024.onehot(sampled_path))[36, :];
		color=:black, markersize=10
	)

	Makie.scatter!(ax, 
		-inputs_h_from_v(TransitionPaths2024.Eugenio_RBM_20230419(:global), TransitionPaths2024.onehot(probed_seqs))[38, :],
		 inputs_h_from_v(TransitionPaths2024.Eugenio_RBM_20230419(:global), TransitionPaths2024.onehot(probed_seqs))[36, :];
		color=:red, markersize=7
	)
	
	Makie.arrows!(ax,
		-inputs_h_from_v(TransitionPaths2024.Eugenio_RBM_20230419(:global), TransitionPaths2024.onehot(sampled_path))[38, 1:1] .+ 0.35,
		 inputs_h_from_v(TransitionPaths2024.Eugenio_RBM_20230419(:global), TransitionPaths2024.onehot(sampled_path))[36, 1:1] .- 0.35,
		[-0.15], [0.15]; linewidth=2, color=:gray, arrowsize=10
	)
	
	Makie.xlims!(ax, -4.5, 2)
	Makie.ylims!(ax, -7, 3)

	Makie.resize_to_layout!(fig)
	#Makie.save("/DATA-SSD/cossio/projects/2024/Transition_Paths/TransitionPaths2024.jl/figures/Fig2.pdf", fig)
	fig
end

# ╔═╡ abe10f1a-c2a3-4596-89d4-80252932e3ca
Makie.scatterlines(-inputs_h_from_v(TransitionPaths2024.Eugenio_RBM_20230419(:global), TransitionPaths2024.onehot(unique(TransitionPaths2024.sampled_path_1to1_20240703())))[38, :])

# ╔═╡ f85a0a90-b33a-4533-ac39-6a497a74f685
Makie.scatterlines(-inputs_h_from_v(TransitionPaths2024.Eugenio_RBM_20230419(:global), TransitionPaths2024.onehot(unique(TransitionPaths2024.sampled_path_1to1_20240703())))[36, :])

# ╔═╡ 02a01dae-ac90-4e93-9351-716f9ec8b4f3
unique(TransitionPaths2024.sampled_path_1to1_20240703())[8:9]

# ╔═╡ 81b9f641-4a13-42cc-8d43-c1b2ae00fc2a
let fig = Makie.Figure()
	sampled_path = unique(TransitionPaths2024.sampled_path_1to1_20240703())
	path_response_data = TransitionPaths2024.response_data_path_1to1_20240703()

	C1_responses = df_20241213.C1[indexin(path_response_data.ww_names, df_20241213.ww_id)]
	C2_responses = df_20241213.C2[indexin(path_response_data.ww_names, df_20241213.ww_id)]
	C1_responses_err = df_20241213.eC1[indexin(path_response_data.ww_names, df_20241213.ww_id)]
	C2_responses_err = df_20241213.eC2[indexin(path_response_data.ww_names, df_20241213.ww_id)]
	
	probed_seqs = TransitionPaths2024.Exp_20240703_sequences().sequences[indexin(path_response_data.ww_names, TransitionPaths2024.Exp_20240703_sequences().names)]
	probed_seqs_idx::Vector{Int} = indexin(probed_seqs, sampled_path)

	_dummy_ax = Makie.Axis(fig[2,1:3]; height=10)
	Makie.hidedecorations!(_dummy_ax)
	Makie.hidespines!(_dummy_ax)

	ax = Makie.Axis(fig[3,1]; width=200, height=75, xgridvisible=false, ygridvisible=false, xlabel="step", ylabel=L"\ln P_\mathrm{RBM}(\text{glob.})")
	Makie.scatterlines!(ax, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:global, TransitionPaths2024.onehot(sampled_path)); color=:black)
	Makie.scatter!(ax, probed_seqs_idx, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:global, TransitionPaths2024.onehot(probed_seqs)); color=:red, markersize=10)
	#Makie.hidexdecorations!(ax)

	ax_I = Makie.Axis(
		fig[3,2]; width=200, height=75, xgridvisible=false, ygridvisible=false,
		xlabel="step", ylabel=L"\log P_\mathrm{RBM}(\text{I})", #yticklabelcolor=_class_scatter_colors[1], ylabelcolor=_class_scatter_colors[1]
	)
	Makie.scatterlines!(ax_I, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:I, TransitionPaths2024.onehot(sampled_path)); color=:black)
	Makie.scatter!(ax_I, probed_seqs_idx, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:I, TransitionPaths2024.onehot(probed_seqs)); color=:red, markersize=10)
	#Makie.hidexdecorations!(ax)

	ax = Makie.Axis(fig[3,3]; width=200, height=75, xgridvisible=false, ygridvisible=false, xlabel="step", ylabel="Norm. resp. I", xticks=1:3:20)
	Makie.scatterlines!(ax, probed_seqs_idx, C1_responses; color=:black)
	Makie.errorbars!(ax, probed_seqs_idx, C1_responses, C1_responses_err; color=:gray, direction=:y, linewidth=2, whiskerwidth=5)
	Makie.scatter!(ax, probed_seqs_idx, C1_responses; color=:red)
	Makie.hspan!(ax, -5, thresh_significance; color=(:red, 0.5))
	Makie.ylims!(ax, -2, 35)

	ax = Makie.Axis(fig[1,1], width=200, height=200, xlabel=L"I_{1}", ylabel=L"I_{2}", xgridvisible=false, ygridvisible=false) # h_38 -> h_1, h_36 -> h_2
	Makie.vlines!(ax, -1; color=:gray, linestyle=:dash, linewidth=1)
	Makie.hlines!(ax, -3; color=:gray, linestyle=:dash, linewidth=1)
	for (class, color) = zip((:I, :II, :IV), _class_scatter_colors)
		Makie.scatter!(
			ax,
			-inputs_h_from_v(TransitionPaths2024.Eugenio_RBM_20230419(:global), TransitionPaths2024.onehot(TransitionPaths2024.Eugenio_MSA_20230419()))[38, msa_local_rbm_class_prediction .== class],
			 inputs_h_from_v(TransitionPaths2024.Eugenio_RBM_20230419(:global), TransitionPaths2024.onehot(TransitionPaths2024.Eugenio_MSA_20230419()))[36, msa_local_rbm_class_prediction .== class];
			color=(color, 0.5), markersize=4
		)
	end

	Makie.scatterlines!(ax, 
		-inputs_h_from_v(TransitionPaths2024.Eugenio_RBM_20230419(:global), TransitionPaths2024.onehot(sampled_path))[38, :],
		 inputs_h_from_v(TransitionPaths2024.Eugenio_RBM_20230419(:global), TransitionPaths2024.onehot(sampled_path))[36, :];
		color=:black
	)

	Makie.scatter!(ax, 
		-inputs_h_from_v(TransitionPaths2024.Eugenio_RBM_20230419(:global), TransitionPaths2024.onehot(probed_seqs))[38, :],
		 inputs_h_from_v(TransitionPaths2024.Eugenio_RBM_20230419(:global), TransitionPaths2024.onehot(probed_seqs))[36, :];
		color=:red, markersize=10
	)
	
	Makie.arrows!(ax,
		-inputs_h_from_v(TransitionPaths2024.Eugenio_RBM_20230419(:global), TransitionPaths2024.onehot(sampled_path))[38, 1:1] .+ 0.35,
		 inputs_h_from_v(TransitionPaths2024.Eugenio_RBM_20230419(:global), TransitionPaths2024.onehot(sampled_path))[36, 1:1] .- 0.35,
		[-0.15], [0.15]; linewidth=2, color=:gray, arrowsize=10
	)
	
	Makie.xlims!(ax, -4.5, 2)
	Makie.ylims!(ax, -7, 3)

	Makie.resize_to_layout!(fig)
	#Makie.save("/DATA-SSD/cossio/projects/2024/Transition_Paths/TransitionPaths2024.jl/figures/Fig2.pdf", fig)
	fig
end

# ╔═╡ c1a63405-1e5a-4a9a-b7ca-63664460017f
md"## Fig. 3. I -> II"

# ╔═╡ 8de38839-0452-413b-ba6f-f09016cb4e39
let fig = Makie.Figure()
	path_response_data_new = TransitionPaths2024.artifact_20250531_new_tested_sequences_load()
	sampled_path = unique(TransitionPaths2024.sampled_path_1to2rep1_20240703())
	
	probed_seqs = map(LongAA, path_response_data_new.sequence)
	@assert probed_seqs ⊆ sampled_path
	probed_seqs_idx::Vector{Int} = indexin(probed_seqs, sampled_path)
	
	shuffled_path_response_data = TransitionPaths2024.response_data_path_1to2rev_20240703()
	shuffled_path_ww_names = shuffled_path_response_data.ww_names
	shuffled_path_seqs = TransitionPaths2024.Exp_20240703_sequences().sequences[indexin(shuffled_path_ww_names, TransitionPaths2024.Exp_20240703_sequences().names)]

	path_response_data_old = TransitionPaths2024.response_data_path_1to2rep1_20240703()
	probed_seqs_old = TransitionPaths2024.Exp_20240703_sequences().sequences[indexin(path_response_data_old.ww_names, TransitionPaths2024.Exp_20240703_sequences().names)]
	probed_seqs_idx_for_shuffled::Vector{Int} = indexin(probed_seqs_old, sampled_path)

	C1_responses = path_response_data_new.C1_response
	C2_responses = path_response_data_new.C2_response
	C1_responses_err = path_response_data_new.C1_response_err
	C2_responses_err = path_response_data_new.C2_response_err

	C1_responses_shuffled = df_20241213.C1[indexin(shuffled_path_ww_names, df_20241213.ww_id)]
	C2_responses_shuffled = df_20241213.C2[indexin(shuffled_path_ww_names, df_20241213.ww_id)]
	C1_responses_err_shuffled = df_20241213.eC1[indexin(shuffled_path_ww_names, df_20241213.ww_id)]
	C2_responses_err_shuffled = df_20241213.eC2[indexin(shuffled_path_ww_names, df_20241213.ww_id)]

	ax = Makie.Axis(fig[1,1][1,1]; width=300, height=120, xgridvisible=false, ygridvisible=false, ylabel=L"\log P_\mathrm{RBM}(\text{global})")
	Makie.scatterlines!(ax, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:global, TransitionPaths2024.onehot(sampled_path)); color=:black, label="Designed")
	Makie.scatter!(ax, probed_seqs_idx, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:global, TransitionPaths2024.onehot(probed_seqs)); color=:red, markersize=10)
	Makie.lines!(ax, probed_seqs_idx_for_shuffled, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:global, TransitionPaths2024.onehot(shuffled_path_seqs)); color=:brown, linestyle=:dash, label="Shuffled")
	Makie.scatter!(ax, probed_seqs_idx_for_shuffled, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:global, TransitionPaths2024.onehot(shuffled_path_seqs)); color=:brown)
	Makie.hidexdecorations!(ax)
	Makie.axislegend(ax; position=:lb, framevisible=false)

	ax = Makie.Axis(fig[1,1][2,1]; width=300, height=120, xgridvisible=false, ygridvisible=false, ylabel=L"\log P_\mathrm{RBM}(\text{I or II})")
	Makie.lines!(ax,  TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:I,  TransitionPaths2024.onehot(sampled_path)); color=_class_scatter_colors[1], label="I")
	Makie.lines!(ax, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:II, TransitionPaths2024.onehot(sampled_path)); color=_class_scatter_colors[2], label="II")
	Makie.scatter!(ax,  TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:I,  TransitionPaths2024.onehot(sampled_path)); color=_class_scatter_colors[1])
	Makie.scatter!(ax, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:II, TransitionPaths2024.onehot(sampled_path)); color=_class_scatter_colors[2])
	Makie.scatter!(ax,  probed_seqs_idx, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:I,  TransitionPaths2024.onehot(probed_seqs)); color=:red, markersize=10)
	Makie.scatter!(ax, probed_seqs_idx, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:II, TransitionPaths2024.onehot(probed_seqs)); color=:red, markersize=10)
	Makie.lines!(ax,  probed_seqs_idx_for_shuffled, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:I,  TransitionPaths2024.onehot(shuffled_path_seqs)); color=_class_scatter_colors[1], linestyle=:dash)
	Makie.lines!(ax, probed_seqs_idx_for_shuffled, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:II, TransitionPaths2024.onehot(shuffled_path_seqs)); color=_class_scatter_colors[2], linestyle=:dash)
	Makie.axislegend(ax; position=:lb, framevisible=false)

	# ax = Makie.Axis(fig[1,2][4,1]; width=300, height=30, xgridvisible=false, ygridvisible=false, yscale=log10, xticks=0:5:50)
	# Makie.hidedecorations!(ax)
	# Makie.hidespines!(ax)
	
	ax = Makie.Axis(fig[1,1][5,1]; width=300, height=120, xgridvisible=false, ygridvisible=false, ylabel="I/II rel. resp.", yscale=log10, xticks=0:5:50)
	Makie.hspan!(ax, 1, 1e6; color=(:cyan, 0.3))
	Makie.hspan!(ax, 1e-6, 1; color=(:darkorange, 0.3))
	Makie.scatterlines!(ax, probed_seqs_idx, C1_responses ./ C2_responses; color=:black)
	Makie.scatter!(ax, probed_seqs_idx, C1_responses ./ C2_responses; color=:red, markersize=10)
	Makie.ylims!(ax, 1e-3, 1e3)

	
	ax = Makie.Axis(fig[1,2][1,1], width=200, height=200, xlabel=L"I_{1}", ylabel=L"I_{2}", xgridvisible=false, ygridvisible=false) # h_38 -> h_1, h_36 -> h_2
	Makie.vlines!(ax, -1; color=:gray, linestyle=:dash, linewidth=1)
	Makie.hlines!(ax, -3; color=:gray, linestyle=:dash, linewidth=1)
	for (class, color) = zip((:I, :II, :IV), _class_scatter_colors)
		Makie.scatter!(
			ax,
			-inputs_h_from_v(TransitionPaths2024.Eugenio_RBM_20230419(:global), TransitionPaths2024.onehot(TransitionPaths2024.Eugenio_MSA_20230419()))[38, msa_local_rbm_class_prediction .== class],
			 inputs_h_from_v(TransitionPaths2024.Eugenio_RBM_20230419(:global), TransitionPaths2024.onehot(TransitionPaths2024.Eugenio_MSA_20230419()))[36, msa_local_rbm_class_prediction .== class];
			color=(color, 0.5), markersize=4
		)
	end

	Makie.scatterlines!(ax, 
		-inputs_h_from_v(TransitionPaths2024.Eugenio_RBM_20230419(:global), TransitionPaths2024.onehot(sampled_path))[38, :],
		 inputs_h_from_v(TransitionPaths2024.Eugenio_RBM_20230419(:global), TransitionPaths2024.onehot(sampled_path))[36, :];
		color=:black
	)

	Makie.scatter!(ax, 
		-inputs_h_from_v(TransitionPaths2024.Eugenio_RBM_20230419(:global), TransitionPaths2024.onehot(probed_seqs))[38, :],
		 inputs_h_from_v(TransitionPaths2024.Eugenio_RBM_20230419(:global), TransitionPaths2024.onehot(probed_seqs))[36, :];
		color=:red, markersize=10
	)
	
	Makie.arrows!(ax,
		-inputs_h_from_v(TransitionPaths2024.Eugenio_RBM_20230419(:global), TransitionPaths2024.onehot(sampled_path))[38, 1:1] .+ 0.35,
		 inputs_h_from_v(TransitionPaths2024.Eugenio_RBM_20230419(:global), TransitionPaths2024.onehot(sampled_path))[36, 1:1] .- 0.35,
		[-0.15], [0.15]; linewidth=2, color=:gray, arrowsize=10
	)
	
	Makie.xlims!(ax, -4.5, 2)
	Makie.ylims!(ax, -7, 3)

	ax = Makie.Axis(fig[1,2][2,1]; width=200, height=200, xgridvisible=false, ygridvisible=false, xlabel="Norm. response I", ylabel="Norm. response II")
	Makie.band!(ax, [-10, thresh_significance], [-10, -10], [thresh_significance, thresh_significance]; color=(:red, 0.3)) # neither specificity is significant region
	Makie.band!(ax, [-10, thresh_significance], [thresh_significance, thresh_significance], [100, 100]; color=(:darkorange, 0.3)) # specific to II
	Makie.band!(ax, [thresh_significance, 100], [-100, -100], [thresh_significance, thresh_significance]; color=(:cyan, 0.3)) # specific to II
	Makie.scatterlines!(ax, C1_responses, C2_responses; color=:black, linewidth=1)
	Makie.errorbars!(ax, C1_responses, C2_responses, C1_responses_err; color=:gray, direction=:x, linewidth=2, whiskerwidth=5)
	Makie.errorbars!(ax, C1_responses, C2_responses, C2_responses_err; color=:gray, direction=:y, linewidth=2, whiskerwidth=5)
	Makie.scatter!(ax, C1_responses, C2_responses; color=:red)
	Makie.arrows!(ax, C1_responses[1:1] .+ 0.2, C2_responses[1:1] .- 0.15, [-0.07], [0.07]; linewidth=2, color=:gray, arrowsize=10)
	Makie.xlims!(ax, -0.5, 13)
	Makie.ylims!(ax, -0.1, 4)

	ax_inset = Makie.Axis(
		fig[1,2][2,1], width=75, height=75, halign=:right, valign=:top, alignmode=Makie.Mixed(; right=10, top=10), xgridvisible=false, ygridvisible=false, title="Shuffled", bottomspinecolor=:brown, leftspinecolor=:brown, titlecolor=:brown, topspinevisible=false, rightspinevisible=false, titlegap=-15
	)
	Makie.band!(ax_inset, [-10, thresh_significance], [-10, -10], [thresh_significance, thresh_significance]; color=(:red, 0.3)) # neither specificity is significant region
	Makie.band!(ax_inset, [-10, thresh_significance], [thresh_significance, thresh_significance], [100, 100]; color=(:darkorange, 0.3)) # specific to II
	Makie.band!(ax_inset, [thresh_significance, 100], [-100, -100], [thresh_significance, thresh_significance]; color=(:cyan, 0.3)) # specific to II
	Makie.scatterlines!(ax_inset, C1_responses_shuffled, C2_responses_shuffled; color=:black, linewidth=1)
	Makie.errorbars!(ax_inset, C1_responses_shuffled, C2_responses_shuffled, C1_responses_err_shuffled; color=:gray, direction=:x, linewidth=2, whiskerwidth=5)
	Makie.errorbars!(ax_inset, C1_responses_shuffled, C2_responses_shuffled, C2_responses_err_shuffled; color=:gray, direction=:y, linewidth=2, whiskerwidth=5)
	Makie.scatter!(ax_inset, C1_responses_shuffled, C2_responses_shuffled; color=:red)
	Makie.xlims!(ax_inset, -0.5, 7)
	Makie.ylims!(ax_inset, -0.1, 1.2)
	Makie.translate!(ax_inset.scene, 0, 0, 10)
	
	Makie.resize_to_layout!(fig)
	#Makie.save(Base.homedir() * "/projects/2024/Transition_Paths/TransitionPaths2024.jl/figures/Fig3.pdf", fig)
	fig
end

# ╔═╡ b55b9ccd-2f0c-4d66-b929-68a17361dee4
let fig = Makie.Figure()
	path_response_data_new = TransitionPaths2024.artifact_20250531_new_tested_sequences_load()
	sampled_path = unique(TransitionPaths2024.sampled_path_1to2rep1_20240703())
	
	probed_seqs = map(LongAA, path_response_data_new.sequence)
	@assert probed_seqs ⊆ sampled_path
	probed_seqs_idx::Vector{Int} = indexin(probed_seqs, sampled_path)
	
	shuffled_path_response_data = TransitionPaths2024.response_data_path_1to2rev_20240703()
	shuffled_path_ww_names = shuffled_path_response_data.ww_names
	shuffled_path_seqs = TransitionPaths2024.Exp_20240703_sequences().sequences[indexin(shuffled_path_ww_names, TransitionPaths2024.Exp_20240703_sequences().names)]

	path_response_data_old = TransitionPaths2024.response_data_path_1to2rep1_20240703()
	probed_seqs_old = TransitionPaths2024.Exp_20240703_sequences().sequences[indexin(path_response_data_old.ww_names, TransitionPaths2024.Exp_20240703_sequences().names)]
	probed_seqs_idx_for_shuffled::Vector{Int} = indexin(probed_seqs_old, sampled_path)

	C1_responses = path_response_data_new.C1_response
	C2_responses = path_response_data_new.C2_response
	C1_responses_err = path_response_data_new.C1_response_err
	C2_responses_err = path_response_data_new.C2_response_err

	C1_responses_shuffled = df_20241213.C1[indexin(shuffled_path_ww_names, df_20241213.ww_id)]
	C2_responses_shuffled = df_20241213.C2[indexin(shuffled_path_ww_names, df_20241213.ww_id)]
	C1_responses_err_shuffled = df_20241213.eC1[indexin(shuffled_path_ww_names, df_20241213.ww_id)]
	C2_responses_err_shuffled = df_20241213.eC2[indexin(shuffled_path_ww_names, df_20241213.ww_id)]

	ax = Makie.Axis(fig[1,1]; width=500, height=500, xgridvisible=false, ygridvisible=false, xlabel="Norm. response I", ylabel="Norm. response II")
	Makie.band!(ax, [-10, thresh_significance], [-10, -10], [thresh_significance, thresh_significance]; color=(:red, 0.3)) # neither specificity is significant region
	Makie.band!(ax, [-10, thresh_significance], [thresh_significance, thresh_significance], [100, 100]; color=(:darkorange, 0.3)) # specific to II
	Makie.band!(ax, [thresh_significance, 100], [-100, -100], [thresh_significance, thresh_significance]; color=(:cyan, 0.3)) # specific to II
	Makie.scatterlines!(ax, C1_responses, C2_responses; color=:black, linewidth=1)
	Makie.text!(ax, tuple.(C1_responses, C2_responses); text=string.(probed_seqs_idx))
	# Makie.errorbars!(ax, C1_responses, C2_responses, C1_responses_err; color=:gray, direction=:x, linewidth=2, whiskerwidth=5)
	# Makie.errorbars!(ax, C1_responses, C2_responses, C2_responses_err; color=:gray, direction=:y, linewidth=2, whiskerwidth=5)
	Makie.scatter!(ax, C1_responses, C2_responses; color=:red)
	Makie.arrows!(ax, C1_responses[1:1] .+ 0.2, C2_responses[1:1] .- 0.15, [-0.07], [0.07]; linewidth=2, color=:gray, arrowsize=10)
	Makie.xlims!(ax, -0.5, 13)
	Makie.ylims!(ax, -0.1, 4)

	@show probed_seqs_idx

	Makie.resize_to_layout!(fig)
	fig
end

# ╔═╡ 0ca0ac10-2c96-4cc5-b6fe-ea2f142cba76
tuple.([1,2,3], [5,6,7])

# ╔═╡ 449467a6-8cde-4653-9402-557466ffadba
let fig = Makie.Figure()
	path_response_data_new = TransitionPaths2024.artifact_20250531_new_tested_sequences_load()
	sampled_path = unique(TransitionPaths2024.sampled_path_1to2rep1_20240703())
	
	probed_seqs = map(LongAA, path_response_data_new.sequence)
	@assert probed_seqs ⊆ sampled_path
	probed_seqs_idx::Vector{Int} = indexin(probed_seqs, sampled_path)
	
	shuffled_path_response_data = TransitionPaths2024.response_data_path_1to2rev_20240703()
	shuffled_path_ww_names = shuffled_path_response_data.ww_names
	shuffled_path_seqs = TransitionPaths2024.Exp_20240703_sequences().sequences[indexin(shuffled_path_ww_names, TransitionPaths2024.Exp_20240703_sequences().names)]

	path_response_data_old = TransitionPaths2024.response_data_path_1to2rep1_20240703()
	probed_seqs_old = TransitionPaths2024.Exp_20240703_sequences().sequences[indexin(path_response_data_old.ww_names, TransitionPaths2024.Exp_20240703_sequences().names)]
	probed_seqs_idx_for_shuffled::Vector{Int} = indexin(probed_seqs_old, sampled_path)

	C1_responses = path_response_data_new.C1_response
	C2_responses = path_response_data_new.C2_response
	C1_responses_err = path_response_data_new.C1_response_err
	C2_responses_err = path_response_data_new.C2_response_err

	C1_responses_shuffled = df_20241213.C1[indexin(shuffled_path_ww_names, df_20241213.ww_id)]
	C2_responses_shuffled = df_20241213.C2[indexin(shuffled_path_ww_names, df_20241213.ww_id)]
	C1_responses_err_shuffled = df_20241213.eC1[indexin(shuffled_path_ww_names, df_20241213.ww_id)]
	C2_responses_err_shuffled = df_20241213.eC2[indexin(shuffled_path_ww_names, df_20241213.ww_id)]

	ax = Makie.Axis(fig[1,1][1,1]; width=300, height=120, xgridvisible=false, ygridvisible=false, ylabel=L"\log P_\mathrm{RBM}(\text{global})")
	Makie.scatterlines!(ax, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:global, TransitionPaths2024.onehot(sampled_path)); color=:black, label="Designed")
	Makie.scatter!(ax, probed_seqs_idx, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:global, TransitionPaths2024.onehot(probed_seqs)); color=:red, markersize=10)
	Makie.lines!(ax, probed_seqs_idx_for_shuffled, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:global, TransitionPaths2024.onehot(shuffled_path_seqs)); color=:brown, linestyle=:dash, label="Shuffled")
	Makie.scatter!(ax, probed_seqs_idx_for_shuffled, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:global, TransitionPaths2024.onehot(shuffled_path_seqs)); color=:brown)
	Makie.hidexdecorations!(ax)
	Makie.axislegend(ax; position=:lb, framevisible=false)

	ax_I = Makie.Axis(
		fig[1,1][2,1]; width=300, height=120, xgridvisible=false, ygridvisible=false,
		ylabel=L"\log P_\mathrm{RBM}(\text{I or II})", yticklabelcolor=_class_scatter_colors[1], ylabelcolor=_class_scatter_colors[1]
	)
	ax_II = Makie.Axis(
		fig[1,1][2,1]; width=300, height=120, xgridvisible=false, ygridvisible=false, yaxisposition=:right,
		ylabel=L"\log P_\mathrm{RBM}(\text{II})", yticklabelcolor=_class_scatter_colors[2], ylabelcolor=_class_scatter_colors[2],
	)
	Makie.lines!(ax_I,  TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:I,  TransitionPaths2024.onehot(sampled_path)); color=_class_scatter_colors[1])
	Makie.lines!(ax_II, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:II, TransitionPaths2024.onehot(sampled_path)); color=_class_scatter_colors[2])
	Makie.scatter!(ax_I,  TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:I,  TransitionPaths2024.onehot(sampled_path)); color=_class_scatter_colors[1], label="Designed")
	Makie.scatter!(ax_II, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:II, TransitionPaths2024.onehot(sampled_path)); color=_class_scatter_colors[2])
	Makie.scatter!(ax_I,  probed_seqs_idx, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:I,  TransitionPaths2024.onehot(probed_seqs)); color=:red, markersize=10)
	Makie.scatter!(ax_II, probed_seqs_idx, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:II, TransitionPaths2024.onehot(probed_seqs)); color=:red, markersize=10)
	Makie.lines!(ax_I,  probed_seqs_idx_for_shuffled, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:I,  TransitionPaths2024.onehot(shuffled_path_seqs)); color=_class_scatter_colors[1], linestyle=:dash)
	Makie.lines!(ax_II, probed_seqs_idx_for_shuffled, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:II, TransitionPaths2024.onehot(shuffled_path_seqs)); color=_class_scatter_colors[2], linestyle=:dash)
	# Makie.scatter!(ax_I,  probed_seqs_idx, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:I,  TransitionPaths2024.onehot(shuffled_path_seqs)); color=:brown, label="Shuffled")
	# Makie.scatter!(ax_II, probed_seqs_idx, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:II, TransitionPaths2024.onehot(shuffled_path_seqs)); color=:brown)
	Makie.ylims!(ax_I, -130, 0)
	Makie.ylims!(ax_II, -130, 0)
	Makie.hidexdecorations!(ax_I)
	Makie.hidexdecorations!(ax_II)
	Makie.hidespines!(ax_II)

	# ax = Makie.Axis(fig[1,2][4,1]; width=300, height=30, xgridvisible=false, ygridvisible=false, yscale=log10, xticks=0:5:50)
	# Makie.hidedecorations!(ax)
	# Makie.hidespines!(ax)
	
	ax = Makie.Axis(fig[1,1][5,1]; width=300, height=120, xgridvisible=false, ygridvisible=false, ylabel="I/II rel. resp.", yscale=log10, xticks=0:5:50)
	Makie.hspan!(ax, 1, 1e6; color=(:cyan, 0.3))
	Makie.hspan!(ax, 1e-6, 1; color=(:darkorange, 0.3))
	Makie.scatterlines!(ax, probed_seqs_idx, C1_responses ./ C2_responses; color=:black)
	Makie.scatter!(ax, probed_seqs_idx, C1_responses ./ C2_responses; color=:red, markersize=10)
	Makie.ylims!(ax, 1e-3, 1e3)

	
	ax = Makie.Axis(fig[1,2][1,1], width=200, height=200, xlabel=L"I_{38}", ylabel=L"I_{36}", xgridvisible=false, ygridvisible=false)
	Makie.vlines!(ax, -1; color=:gray, linestyle=:dash, linewidth=1)
	Makie.hlines!(ax, -3; color=:gray, linestyle=:dash, linewidth=1)
	for (class, color) = zip((:I, :II, :IV), _class_scatter_colors)
		Makie.scatter!(
			ax,
			-inputs_h_from_v(TransitionPaths2024.Eugenio_RBM_20230419(:global), TransitionPaths2024.onehot(TransitionPaths2024.Eugenio_MSA_20230419()))[38, msa_local_rbm_class_prediction .== class],
			 inputs_h_from_v(TransitionPaths2024.Eugenio_RBM_20230419(:global), TransitionPaths2024.onehot(TransitionPaths2024.Eugenio_MSA_20230419()))[36, msa_local_rbm_class_prediction .== class];
			color=(color, 0.5), markersize=4
		)
	end

	Makie.scatterlines!(ax, 
		-inputs_h_from_v(TransitionPaths2024.Eugenio_RBM_20230419(:global), TransitionPaths2024.onehot(sampled_path))[38, :],
		 inputs_h_from_v(TransitionPaths2024.Eugenio_RBM_20230419(:global), TransitionPaths2024.onehot(sampled_path))[36, :];
		color=:black
	)

	Makie.scatter!(ax, 
		-inputs_h_from_v(TransitionPaths2024.Eugenio_RBM_20230419(:global), TransitionPaths2024.onehot(probed_seqs))[38, :],
		 inputs_h_from_v(TransitionPaths2024.Eugenio_RBM_20230419(:global), TransitionPaths2024.onehot(probed_seqs))[36, :];
		color=:red, markersize=10, #marker=:x
	)
	
	Makie.arrows!(ax,
		-inputs_h_from_v(TransitionPaths2024.Eugenio_RBM_20230419(:global), TransitionPaths2024.onehot(sampled_path))[38, 1:1] .+ 0.35,
		 inputs_h_from_v(TransitionPaths2024.Eugenio_RBM_20230419(:global), TransitionPaths2024.onehot(sampled_path))[36, 1:1] .- 0.35,
		[-0.15], [0.15]; linewidth=2, color=:gray, arrowsize=10
	)
	
	Makie.xlims!(ax, -4.5, 2)
	Makie.ylims!(ax, -7, 3)

	ax = Makie.Axis(fig[1,2][2,1]; width=200, height=200, xgridvisible=false, ygridvisible=false, xlabel="Norm. response I", ylabel="Norm. response II")
	Makie.band!(ax, [-10, thresh_significance], [-10, -10], [thresh_significance, thresh_significance]; color=(:red, 0.3)) # neither specificity is significant region
	Makie.band!(ax, [-10, thresh_significance], [thresh_significance, thresh_significance], [100, 100]; color=(:darkorange, 0.3)) # specific to II
	Makie.band!(ax, [thresh_significance, 100], [-100, -100], [thresh_significance, thresh_significance]; color=(:cyan, 0.3)) # specific to II
	Makie.scatterlines!(ax, C1_responses, C2_responses; color=:black, linewidth=1)
	Makie.errorbars!(ax, C1_responses, C2_responses, C1_responses_err; color=:gray, direction=:x, linewidth=2, whiskerwidth=5)
	Makie.errorbars!(ax, C1_responses, C2_responses, C2_responses_err; color=:gray, direction=:y, linewidth=2, whiskerwidth=5)
	Makie.scatter!(ax, C1_responses, C2_responses; color=:red)
	Makie.arrows!(ax, C1_responses[1:1] .+ 0.2, C2_responses[1:1] .- 0.15, [-0.07], [0.07]; linewidth=2, color=:gray, arrowsize=10)
	Makie.xlims!(ax, -0.5, 13)
	Makie.ylims!(ax, -0.1, 4)

	ax_inset = Makie.Axis(
		fig[1,2][2,1], width=75, height=75, halign=:right, valign=:top, alignmode=Makie.Mixed(; right=10, top=10), xgridvisible=false, ygridvisible=false, title="Shuffled", bottomspinecolor=:brown, leftspinecolor=:brown, titlecolor=:brown, topspinevisible=false, rightspinevisible=false, titlegap=-15
	)
	Makie.band!(ax_inset, [-10, thresh_significance], [-10, -10], [thresh_significance, thresh_significance]; color=(:red, 0.3)) # neither specificity is significant region
	Makie.band!(ax_inset, [-10, thresh_significance], [thresh_significance, thresh_significance], [100, 100]; color=(:darkorange, 0.3)) # specific to II
	Makie.band!(ax_inset, [thresh_significance, 100], [-100, -100], [thresh_significance, thresh_significance]; color=(:cyan, 0.3)) # specific to II
	Makie.scatterlines!(ax_inset, C1_responses_shuffled, C2_responses_shuffled; color=:black, linewidth=1)
	Makie.errorbars!(ax_inset, C1_responses_shuffled, C2_responses_shuffled, C1_responses_err_shuffled; color=:gray, direction=:x, linewidth=2, whiskerwidth=5)
	Makie.errorbars!(ax_inset, C1_responses_shuffled, C2_responses_shuffled, C2_responses_err_shuffled; color=:gray, direction=:y, linewidth=2, whiskerwidth=5)
	Makie.scatter!(ax_inset, C1_responses_shuffled, C2_responses_shuffled; color=:red)
	Makie.xlims!(ax_inset, -0.5, 7)
	Makie.ylims!(ax_inset, -0.1, 1.2)
	Makie.translate!(ax_inset.scene, 0, 0, 10)
	
	Makie.resize_to_layout!(fig)
	#Makie.save(Base.homedir() * "/projects/2024/Transition_Paths/TransitionPaths2024.jl/figures/Fig3.pdf", fig)
	fig
end

# ╔═╡ 6f35366a-b6c5-4b05-a236-4b13b0f3b546
md"## Fig. 4. I -> IV"

# ╔═╡ 606f2eea-b7f0-47b2-93cc-7fa482f342fe
let fig = Makie.Figure()
	sampled_path = unique(TransitionPaths2024.sampled_path_1to4_20240703())
	path_response_data = TransitionPaths2024.response_data_path_1to4_20240703()
	probed_seqs = TransitionPaths2024.Exp_20240703_sequences().sequences[indexin(path_response_data.ww_names, TransitionPaths2024.Exp_20240703_sequences().names)]
	probed_seqs_idx::Vector{Int} = indexin(probed_seqs, sampled_path)

	shuffled_path_response_data = TransitionPaths2024.response_data_path_1to4rev_20240703()
	shuffled_path_ww_names = shuffled_path_response_data.ww_names
	shuffled_path_seqs = TransitionPaths2024.Exp_20240703_sequences().sequences[indexin(shuffled_path_ww_names, TransitionPaths2024.Exp_20240703_sequences().names)]
	shuffled_probed_seqs_idx = probed_seqs_idx[[1:6; 9; 11]] # done by hand because there were some sequences added to experiments

	path_response_data_names = replace(path_response_data.ww_names,
		"WW147" => "WW44/147", "WW148" => "WW148/58", "WW149" => "WW9/149"
	)
	
	C1_responses = df_20241213.C1[indexin(path_response_data_names, df_20241213.ww_id)]
	C4_responses = df_20241213.C4[indexin(path_response_data_names, df_20241213.ww_id)]
	C1_responses_err = df_20241213.eC1[indexin(path_response_data_names, df_20241213.ww_id)]
	C4_responses_err = df_20241213.eC4[indexin(path_response_data_names, df_20241213.ww_id)]

	C1_responses_shuffled = df_20241213.C1[indexin(shuffled_path_ww_names, df_20241213.ww_id)]
	C4_responses_shuffled = df_20241213.C4[indexin(shuffled_path_ww_names, df_20241213.ww_id)]
	C1_responses_err_shuffled = df_20241213.eC1[indexin(shuffled_path_ww_names, df_20241213.ww_id)]
	C4_responses_err_shuffled = df_20241213.eC4[indexin(shuffled_path_ww_names, df_20241213.ww_id)]
	
	rel_response = C1_responses ./ C4_responses
	
	ax = Makie.Axis(fig[1,1][1,1]; width=300, height=100, xgridvisible=false, ygridvisible=false, ylabel=L"\log P_\mathrm{RBM}(\text{global})")
	Makie.scatterlines!(ax, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:global, TransitionPaths2024.onehot(sampled_path)); color=:black, label="Designed")
	Makie.scatter!(ax, probed_seqs_idx, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:global, TransitionPaths2024.onehot(probed_seqs)); color=:red, markersize=10)
	Makie.scatterlines!(ax, shuffled_probed_seqs_idx, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:global, TransitionPaths2024.onehot(shuffled_path_seqs)); color=:brown, linestyle=:dash, label="Shuffled")
	Makie.axislegend(ax; position=:lb, framevisible=false)
	Makie.hidexdecorations!(ax)

	ax = Makie.Axis(fig[1,1][2,1]; width=300, height=100, xgridvisible=false, ygridvisible=false, ylabel=L"\log P_\mathrm{RBM}(\text{I or IV})")
	Makie.lines!(ax,  TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:I,  TransitionPaths2024.onehot(sampled_path)); color=_class_scatter_colors[1], label="I")
	Makie.lines!(ax, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:IV, TransitionPaths2024.onehot(sampled_path)); color=_class_scatter_colors[3], label="IV")
	Makie.scatter!(ax,  TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:I,  TransitionPaths2024.onehot(sampled_path)); color=_class_scatter_colors[1])
	Makie.scatter!(ax, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:IV, TransitionPaths2024.onehot(sampled_path)); color=_class_scatter_colors[3])
	Makie.scatter!(ax,  probed_seqs_idx, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:I,  TransitionPaths2024.onehot(probed_seqs)); color=:red, markersize=10)
	Makie.scatter!(ax, probed_seqs_idx, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:IV, TransitionPaths2024.onehot(probed_seqs)); color=:red, markersize=10)
	Makie.lines!(ax,  shuffled_probed_seqs_idx, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:I,  TransitionPaths2024.onehot(shuffled_path_seqs)); color=_class_scatter_colors[1], linestyle=:dash)
	Makie.lines!(ax, shuffled_probed_seqs_idx, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:IV, TransitionPaths2024.onehot(shuffled_path_seqs)); color=_class_scatter_colors[3], linestyle=:dash)
	# Makie.scatter!(ax_I,  probed_seqs_idx, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:I,  TransitionPaths2024.onehot(shuffled_path_seqs)); color=:brown, label="Shuffled")
	# Makie.scatter!(ax_II, probed_seqs_idx, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:II, TransitionPaths2024.onehot(shuffled_path_seqs)); color=:brown)
	Makie.axislegend(ax; position=:lb, framevisible=false)

	ax = Makie.Axis(fig[1,1][3,1]; width=300, height=100, xgridvisible=false, ygridvisible=false, ylabel=L"\log P_\mathrm{RBM}(\text{I or II})")
	Makie.lines!(ax,  TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:I,  TransitionPaths2024.onehot(sampled_path)); color=_class_scatter_colors[1], label="I")
	Makie.lines!(ax, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:II, TransitionPaths2024.onehot(sampled_path)); color=_class_scatter_colors[2], label="II")
	Makie.scatter!(ax,  TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:I,  TransitionPaths2024.onehot(sampled_path)); color=_class_scatter_colors[1])
	Makie.scatter!(ax, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:II, TransitionPaths2024.onehot(sampled_path)); color=_class_scatter_colors[2])
	Makie.scatter!(ax,  probed_seqs_idx, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:I,  TransitionPaths2024.onehot(probed_seqs)); color=:red, markersize=10)
	Makie.scatter!(ax, probed_seqs_idx, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:II, TransitionPaths2024.onehot(probed_seqs)); color=:red, markersize=10)
	Makie.lines!(ax,  shuffled_probed_seqs_idx, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:I,  TransitionPaths2024.onehot(shuffled_path_seqs)); color=_class_scatter_colors[1], linestyle=:dash)
	Makie.lines!(ax, shuffled_probed_seqs_idx, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:II, TransitionPaths2024.onehot(shuffled_path_seqs)); color=_class_scatter_colors[2], linestyle=:dash)
	Makie.axislegend(ax; position=:lc, framevisible=false)
	
	ax = Makie.Axis(fig[1,1][4,1]; width=300, height=120, xgridvisible=false, ygridvisible=false, ylabel="I/IV rel. resp.", yscale=log10, xticks=0:5:50)
	Makie.hspan!(ax, 1, 1e6; color=(:cyan, 0.3))
	Makie.hspan!(ax, 1e-6, 1; color=(:green, 0.3))
	Makie.scatterlines!(ax, [i for (i,x) = zip(probed_seqs_idx, rel_response) if x > 0], [x for x = rel_response if x > 0]; color=:black)
	Makie.scatter!(ax, [i for (i,x) = zip(probed_seqs_idx, rel_response) if x > 0], [x for x = rel_response if x > 0]; color=:red, markersize=10)
	Makie.ylims!(ax, 1e-3, 5e3)


	ax = Makie.Axis(fig[1,2][1,1], width=200, height=200, xlabel=L"I_{38}", ylabel=L"I_{36}", xgridvisible=false, ygridvisible=false)
	Makie.vlines!(ax, -1; color=:gray, linestyle=:dash, linewidth=1)
	Makie.hlines!(ax, -3; color=:gray, linestyle=:dash, linewidth=1)
	for (class, color) = zip((:I, :II, :IV), _class_scatter_colors)
		Makie.scatter!(
			ax,
			-inputs_h_from_v(TransitionPaths2024.Eugenio_RBM_20230419(:global), TransitionPaths2024.onehot(TransitionPaths2024.Eugenio_MSA_20230419()))[38, msa_local_rbm_class_prediction .== class],
			 inputs_h_from_v(TransitionPaths2024.Eugenio_RBM_20230419(:global), TransitionPaths2024.onehot(TransitionPaths2024.Eugenio_MSA_20230419()))[36, msa_local_rbm_class_prediction .== class];
			color=(color, 0.5), markersize=4
		)
	end

	Makie.scatterlines!(ax, 
		-inputs_h_from_v(TransitionPaths2024.Eugenio_RBM_20230419(:global), TransitionPaths2024.onehot(sampled_path))[38, :],
		 inputs_h_from_v(TransitionPaths2024.Eugenio_RBM_20230419(:global), TransitionPaths2024.onehot(sampled_path))[36, :];
		color=:black
	)

	Makie.scatter!(ax, 
		-inputs_h_from_v(TransitionPaths2024.Eugenio_RBM_20230419(:global), TransitionPaths2024.onehot(probed_seqs))[38, :],
		 inputs_h_from_v(TransitionPaths2024.Eugenio_RBM_20230419(:global), TransitionPaths2024.onehot(probed_seqs))[36, :];
		color=:red, markersize=10
	)
	
	Makie.arrows!(ax,
		-inputs_h_from_v(TransitionPaths2024.Eugenio_RBM_20230419(:global), TransitionPaths2024.onehot(sampled_path))[38, 1:1] .+ 0.35,
		 inputs_h_from_v(TransitionPaths2024.Eugenio_RBM_20230419(:global), TransitionPaths2024.onehot(sampled_path))[36, 1:1] .- 0.35,
		[-0.15], [0.15]; linewidth=2, color=:gray, arrowsize=10
	)
	
	Makie.xlims!(ax, -4.5, 2)
	Makie.ylims!(ax, -7, 3)

	ax = Makie.Axis(fig[1,2][2,1]; width=200, height=200, xgridvisible=false, ygridvisible=false, xlabel="Norm. response I", ylabel="Norm. response IV")
	Makie.band!(ax, [-10, thresh_significance], [-10, -10], [thresh_significance, thresh_significance]; color=(:red, 0.3)) # neither specificity is significant region
	Makie.band!(ax, [-10, thresh_significance], [thresh_significance, thresh_significance], [100, 100]; color=(:green, 0.3)) # specific to IV
	Makie.band!(ax, [thresh_significance, 100], [-100, -100], [thresh_significance, thresh_significance]; color=(:cyan, 0.3)) # specific to I
	Makie.scatterlines!(ax, C1_responses, C4_responses; color=:black, linewidth=1)
	Makie.errorbars!(ax, C1_responses, C4_responses, C1_responses_err; color=:gray, direction=:x, linewidth=2, whiskerwidth=5)
	Makie.errorbars!(ax, C1_responses, C4_responses, C4_responses_err; color=:gray, direction=:y, linewidth=2, whiskerwidth=5)
	Makie.scatter!(ax, C1_responses, C4_responses; color=:red)
	Makie.arrows!(ax, C1_responses[1:1] .+ 0.2, C4_responses[1:1] .- 0.15, [-0.07], [0.07]; linewidth=2, color=:gray, arrowsize=10)
	Makie.xlims!(ax, -0.5, 12)
	Makie.ylims!(ax, -0.25, 2)

	ax_inset = Makie.Axis(
		fig[1,2][2,1], width=75, height=75, halign=:right, valign=:top, alignmode=Makie.Mixed(; right=10, top=10), xgridvisible=false, ygridvisible=false, title="Shuffled", bottomspinecolor=:brown, leftspinecolor=:brown, titlecolor=:brown, topspinevisible=false, rightspinevisible=false, titlegap=-15
	)
	Makie.band!(ax_inset, [-10, thresh_significance], [-10, -10], [thresh_significance, thresh_significance]; color=(:red, 0.3)) # neither specificity is significant region
	Makie.band!(ax_inset, [-10, thresh_significance], [thresh_significance, thresh_significance], [100, 100]; color=(:green, 0.3)) # specific to IV
	Makie.band!(ax_inset, [thresh_significance, 100], [-100, -100], [thresh_significance, thresh_significance]; color=(:cyan, 0.3)) # specific to I
	Makie.scatterlines!(ax_inset, C1_responses_shuffled, C4_responses_shuffled; color=:black, linewidth=1)
	Makie.errorbars!(ax_inset, C1_responses_shuffled, C4_responses_shuffled, C4_responses_err_shuffled; color=:gray, direction=:x, linewidth=2, whiskerwidth=5)
	Makie.errorbars!(ax_inset, C1_responses_shuffled, C4_responses_shuffled, C4_responses_err_shuffled; color=:gray, direction=:y, linewidth=2, whiskerwidth=5)
	Makie.scatter!(ax_inset, C1_responses_shuffled, C4_responses_shuffled; color=:red)
	Makie.xlims!(ax_inset, -0.5, 5)
	Makie.ylims!(ax_inset, -0.2, 2)
	Makie.translate!(ax_inset.scene, 0, 0, 10)	

	Makie.resize_to_layout!(fig)
#	Makie.save("/Users/jfdcd/projects/2024/Transition_Paths/TransitionPaths2024.jl/figures/Fig4.pdf", fig)
	fig
end

# ╔═╡ 71f7257d-6998-4833-98c7-5513e5b23b93
let fig = Makie.Figure()
	sampled_path = unique(TransitionPaths2024.sampled_path_1to4_20240703())
	path_response_data = TransitionPaths2024.response_data_path_1to4_20240703()
	probed_seqs = TransitionPaths2024.Exp_20240703_sequences().sequences[indexin(path_response_data.ww_names, TransitionPaths2024.Exp_20240703_sequences().names)]
	probed_seqs_idx::Vector{Int} = indexin(probed_seqs, sampled_path)

	shuffled_path_response_data = TransitionPaths2024.response_data_path_1to4rev_20240703()
	shuffled_path_ww_names = shuffled_path_response_data.ww_names
	shuffled_path_seqs = TransitionPaths2024.Exp_20240703_sequences().sequences[indexin(shuffled_path_ww_names, TransitionPaths2024.Exp_20240703_sequences().names)]
	shuffled_probed_seqs_idx = probed_seqs_idx[[1:6; 9; 11]] # done by hand because there were some sequences added to experiments

	path_response_data_names = replace(path_response_data.ww_names,
		"WW147" => "WW44/147", "WW148" => "WW148/58", "WW149" => "WW9/149"
	)
	
	C1_responses = df_20241213.C1[indexin(path_response_data_names, df_20241213.ww_id)]
	C4_responses = df_20241213.C4[indexin(path_response_data_names, df_20241213.ww_id)]
	C1_responses_err = df_20241213.eC1[indexin(path_response_data_names, df_20241213.ww_id)]
	C4_responses_err = df_20241213.eC4[indexin(path_response_data_names, df_20241213.ww_id)]

	C1_responses_shuffled = df_20241213.C1[indexin(shuffled_path_ww_names, df_20241213.ww_id)]
	C4_responses_shuffled = df_20241213.C4[indexin(shuffled_path_ww_names, df_20241213.ww_id)]
	C1_responses_err_shuffled = df_20241213.eC1[indexin(shuffled_path_ww_names, df_20241213.ww_id)]
	C4_responses_err_shuffled = df_20241213.eC4[indexin(shuffled_path_ww_names, df_20241213.ww_id)]
	
	rel_response = C1_responses ./ C4_responses
	
	ax = Makie.Axis(fig[1,1][1,1]; width=300, height=120, xgridvisible=false, ygridvisible=false, ylabel=L"\log P_\mathrm{RBM}(\text{global})")
	Makie.scatterlines!(ax, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:global, TransitionPaths2024.onehot(sampled_path)); color=:black, label="Designed")
	Makie.scatter!(ax, probed_seqs_idx, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:global, TransitionPaths2024.onehot(probed_seqs)); color=:red, markersize=10)
	Makie.scatterlines!(ax, shuffled_probed_seqs_idx, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:global, TransitionPaths2024.onehot(shuffled_path_seqs)); color=:brown, linestyle=:dash, label="Shuffled")
	Makie.axislegend(ax; position=:lb, framevisible=false)
	Makie.hidexdecorations!(ax)

	ax_I = Makie.Axis(
		fig[1,1][2,1]; width=300, height=120, xgridvisible=false, ygridvisible=false,
		ylabel=L"\log P_\mathrm{RBM}(\text{I})", yticklabelcolor=_class_scatter_colors[1], ylabelcolor=_class_scatter_colors[1]
	)
	ax_IV = Makie.Axis(
		fig[1,1][2,1]; width=300, height=120, xgridvisible=false, ygridvisible=false, yaxisposition=:right,
		ylabel=L"\log P_\mathrm{RBM}(\text{IV})", yticklabelcolor=_class_scatter_colors[3], ylabelcolor=_class_scatter_colors[3],
	)
	Makie.lines!(ax_I,  TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:I,  TransitionPaths2024.onehot(sampled_path)); color=_class_scatter_colors[1])
	Makie.lines!(ax_IV, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:II, TransitionPaths2024.onehot(sampled_path)); color=_class_scatter_colors[3])
	Makie.scatter!(ax_I,  TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:I,  TransitionPaths2024.onehot(sampled_path)); color=_class_scatter_colors[1], label="Designed")
	Makie.scatter!(ax_IV, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:II, TransitionPaths2024.onehot(sampled_path)); color=_class_scatter_colors[3])
	Makie.scatter!(ax_I,  probed_seqs_idx, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:I,  TransitionPaths2024.onehot(probed_seqs)); color=:red, markersize=10)
	Makie.scatter!(ax_IV, probed_seqs_idx, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:II, TransitionPaths2024.onehot(probed_seqs)); color=:red, markersize=10)
	Makie.lines!(ax_I,  shuffled_probed_seqs_idx, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:I,  TransitionPaths2024.onehot(shuffled_path_seqs)); color=_class_scatter_colors[1], linestyle=:dash)
	Makie.lines!(ax_IV, shuffled_probed_seqs_idx, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:II, TransitionPaths2024.onehot(shuffled_path_seqs)); color=_class_scatter_colors[3], linestyle=:dash)
	# Makie.scatter!(ax_I,  probed_seqs_idx, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:I,  TransitionPaths2024.onehot(shuffled_path_seqs)); color=:brown, label="Shuffled")
	# Makie.scatter!(ax_II, probed_seqs_idx, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:II, TransitionPaths2024.onehot(shuffled_path_seqs)); color=:brown)
	Makie.hidexdecorations!(ax_I)
	Makie.hidexdecorations!(ax_IV)
	Makie.hidespines!(ax_IV)

	
	# ax = Makie.Axis(fig[1,1][2,1]; width=300, height=100, xgridvisible=false, ygridvisible=false, ylabel=L"\log P_\mathrm{RBM}(\text{I})")
	# Makie.scatterlines!(ax, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:I, TransitionPaths2024.onehot(sampled_path)); color=:black, label="Designed")
	# Makie.scatter!(ax, probed_seqs_idx, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:I, TransitionPaths2024.onehot(probed_seqs)); color=:red, markersize=10)
	# Makie.scatterlines!(ax, shuffled_probed_seqs_idx, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:I, TransitionPaths2024.onehot(shuffled_path_seqs)); color=:brown, linestyle=:dash, label="Shuffled")
	# Makie.hidexdecorations!(ax)
	# Makie.axislegend(ax; position=:lb, framevisible=false)

	# ax = Makie.Axis(fig[1,1][3,1]; width=300, height=100, xgridvisible=false, ygridvisible=false, ylabel=L"\log P_\mathrm{RBM}(\text{IV})")
	# Makie.scatterlines!(ax, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:II, TransitionPaths2024.onehot(sampled_path)); color=:black)
	# Makie.scatter!(ax, probed_seqs_idx, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:II, TransitionPaths2024.onehot(probed_seqs)); color=:red, markersize=10)
	# Makie.scatterlines!(ax, shuffled_probed_seqs_idx, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:II, TransitionPaths2024.onehot(shuffled_path_seqs)); color=:brown, linestyle=:dash, label="Shuffled")
	# Makie.hidexdecorations!(ax)

	# ax = Makie.Axis(fig[1,1][4,1]; width=300, height=30, xgridvisible=false, ygridvisible=false, yscale=log10, xticks=0:5:50)
	# Makie.hidedecorations!(ax)
	# Makie.hidespines!(ax)

	ax = Makie.Axis(fig[1,1][5,1]; width=300, height=120, xgridvisible=false, ygridvisible=false, ylabel="I/IV rel. resp.", yscale=log10, xticks=0:5:50)
	Makie.hspan!(ax, 1, 1e6; color=(:cyan, 0.3))
	Makie.hspan!(ax, 1e-6, 1; color=(:green, 0.3))
	Makie.scatterlines!(ax, [i for (i,x) = zip(probed_seqs_idx, rel_response) if x > 0], [x for x = rel_response if x > 0]; color=:black)
	Makie.scatter!(ax, [i for (i,x) = zip(probed_seqs_idx, rel_response) if x > 0], [x for x = rel_response if x > 0]; color=:red, markersize=10)
	Makie.ylims!(ax, 1e-3, 5e3)


	ax = Makie.Axis(fig[1,2][1,1], width=200, height=200, xlabel=L"I_{38}", ylabel=L"I_{36}", xgridvisible=false, ygridvisible=false)
	Makie.vlines!(ax, -1; color=:gray, linestyle=:dash, linewidth=1)
	Makie.hlines!(ax, -3; color=:gray, linestyle=:dash, linewidth=1)
	for (class, color) = zip((:I, :II, :IV), _class_scatter_colors)
		Makie.scatter!(
			ax,
			-inputs_h_from_v(TransitionPaths2024.Eugenio_RBM_20230419(:global), TransitionPaths2024.onehot(TransitionPaths2024.Eugenio_MSA_20230419()))[38, msa_local_rbm_class_prediction .== class],
			 inputs_h_from_v(TransitionPaths2024.Eugenio_RBM_20230419(:global), TransitionPaths2024.onehot(TransitionPaths2024.Eugenio_MSA_20230419()))[36, msa_local_rbm_class_prediction .== class];
			color=(color, 0.5), markersize=4
		)
	end

	Makie.scatterlines!(ax, 
		-inputs_h_from_v(TransitionPaths2024.Eugenio_RBM_20230419(:global), TransitionPaths2024.onehot(sampled_path))[38, :],
		 inputs_h_from_v(TransitionPaths2024.Eugenio_RBM_20230419(:global), TransitionPaths2024.onehot(sampled_path))[36, :];
		color=:black
	)

	Makie.scatter!(ax, 
		-inputs_h_from_v(TransitionPaths2024.Eugenio_RBM_20230419(:global), TransitionPaths2024.onehot(probed_seqs))[38, :],
		 inputs_h_from_v(TransitionPaths2024.Eugenio_RBM_20230419(:global), TransitionPaths2024.onehot(probed_seqs))[36, :];
		color=:red, markersize=10
	)
	
	Makie.arrows!(ax,
		-inputs_h_from_v(TransitionPaths2024.Eugenio_RBM_20230419(:global), TransitionPaths2024.onehot(sampled_path))[38, 1:1] .+ 0.35,
		 inputs_h_from_v(TransitionPaths2024.Eugenio_RBM_20230419(:global), TransitionPaths2024.onehot(sampled_path))[36, 1:1] .- 0.35,
		[-0.15], [0.15]; linewidth=2, color=:gray, arrowsize=10
	)
	
	Makie.xlims!(ax, -4.5, 2)
	Makie.ylims!(ax, -7, 3)

	ax = Makie.Axis(fig[1,2][2,1]; width=200, height=200, xgridvisible=false, ygridvisible=false, xlabel="Norm. response I", ylabel="Norm. response IV")
	Makie.band!(ax, [-10, thresh_significance], [-10, -10], [thresh_significance, thresh_significance]; color=(:red, 0.3)) # neither specificity is significant region
	Makie.band!(ax, [-10, thresh_significance], [thresh_significance, thresh_significance], [100, 100]; color=(:green, 0.3)) # specific to IV
	Makie.band!(ax, [thresh_significance, 100], [-100, -100], [thresh_significance, thresh_significance]; color=(:cyan, 0.3)) # specific to I
	Makie.scatterlines!(ax, C1_responses, C4_responses; color=:black, linewidth=1)
	Makie.errorbars!(ax, C1_responses, C4_responses, C1_responses_err; color=:gray, direction=:x, linewidth=2, whiskerwidth=5)
	Makie.errorbars!(ax, C1_responses, C4_responses, C4_responses_err; color=:gray, direction=:y, linewidth=2, whiskerwidth=5)
	Makie.scatter!(ax, C1_responses, C4_responses; color=:red)
	Makie.arrows!(ax, C1_responses[1:1] .+ 0.2, C4_responses[1:1] .- 0.15, [-0.07], [0.07]; linewidth=2, color=:gray, arrowsize=10)
	Makie.xlims!(ax, -0.5, 12)
	Makie.ylims!(ax, -0.25, 2)

	ax_inset = Makie.Axis(
		fig[1,2][2,1], width=75, height=75, halign=:right, valign=:top, alignmode=Makie.Mixed(; right=10, top=10), xgridvisible=false, ygridvisible=false, title="Shuffled", bottomspinecolor=:brown, leftspinecolor=:brown, titlecolor=:brown, topspinevisible=false, rightspinevisible=false, titlegap=-15
	)
	Makie.band!(ax_inset, [-10, thresh_significance], [-10, -10], [thresh_significance, thresh_significance]; color=(:red, 0.3)) # neither specificity is significant region
	Makie.band!(ax_inset, [-10, thresh_significance], [thresh_significance, thresh_significance], [100, 100]; color=(:green, 0.3)) # specific to IV
	Makie.band!(ax_inset, [thresh_significance, 100], [-100, -100], [thresh_significance, thresh_significance]; color=(:cyan, 0.3)) # specific to I
	Makie.scatterlines!(ax_inset, C1_responses_shuffled, C4_responses_shuffled; color=:black, linewidth=1)
	Makie.errorbars!(ax_inset, C1_responses_shuffled, C4_responses_shuffled, C4_responses_err_shuffled; color=:gray, direction=:x, linewidth=2, whiskerwidth=5)
	Makie.errorbars!(ax_inset, C1_responses_shuffled, C4_responses_shuffled, C4_responses_err_shuffled; color=:gray, direction=:y, linewidth=2, whiskerwidth=5)
	Makie.scatter!(ax_inset, C1_responses_shuffled, C4_responses_shuffled; color=:red)
	Makie.xlims!(ax_inset, -0.5, 5)
	Makie.ylims!(ax_inset, -0.2, 2)
	Makie.translate!(ax_inset.scene, 0, 0, 10)	

	Makie.resize_to_layout!(fig)
	#Makie.save("/DATA-SSD/cossio/projects/2024/Transition_Paths/TransitionPaths2024.jl/figures/Fig4.pdf", fig)
	fig
end

# ╔═╡ 9a991c6a-8479-4e3d-90f4-c88f7d298c99
md"## Fig. 5. Correlations"

# ╔═╡ e5fb8b94-2b1a-4544-88e1-3636bd4798aa
findall(TransitionPaths2024.artifact_20250531_new_tested_sequences_load().seq_name .∈ Ref(["Seq84", "Seq85", "Seq86", "Seq87"]))

# ╔═╡ 0ae8f8b8-225d-4022-8971-4fb187c726d3
TransitionPaths2024.artifact_20250531_new_tested_sequences_load().sequence[6:9] ∩ map(string, response_data_sequences_new)

# ╔═╡ 0ebe69bb-00c4-4381-ab36-770b27ae356e
TransitionPaths2024.artifact_20250531_new_tested_sequences_load().sequence

# ╔═╡ d4216956-091a-4780-94ca-6077ee29a1e6
	# path_response_data_new = TransitionPaths2024.artifact_20250531_new_tested_sequences_load()
	# sampled_path = unique(TransitionPaths2024.sampled_path_1to2rep1_20240703())
	
	# probed_seqs = map(LongAA, path_response_data_new.sequence)
	# @assert probed_seqs ⊆ sampled_path
	# probed_seqs_idx::Vector{Int} = indexin(probed_seqs, sampled_path)
	
	# shuffled_path_response_data = TransitionPaths2024.response_data_path_1to2rev_20240703()
	# shuffled_path_ww_names = shuffled_path_response_data.ww_names
	# shuffled_path_seqs = TransitionPaths2024.Exp_20240703_sequences().sequences[indexin(shuffled_path_ww_names, TransitionPaths2024.Exp_20240703_sequences().names)]

	# path_response_data_old = TransitionPaths2024.response_data_path_1to2rep1_20240703()
	# probed_seqs_old = TransitionPaths2024.Exp_20240703_sequences().sequences[indexin(path_response_data.ww_names, TransitionPaths2024.Exp_20240703_sequences().names)]
	# probed_seqs_idx_for_shuffled::Vector{Int} = indexin(probed_seqs_old, sampled_path)


# ╔═╡ bb6e6dc0-9e56-4c5f-80b6-89637573b4eb
let fig = Makie.Figure()
	# add new tested sequences part of Path I -> II
	response_data_sequences_new_new = [response_data_sequences_new; LongAA.(TransitionPaths2024.artifact_20250531_new_tested_sequences_load().sequence[6:9])]
	
	_sz = 150

	pred_class = TransitionPaths2024.Eugenio_Predict_Class_from_Inputs_20240715(TransitionPaths2024.onehot(response_data_sequences_new_new))
	likelihood_I  = TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:I,  TransitionPaths2024.onehot(response_data_sequences_new_new))
	likelihood_II = TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:II, TransitionPaths2024.onehot(response_data_sequences_new_new))
	likelihood_IV = TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:IV, TransitionPaths2024.onehot(response_data_sequences_new_new))
	response_I  = [ismissing(r) || r ≤ 0 ? NaN : r for r = [df_20241213.C1; TransitionPaths2024.artifact_20250531_new_tested_sequences_load().C1_response[6:9]]]
	response_II = [ismissing(r) || r ≤ 0 ? NaN : r for r = [df_20241213.C2; TransitionPaths2024.artifact_20250531_new_tested_sequences_load().C2_response[6:9]]]
	response_IV = [ismissing(r) || r ≤ 0 ? NaN : r for r = [df_20241213.C4; TransitionPaths2024.artifact_20250531_new_tested_sequences_load().C4_response[6:9]]]
	
	ax_I  = Makie.Axis(fig[1,1][1,1], width=_sz, height=_sz, xlabel=L"\ln P_\mathrm{rbm}(\text{local I})",  ylabel="Class I. norm. resp.",  yscale=log10, xgridvisible=false, ygridvisible=false, xticks=-100:50:0)
	ax_II = Makie.Axis(fig[1,1][1,2], width=_sz, height=_sz, xlabel=L"\ln P_\mathrm{rbm}(\text{local II})", ylabel="Class II. norm. resp.", yscale=log10, xgridvisible=false, ygridvisible=false, xticks=-80:20:0)
	ax_IV = Makie.Axis(fig[1,1][1,3], width=_sz, height=_sz, xlabel=L"\ln P_\mathrm{rbm}(\text{local IV})", ylabel="Class IV. norm. resp.", yscale=log10, xgridvisible=false, ygridvisible=false, xticks=-100:25:0)

	for (class, color) = zip((:I, :II, :IV), _class_scatter_colors)
		Makie.scatter!(ax_I,  likelihood_I[pred_class .== class],  response_I[pred_class .== class]; color)
		Makie.scatter!(ax_II, likelihood_II[pred_class .== class], response_II[pred_class .== class]; color)
		Makie.scatter!(ax_IV, likelihood_IV[pred_class .== class], response_IV[pred_class .== class]; color)
	end

	labels = ["I", "II", "IV"]
	elements = [Makie.MarkerElement(; marker=:circle, color = _class_scatter_colors[j]) for j = 1:3]
	Makie.Legend(fig[1,1][1,4], elements, labels; framevisible=false)

	# Correlations barplot
	responses = Dict(
		:I => [df_20241213.C1; TransitionPaths2024.artifact_20250531_new_tested_sequences_load().C1_response[6:9]],
		:II => [df_20241213.C2; TransitionPaths2024.artifact_20250531_new_tested_sequences_load().C2_response[6:9]],
		:IV => [df_20241213.C4; TransitionPaths2024.artifact_20250531_new_tested_sequences_load().C4_response[6:9]]
	)

	cat = [i for (i, target) = enumerate((:I, :II, :IV)) for (j, rbm) = enumerate((:global, :I, :II, :IV))]
	grp = [j for (i, target) = enumerate((:I, :II, :IV)) for (j, rbm) = enumerate((:global, :I, :II, :IV))]
				
	# correlations = [
	# 	nancor(TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(rbm, TransitionPaths2024.onehot(response_data_sequences_new)), replace(responses[target], missing => NaN))
	# 	for target = (:I, :II, :IV) for rbm = (:global, :I, :II, :IV)
	# ]

	response_data_sequences_new_uniq = unique(response_data_sequences_new_new)
	responses_uniq = Dict(target => last.(unique(first, zip(response_data_sequences_new_new, responses[target]))) for target = (:I, :II, :IV))
	
	correlations = [
		nancor(
			TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(rbm, TransitionPaths2024.onehot(response_data_sequences_new_new)),
			[log(ifelse(ismissing(r) || isnan(r) || r ≤ 0, NaN, r)) for r = responses[target]]
			#log.(max.(replace(responses[target], missing => NaN), 1e-10))
		)
		for target = (:I, :II, :IV) for rbm = (:global, :I, :II, :IV)
	]
			
	correlations_uniq = [
		nancor(
			TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(rbm, TransitionPaths2024.onehot(response_data_sequences_new_uniq)),
			[log(ifelse(ismissing(r) || isnan(r) || r ≤ 0, NaN, r)) for r = responses_uniq[target]]
		)
		for target = (:I, :II, :IV) for rbm = (:global, :I, :II, :IV)
	]

	
	colors = [:gray, :cyan, :orange, :green]
		
	ax = Makie.Axis(fig[2,1][1,2], width=200, height=_sz, xticks=(1:3, ["I", "II", "IV"]), ylabel="Exp. resp. corr.", xgridvisible=false, ygridvisible=false)
	Makie.barplot!(ax, cat, correlations_uniq; dodge=grp, color=colors[grp], dodge_gap=0.05, width=0.9)
	Makie.hlines!(ax, 0.0; color=:black, linestyle=:dash)
	Makie.ylims!(-0.4, 0.8)

	labels = ["RBM(glob.)", "RBM(I)", "RBM(II)", "RBM(IV)"]
	elements = [Makie.PolyElement(; polycolor = colors[j]) for j = 1:4]
	Makie.Legend(fig[2,1][1,3], elements, labels; framevisible=false)

	# response_data_v2_sequences = TransitionPaths2024.Exp_20240703_sequences().sequences[
	# 	indexin(TransitionPaths2024.response_data_all_v2_20240708().name, TransitionPaths2024.Exp_20240703_sequences().names)
	# ]
	# rbm_log_likelihoods_v2 = TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:global, TransitionPaths2024.onehot(response_data_v2_sequences))
	
	ax = Makie.Axis(fig[2,1][1,1], width=250, height=_sz, xlabel=L"P_{\mathrm{RBM}}(\text{glob.})", ylabel="Frequency", xgridvisible=false, ygridvisible=false)
	bins = -90:10:10
	Makie.hist!(ax, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:global, TransitionPaths2024.onehot(TransitionPaths2024.Eugenio_MSA_20230419())); normalization=:pdf, bins=-90:2:10, color=(:green, 0.3), label="MSA")
	# Makie.stephist!(ax, rbm_log_likelihoods_v2[TransitionPaths2024.response_data_all_v2_20240708().type .== "Natural"], color=:green, label="Natural", normalization=:pdf, bins=-90:5:10, linewidth=4)
	# Makie.stephist!(ax, rbm_log_likelihoods_v2[TransitionPaths2024.response_data_all_v2_20240708().type .== "Designed"], color=:blue, label="Designed", normalization=:pdf, bins=-90:5:10, linewidth=4)
	# Makie.stephist!(ax, rbm_log_likelihoods_v2[TransitionPaths2024.response_data_all_v2_20240708().type .== "Shuffled"], color=:red, label="Shuffled", normalization=:pdf, bins=-90:5:10, linewidth=4)
	Makie.stephist!(ax, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:global, TransitionPaths2024.onehot(unique(df_20241213.aligned_sequences[df_20241213.group .== "Wild-types"]))); color=:green, label="Natural", normalization=:pdf, bins, linewidth=4)
	Makie.stephist!(ax, 
		TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:global, TransitionPaths2024.onehot(
			unique(df_20241213.aligned_sequences[df_20241213.group .∈ Ref(("1st", "2nd", "3rd", "4th") .* "_batch")])));
		color=:blue, label="Designed", normalization=:pdf, bins, linewidth=4)
	Makie.stephist!(ax, 
		TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:global, TransitionPaths2024.onehot(
			unique(df_20241213.aligned_sequences[df_20241213.group .== "Scrambled_paths"])));
		color=:red, label="Shuffled", normalization=:pdf, bins, linewidth=2, linestyle=(:dash, 1.2))
	Makie.xlims!(ax, -90, 0)
	Makie.ylims!(ax, -0.001, 0.1)
	Makie.axislegend(ax; position=:lt, framevisible=false)

	Makie.resize_to_layout!(fig)
	#Makie.save(homedir() * "/projects/2024/Transition_Paths/TransitionPaths2024.jl/figures/Fig5.pdf", fig)
	fig
end

# ╔═╡ 10489d92-333c-41e6-a3b4-5f29e3a3689a
let fig = Makie.Figure()
	_sz = 150

	pred_class = TransitionPaths2024.Eugenio_Predict_Class_from_Inputs_20240715(TransitionPaths2024.onehot(response_data_sequences_new))
	likelihood_I  = TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:I,  TransitionPaths2024.onehot(response_data_sequences_new))
	likelihood_II = TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:II, TransitionPaths2024.onehot(response_data_sequences_new))
	likelihood_IV = TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:IV, TransitionPaths2024.onehot(response_data_sequences_new))
	response_I  = [ismissing(r) || r ≤ 0 ? NaN : r for r = df_20241213.C1]
	response_II = [ismissing(r) || r ≤ 0 ? NaN : r for r = df_20241213.C2]
	response_IV = [ismissing(r) || r ≤ 0 ? NaN : r for r = df_20241213.C4]
	
	ax_I  = Makie.Axis(fig[1,1][1,1], width=_sz, height=_sz, xlabel=L"\ln P_\mathrm{rbm}(\text{local I})",  ylabel="Class I. norm. resp.",  yscale=log10, xgridvisible=false, ygridvisible=false, xticks=-100:50:0)
	ax_II = Makie.Axis(fig[1,1][1,2], width=_sz, height=_sz, xlabel=L"\ln P_\mathrm{rbm}(\text{local II})", ylabel="Class II. norm. resp.", yscale=log10, xgridvisible=false, ygridvisible=false, xticks=-80:20:0)
	ax_IV = Makie.Axis(fig[1,1][1,3], width=_sz, height=_sz, xlabel=L"\ln P_\mathrm{rbm}(\text{local IV})", ylabel="Class IV. norm. resp.", yscale=log10, xgridvisible=false, ygridvisible=false, xticks=-100:25:0)

	for (class, color) = zip((:I, :II, :IV), _class_scatter_colors)
		Makie.scatter!(ax_I,  likelihood_I[pred_class .== class],  response_I[pred_class .== class]; color)
		Makie.scatter!(ax_II, likelihood_II[pred_class .== class], response_II[pred_class .== class]; color)
		Makie.scatter!(ax_IV, likelihood_IV[pred_class .== class], response_IV[pred_class .== class]; color)
	end

	labels = ["I", "II", "IV"]
	elements = [Makie.MarkerElement(; marker=:circle, color = _class_scatter_colors[j]) for j = 1:3]
	Makie.Legend(fig[1,1][1,4], elements, labels; framevisible=false)

	# Correlations barplot
	responses = Dict(:I => df_20241213.C1, :II => df_20241213.C2, :IV => df_20241213.C4)

	cat = [i for (i, target) = enumerate((:I, :II, :IV)) for (j, rbm) = enumerate((:global, :I, :II, :IV))]
	grp = [j for (i, target) = enumerate((:I, :II, :IV)) for (j, rbm) = enumerate((:global, :I, :II, :IV))]
				
	# correlations = [
	# 	nancor(TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(rbm, TransitionPaths2024.onehot(response_data_sequences_new)), replace(responses[target], missing => NaN))
	# 	for target = (:I, :II, :IV) for rbm = (:global, :I, :II, :IV)
	# ]

	response_data_sequences_new_uniq = unique(response_data_sequences_new)
	responses_uniq = Dict(target => last.(unique(first, zip(response_data_sequences_new, responses[target]))) for target = (:I, :II, :IV))
	
	correlations = [
		nancor(
			TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(rbm, TransitionPaths2024.onehot(response_data_sequences_new)),
			[log(ifelse(ismissing(r) || isnan(r) || r ≤ 0, NaN, r)) for r = responses[target]]
			#log.(max.(replace(responses[target], missing => NaN), 1e-10))
		)
		for target = (:I, :II, :IV) for rbm = (:global, :I, :II, :IV)
	]
			
	correlations_uniq = [
		nancor(
			TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(rbm, TransitionPaths2024.onehot(response_data_sequences_new_uniq)),
			[log(ifelse(ismissing(r) || isnan(r) || r ≤ 0, NaN, r)) for r = responses_uniq[target]]
		)
		for target = (:I, :II, :IV) for rbm = (:global, :I, :II, :IV)
	]

	
	colors = [:gray, :cyan, :orange, :green]
		
	ax = Makie.Axis(fig[2,1][1,2], width=200, height=_sz, xticks=(1:3, ["I", "II", "IV"]), ylabel="Exp. resp. corr.", xgridvisible=false, ygridvisible=false)
	Makie.barplot!(ax, cat, correlations_uniq; dodge=grp, color=colors[grp], dodge_gap=0.05, width=0.9)
	Makie.hlines!(ax, 0.0; color=:black, linestyle=:dash)
	Makie.ylims!(-0.4, 0.8)

	labels = ["RBM(glob.)", "RBM(I)", "RBM(II)", "RBM(IV)"]
	elements = [Makie.PolyElement(; polycolor = colors[j]) for j = 1:4]
	Makie.Legend(fig[2,1][1,3], elements, labels; framevisible=false)

	# response_data_v2_sequences = TransitionPaths2024.Exp_20240703_sequences().sequences[
	# 	indexin(TransitionPaths2024.response_data_all_v2_20240708().name, TransitionPaths2024.Exp_20240703_sequences().names)
	# ]
	# rbm_log_likelihoods_v2 = TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:global, TransitionPaths2024.onehot(response_data_v2_sequences))
	
	ax = Makie.Axis(fig[2,1][1,1], width=250, height=_sz, xlabel=L"P_{\mathrm{RBM}}(\text{glob.})", ylabel="Frequency", xgridvisible=false, ygridvisible=false)
	bins = -90:10:10
	Makie.hist!(ax, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:global, TransitionPaths2024.onehot(TransitionPaths2024.Eugenio_MSA_20230419())); normalization=:pdf, bins=-90:2:10, color=(:green, 0.3), label="MSA")
	# Makie.stephist!(ax, rbm_log_likelihoods_v2[TransitionPaths2024.response_data_all_v2_20240708().type .== "Natural"], color=:green, label="Natural", normalization=:pdf, bins=-90:5:10, linewidth=4)
	# Makie.stephist!(ax, rbm_log_likelihoods_v2[TransitionPaths2024.response_data_all_v2_20240708().type .== "Designed"], color=:blue, label="Designed", normalization=:pdf, bins=-90:5:10, linewidth=4)
	# Makie.stephist!(ax, rbm_log_likelihoods_v2[TransitionPaths2024.response_data_all_v2_20240708().type .== "Shuffled"], color=:red, label="Shuffled", normalization=:pdf, bins=-90:5:10, linewidth=4)
	Makie.stephist!(ax, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:global, TransitionPaths2024.onehot(unique(df_20241213.aligned_sequences[df_20241213.group .== "Wild-types"]))); color=:green, label="Natural", normalization=:pdf, bins, linewidth=4)
	Makie.stephist!(ax, 
		TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:global, TransitionPaths2024.onehot(
			unique(df_20241213.aligned_sequences[df_20241213.group .∈ Ref(("1st", "2nd", "3rd", "4th") .* "_batch")])));
		color=:blue, label="Designed", normalization=:pdf, bins, linewidth=4)
	Makie.stephist!(ax, 
		TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:global, TransitionPaths2024.onehot(
			unique(df_20241213.aligned_sequences[df_20241213.group .== "Scrambled_paths"])));
		color=:red, label="Shuffled", normalization=:pdf, bins, linewidth=2, linestyle=(:dash, 1.2))
	Makie.xlims!(ax, -90, 0)
	Makie.ylims!(ax, -0.001, 0.1)
	Makie.axislegend(ax; position=:lt, framevisible=false)

	Makie.resize_to_layout!(fig)
	#Makie.save("/DATA-SSD/cossio/projects/2024/Transition_Paths/TransitionPaths2024.jl/figures/Fig5.pdf", fig)
	fig
end

# ╔═╡ cf6d18bc-90da-4576-a8a4-1ff412124bdb
md"### with Spearman"

# ╔═╡ 92da5b0e-75e1-4aaa-911d-f375a36174f4
let fig = Makie.Figure()
	_sz = 150

	pred_class = TransitionPaths2024.Eugenio_Predict_Class_from_Inputs_20240715(TransitionPaths2024.onehot(response_data_sequences_new))
	likelihood_I  = TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:I,  TransitionPaths2024.onehot(response_data_sequences_new))
	likelihood_II = TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:II, TransitionPaths2024.onehot(response_data_sequences_new))
	likelihood_IV = TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:IV, TransitionPaths2024.onehot(response_data_sequences_new))
	response_I  = [ismissing(r) || r ≤ 0 ? NaN : r for r = df_20241213.C1]
	response_II = [ismissing(r) || r ≤ 0 ? NaN : r for r = df_20241213.C2]
	response_IV = [ismissing(r) || r ≤ 0 ? NaN : r for r = df_20241213.C4]
	
	ax_I  = Makie.Axis(fig[1,1][1,1], width=_sz, height=_sz, xlabel=L"\ln P_\mathrm{rbm}(\text{local I})",  ylabel="Class I. norm. resp.",  yscale=log10, xgridvisible=false, ygridvisible=false, xticks=-100:50:0)
	ax_II = Makie.Axis(fig[1,1][1,2], width=_sz, height=_sz, xlabel=L"\ln P_\mathrm{rbm}(\text{local II})", ylabel="Class II. norm. resp.", yscale=log10, xgridvisible=false, ygridvisible=false, xticks=-80:20:0)
	ax_IV = Makie.Axis(fig[1,1][1,3], width=_sz, height=_sz, xlabel=L"\ln P_\mathrm{rbm}(\text{local IV})", ylabel="Class IV. norm. resp.", yscale=log10, xgridvisible=false, ygridvisible=false, xticks=-100:25:0)

	for (class, color) = zip((:I, :II, :IV), _class_scatter_colors)
		Makie.scatter!(ax_I,  likelihood_I[pred_class .== class],  response_I[pred_class .== class]; color)
		Makie.scatter!(ax_II, likelihood_II[pred_class .== class], response_II[pred_class .== class]; color)
		Makie.scatter!(ax_IV, likelihood_IV[pred_class .== class], response_IV[pred_class .== class]; color)
	end

	labels = ["I", "II", "IV"]
	elements = [Makie.MarkerElement(; marker=:circle, color = _class_scatter_colors[j]) for j = 1:3]
	Makie.Legend(fig[1,1][1,4], elements, labels; framevisible=false)

	# Correlations barplot
	responses = Dict(:I => df_20241213.C1, :II => df_20241213.C2, :IV => df_20241213.C4)

	cat = [i for (i, target) = enumerate((:I, :II, :IV)) for (j, rbm) = enumerate((:global, :I, :II, :IV))]
	grp = [j for (i, target) = enumerate((:I, :II, :IV)) for (j, rbm) = enumerate((:global, :I, :II, :IV))]
				
	# correlations = [
	# 	nancor(TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(rbm, TransitionPaths2024.onehot(response_data_sequences_new)), replace(responses[target], missing => NaN))
	# 	for target = (:I, :II, :IV) for rbm = (:global, :I, :II, :IV)
	# ]

	response_data_sequences_new_uniq = unique(response_data_sequences_new)
	responses_uniq = Dict(target => last.(unique(first, zip(response_data_sequences_new, responses[target]))) for target = (:I, :II, :IV))

	__idx = Dict((target, rbm) => intersect(
		findall(isfinite, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(rbm, TransitionPaths2024.onehot(response_data_sequences_new_uniq))),
		findall(isfinite, [log(ifelse(ismissing(r) || isnan(r) || r ≤ 0, NaN, r)) for r = responses_uniq[target]])
	) for target = (:I, :II, :IV) for rbm = (:global, :I, :II, :IV))
	
	correlations_uniq = [
		corspearman(
			TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(rbm, TransitionPaths2024.onehot(response_data_sequences_new_uniq))[__idx[target, rbm]],
			[log(ifelse(ismissing(r) || isnan(r) || r ≤ 0, NaN, r)) for r = responses_uniq[target]][__idx[target, rbm]]
		)
		for target = (:I, :II, :IV) for rbm = (:global, :I, :II, :IV)
	]

	
	colors = [:gray, :cyan, :orange, :green]
		
	ax = Makie.Axis(fig[2,1][1,2], width=200, height=_sz, xticks=(1:3, ["I", "II", "IV"]), ylabel="Exp. resp. corr.", xgridvisible=false, ygridvisible=false)
	Makie.barplot!(ax, cat, correlations_uniq; dodge=grp, color=colors[grp], dodge_gap=0.05, width=0.9)
	Makie.hlines!(ax, 0.0; color=:black, linestyle=:dash)
	Makie.ylims!(-0.4, 0.9)

	labels = ["RBM(glob.)", "RBM(I)", "RBM(II)", "RBM(IV)"]
	elements = [Makie.PolyElement(; polycolor = colors[j]) for j = 1:4]
	Makie.Legend(fig[2,1][1,3], elements, labels; framevisible=false)

	# response_data_v2_sequences = TransitionPaths2024.Exp_20240703_sequences().sequences[
	# 	indexin(TransitionPaths2024.response_data_all_v2_20240708().name, TransitionPaths2024.Exp_20240703_sequences().names)
	# ]
	# rbm_log_likelihoods_v2 = TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:global, TransitionPaths2024.onehot(response_data_v2_sequences))
	
	ax = Makie.Axis(fig[2,1][1,1], width=250, height=_sz, xlabel=L"P_{\mathrm{RBM}}(\text{glob.})", ylabel="Frequency", xgridvisible=false, ygridvisible=false)
	bins = -90:10:10
	Makie.hist!(ax, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:global, TransitionPaths2024.onehot(TransitionPaths2024.Eugenio_MSA_20230419())); normalization=:pdf, bins=-90:2:10, color=(:green, 0.3), label="MSA")
	# Makie.stephist!(ax, rbm_log_likelihoods_v2[TransitionPaths2024.response_data_all_v2_20240708().type .== "Natural"], color=:green, label="Natural", normalization=:pdf, bins=-90:5:10, linewidth=4)
	# Makie.stephist!(ax, rbm_log_likelihoods_v2[TransitionPaths2024.response_data_all_v2_20240708().type .== "Designed"], color=:blue, label="Designed", normalization=:pdf, bins=-90:5:10, linewidth=4)
	# Makie.stephist!(ax, rbm_log_likelihoods_v2[TransitionPaths2024.response_data_all_v2_20240708().type .== "Shuffled"], color=:red, label="Shuffled", normalization=:pdf, bins=-90:5:10, linewidth=4)
	Makie.stephist!(ax, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:global, TransitionPaths2024.onehot(unique(df_20241213.aligned_sequences[df_20241213.group .== "Wild-types"]))); color=:green, label="Natural", normalization=:pdf, bins, linewidth=4)
	Makie.stephist!(ax, 
		TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:global, TransitionPaths2024.onehot(
			unique(df_20241213.aligned_sequences[df_20241213.group .∈ Ref(("1st", "2nd", "3rd", "4th") .* "_batch")])));
		color=:blue, label="Designed", normalization=:pdf, bins, linewidth=4)
	Makie.stephist!(ax, 
		TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:global, TransitionPaths2024.onehot(
			unique(df_20241213.aligned_sequences[df_20241213.group .== "Scrambled_paths"])));
		color=:red, label="Shuffled", normalization=:pdf, bins, linewidth=2, linestyle=(:dash, 1.2))
	Makie.xlims!(ax, -90, 0)
	Makie.ylims!(ax, -0.001, 0.1)
	Makie.axislegend(ax; position=:lt, framevisible=false)

	Makie.resize_to_layout!(fig)
	#Makie.save("/DATA-SSD/cossio/projects/2024/Transition_Paths/TransitionPaths2024.jl/figures/Fig5.pdf", fig)
	fig
end

# ╔═╡ 191d062a-6da4-46f2-8dfe-e70e35adcfd5
md"# Supplementary figures"

# ╔═╡ 6887c9cb-f9be-43a7-82a4-654fb6bc008f
df_20241213

# ╔═╡ 776e0bf2-2fbb-4377-9ceb-08bbab82f4a6
indexin(["Seq01", "Seq28", "Seq29", "Seq30", "Seq31", "Seq32", "Seq33", "Seq34", "Seq35", "Seq27"], df_20241213.seq_id)

# ╔═╡ 19fa19d9-74b5-44e5-8fdc-dd5f8d1bf8c9
md"## Supp. I -> IV batch 1"

# ╔═╡ 76b90e10-3325-49f8-95c0-9d48a2cbdbb7
let fig = Makie.Figure()
	sampled_path = unique(TransitionPaths2024.sampled_path_1to4batch1_20240703())
	probed_seq_ids = ["Seq01", "Seq28", "Seq29", "Seq30", "Seq31", "Seq32", "Seq33", "Seq34", "Seq35", "Seq27"]
	probed_seqs = df_20241213.aligned_sequences[indexin(probed_seq_ids, df_20241213.seq_id)]
	probed_seqs_idx::Vector{Int} = indexin(probed_seqs, sampled_path)

	C1_responses = df_20241213.C1[indexin(probed_seq_ids, df_20241213.seq_id)]
	C4_responses = df_20241213.C4[indexin(probed_seq_ids, df_20241213.seq_id)]
	C1_responses_err = df_20241213.eC1[indexin(probed_seq_ids, df_20241213.seq_id)]
	C4_responses_err = df_20241213.eC4[indexin(probed_seq_ids, df_20241213.seq_id)]
	rel_response = C1_responses ./ C4_responses
	@show rel_response
	
	ax = Makie.Axis(fig[1,1][1,1]; width=300, height=120, xgridvisible=false, ygridvisible=false, ylabel=L"\log P_\mathrm{RBM}(\text{global})")
	Makie.scatterlines!(ax, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:global, TransitionPaths2024.onehot(sampled_path)); color=:black, label="Designed")
	Makie.scatter!(ax, probed_seqs_idx, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:global, TransitionPaths2024.onehot(probed_seqs)); color=:red, markersize=10)
	Makie.hidexdecorations!(ax)

	ax_I = Makie.Axis(
		fig[1,1][2,1]; width=300, height=120, xgridvisible=false, ygridvisible=false,
		ylabel=L"\log P_\mathrm{RBM}(\text{I})", yticklabelcolor=_class_scatter_colors[1], ylabelcolor=_class_scatter_colors[1]
	)
	ax_IV = Makie.Axis(
		fig[1,1][2,1]; width=300, height=120, xgridvisible=false, ygridvisible=false, yaxisposition=:right,
		ylabel=L"\log P_\mathrm{RBM}(\text{IV})", yticklabelcolor=_class_scatter_colors[3], ylabelcolor=_class_scatter_colors[3],
	)
	Makie.lines!(ax_I,  TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:I,  TransitionPaths2024.onehot(sampled_path)); color=_class_scatter_colors[1])
	Makie.lines!(ax_IV, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:IV, TransitionPaths2024.onehot(sampled_path)); color=_class_scatter_colors[3])
	Makie.scatter!(ax_I,  TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:I,  TransitionPaths2024.onehot(sampled_path)); color=_class_scatter_colors[1], label="Designed")
	Makie.scatter!(ax_IV, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:IV, TransitionPaths2024.onehot(sampled_path)); color=_class_scatter_colors[3])
	Makie.scatter!(ax_I,  probed_seqs_idx, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:I,  TransitionPaths2024.onehot(probed_seqs)); color=:red, markersize=10)
	Makie.scatter!(ax_IV, probed_seqs_idx, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:IV, TransitionPaths2024.onehot(probed_seqs)); color=:red, markersize=10)
	Makie.hidexdecorations!(ax_I)
	Makie.hidexdecorations!(ax_IV)
	Makie.hidespines!(ax_IV)

	# ax = Makie.Axis(fig[1,1][5,1]; width=300, height=120, xgridvisible=false, ygridvisible=false, ylabel="I/IV rel. resp.", xticks=0:5:50)
	# Makie.hspan!(ax, 1, 1e6; color=(:cyan, 0.3))
	# Makie.hspan!(ax, 1e-6, 1; color=(:green, 0.3))
	# Makie.scatterlines!(ax, probed_seqs_idx, rel_response; color=:black)
	# #Makie.scatter!(ax, [i for (i,x) = zip(probed_seqs_idx, rel_response) if 0 < x < Inf], [x for x = rel_response if 0 < x < Inf]; color=:red, markersize=10)
	# Makie.ylims!(ax, -1000, 100)

	ax = Makie.Axis(fig[1,2][1,1], width=200, height=200, xlabel=L"I_{38}", ylabel=L"I_{36}", xgridvisible=false, ygridvisible=false)
	Makie.vlines!(ax, -1; color=:gray, linestyle=:dash, linewidth=1)
	Makie.hlines!(ax, -3; color=:gray, linestyle=:dash, linewidth=1)
	for (class, color) = zip((:I, :II, :IV), _class_scatter_colors)
		Makie.scatter!(
			ax,
			-inputs_h_from_v(TransitionPaths2024.Eugenio_RBM_20230419(:global), TransitionPaths2024.onehot(TransitionPaths2024.Eugenio_MSA_20230419()))[38, msa_local_rbm_class_prediction .== class],
			 inputs_h_from_v(TransitionPaths2024.Eugenio_RBM_20230419(:global), TransitionPaths2024.onehot(TransitionPaths2024.Eugenio_MSA_20230419()))[36, msa_local_rbm_class_prediction .== class];
			color=(color, 0.5), markersize=4
		)
	end

	Makie.scatterlines!(ax, 
		-inputs_h_from_v(TransitionPaths2024.Eugenio_RBM_20230419(:global), TransitionPaths2024.onehot(sampled_path))[38, :],
		 inputs_h_from_v(TransitionPaths2024.Eugenio_RBM_20230419(:global), TransitionPaths2024.onehot(sampled_path))[36, :];
		color=:black
	)

	Makie.scatter!(ax, 
		-inputs_h_from_v(TransitionPaths2024.Eugenio_RBM_20230419(:global), TransitionPaths2024.onehot(probed_seqs))[38, :],
		 inputs_h_from_v(TransitionPaths2024.Eugenio_RBM_20230419(:global), TransitionPaths2024.onehot(probed_seqs))[36, :];
		color=:red, markersize=10
	)
	
	Makie.arrows!(ax,
		-inputs_h_from_v(TransitionPaths2024.Eugenio_RBM_20230419(:global), TransitionPaths2024.onehot(sampled_path))[38, 1:1] .+ 0.35,
		 inputs_h_from_v(TransitionPaths2024.Eugenio_RBM_20230419(:global), TransitionPaths2024.onehot(sampled_path))[36, 1:1] .- 0.35,
		[-0.15], [0.15]; linewidth=2, color=:gray, arrowsize=10
	)
	
	Makie.xlims!(ax, -4.5, 2)
	Makie.ylims!(ax, -7, 3)

	ax = Makie.Axis(fig[1,2][2,1]; width=200, height=200, xgridvisible=false, ygridvisible=false, xlabel="Norm. response I", ylabel="Norm. response IV")
	Makie.band!(ax, [-10, thresh_significance], [-10, -10], [thresh_significance, thresh_significance]; color=(:red, 0.3)) # neither specificity is significant region
	Makie.band!(ax, [-10, thresh_significance], [thresh_significance, thresh_significance], [100, 100]; color=(:green, 0.3)) # specific to IV
	Makie.band!(ax, [thresh_significance, 100], [-100, -100], [thresh_significance, thresh_significance]; color=(:cyan, 0.3)) # specific to I
	Makie.scatterlines!(ax, C1_responses, C4_responses; color=:black, linewidth=1)
	Makie.errorbars!(ax, C1_responses, C4_responses, C1_responses_err; color=:gray, direction=:x, linewidth=2, whiskerwidth=5)
	Makie.errorbars!(ax, C1_responses, C4_responses, C4_responses_err; color=:gray, direction=:y, linewidth=2, whiskerwidth=5)
	Makie.scatter!(ax, C1_responses, C4_responses; color=:red)
	Makie.arrows!(ax, C1_responses[1:1] .+ 0.2, C4_responses[1:1] .- 0.15, [-0.07], [0.07]; linewidth=2, color=:gray, arrowsize=10)
	Makie.xlims!(ax, -0.5, 12)
	Makie.ylims!(ax, -0.25, 2)

	Makie.resize_to_layout!(fig)
	fig
end

# ╔═╡ 6184d501-0f4e-4d22-8eb0-88a07b0efab8
md"## Supp. I -> IV batch 2"

# ╔═╡ d0cc1b76-9705-47ff-8c08-3402531002f2
let fig = Makie.Figure()
	sampled_path = unique(TransitionPaths2024.sampled_path_1to4batch2_20240703())
	probed_seq_ids = ["Seq11", "Seq43", "Seq44", "Seq45", "Seq46", "Seq47", "Seq34", "Seq35", "Seq48", "Seq27"]
	probed_seqs = df_20241213.aligned_sequences[indexin(probed_seq_ids, df_20241213.seq_id)]
	probed_seqs_idx::Vector{Int} = indexin(probed_seqs, sampled_path)

	C1_responses = df_20241213.C1[indexin(probed_seq_ids, df_20241213.seq_id)]
	C4_responses = df_20241213.C4[indexin(probed_seq_ids, df_20241213.seq_id)]
	C1_responses_err = df_20241213.eC1[indexin(probed_seq_ids, df_20241213.seq_id)]
	C4_responses_err = df_20241213.eC4[indexin(probed_seq_ids, df_20241213.seq_id)]
	rel_response = C1_responses ./ C4_responses
	@show rel_response
	
	ax = Makie.Axis(fig[1,1][1,1]; width=300, height=120, xgridvisible=false, ygridvisible=false, ylabel=L"\log P_\mathrm{RBM}(\text{global})")
	Makie.scatterlines!(ax, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:global, TransitionPaths2024.onehot(sampled_path)); color=:black, label="Designed")
	Makie.scatter!(ax, probed_seqs_idx, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:global, TransitionPaths2024.onehot(probed_seqs)); color=:red, markersize=10)
	Makie.hidexdecorations!(ax)

	ax_I = Makie.Axis(
		fig[1,1][2,1]; width=300, height=120, xgridvisible=false, ygridvisible=false,
		ylabel=L"\log P_\mathrm{RBM}(\text{I})", yticklabelcolor=_class_scatter_colors[1], ylabelcolor=_class_scatter_colors[1]
	)
	ax_IV = Makie.Axis(
		fig[1,1][2,1]; width=300, height=120, xgridvisible=false, ygridvisible=false, yaxisposition=:right,
		ylabel=L"\log P_\mathrm{RBM}(\text{IV})", yticklabelcolor=_class_scatter_colors[3], ylabelcolor=_class_scatter_colors[3],
	)
	Makie.lines!(ax_I,  TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:I,  TransitionPaths2024.onehot(sampled_path)); color=_class_scatter_colors[1])
	Makie.lines!(ax_IV, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:IV, TransitionPaths2024.onehot(sampled_path)); color=_class_scatter_colors[3])
	Makie.scatter!(ax_I,  TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:I,  TransitionPaths2024.onehot(sampled_path)); color=_class_scatter_colors[1], label="Designed")
	Makie.scatter!(ax_IV, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:IV, TransitionPaths2024.onehot(sampled_path)); color=_class_scatter_colors[3])
	Makie.scatter!(ax_I,  probed_seqs_idx, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:I,  TransitionPaths2024.onehot(probed_seqs)); color=:red, markersize=10)
	Makie.scatter!(ax_IV, probed_seqs_idx, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:IV, TransitionPaths2024.onehot(probed_seqs)); color=:red, markersize=10)
	Makie.hidexdecorations!(ax_I)
	Makie.hidexdecorations!(ax_IV)
	Makie.hidespines!(ax_IV)

	# ax = Makie.Axis(fig[1,1][5,1]; width=300, height=120, xgridvisible=false, ygridvisible=false, ylabel="I/IV rel. resp.", yscale=log10, xticks=0:5:50)
	# Makie.hspan!(ax, 1, 1e6; color=(:cyan, 0.3))
	# Makie.hspan!(ax, 1e-6, 1; color=(:green, 0.3))
	# Makie.scatterlines!(ax, [i for (i,x) = zip(probed_seqs_idx, rel_response) if 0 < x < Inf], [x for x = rel_response if 0 < x < Inf]; color=:black)
	# Makie.scatter!(ax, [i for (i,x) = zip(probed_seqs_idx, rel_response) if 0 < x < Inf], [x for x = rel_response if 0 < x < Inf]; color=:red, markersize=10)
	# Makie.ylims!(ax, 1e-3, 5e3)

	ax = Makie.Axis(fig[1,2][1,1], width=200, height=200, xlabel=L"I_{1}", ylabel=L"I_{2}", xgridvisible=false, ygridvisible=false) # h_38 -> h_1, h_36 -> h_2
	Makie.vlines!(ax, -1; color=:gray, linestyle=:dash, linewidth=1)
	Makie.hlines!(ax, -3; color=:gray, linestyle=:dash, linewidth=1)
	for (class, color) = zip((:I, :II, :IV), _class_scatter_colors)
		Makie.scatter!(
			ax,
			-inputs_h_from_v(TransitionPaths2024.Eugenio_RBM_20230419(:global), TransitionPaths2024.onehot(TransitionPaths2024.Eugenio_MSA_20230419()))[38, msa_local_rbm_class_prediction .== class],
			 inputs_h_from_v(TransitionPaths2024.Eugenio_RBM_20230419(:global), TransitionPaths2024.onehot(TransitionPaths2024.Eugenio_MSA_20230419()))[36, msa_local_rbm_class_prediction .== class];
			color=(color, 0.5), markersize=4
		)
	end

	Makie.scatterlines!(ax, 
		-inputs_h_from_v(TransitionPaths2024.Eugenio_RBM_20230419(:global), TransitionPaths2024.onehot(sampled_path))[38, :],
		 inputs_h_from_v(TransitionPaths2024.Eugenio_RBM_20230419(:global), TransitionPaths2024.onehot(sampled_path))[36, :];
		color=:black
	)

	Makie.scatter!(ax, 
		-inputs_h_from_v(TransitionPaths2024.Eugenio_RBM_20230419(:global), TransitionPaths2024.onehot(probed_seqs))[38, :],
		 inputs_h_from_v(TransitionPaths2024.Eugenio_RBM_20230419(:global), TransitionPaths2024.onehot(probed_seqs))[36, :];
		color=:red, markersize=10
	)
	
	Makie.arrows!(ax,
		-inputs_h_from_v(TransitionPaths2024.Eugenio_RBM_20230419(:global), TransitionPaths2024.onehot(sampled_path))[38, 1:1] .+ 0.35,
		 inputs_h_from_v(TransitionPaths2024.Eugenio_RBM_20230419(:global), TransitionPaths2024.onehot(sampled_path))[36, 1:1] .- 0.35,
		[-0.15], [0.15]; linewidth=2, color=:gray, arrowsize=10
	)
	
	Makie.xlims!(ax, -4.5, 2)
	Makie.ylims!(ax, -7, 3)

	ax = Makie.Axis(fig[1,2][2,1]; width=200, height=200, xgridvisible=false, ygridvisible=false, xlabel="Norm. response I", ylabel="Norm. response IV")
	Makie.band!(ax, [-10, thresh_significance], [-10, -10], [thresh_significance, thresh_significance]; color=(:red, 0.3)) # neither specificity is significant region
	Makie.band!(ax, [-10, thresh_significance], [thresh_significance, thresh_significance], [100, 100]; color=(:green, 0.3)) # specific to IV
	Makie.band!(ax, [thresh_significance, 100], [-100, -100], [thresh_significance, thresh_significance]; color=(:cyan, 0.3)) # specific to I
	Makie.scatterlines!(ax, C1_responses, C4_responses; color=:black, linewidth=1)
	Makie.errorbars!(ax, C1_responses, C4_responses, C1_responses_err; color=:gray, direction=:x, linewidth=2, whiskerwidth=5)
	Makie.errorbars!(ax, C1_responses, C4_responses, C4_responses_err; color=:gray, direction=:y, linewidth=2, whiskerwidth=5)
	Makie.scatter!(ax, C1_responses, C4_responses; color=:red)
	Makie.arrows!(ax, C1_responses[1:1] .+ 0.2, C4_responses[1:1] .- 0.15, [-0.07], [0.07]; linewidth=2, color=:gray, arrowsize=10)
	Makie.xlims!(ax, -0.5, 12)
	Makie.ylims!(ax, -0.25, 2)

	Makie.resize_to_layout!(fig)
	fig
end

# ╔═╡ 4354e416-143a-45ca-916e-6304accf1bdf
md"## Supp. I -> IV direct (Batch 3)"

# ╔═╡ 9734267d-e44f-4f1d-9c7b-3464ed5a2677
let fig = Makie.Figure()
	sampled_path = unique(TransitionPaths2024.sampled_path_1to4direct_20240703())
	probed_seq_ids = ["Seq02", "Seq55", "Seq56", "Seq57", "Seq58", "Seq53", "Seq35", "Seq48", "Seq27"]
	probed_seqs = df_20241213.aligned_sequences[indexin(probed_seq_ids, df_20241213.seq_id)]
	probed_seqs_idx::Vector{Int} = indexin(probed_seqs, sampled_path)

	C1_responses = df_20241213.C1[indexin(probed_seq_ids, df_20241213.seq_id)]
	C4_responses = df_20241213.C4[indexin(probed_seq_ids, df_20241213.seq_id)]
	C1_responses_err = df_20241213.eC1[indexin(probed_seq_ids, df_20241213.seq_id)]
	C4_responses_err = df_20241213.eC4[indexin(probed_seq_ids, df_20241213.seq_id)]
	rel_response = C1_responses ./ C4_responses
	
	ax = Makie.Axis(fig[1,1][1,1]; width=300, height=120, xgridvisible=false, ygridvisible=false, ylabel=L"\log P_\mathrm{RBM}(\text{global})")
	Makie.scatterlines!(ax, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:global, TransitionPaths2024.onehot(sampled_path)); color=:black, label="Designed")
	Makie.scatter!(ax, probed_seqs_idx, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:global, TransitionPaths2024.onehot(probed_seqs)); color=:red, markersize=10)
	Makie.hidexdecorations!(ax)

	ax_I = Makie.Axis(
		fig[1,1][2,1]; width=300, height=120, xgridvisible=false, ygridvisible=false,
		ylabel=L"\log P_\mathrm{RBM}(\text{I})", yticklabelcolor=_class_scatter_colors[1], ylabelcolor=_class_scatter_colors[1]
	)
	ax_IV = Makie.Axis(
		fig[1,1][2,1]; width=300, height=120, xgridvisible=false, ygridvisible=false, yaxisposition=:right,
		ylabel=L"\log P_\mathrm{RBM}(\text{IV})", yticklabelcolor=_class_scatter_colors[3], ylabelcolor=_class_scatter_colors[3],
	)
	Makie.lines!(ax_I,  TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:I,  TransitionPaths2024.onehot(sampled_path)); color=_class_scatter_colors[1])
	Makie.lines!(ax_IV, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:IV, TransitionPaths2024.onehot(sampled_path)); color=_class_scatter_colors[3])
	Makie.scatter!(ax_I,  TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:I,  TransitionPaths2024.onehot(sampled_path)); color=_class_scatter_colors[1], label="Designed")
	Makie.scatter!(ax_IV, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:IV, TransitionPaths2024.onehot(sampled_path)); color=_class_scatter_colors[3])
	Makie.scatter!(ax_I,  probed_seqs_idx, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:I,  TransitionPaths2024.onehot(probed_seqs)); color=:red, markersize=10)
	Makie.scatter!(ax_IV, probed_seqs_idx, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:IV, TransitionPaths2024.onehot(probed_seqs)); color=:red, markersize=10)
	Makie.hidexdecorations!(ax_I)
	Makie.hidexdecorations!(ax_IV)
	Makie.hidespines!(ax_IV)

	# ax = Makie.Axis(fig[1,1][5,1]; width=300, height=120, xgridvisible=false, ygridvisible=false, ylabel="I/IV rel. resp.", yscale=log10, xticks=0:5:50)
	# Makie.hspan!(ax, 1, 1e6; color=(:cyan, 0.3))
	# Makie.hspan!(ax, 1e-6, 1; color=(:green, 0.3))
	# Makie.scatterlines!(ax, [i for (i,x) = zip(probed_seqs_idx, rel_response) if 0 < x < Inf], [x for x = rel_response if 0 < x < Inf]; color=:black)
	# Makie.scatter!(ax, [i for (i,x) = zip(probed_seqs_idx, rel_response) if 0 < x < Inf], [x for x = rel_response if 0 < x < Inf]; color=:red, markersize=10)
	# Makie.ylims!(ax, 1e-3, 5e3)

	ax = Makie.Axis(fig[1,2][1,1], width=200, height=200, xlabel=L"I_{1}", ylabel=L"I_{2}", xgridvisible=false, ygridvisible=false) # h_38 -> h_1, h_36 -> h_2
	Makie.vlines!(ax, -1; color=:gray, linestyle=:dash, linewidth=1)
	Makie.hlines!(ax, -3; color=:gray, linestyle=:dash, linewidth=1)
	for (class, color) = zip((:I, :II, :IV), _class_scatter_colors)
		Makie.scatter!(
			ax,
			-inputs_h_from_v(TransitionPaths2024.Eugenio_RBM_20230419(:global), TransitionPaths2024.onehot(TransitionPaths2024.Eugenio_MSA_20230419()))[38, msa_local_rbm_class_prediction .== class],
			 inputs_h_from_v(TransitionPaths2024.Eugenio_RBM_20230419(:global), TransitionPaths2024.onehot(TransitionPaths2024.Eugenio_MSA_20230419()))[36, msa_local_rbm_class_prediction .== class];
			color=(color, 0.5), markersize=4
		)
	end

	Makie.scatterlines!(ax, 
		-inputs_h_from_v(TransitionPaths2024.Eugenio_RBM_20230419(:global), TransitionPaths2024.onehot(sampled_path))[38, :],
		 inputs_h_from_v(TransitionPaths2024.Eugenio_RBM_20230419(:global), TransitionPaths2024.onehot(sampled_path))[36, :];
		color=:black
	)

	Makie.scatter!(ax, 
		-inputs_h_from_v(TransitionPaths2024.Eugenio_RBM_20230419(:global), TransitionPaths2024.onehot(probed_seqs))[38, :],
		 inputs_h_from_v(TransitionPaths2024.Eugenio_RBM_20230419(:global), TransitionPaths2024.onehot(probed_seqs))[36, :];
		color=:red, markersize=10
	)
	
	Makie.arrows!(ax,
		-inputs_h_from_v(TransitionPaths2024.Eugenio_RBM_20230419(:global), TransitionPaths2024.onehot(sampled_path))[38, 1:1] .+ 0.35,
		 inputs_h_from_v(TransitionPaths2024.Eugenio_RBM_20230419(:global), TransitionPaths2024.onehot(sampled_path))[36, 1:1] .- 0.35,
		[-0.15], [0.15]; linewidth=2, color=:gray, arrowsize=10
	)
	
	Makie.xlims!(ax, -4.5, 2)
	Makie.ylims!(ax, -7, 3)

	ax = Makie.Axis(fig[1,2][2,1]; width=200, height=200, xgridvisible=false, ygridvisible=false, xlabel="Norm. response I", ylabel="Norm. response IV")
	Makie.band!(ax, [-10, thresh_significance], [-10, -10], [thresh_significance, thresh_significance]; color=(:red, 0.3)) # neither specificity is significant region
	Makie.band!(ax, [-10, thresh_significance], [thresh_significance, thresh_significance], [100, 100]; color=(:green, 0.3)) # specific to IV
	Makie.band!(ax, [thresh_significance, 100], [-100, -100], [thresh_significance, thresh_significance]; color=(:cyan, 0.3)) # specific to I
	Makie.scatterlines!(ax, C1_responses, C4_responses; color=:black, linewidth=1)
	Makie.errorbars!(ax, C1_responses, C4_responses, C1_responses_err; color=:gray, direction=:x, linewidth=2, whiskerwidth=5)
	Makie.errorbars!(ax, C1_responses, C4_responses, C4_responses_err; color=:gray, direction=:y, linewidth=2, whiskerwidth=5)
	Makie.scatter!(ax, C1_responses, C4_responses; color=:red)
	Makie.arrows!(ax, C1_responses[1:1] .+ 0.2, C4_responses[1:1] .- 0.15, [-0.07], [0.07]; linewidth=2, color=:gray, arrowsize=10)
	Makie.xlims!(ax, -0.5, 12)
	Makie.ylims!(ax, -0.25, 2)

	Makie.resize_to_layout!(fig)
	fig
end

# ╔═╡ bf689107-5744-4535-b455-f0678071b29e
md"## Supp. I -> II (Batch 4)"

# ╔═╡ 3e59d567-6300-4e01-a3c9-7226acd6174d
let fig = Makie.Figure()
	sampled_path = unique(TransitionPaths2024.sampled_path_1to2rep2_20240703())
	probed_seq_ids = ["Seq02", "Seq59", "Seq65", "Seq66", "Seq67", "Seq68", "Seq69", "Seq64"]
	probed_seqs = df_20241213.aligned_sequences[indexin(probed_seq_ids, df_20241213.seq_id)]
	probed_seqs_idx::Vector{Int} = indexin(probed_seqs, sampled_path)

	C1_responses = df_20241213.C1[indexin(probed_seq_ids, df_20241213.seq_id)]
	C2_responses = df_20241213.C2[indexin(probed_seq_ids, df_20241213.seq_id)]
	C1_responses_err = df_20241213.eC1[indexin(probed_seq_ids, df_20241213.seq_id)]
	C2_responses_err = df_20241213.eC2[indexin(probed_seq_ids, df_20241213.seq_id)]
	rel_response = C1_responses ./ C2_responses
	
	ax = Makie.Axis(fig[1,1][1,1]; width=300, height=120, xgridvisible=false, ygridvisible=false, ylabel=L"\log P_\mathrm{RBM}(\text{global})")
	Makie.scatterlines!(ax, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:global, TransitionPaths2024.onehot(sampled_path)); color=:black, label="Designed")
	Makie.scatter!(ax, probed_seqs_idx, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:global, TransitionPaths2024.onehot(probed_seqs)); color=:red, markersize=10)
	Makie.hidexdecorations!(ax)

	ax_I = Makie.Axis(
		fig[1,1][2,1]; width=300, height=120, xgridvisible=false, ygridvisible=false,
		ylabel=L"\log P_\mathrm{RBM}(\text{I})", yticklabelcolor=_class_scatter_colors[1], ylabelcolor=_class_scatter_colors[1]
	)
	ax_II = Makie.Axis(
		fig[1,1][2,1]; width=300, height=120, xgridvisible=false, ygridvisible=false, yaxisposition=:right,
		ylabel=L"\log P_\mathrm{RBM}(\text{II})", yticklabelcolor=_class_scatter_colors[2], ylabelcolor=_class_scatter_colors[2],
	)
	Makie.lines!(ax_I,  TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:I,  TransitionPaths2024.onehot(sampled_path)); color=_class_scatter_colors[1])
	Makie.lines!(ax_II, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:II, TransitionPaths2024.onehot(sampled_path)); color=_class_scatter_colors[2])
	Makie.scatter!(ax_I,  TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:I,  TransitionPaths2024.onehot(sampled_path)); color=_class_scatter_colors[1], label="Designed")
	Makie.scatter!(ax_II, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:II, TransitionPaths2024.onehot(sampled_path)); color=_class_scatter_colors[2])
	Makie.scatter!(ax_I,  probed_seqs_idx, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:I,  TransitionPaths2024.onehot(probed_seqs)); color=:red, markersize=10)
	Makie.scatter!(ax_II, probed_seqs_idx, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:II, TransitionPaths2024.onehot(probed_seqs)); color=:red, markersize=10)
	Makie.hidexdecorations!(ax_I)
	Makie.hidexdecorations!(ax_II)
	Makie.hidespines!(ax_II)

	# ax = Makie.Axis(fig[1,1][5,1]; width=300, height=120, xgridvisible=false, ygridvisible=false, ylabel="I/II rel. resp.", yscale=log10, xticks=0:5:50)
	# Makie.hspan!(ax, 1, 1e6; color=(:cyan, 0.3))
	# Makie.hspan!(ax, 1e-6, 1; color=(:green, 0.3))
	# Makie.scatterlines!(ax, [i for (i,x) = zip(probed_seqs_idx, rel_response) if 0 < x < Inf], [x for x = rel_response if 0 < x < Inf]; color=:black)
	# Makie.scatter!(ax, [i for (i,x) = zip(probed_seqs_idx, rel_response) if 0 < x < Inf], [x for x = rel_response if 0 < x < Inf]; color=:red, markersize=10)
	# Makie.ylims!(ax, 1e-3, 5e3)

	ax = Makie.Axis(fig[1,2][1,1], width=200, height=200, xlabel=L"I_{1}", ylabel=L"I_{2}", xgridvisible=false, ygridvisible=false) # h_38 -> h_1, h_36 -> h_2
	Makie.vlines!(ax, -1; color=:gray, linestyle=:dash, linewidth=1)
	Makie.hlines!(ax, -3; color=:gray, linestyle=:dash, linewidth=1)
	for (class, color) = zip((:I, :II, :IV), _class_scatter_colors)
		Makie.scatter!(
			ax,
			-inputs_h_from_v(TransitionPaths2024.Eugenio_RBM_20230419(:global), TransitionPaths2024.onehot(TransitionPaths2024.Eugenio_MSA_20230419()))[38, msa_local_rbm_class_prediction .== class],
			 inputs_h_from_v(TransitionPaths2024.Eugenio_RBM_20230419(:global), TransitionPaths2024.onehot(TransitionPaths2024.Eugenio_MSA_20230419()))[36, msa_local_rbm_class_prediction .== class];
			color=(color, 0.5), markersize=4
		)
	end

	Makie.scatterlines!(ax, 
		-inputs_h_from_v(TransitionPaths2024.Eugenio_RBM_20230419(:global), TransitionPaths2024.onehot(sampled_path))[38, :],
		 inputs_h_from_v(TransitionPaths2024.Eugenio_RBM_20230419(:global), TransitionPaths2024.onehot(sampled_path))[36, :];
		color=:black
	)

	Makie.scatter!(ax, 
		-inputs_h_from_v(TransitionPaths2024.Eugenio_RBM_20230419(:global), TransitionPaths2024.onehot(probed_seqs))[38, :],
		 inputs_h_from_v(TransitionPaths2024.Eugenio_RBM_20230419(:global), TransitionPaths2024.onehot(probed_seqs))[36, :];
		color=:red, markersize=10
	)
	
	Makie.arrows2d!(ax,
		-inputs_h_from_v(TransitionPaths2024.Eugenio_RBM_20230419(:global), TransitionPaths2024.onehot(sampled_path))[38, 1:1] .+ 0.4,
		 inputs_h_from_v(TransitionPaths2024.Eugenio_RBM_20230419(:global), TransitionPaths2024.onehot(sampled_path))[36, 1:1] .- 0.4,
		[-0.3], [0.3]; color=:gray, tipwidth=20, shaftwidth=6
	)
	
	Makie.xlims!(ax, -4.5, 2)
	Makie.ylims!(ax, -7, 3)

	ax = Makie.Axis(fig[1,2][2,1]; width=200, height=200, xgridvisible=false, ygridvisible=false, xlabel="Norm. response I", ylabel="Norm. response II")
	Makie.band!(ax, [-10, thresh_significance], [-10, -10], [thresh_significance, thresh_significance]; color=(:red, 0.3)) # neither specificity is significant region
	Makie.band!(ax, [-10, thresh_significance], [thresh_significance, thresh_significance], [100, 100]; color=(:green, 0.3)) # specific to IV
	Makie.band!(ax, [thresh_significance, 100], [-100, -100], [thresh_significance, thresh_significance]; color=(:cyan, 0.3)) # specific to I
	Makie.scatterlines!(ax, C1_responses, C2_responses; color=:black, linewidth=1)
	Makie.errorbars!(ax, C1_responses, C2_responses, C1_responses_err; color=:gray, direction=:x, linewidth=2, whiskerwidth=5)
	Makie.errorbars!(ax, C1_responses, C2_responses, C2_responses_err; color=:gray, direction=:y, linewidth=2, whiskerwidth=5)
	Makie.scatter!(ax, C1_responses, C2_responses; color=:red)
	Makie.arrows!(ax, C1_responses[1:1] .+ 0.2, C2_responses[1:1] .- 0.15, [-0.07], [0.07]; linewidth=2, color=:gray, arrowsize=10)
	Makie.xlims!(ax, -0.5, 12)
	Makie.ylims!(ax, -0.25, 2)

	Makie.resize_to_layout!(fig)
	fig
end

# ╔═╡ 4a19b62a-8644-4ad1-a516-8b04bd84464a
md"## Suppl. path I->I"

# ╔═╡ ef668901-d457-42fe-bcb9-0dd53bf04e28
let fig = Makie.Figure()

	sampled_path = unique(TransitionPaths2024.sampled_path_1to1_20240703())
	path_response_data = TransitionPaths2024.response_data_path_1to1_20240703()

	C1_responses = df_20241213.C1[indexin(path_response_data.ww_names, df_20241213.ww_id)]
	C2_responses = df_20241213.C2[indexin(path_response_data.ww_names, df_20241213.ww_id)]
	C1_responses_err = df_20241213.eC1[indexin(path_response_data.ww_names, df_20241213.ww_id)]
	C2_responses_err = df_20241213.eC2[indexin(path_response_data.ww_names, df_20241213.ww_id)]
	
	probed_seqs = TransitionPaths2024.Exp_20240703_sequences().sequences[indexin(path_response_data.ww_names, TransitionPaths2024.Exp_20240703_sequences().names)]
	probed_seqs_idx::Vector{Int} = indexin(probed_seqs, sampled_path)

	ax = Makie.Axis(fig[1,1][1,1]; width=175, height=75, xgridvisible=false, ygridvisible=false, ylabel=L"\ln P_\mathrm{RBM}(\text{glob.})")
	Makie.scatterlines!(ax, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:global, TransitionPaths2024.onehot(sampled_path)); color=:black)
	Makie.scatter!(ax, probed_seqs_idx, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:global, TransitionPaths2024.onehot(probed_seqs)); color=:red, markersize=10)
	Makie.hidexdecorations!(ax)

	ax_I = Makie.Axis(
		fig[1,1][2,1]; width=175, height=75, xgridvisible=false, ygridvisible=false,
		ylabel=L"\log P_\mathrm{RBM}(\text{I})", yticklabelcolor=_class_scatter_colors[1], ylabelcolor=_class_scatter_colors[1]
	)
	ax_II = Makie.Axis(
		fig[1,1][2,1]; width=175, height=75, xgridvisible=false, ygridvisible=false, yaxisposition=:right,
		ylabel=L"\log P_\mathrm{RBM}(\text{II})", yticklabelcolor=_class_scatter_colors[2], ylabelcolor=_class_scatter_colors[2],
	)
	Makie.lines!(ax_II, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:II, TransitionPaths2024.onehot(sampled_path)); color=_class_scatter_colors[2], linewidth=0.75)
	Makie.scatterlines!(ax_I, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:I, TransitionPaths2024.onehot(sampled_path)); color=_class_scatter_colors[1])
	Makie.scatter!(ax_I, probed_seqs_idx, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:I, TransitionPaths2024.onehot(probed_seqs)); color=:red, markersize=10)
	Makie.hidexdecorations!(ax)

	ax = Makie.Axis(fig[1,1][3,1]; width=175, height=75, xgridvisible=false, ygridvisible=false, xlabel="step", ylabel="Norm. resp. I", xticks=1:3:20)
	# Makie.scatterlines!(ax, probed_seqs_idx, path_response_data.C1_response; color=:black)
	# Makie.errorbars!(ax, probed_seqs_idx, path_response_data.C1_response, path_response_data.C1_response_err; color=:gray, direction=:y, linewidth=2, whiskerwidth=5)
	# Makie.scatter!(ax, probed_seqs_idx, path_response_data.C1_response; color=:red)

	Makie.scatterlines!(ax, probed_seqs_idx, C1_responses; color=:black)
	Makie.errorbars!(ax, probed_seqs_idx, C1_responses, C1_responses_err; color=:gray, direction=:y, linewidth=2, whiskerwidth=5)
	Makie.scatter!(ax, probed_seqs_idx, C1_responses; color=:red)

	Makie.hspan!(ax, -5, thresh_significance; color=(:red, 0.5))
	Makie.ylims!(ax, -2, 35)

	ax = Makie.Axis(fig[1,2], width=200, height=200, xlabel=L"I_{1}", ylabel=L"I_{2}", xgridvisible=false, ygridvisible=false) # h_38 -> h_1, h_36 -> h_2
	Makie.vlines!(ax, -1; color=:gray, linestyle=:dash, linewidth=1)
	Makie.hlines!(ax, -3; color=:gray, linestyle=:dash, linewidth=1)
	for (class, color) = zip((:I, :II, :IV), _class_scatter_colors)
		Makie.scatter!(
			ax,
			-inputs_h_from_v(TransitionPaths2024.Eugenio_RBM_20230419(:global), TransitionPaths2024.onehot(TransitionPaths2024.Eugenio_MSA_20230419()))[38, msa_local_rbm_class_prediction .== class],
			 inputs_h_from_v(TransitionPaths2024.Eugenio_RBM_20230419(:global), TransitionPaths2024.onehot(TransitionPaths2024.Eugenio_MSA_20230419()))[36, msa_local_rbm_class_prediction .== class];
			color=(color, 0.5), markersize=4
		)
	end

	Makie.scatterlines!(ax, 
		-inputs_h_from_v(TransitionPaths2024.Eugenio_RBM_20230419(:global), TransitionPaths2024.onehot(sampled_path))[38, :],
		 inputs_h_from_v(TransitionPaths2024.Eugenio_RBM_20230419(:global), TransitionPaths2024.onehot(sampled_path))[36, :];
		color=:black
	)

	Makie.scatter!(ax, 
		-inputs_h_from_v(TransitionPaths2024.Eugenio_RBM_20230419(:global), TransitionPaths2024.onehot(probed_seqs))[38, :],
		 inputs_h_from_v(TransitionPaths2024.Eugenio_RBM_20230419(:global), TransitionPaths2024.onehot(probed_seqs))[36, :];
		color=:red, markersize=10
	)
	
	Makie.arrows2d!(ax,
		-inputs_h_from_v(TransitionPaths2024.Eugenio_RBM_20230419(:global), TransitionPaths2024.onehot(sampled_path))[38, 1:1] .+ 0.35,
		 inputs_h_from_v(TransitionPaths2024.Eugenio_RBM_20230419(:global), TransitionPaths2024.onehot(sampled_path))[36, 1:1] .- 0.35,
		[-0.2], [0.2]; color=:gray, tipwidth=20, shaftwidth=6
	)
	
	Makie.xlims!(ax, -4.5, 2)
	Makie.ylims!(ax, -7, 3)

	Makie.resize_to_layout!(fig)
	#Makie.save("/DATA-SSD/cossio/projects/2024/Transition_Paths/TransitionPaths2024.jl/figures/Fig2.pdf", fig)
	fig
end

# ╔═╡ d2772e25-a731-433a-b10c-a92901172ea1
md"## Paths alignment figures"

# ╔═╡ 38d65a29-446d-4e04-b876-b039baa00369
TransitionPaths2024.aligned_path_plot(map(string, TransitionPaths2024.sampled_path_1to1_20240703()))

# ╔═╡ 340db97a-f8e3-4523-8101-d9cf6bfb549d
TransitionPaths2024.aligned_path_plot(map(string, TransitionPaths2024.sampled_path_1to2rep1_20240703()))

# ╔═╡ 735b7450-6722-46e8-a1b9-3b963d4324e2
TransitionPaths2024.aligned_path_plot(map(string, TransitionPaths2024.sampled_path_1to4_20240703()))

# ╔═╡ 3a544f67-efe6-47c3-b522-d7f3a2049b58
md"## Distance histograms"

# ╔═╡ 33cae588-2a3d-4755-a856-517ca4593384
function hamming(seq1::BitMatrix, seq2::BitMatrix)
    @assert size(seq1) == size(seq2)
    return sum(any(seq1 .!= seq2, dims=1))
end

# ╔═╡ 2b5cdff8-7b1d-4f4a-8b13-bc7bb5a73071
function hamming(s1::BitArray{3}, s2::BitArray{3})
    @assert size(s1, 1) == size(s2, 1)
    @assert size(s1, 2) == size(s2, 2)
    return [hamming(s1[:,:,n], s2[:,:,m]) for n in axes(s1,3), m in axes(s2,3)]
end

# ╔═╡ 7742b870-4495-498c-b430-f83635a03536
function sample_rbm(rbm)
	sampled_f = zeros(10_000)
	sampled_v = sample_from_inputs(rbm.visible, zeros(size(rbm.visible)..., 1000))
	sampled_f[1] = mean(free_energy(rbm, sampled_v))
	@progress for t in 2:length(sampled_f)
	    sampled_v .= sample_v_from_v(rbm, sampled_v)
	    sampled_f[t] = mean(free_energy(rbm, sampled_v))
	end
	return sampled_f, sampled_v
end

# ╔═╡ 54cab97b-82b1-4c07-b480-f2477c4a8ce9
equilibrated_samples_from_rbm_global_free_energy, equilibrated_samples_from_rbm_global = sample_rbm(TransitionPaths2024.Eugenio_RBM_20230419(:global));

# ╔═╡ f70aa4f9-2b77-48f3-95a3-ee9e0a4054d5
let fig = Makie.Figure()
	ax = Makie.Axis(fig[1,1], width=300, height=200)
	Makie.lines!(ax, equilibrated_samples_from_rbm_global_free_energy)
	Makie.resize_to_layout!(fig)
	fig
end

# ╔═╡ ba030665-6baf-4fd1-8847-8df4fa81887b
all_sampled_path_unique_sequences = TransitionPaths2024.onehot(unique([
	TransitionPaths2024.sampled_path_1to1_20240703();
	TransitionPaths2024.sampled_path_1to2rep1_20240703();
	TransitionPaths2024.sampled_path_1to4_20240703()
]));

# ╔═╡ d29bbfc0-d0b7-4188-84b5-24e76ea96b4b
path_natural_extremeties = [
	TransitionPaths2024.sampled_path_1to1_20240703()[[1,end]];
	TransitionPaths2024.sampled_path_1to2rep1_20240703()[[1,end]];
	TransitionPaths2024.sampled_path_1to4_20240703()[[1,end]];
];

# ╔═╡ f38bfcd8-c5fb-4b66-8a36-779c3ee41981
natural_to_eq_path_samples_distances = hamming(
	TransitionPaths2024.onehot(setdiff(TransitionPaths2024.Eugenio_MSA_20230419(), path_natural_extremeties)), all_sampled_path_unique_sequences
);

# ╔═╡ 4ac68e13-6137-45e0-a118-c36f14794318
natural_to_eq_samples_distances = hamming(TransitionPaths2024.onehot(unique(TransitionPaths2024.Eugenio_MSA_20230419())), equilibrated_samples_from_rbm_global);

# ╔═╡ e3744752-d14d-4385-b39a-06fd0a6b25e6
let fig = Makie.Figure()
	ax = Makie.Axis(fig[1,1]; width=300, height=250, xlabel="Hamming distance to closest natural", xticks=0:5:20, ylabel="Frequency")
	Makie.hist!(ax, vec(minimum(natural_to_eq_samples_distances; dims=1)); normalization=:pdf, bins=0:17, label="equilibrium samples")
	Makie.hist!(ax, vec(minimum(natural_to_eq_path_samples_distances; dims=1)); normalization=:pdf, bins=0:17, label="sampled paths")
	Makie.axislegend(ax)
	Makie.resize_to_layout!(fig)
	fig
end

# ╔═╡ e7ae3d8f-998e-4422-b235-f57a1308ed6d
md"# Sequene logos"

# ╔═╡ ab5b254b-43f9-474f-8496-8d708e9d2e82
TransitionPaths2024.seqlogo_entropic(dropdims(mean(TransitionPaths2024.onehot(TransitionPaths2024.Eugenio_MSA_20230419()); dims=3); dims=3)).fig

# ╔═╡ 3b596356-3652-42c8-9dcc-1e744225e2bb
TransitionPaths2024.seqlogo_entropic(dropdims(mean(TransitionPaths2024.onehot(TransitionPaths2024.Eugenio_MSA_20230419()[msa_local_rbm_class_prediction .== :I]); dims=3); dims=3)).fig

# ╔═╡ c55755eb-8fdf-4cdb-9159-2f5fd3ca8dc9
TransitionPaths2024.seqlogo_entropic(dropdims(mean(TransitionPaths2024.onehot(TransitionPaths2024.Eugenio_MSA_20230419()[msa_local_rbm_class_prediction .== :II]); dims=3); dims=3)).fig

# ╔═╡ 2e808e48-8167-451b-a175-535d194292a1
TransitionPaths2024.seqlogo_entropic(dropdims(mean(TransitionPaths2024.onehot(TransitionPaths2024.Eugenio_MSA_20230419()[msa_local_rbm_class_prediction .== :IV]); dims=3); dims=3)).fig

# ╔═╡ df92a6e3-f2ea-4bde-a6b5-e620632ecc26
md"## Weights determining specificity"

# ╔═╡ 6c0cf19b-79e3-4456-aa25-6d3acbc720bb
begin
	#color_scheme = TransitionPaths2024.sequence_logo_color_scheme()
	AAs = collect("ACDEFGHIKLMNPQRSTVWY⊟")
	color_scheme = Logomaker.color_scheme(
        'C' => "green",
        (a => "orange" for a = ('F', 'W', 'Y'))...,
        (a => "purple" for a = ('Q', 'N', 'S', 'T'))...,
        (a => "black" for a = ('V', 'L', 'I', 'M'))...,
        (a => "blue" for a = ('K', 'R', 'H'))...,
        (a => "red" for a = ('D', 'E'))...,
        (a => "grey" for a = ('A', 'P', 'G'))...,
        '⊟' => "black"
    )
	
	fig, ax = PythonPlot.subplots(2, 1, figsize=[6,3])
	
	logo = Logomaker.Logo(-Float64.(TransitionPaths2024.Eugenio_RBM_20230419(:global).w[:,:,38]), AAs; ax=ax[0], color_scheme, flip_below=false)
	logo.ax.set_ylabel("W h.u. #38")
	#logo.ax.set_xlabel("site")

	logo = Logomaker.Logo(+Float64.(TransitionPaths2024.Eugenio_RBM_20230419(:global).w[:,:,36]), AAs; ax=ax[1], color_scheme, flip_below=false)
	logo.ax.set_ylabel("W h.u. #36")
	logo.ax.set_xlabel("site")

	fig.tight_layout()
	
	fig
end

# ╔═╡ 46bf2dc8-1fd8-4a07-8318-ebfd7146f608
let fig = Makie.Figure()
	w2 = sum(abs2, TransitionPaths2024.Eugenio_RBM_20230419(:global).w[:,:,38]; dims=1)
	ax = Makie.Axis(fig[1,1]; width=600, height=10)
	Makie.heatmap!(ax, w2'; colormap=:tempo)
	Makie.hidedecorations!(ax)
	Makie.hidespines!(ax)
	Makie.resize_to_layout!(fig)
	fig
end

# ╔═╡ 1b6ea41a-c1fc-4270-84b7-c8bc35546898
let fig = Makie.Figure()
	w2 = sum(abs2, TransitionPaths2024.Eugenio_RBM_20230419(:global).w[:,:,36]; dims=1)
	ax = Makie.Axis(fig[1,1]; width=600, height=10)
	Makie.heatmap!(ax, w2'; colormap=:tempo)
	Makie.hidedecorations!(ax)
	Makie.hidespines!(ax)	
	Makie.resize_to_layout!(fig)
	fig
end

# ╔═╡ ee72ebf0-d19a-4d7d-856d-6a34d39c6f64
md"# Scores per class"

# ╔═╡ d66e29e4-9304-484d-8bc9-64d3c372420f
let fig = Makie.Figure()
	bins = -65:-15
	width = 500
	height = 100
	ylabel = "frequency"
	ax = Makie.Axis(fig[1,1]; width, height, xgridvisible=false, ygridvisible=false, ylabel, title="Full MSA")
	Makie.hist!(ax, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:global, TransitionPaths2024.onehot(TransitionPaths2024.Eugenio_MSA_20230419())); normalization=:pdf, bins)
	ax = Makie.Axis(fig[2,1]; width, height, xgridvisible=false, ygridvisible=false, ylabel, title="Class I")
	Makie.hist!(ax, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:global, TransitionPaths2024.onehot(TransitionPaths2024.Eugenio_MSA_20230419()[msa_local_rbm_class_prediction .== :I])); normalization=:pdf, bins)
	ax = Makie.Axis(fig[3,1]; width, height, xgridvisible=false, ygridvisible=false, ylabel, title="Class II")
	Makie.hist!(ax, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:global, TransitionPaths2024.onehot(TransitionPaths2024.Eugenio_MSA_20230419()[msa_local_rbm_class_prediction .== :II])); normalization=:pdf, bins)
	ax = Makie.Axis(fig[4,1]; width, height, xgridvisible=false, ygridvisible=false, ylabel, title="Class IV", xlabel=L"\log(P_\mathrm{rbm}(\text{global}))")
	Makie.hist!(ax, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:global, TransitionPaths2024.onehot(TransitionPaths2024.Eugenio_MSA_20230419()[msa_local_rbm_class_prediction .== :IV])); normalization=:pdf, bins)
	Makie.resize_to_layout!(fig)
	fig
end

# ╔═╡ ddc3a193-e723-48b7-b751-d5d55a1c668d
md"# Contact maps"

# ╔═╡ d280337a-4d54-4037-ab05-87eee053be30
# gives position in full sequence corresponding to a MSA match column
function align_map_msa_to_full(aln_seq_with_inserts::AbstractString)
	aln_map = Int[]
	position_in_msa = 0
	position_in_full = 0
	for (i,c) = enumerate(aln_seq_with_inserts)
		if islowercase(c)
			position_in_full += 1
		elseif c == '-'
			position_in_msa += 1
			push!(aln_map, 0)
		elseif isuppercase(c)
			position_in_full += 1
			position_in_msa += 1
			push!(aln_map, position_in_full)
		else
			error()
		end
	end
	return aln_map
end

# ╔═╡ ca227bca-6e39-439a-a0b1-314905427d2f
temp_pdb_dir = mktempdir()

# ╔═╡ f8baf3e0-d13a-4e88-8d45-0db5b27ec77a
pdb_1ywi = BioStructures.retrievepdb("1ywi"; dir=temp_pdb_dir)

# ╔═╡ b54c6ad1-6b53-4d67-89b6-0f3e2b928887
pdb_2ltw = BioStructures.retrievepdb("2ltw"; dir=temp_pdb_dir)

# ╔═╡ 3e7fb491-4697-4d83-ac43-8ec002ad801b
pdb_1i8g = BioStructures.retrievepdb("1I8G"; dir=temp_pdb_dir)

# ╔═╡ 2103df73-a592-4dda-bdac-df585bfdb3a7
full_seq_1ywi = LongAA(pdb_1ywi['A'])

# ╔═╡ 99d18896-729b-44f8-b594-b140811141f0
full_seq_2ltw = LongAA(pdb_2ltw['A'])

# ╔═╡ c8059e51-6122-430a-98d6-128d7f7aba21
full_seq_1i8g = LongAA(pdb_1i8g['B'])

# ╔═╡ c20b9b9e-45b9-47a2-8d99-4cecec1452eb
aln_seq_1ywi_with_inserts, aln_seq_2ltw_with_inserts, aln_seq_1i8g_with_inserts = map([full_seq_1ywi, full_seq_2ltw, full_seq_1i8g]) do seq
	return only(TransitionPaths2024.align_to_PF00397([seq]; remove_inserts=false))
end

# ╔═╡ f1180d4d-da49-44a4-b8f2-4f5a83a5c425
dist_map_1ywi = BioStructures.DistanceMap(pdb_1ywi['A'])

# ╔═╡ d800b387-dc5f-4e96-839a-72f18d949a55
dist_map_2ltw = BioStructures.DistanceMap(pdb_2ltw['A'])

# ╔═╡ 508812bf-8d45-43d2-9e9f-70d6485b5f2d
dist_map_1i8g = BioStructures.DistanceMap(pdb_1i8g['B'])

# ╔═╡ 06c406c0-d8c4-4ca1-8357-fc04b5f4b281
aln_map_1ywi = align_map_msa_to_full(aln_seq_1ywi_with_inserts)

# ╔═╡ f2d74f94-0e12-40dc-b834-55af6fb5076f
aln_map_2ltw = align_map_msa_to_full(aln_seq_2ltw_with_inserts)

# ╔═╡ a7f991c7-f088-402c-951c-74bf7481fb8d
aln_map_1i8g = align_map_msa_to_full(aln_seq_1i8g_with_inserts)

# ╔═╡ ec0eb26b-8818-4197-a75a-bd694728bca0
md"## Predicted contacts"

# ╔═╡ d5f3b04c-87e5-4a8b-9ed3-dce071b9d093
pdb_contact_distance_threshold = 8 # in Angstroms

# ╔═╡ 8837f58a-73b8-4e6e-b0a7-c7b8faba94be
begin
	# APC-corrected contact map
	Cmap = TransitionPaths2024.effective_contacts(TransitionPaths2024.Eugenio_RBM_20230419(:global), TransitionPaths2024.onehot(TransitionPaths2024.Eugenio_MSA_20230419()))
	
	# zero-out diagonal
	for i in axes(Cmap, 1)
	    Cmap[i,i] = 0
	end
end

# ╔═╡ 3bd5e451-753d-4fe8-a8f3-ac106cf762c8
let fig = Makie.Figure()
	pdb_map = [iszero(i) || iszero(j) ? NaN : dist_map_1ywi[i,j] for i = aln_map_1ywi, j = aln_map_1ywi] .< pdb_contact_distance_threshold

	ax = Makie.Axis(fig[1,1], yreversed=true, width=400, height=400, xticks=0:5:size(Cmap,1), yticks=0:5:size(Cmap,2), xlabel="sequence position", ylabel="sequence position")
	#hm = Makie.heatmap!(ax, abs.(plotted_matrix), colormap=["white", "white", "gray", "dimgray", "black", "black"])
	
	plotted_matrix = [i < j ? Cmap[i,j] : NaN for i = axes(Cmap, 1), j = axes(Cmap, 2)]
	hm = Makie.heatmap!(ax, plotted_matrix; colormap=:balance)
	Makie.hidespines!(ax, :t, :r)
	Makie.Colorbar(fig[1,2], hm, height=200, label="Epistasis score")
	plotted_matrix = [i < j ? NaN : pdb_map[i,j] for i = axes(Cmap, 1), j = axes(Cmap, 2)]
	hm = Makie.heatmap!(ax, plotted_matrix; colormap=:reds)
	Makie.hidespines!(ax, :t, :r)
	Makie.resize_to_layout!(fig)
	fig
end

# ╔═╡ c9794e23-393e-4ea5-9e67-41ab7442470a
let fig = Makie.Figure()
	pdb_map = [iszero(i) || iszero(j) ? NaN : dist_map_2ltw[i,j] for i = aln_map_2ltw, j = aln_map_2ltw] .< pdb_contact_distance_threshold
	
	ax = Makie.Axis(fig[1,1], yreversed=true, width=400, height=400, xticks=0:5:size(Cmap,1), yticks=0:5:size(Cmap,2), xlabel="sequence position", ylabel="sequence position")
	#hm = Makie.heatmap!(ax, abs.(plotted_matrix), colormap=["white", "white", "gray", "dimgray", "black", "black"])
	
	plotted_matrix = [i < j ? NaN : pdb_map[i,j] for i = axes(Cmap, 1), j = axes(Cmap, 2)]
	hm = Makie.heatmap!(ax, plotted_matrix; colormap=:greens)
	
	plotted_matrix = [i < j ? Cmap[i,j] : NaN for i = axes(Cmap, 1), j = axes(Cmap, 2)]
	hm = Makie.heatmap!(ax, plotted_matrix; colormap=:balance)
	
	Makie.hidespines!(ax, :t, :r)
	Makie.Colorbar(fig[1,2], hm, height=200, label="Epistasis score")
	Makie.resize_to_layout!(fig)
	fig
end

# ╔═╡ a7762c4a-3cf9-4f21-982f-dcadaa4c417c
let fig = Makie.Figure()
	pdb_map = [iszero(i) || iszero(j) ? NaN : dist_map_1i8g[i,j] for i = aln_map_1i8g, j = aln_map_1i8g] .< pdb_contact_distance_thre shold

	
	ax = Makie.Axis(fig[1,1], yreversed=true, width=400, height=400, xticks=0:5:size(Cmap,1), yticks=0:5:size(Cmap,2), xlabel="sequence position", ylabel="sequence position")
	#hm = Makie.heatmap!(ax, abs.(plotted_matrix), colormap=["white", "white", "gray", "dimgray", "black", "black"])
	
	plotted_matrix = [i < j ? NaN : pdb_map[i,j] for i = axes(Cmap, 1), j = axes(Cmap, 2)]
	hm = Makie.heatmap!(ax, plotted_matrix; colormap=:greens)
	
	plotted_matrix = [i < j ? Cmap[i,j] : NaN for i = axes(Cmap, 1), j = axes(Cmap, 2)]
	hm = Makie.heatmap!(ax, plotted_matrix; colormap=:balance)
	
	Makie.hidespines!(ax, :t, :r)
	Makie.Colorbar(fig[1,2], hm, height=200, label="Epistasis score")
	Makie.resize_to_layout!(fig)
	fig
end

# ╔═╡ 4e5cfd73-e7d9-4a5d-8bc5-3fd8bea45a11
md"# Wild-types"

# ╔═╡ 3376c2a3-f407-4c69-b391-4958b46658ca
let fig = Makie.Figure()
	ax = Makie.Axis(fig[1:2,1], width=250, height=250, xlabel=L"I_{1}", ylabel=L"I_{2}", xgridvisible=false, ygridvisible=false) # h_38 -> h_1, h_36 -> h_2
	Makie.vlines!(ax, -1; color=:gray, linestyle=:dash, linewidth=1)
	Makie.hlines!(ax, -3; color=:gray, linestyle=:dash, linewidth=1)
	for (class, color) = zip((:I, :II, :IV), _class_scatter_colors)
		Makie.scatter!(
			ax,
			-inputs_h_from_v(TransitionPaths2024.Eugenio_RBM_20230419(:global), TransitionPaths2024.onehot(TransitionPaths2024.Eugenio_MSA_20230419()))[38, msa_local_rbm_class_prediction .== class],
			 inputs_h_from_v(TransitionPaths2024.Eugenio_RBM_20230419(:global), TransitionPaths2024.onehot(TransitionPaths2024.Eugenio_MSA_20230419()))[36, msa_local_rbm_class_prediction .== class];
			color=(color, 0.5), markersize=4
		)
	end
	Makie.xlims!(ax, -4.5, 2)
	Makie.ylims!(ax, -7, 3)

	response_data_v2_sequences = TransitionPaths2024.Exp_20240703_sequences().sequences[
		indexin(TransitionPaths2024.response_data_all_v2_20240708().name, TransitionPaths2024.Exp_20240703_sequences().names)
	]

	Makie.scatter!(
		ax,
		-inputs_h_from_v(TransitionPaths2024.Eugenio_RBM_20230419(:global), TransitionPaths2024.onehot(df_20241213.aligned_sequences[df_20241213.group .== "Wild-types"]))[38, 1:11],
		 inputs_h_from_v(TransitionPaths2024.Eugenio_RBM_20230419(:global), TransitionPaths2024.onehot(df_20241213.aligned_sequences[df_20241213.group .== "Wild-types"]))[36, 1:11];
		color=:blue, markersize=20, marker=:cross, label="I"
	)

	Makie.scatter!(
		ax,
		-inputs_h_from_v(TransitionPaths2024.Eugenio_RBM_20230419(:global), TransitionPaths2024.onehot(df_20241213.aligned_sequences[df_20241213.group .== "Wild-types"]))[38, 12:26],
		 inputs_h_from_v(TransitionPaths2024.Eugenio_RBM_20230419(:global), TransitionPaths2024.onehot(df_20241213.aligned_sequences[df_20241213.group .== "Wild-types"]))[36, 12:26];
		color=:red, markersize=20, marker=:cross, label="II"
	)

	Makie.scatter!(
		ax,
		-inputs_h_from_v(TransitionPaths2024.Eugenio_RBM_20230419(:global), TransitionPaths2024.onehot(df_20241213.aligned_sequences[df_20241213.group .== "Wild-types"]))[38, 27:27],
		 inputs_h_from_v(TransitionPaths2024.Eugenio_RBM_20230419(:global), TransitionPaths2024.onehot(df_20241213.aligned_sequences[df_20241213.group .== "Wild-types"]))[36, 27:27];
		color=:black, markersize=20, marker=:cross, label="IV"
	)

	Makie.axislegend(ax; framevisible=false, position=:lb)


	_ylim = (-190, 10)
	
	ax = Makie.Axis(fig[1,2], width=400, height=70, xgridvisible=false, ygridvisible=false)
	Makie.scatter!(ax, 1:11, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:global, TransitionPaths2024.onehot(df_20241213.aligned_sequences[df_20241213.group .== "Wild-types"][1:11])), color=:blue, label="I", marker=:cross, markersize=10)
	Makie.scatter!(ax, 12:26, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:global, TransitionPaths2024.onehot(df_20241213.aligned_sequences[df_20241213.group .== "Wild-types"][12:26])), color=:red, label="II", marker=:cross, markersize=10)
	Makie.scatter!(ax, 27:27, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:global, TransitionPaths2024.onehot(df_20241213.aligned_sequences[df_20241213.group .== "Wild-types"][27:27])), color=:black, label="IV", marker=:cross, markersize=10)
	Makie.ylims!(ax, -100, 0)
	#Makie.axislegend(ax; framevisible=false, position=:lb)

	ax = Makie.Axis(fig[2,2], width=400, height=150, xgridvisible=false, ygridvisible=false)
	Makie.lines!(ax, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:I, TransitionPaths2024.onehot(df_20241213.aligned_sequences[df_20241213.group .== "Wild-types"])), color=:cyan, linewidth=4)
	Makie.scatter!(ax, 1:11, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:I, TransitionPaths2024.onehot(df_20241213.aligned_sequences[df_20241213.group .== "Wild-types"][1:11])), color=:blue, marker=:cross, markersize=10)
	Makie.scatter!(ax, 12:26, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:I, TransitionPaths2024.onehot(df_20241213.aligned_sequences[df_20241213.group .== "Wild-types"][12:26])), color=:red, marker=:cross, markersize=10)
	Makie.scatter!(ax, 27:27, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:I, TransitionPaths2024.onehot(df_20241213.aligned_sequences[df_20241213.group .== "Wild-types"][27:27])), color=:black, marker=:cross, markersize=10)
	Makie.ylims!(ax, _ylim...)
	
	#ax = Makie.Axis(fig[3,1], width=500, height=75, xgridvisible=false, ygridvisible=false)
	Makie.lines!(ax, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:II, TransitionPaths2024.onehot(df_20241213.aligned_sequences[df_20241213.group .== "Wild-types"])), color=:orange, linewidth=4)
	Makie.scatter!(ax, 1:11, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:II, TransitionPaths2024.onehot(df_20241213.aligned_sequences[df_20241213.group .== "Wild-types"][1:11])), color=:blue, marker=:cross, markersize=10)
	Makie.scatter!(ax, 12:26, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:II, TransitionPaths2024.onehot(df_20241213.aligned_sequences[df_20241213.group .== "Wild-types"][12:26])), color=:red, marker=:cross, markersize=10)
	Makie.scatter!(ax, 27:27, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:II, TransitionPaths2024.onehot(df_20241213.aligned_sequences[df_20241213.group .== "Wild-types"][27:27])), color=:black, marker=:cross, markersize=10)
	Makie.ylims!(ax, _ylim...)
	
	#ax = Makie.Axis(fig[4,1], width=500, height=75, xgridvisible=false, ygridvisible=false)
	Makie.lines!(ax, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:IV, TransitionPaths2024.onehot(df_20241213.aligned_sequences[df_20241213.group .== "Wild-types"])), color=:green, linewidth=4)
	Makie.scatter!(ax, 1:11, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:IV, TransitionPaths2024.onehot(df_20241213.aligned_sequences[df_20241213.group .== "Wild-types"][1:11])), color=:blue, marker=:cross, markersize=10)
	Makie.scatter!(ax, 12:26, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:IV, TransitionPaths2024.onehot(df_20241213.aligned_sequences[df_20241213.group .== "Wild-types"][12:26])), color=:red, marker=:cross, markersize=10)
	Makie.scatter!(ax, 27:27, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:IV, TransitionPaths2024.onehot(df_20241213.aligned_sequences[df_20241213.group .== "Wild-types"][27:27])), color=:black, marker=:cross, markersize=10)
	Makie.ylims!(ax, _ylim...)

	Makie.resize_to_layout!(fig)
	fig
end

# ╔═╡ 2a23e69c-0169-4e54-9f95-da5fbf97a657
let fig = Makie.Figure()
	_ylim = (-190, 10)
	
	ax = Makie.Axis(fig[1,1], width=500, height=75, xgridvisible=false, ygridvisible=false)
	Makie.scatter!(ax, 1:11, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:global, TransitionPaths2024.onehot(df_20241213.aligned_sequences[df_20241213.group .== "Wild-types"][1:11])), color=:blue, label="I", marker=:cross, markersize=10)
	Makie.scatter!(ax, 12:26, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:global, TransitionPaths2024.onehot(df_20241213.aligned_sequences[df_20241213.group .== "Wild-types"][12:26])), color=:red, label="II", marker=:cross, markersize=10)
	Makie.scatter!(ax, 27:27, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:global, TransitionPaths2024.onehot(df_20241213.aligned_sequences[df_20241213.group .== "Wild-types"][27:27])), color=:black, label="IV", marker=:cross, markersize=10)
	Makie.ylims!(ax, -100, 0)
	#Makie.axislegend(ax; framevisible=false, position=:lb)

	ax = Makie.Axis(fig[2,1], width=500, height=200, xgridvisible=false, ygridvisible=false)
	Makie.lines!(ax, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:I, TransitionPaths2024.onehot(df_20241213.aligned_sequences[df_20241213.group .== "Wild-types"])), color=:cyan, linewidth=4)
	Makie.scatter!(ax, 1:11, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:I, TransitionPaths2024.onehot(df_20241213.aligned_sequences[df_20241213.group .== "Wild-types"][1:11])), color=:blue, marker=:cross, markersize=10)
	Makie.scatter!(ax, 12:26, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:I, TransitionPaths2024.onehot(df_20241213.aligned_sequences[df_20241213.group .== "Wild-types"][12:26])), color=:red, marker=:cross, markersize=10)
	Makie.scatter!(ax, 27:27, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:I, TransitionPaths2024.onehot(df_20241213.aligned_sequences[df_20241213.group .== "Wild-types"][27:27])), color=:black, marker=:cross, markersize=10)
	Makie.ylims!(ax, _ylim...)
	
	#ax = Makie.Axis(fig[3,1], width=500, height=75, xgridvisible=false, ygridvisible=false)
	Makie.lines!(ax, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:II, TransitionPaths2024.onehot(df_20241213.aligned_sequences[df_20241213.group .== "Wild-types"])), color=:orange, linewidth=4)
	Makie.scatter!(ax, 1:11, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:II, TransitionPaths2024.onehot(df_20241213.aligned_sequences[df_20241213.group .== "Wild-types"][1:11])), color=:blue, marker=:cross, markersize=10)
	Makie.scatter!(ax, 12:26, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:II, TransitionPaths2024.onehot(df_20241213.aligned_sequences[df_20241213.group .== "Wild-types"][12:26])), color=:red, marker=:cross, markersize=10)
	Makie.scatter!(ax, 27:27, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:II, TransitionPaths2024.onehot(df_20241213.aligned_sequences[df_20241213.group .== "Wild-types"][27:27])), color=:black, marker=:cross, markersize=10)
	Makie.ylims!(ax, _ylim...)
	
	#ax = Makie.Axis(fig[4,1], width=500, height=75, xgridvisible=false, ygridvisible=false)
	Makie.lines!(ax, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:IV, TransitionPaths2024.onehot(df_20241213.aligned_sequences[df_20241213.group .== "Wild-types"])), color=:green, linewidth=4)
	Makie.scatter!(ax, 1:11, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:IV, TransitionPaths2024.onehot(df_20241213.aligned_sequences[df_20241213.group .== "Wild-types"][1:11])), color=:blue, marker=:cross, markersize=10)
	Makie.scatter!(ax, 12:26, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:IV, TransitionPaths2024.onehot(df_20241213.aligned_sequences[df_20241213.group .== "Wild-types"][12:26])), color=:red, marker=:cross, markersize=10)
	Makie.scatter!(ax, 27:27, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:IV, TransitionPaths2024.onehot(df_20241213.aligned_sequences[df_20241213.group .== "Wild-types"][27:27])), color=:black, marker=:cross, markersize=10)
	Makie.ylims!(ax, _ylim...)
	
	Makie.resize_to_layout!(fig)
	fig
end

# ╔═╡ 24493642-3e88-4312-ad51-4345ef6aafd8
let fig = Makie.Figure()
	sampled_path = unique(TransitionPaths2024.sampled_path_1to1_20240703())

	ax = Makie.Axis(fig[1,2][1,1]; width=250, height=75, xgridvisible=false, ygridvisible=false, ylabel=L"\ln P_\mathrm{RBM}(\text{glob.})")
	Makie.scatterlines!(ax, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:global, TransitionPaths2024.onehot(sampled_path)); color=:black)
	Makie.hidexdecorations!(ax)

	ax = Makie.Axis(
		fig[1,2][2,1]; width=250, height=75, xgridvisible=false, ygridvisible=false,
		ylabel=L"\log P_\mathrm{RBM}(\text{loc.})", #yticklabelcolor=_class_scatter_colors[1], ylabelcolor=_class_scatter_colors[1]
	)
	Makie.scatterlines!(ax, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:I, TransitionPaths2024.onehot(sampled_path)); color=:cyan)
	Makie.scatterlines!(ax, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:II, TransitionPaths2024.onehot(sampled_path)); color=:orange)
	Makie.scatterlines!(ax, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:IV, TransitionPaths2024.onehot(sampled_path)); color=:green)
	Makie.hidexdecorations!(ax)


	ax = Makie.Axis(fig[1,1], width=200, height=200, xlabel=L"I_{1}", ylabel=L"I_{2}", xgridvisible=false, ygridvisible=false) # h_38 -> h_1, h_36 -> h_2
	Makie.vlines!(ax, -1; color=:gray, linestyle=:dash, linewidth=1)
	Makie.hlines!(ax, -3; color=:gray, linestyle=:dash, linewidth=1)
	for (class, color) = zip((:I, :II, :IV), _class_scatter_colors)
		Makie.scatter!(
			ax,
			-inputs_h_from_v(TransitionPaths2024.Eugenio_RBM_20230419(:global), TransitionPaths2024.onehot(TransitionPaths2024.Eugenio_MSA_20230419()))[38, msa_local_rbm_class_prediction .== class],
			 inputs_h_from_v(TransitionPaths2024.Eugenio_RBM_20230419(:global), TransitionPaths2024.onehot(TransitionPaths2024.Eugenio_MSA_20230419()))[36, msa_local_rbm_class_prediction .== class];
			color=(color, 0.5), markersize=4
		)
	end

	Makie.scatterlines!(ax, 
		-inputs_h_from_v(TransitionPaths2024.Eugenio_RBM_20230419(:global), TransitionPaths2024.onehot(sampled_path))[38, :],
		 inputs_h_from_v(TransitionPaths2024.Eugenio_RBM_20230419(:global), TransitionPaths2024.onehot(sampled_path))[36, :];
		color=:black
	)

	Makie.arrows!(ax,
		-inputs_h_from_v(TransitionPaths2024.Eugenio_RBM_20230419(:global), TransitionPaths2024.onehot(sampled_path))[38, 1:1] .+ 0.35,
		 inputs_h_from_v(TransitionPaths2024.Eugenio_RBM_20230419(:global), TransitionPaths2024.onehot(sampled_path))[36, 1:1] .- 0.35,
		[-0.15], [0.15]; linewidth=2, color=:gray, arrowsize=10
	)
	
	Makie.xlims!(ax, -4.5, 2)
	Makie.ylims!(ax, -7, 3)

	Makie.resize_to_layout!(fig)
	#Makie.save("/DATA-SSD/cossio/projects/2024/Transition_Paths/TransitionPaths2024.jl/figures/Fig2.pdf", fig)
	fig
end

# ╔═╡ 2fecd867-e4d7-4267-b12e-52fc02cea62d
let fig = Makie.Figure()
	sampled_path_1 = unique(TransitionPaths2024.sampled_path_1to2rep1_20240703())
	sampled_path_2 = unique(TransitionPaths2024.sampled_path_1to2rep2_20240703())

	ax = Makie.Axis(fig[1,2][1,1]; width=250, height=75, xgridvisible=false, ygridvisible=false, ylabel=L"\ln P_\mathrm{RBM}(\text{glob.})")
	Makie.lines!(ax, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:global, TransitionPaths2024.onehot(sampled_path_1)); color=:black)
	Makie.lines!(ax, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:global, TransitionPaths2024.onehot(sampled_path_2)); color=:black, linestyle=:dash)
	Makie.hidexdecorations!(ax)

	ax = Makie.Axis(
		fig[1,2][2,1]; width=250, height=75, xgridvisible=false, ygridvisible=false,
		ylabel=L"\log P_\mathrm{RBM}(\text{loc.})", #yticklabelcolor=_class_scatter_colors[1], ylabelcolor=_class_scatter_colors[1]
	)
	Makie.lines!(ax, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:I, TransitionPaths2024.onehot(sampled_path_1)); color=:cyan)
	Makie.lines!(ax, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:I, TransitionPaths2024.onehot(sampled_path_2)); color=:cyan, linestyle=:dash)
	Makie.lines!(ax, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:II, TransitionPaths2024.onehot(sampled_path_1)); color=:orange)
	Makie.lines!(ax, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:II, TransitionPaths2024.onehot(sampled_path_2)); color=:orange, linestyle=:dash)
	Makie.lines!(ax, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:IV, TransitionPaths2024.onehot(sampled_path_1)); color=:green)
	Makie.lines!(ax, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:IV, TransitionPaths2024.onehot(sampled_path_2)); color=:green, linestyle=:dash)
	Makie.hidexdecorations!(ax)

	ax = Makie.Axis(fig[1,1], width=200, height=200, xlabel=L"I_{1}", ylabel=L"I_{2}", xgridvisible=false, ygridvisible=false) # h_38 -> h_1, h_36 -> h_2
	Makie.vlines!(ax, -1; color=:gray, linestyle=:dash, linewidth=1)
	Makie.hlines!(ax, -3; color=:gray, linestyle=:dash, linewidth=1)
	for (class, color) = zip((:I, :II, :IV), _class_scatter_colors)
		Makie.scatter!(
			ax,
			-inputs_h_from_v(TransitionPaths2024.Eugenio_RBM_20230419(:global), TransitionPaths2024.onehot(TransitionPaths2024.Eugenio_MSA_20230419()))[38, msa_local_rbm_class_prediction .== class],
			 inputs_h_from_v(TransitionPaths2024.Eugenio_RBM_20230419(:global), TransitionPaths2024.onehot(TransitionPaths2024.Eugenio_MSA_20230419()))[36, msa_local_rbm_class_prediction .== class];
			color=(color, 0.5), markersize=4
		)
	end

	Makie.lines!(ax, 
		-inputs_h_from_v(TransitionPaths2024.Eugenio_RBM_20230419(:global), TransitionPaths2024.onehot(sampled_path_1))[38, :],
		 inputs_h_from_v(TransitionPaths2024.Eugenio_RBM_20230419(:global), TransitionPaths2024.onehot(sampled_path_1))[36, :];
		color=:black, linewidth=2
	)
	Makie.lines!(ax, 
		-inputs_h_from_v(TransitionPaths2024.Eugenio_RBM_20230419(:global), TransitionPaths2024.onehot(sampled_path_2))[38, :],
		 inputs_h_from_v(TransitionPaths2024.Eugenio_RBM_20230419(:global), TransitionPaths2024.onehot(sampled_path_2))[36, :];
		color=:black, linestyle=:dash, linewidth=2
	)
		
	Makie.xlims!(ax, -4.5, 2)
	Makie.ylims!(ax, -7, 3)

	Makie.resize_to_layout!(fig)
	#Makie.save("/DATA-SSD/cossio/projects/2024/Transition_Paths/TransitionPaths2024.jl/figures/Fig2.pdf", fig)
	fig
end

# ╔═╡ ef351be2-07bb-4f2a-9d0e-80518788c890
let fig = Makie.Figure()
	sampled_path_main = unique(TransitionPaths2024.sampled_path_1to4_20240703())
	sampled_path_batch1 = unique(TransitionPaths2024.sampled_path_1to4batch1_20240703())
	sampled_path_batch2 = unique(TransitionPaths2024.sampled_path_1to4batch2_20240703())
	sampled_path_direct = unique(TransitionPaths2024.sampled_path_1to4direct_20240703())

	ax = Makie.Axis(fig[1,2][1,1]; width=250, height=75, xgridvisible=false, ygridvisible=false, ylabel=L"\ln P_\mathrm{RBM}(\text{glob.})")
	Makie.lines!(ax, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:global, TransitionPaths2024.onehot(sampled_path_main)); color=:black)
	Makie.lines!(ax, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:global, TransitionPaths2024.onehot(sampled_path_direct)); color=:black, linestyle=:dash)
	Makie.lines!(ax, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:global, TransitionPaths2024.onehot(sampled_path_batch1)); color=:black, linestyle=:dot)
	Makie.lines!(ax, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:global, TransitionPaths2024.onehot(sampled_path_batch2)); color=:black, linestyle=:dashdot)
	Makie.hidexdecorations!(ax)

	ax = Makie.Axis(
		fig[1,2][2,1]; width=250, height=75, xgridvisible=false, ygridvisible=false,
		ylabel=L"\log P_\mathrm{RBM}(\text{loc.})", #yticklabelcolor=_class_scatter_colors[1], ylabelcolor=_class_scatter_colors[1]
	)
	for (model, color) = zip((:I, :II, :IV), (:cyan, :orange, :green))
		Makie.lines!(ax, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(model, TransitionPaths2024.onehot(sampled_path_main)); color)
		Makie.lines!(ax, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(model, TransitionPaths2024.onehot(sampled_path_direct)); color, linestyle=:dash)
		Makie.lines!(ax, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(model, TransitionPaths2024.onehot(sampled_path_batch1)); color, linestyle=:dot)
		Makie.lines!(ax, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(model, TransitionPaths2024.onehot(sampled_path_batch2)); color, linestyle=:dashdot)
	end
	Makie.hidexdecorations!(ax)

	ax = Makie.Axis(fig[1,1], width=200, height=200, xlabel=L"I_{1}", ylabel=L"I_{2}", xgridvisible=false, ygridvisible=false) # h_38 -> h_1, h_36 -> h_2
	Makie.vlines!(ax, -1; color=:gray, linestyle=:dash, linewidth=1)
	Makie.hlines!(ax, -3; color=:gray, linestyle=:dash, linewidth=1)
	for (class, color) = zip((:I, :II, :IV), _class_scatter_colors)
		Makie.scatter!(
			ax,
			-inputs_h_from_v(TransitionPaths2024.Eugenio_RBM_20230419(:global), TransitionPaths2024.onehot(TransitionPaths2024.Eugenio_MSA_20230419()))[38, msa_local_rbm_class_prediction .== class],
			 inputs_h_from_v(TransitionPaths2024.Eugenio_RBM_20230419(:global), TransitionPaths2024.onehot(TransitionPaths2024.Eugenio_MSA_20230419()))[36, msa_local_rbm_class_prediction .== class];
			color=(color, 0.5), markersize=4
		)
	end

	for path = (sampled_path_main, sampled_path_direct, sampled_path_batch1, sampled_path_batch2)
		Makie.lines!(ax, 
			-inputs_h_from_v(TransitionPaths2024.Eugenio_RBM_20230419(:global), TransitionPaths2024.onehot(path))[38, :],
			 inputs_h_from_v(TransitionPaths2024.Eugenio_RBM_20230419(:global), TransitionPaths2024.onehot(path))[36, :];
			color=:black, linewidth=1
		)
	end
		
	Makie.xlims!(ax, -4.5, 2)
	Makie.ylims!(ax, -7, 3)

	Makie.resize_to_layout!(fig)
	#Makie.save("/DATA-SSD/cossio/projects/2024/Transition_Paths/TransitionPaths2024.jl/figures/Fig2.pdf", fig)
	fig
end

# ╔═╡ 3ed7f78d-e40a-4c8d-b2db-c04451663b0c
md"# RBM IV and epistasis"

# ╔═╡ f29f4108-aa2a-420c-97cc-34a472f0f62d
# ╠═╡ disabled = true
#=╠═╡
single_site_mutation_costs_per_rbm = Dict(
	rbm => [
		begin
			seqs = TransitionPaths2024.onehot(TransitionPaths2024.Eugenio_MSA_20230419())
			
			seqs[:,i,:] .= 0
			seqs[a,i,:] .= 1
			ll0 = TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:I, seqs)
	
			seqs[:,i,:] .= 0
			seqs[b,i,:] .= 1
			ll1 = TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:I, seqs)
			
			mean(ll1 - ll0)
		end for a = 1:21, b = 1:21, i = 1:31
	]
	for rbm = (:global, :I, :II, :IV)
)
  ╠═╡ =#

# ╔═╡ 516f1104-1f8d-41a8-a73b-488526413a2b
#=╠═╡
single_site_mutation_costs_per_rbm_mean = Dict(
	rbm => dropdims(mean(single_site_mutation_costs_per_rbm[rbm]; dims=1); dims=1)
	for rbm = (:global, :I, :II, :IV)
)
  ╠═╡ =#

# ╔═╡ 12ec5e4d-8db0-4c65-880a-e0c7b36ce9bc
#=╠═╡
let fig = Makie.Figure()
	for (n, rbm) = enumerate((:global, :I, :II, :IV))
		seqs = TransitionPaths2024.onehot(TransitionPaths2024.Eugenio_MSA_20230419())
		lls_full = TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(rbm, seqs)
		lls_inde = vec(single_site_mutation_costs_per_rbm_mean[rbm])' * reshape(seqs, :, size(seqs)[end])
	
		ax = Makie.Axis(fig[1,n]; with=300, height=300)
		Makie.scatter!(ax, lls_full .- mean(lls_full), lls_inde .- mean(lls_inde))
	end
	Makie.resize_to_layout!(fig)
	fig
end
  ╠═╡ =#

# ╔═╡ 8fc31f9d-de38-4d8a-90ac-53a716d55da9
# ╠═╡ disabled = true
#=╠═╡
single_site_IV_mutation_costs = [
	begin
		seqs = TransitionPaths2024.onehot(TransitionPaths2024.Eugenio_MSA_20230419())
		
		seqs[:,i,:] .= 0
		seqs[a,i,:] .= 1
		ll0 = TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:IV, seqs)

		seqs[:,i,:] .= 0
		seqs[b,i,:] .= 1
		ll1 = TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:IV, seqs)
		
		mean(ll1 - ll0)
	end for a = 1:21, b = 1:21, i = 1:31
]
  ╠═╡ =#

# ╔═╡ 251a2b2d-3a52-4a65-bf94-2228cf599fc0
#=╠═╡
dropdims(mean(single_site_IV_mutation_costs; dims=1); dims=1)
  ╠═╡ =#

# ╔═╡ 9d150281-a7b0-449a-af71-ee77993502de
let fig = Makie.Figure()
	for (n, rbm) = enumerate((:global, :I, :II, :IV))
		# APC-corrected contact map
		Cmap = TransitionPaths2024.effective_contacts(TransitionPaths2024.Eugenio_RBM_20230419(rbm), TransitionPaths2024.onehot(TransitionPaths2024.Eugenio_MSA_20230419()))
		
		# zero-out diagonal
		for i in axes(Cmap, 1)
		    Cmap[i,i] = 0
		end
	
		plotted_matrix = copy(Cmap)
		for i = axes(plotted_matrix, 1), j = axes(plotted_matrix, 2)
		    if i < j
		        plotted_matrix[i,j] = Cmap[i,j]
		    else
		        plotted_matrix[i,j] = 0
		    end
		end
	
		ax = Makie.Axis(fig[1,n]; yreversed=true, title=string(rbm), width=175, height=175, xticks=0:10:size(Cmap,1), yticks=0:10:size(Cmap,2), xlabel="sequence position", ylabel="sequence position")
		#hm = Makie.heatmap!(ax, abs.(plotted_matrix), colormap=["white", "white", "gray", "dimgray", "black", "black"])
		hm = Makie.heatmap!(ax, plotted_matrix; colormap=:balance, colorrange=(-1.5,5))
		Makie.hidespines!(ax, :t, :r)

		if n == 4
			Makie.Colorbar(fig[1,5], hm, height=200, label="Epistasis score")
		end
	end

	epistasis_norms = [
		begin
			# APC-corrected contact map
			Cmap = TransitionPaths2024.effective_contacts(TransitionPaths2024.Eugenio_RBM_20230419(rbm), TransitionPaths2024.onehot(TransitionPaths2024.Eugenio_MSA_20230419()))
			
			# zero-out diagonal
			for i in axes(Cmap, 1)
			    Cmap[i,i] = 0
			end
	
			norm(Cmap)
		end for rbm = (:global, :I, :II, :IV)
	]

	ax = Makie.Axis(fig[1,6]; width=100, height=175, xticks=(1:4, ["glob.", "I", "II", "IV"]), ylabel="Epistasis norm")
	Makie.barplot!(ax, 1:4, epistasis_norms)
	Makie.ylims!(ax, 0, 33)

	Makie.resize_to_layout!(fig)
	#Makie.save("Figures/Contacts bw.pdf", fig)
	fig
end

# ╔═╡ 538c4027-6cf0-47fc-bbe7-d277df95bf2b
md"# Send paths to PG"

# ╔═╡ 00ecfd29-ca78-41d1-95e2-d0ba3bd601fc
TransitionPaths2024.sampled_path_1to1_20240703()

# ╔═╡ 28c5e03b-b8fa-4896-983f-8a2c5f639ba9
TransitionPaths2024.sampled_path_1to2rep1_20240703()

# ╔═╡ 343cee63-9ebc-4f07-ac33-c8b23e47ff70
md"# Trash"

# ╔═╡ 2aadc194-117e-4100-abdc-755626d6577a
let fig = Makie.Figure()
	_sz = 200

	pred_class = TransitionPaths2024.Eugenio_Predict_Class_from_Inputs_20240715(TransitionPaths2024.onehot(response_data_sequences))
	likelihood_I  = TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:I,  TransitionPaths2024.onehot(response_data_sequences))
	likelihood_II = TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:II, TransitionPaths2024.onehot(response_data_sequences))
	likelihood_IV = TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:IV, TransitionPaths2024.onehot(response_data_sequences))
	response_I  = [r ≤ 0 ? NaN : r for r = all_response_data.response_data_C1]
	response_II = [r ≤ 0 ? NaN : r for r = all_response_data.response_data_C2]
	response_IV = [r ≤ 0 ? NaN : r for r = all_response_data.response_data_C4]
	
	ax_I  = Makie.Axis(fig[1,1][1,1], width=_sz, height=_sz, xlabel=L"\log(P_\mathrm{rbm}(\text{local I}))",  ylabel="Class I. normalized response",  yscale=log10, xgridvisible=false, ygridvisible=false)
	ax_II = Makie.Axis(fig[1,1][1,2], width=_sz, height=_sz, xlabel=L"\log(P_\mathrm{rbm}(\text{local II}))", ylabel="Class II. normalized response", yscale=log10, xgridvisible=false, ygridvisible=false)
	ax_IV = Makie.Axis(fig[1,1][1,3], width=_sz, height=_sz, xlabel=L"\log(P_\mathrm{rbm}(\text{local IV}))", ylabel="Class IV. normalized response", yscale=log10, xgridvisible=false, ygridvisible=false)

	for (class, color) = zip((:I, :II, :IV), _class_scatter_colors)
		Makie.scatter!(ax_I,  likelihood_I[pred_class .== class],  response_I[pred_class .== class]; color)
		Makie.scatter!(ax_II, likelihood_II[pred_class .== class], response_II[pred_class .== class]; color)
		Makie.scatter!(ax_IV, likelihood_IV[pred_class .== class], response_IV[pred_class .== class]; color)
	end

	labels = ["I", "II", "IV"]
	elements = [Makie.MarkerElement(; marker=:circle, color = _class_scatter_colors[j]) for j = 1:3]
	Makie.Legend(fig[1,1][1,4], elements, labels; framed=false)

	# Correlations barplot
	responses = Dict(
		:I => all_response_data.response_data_C1,
		:II => all_response_data.response_data_C2,
		:IV => all_response_data.response_data_C4
	)

	cat = [i for (i, target) = enumerate((:I, :II, :IV)) for (j, rbm) = enumerate((:global, :I, :II, :IV))]
	grp = [j for (i, target) = enumerate((:I, :II, :IV)) for (j, rbm) = enumerate((:global, :I, :II, :IV))]
				
	correlations = [
		cor(TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(rbm, TransitionPaths2024.onehot(response_data_sequences)), responses[target])
		for target = (:I, :II, :IV) for rbm = (:global, :I, :II, :IV)
	]

	colors = [:gray, :cyan, :orange, :green]
		
	ax = Makie.Axis(fig[2,1][1,2], width=250, height=_sz, xticks=(1:3, ["I", "II", "IV"]), ylabel="Experimental response corr.", xgridvisible=false, ygridvisible=false)
	Makie.barplot!(ax, cat, correlations; dodge=grp, color=colors[grp], dodge_gap=0.05, width=0.9)
	Makie.hlines!(ax, 0.0; color=:black, linestyle=:dash)
	Makie.ylims!(-0.4, 0.8)

	labels = ["RBM(glob.)", "RBM(I)", "RBM(II)", "RBM(IV)"]
	elements = [Makie.PolyElement(; polycolor = colors[j]) for j = 1:4]
	Makie.Legend(fig[2,1][1,3], elements, labels; framed=false)

	response_data_v2_sequences = TransitionPaths2024.Exp_20240703_sequences().sequences[
		indexin(TransitionPaths2024.response_data_all_v2_20240708().name, TransitionPaths2024.Exp_20240703_sequences().names)
	]
	rbm_log_likelihoods_v2 = TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:global, TransitionPaths2024.onehot(response_data_v2_sequences))
	
	ax = Makie.Axis(fig[2,1][1,1], width=300, height=_sz, xlabel=L"P_{\mathrm{RBM}}(\text{global})", ylabel="Frequency", xgridvisible=false, ygridvisible=false)
	Makie.hist!(ax, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:global, TransitionPaths2024.onehot(TransitionPaths2024.Eugenio_MSA_20230419())); normalization=:pdf, bins=-90:5:10, color=(:green, 0.3), label="MSA")
	Makie.stephist!(ax, rbm_log_likelihoods_v2[TransitionPaths2024.response_data_all_v2_20240708().type .== "Natural"], color=:green, label="Natural", normalization=:pdf, bins=-90:5:10, linewidth=4)
	Makie.stephist!(ax, rbm_log_likelihoods_v2[TransitionPaths2024.response_data_all_v2_20240708().type .== "Designed"], color=:blue, label="Designed", normalization=:pdf, bins=-90:5:10, linewidth=4)
	Makie.stephist!(ax, rbm_log_likelihoods_v2[TransitionPaths2024.response_data_all_v2_20240708().type .== "Shuffled"], color=:red, label="Shuffled", normalization=:pdf, bins=-90:5:10, linewidth=4)
	Makie.xlims!(ax, -90, 0)
	Makie.ylims!(ax, 0, 0.1)
	Makie.axislegend(ax; position=:lt, framevisible=false)

	Makie.resize_to_layout!(fig)
	#Makie.save("/DATA-SSD/cossio/projects/2024/Transition_Paths/TransitionPaths2024.jl/figures/Fig5.pdf", fig)
	fig
end

# ╔═╡ 971044ac-8df4-4fec-8d4a-5cf20ab2ece1
let fig = Makie.Figure()
	_sz = 200
	
	ax = Makie.Axis(fig[1,1], width=_sz, height=_sz, xlabel=L"\log(P_\mathrm{rbm}(\text{local I}))", ylabel="Class I. normalized response", yscale=log10, xgridvisible=false, ygridvisible=false)
	Makie.scatter!(ax, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:I, TransitionPaths2024.onehot(response_data_sequences)), [r ≤ 0 ? NaN : r for r = all_response_data.response_data_C1]; color=:black)

	ax = Makie.Axis(fig[1,2], width=_sz, height=_sz, xlabel=L"\log(P_\mathrm{rbm}(\text{local II}))", ylabel="Class II. normalized response", yscale=log10, xgridvisible=false, ygridvisible=false)
	Makie.scatter!(ax, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:II, TransitionPaths2024.onehot(response_data_sequences)), [r ≤ 0 ? NaN : r for r = all_response_data.response_data_C2]; color=:black)
	
	ax = Makie.Axis(fig[1,3], width=_sz, height=_sz, xlabel=L"\log(P_\mathrm{rbm}(\text{local IV}))", ylabel="Class IV. normalized response", yscale=log10, xgridvisible=false, ygridvisible=false)
	Makie.scatter!(ax, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:IV, TransitionPaths2024.onehot(response_data_sequences)), [r ≤ 0 ? NaN : r for r = all_response_data.response_data_C4]; color=:black)

	ax = Makie.Axis(fig[2,1], width=_sz, height=_sz, xlabel=L"\log(P_\mathrm{rbm}(\text{Global}))", ylabel="Class I. normalized response", yscale=log10, xgridvisible=false, ygridvisible=false)
	Makie.scatter!(ax, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:global, TransitionPaths2024.onehot(response_data_sequences)), [r ≤ 0 ? NaN : r for r = all_response_data.response_data_C1]; color=:black)

	ax = Makie.Axis(fig[2,2], width=_sz, height=_sz, xlabel=L"\log(P_\mathrm{rbm}(\text{Global}))", ylabel="Class II. normalized response", yscale=log10, xgridvisible=false, ygridvisible=false)
	Makie.scatter!(ax, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:global, TransitionPaths2024.onehot(response_data_sequences)), [r ≤ 0 ? NaN : r for r = all_response_data.response_data_C2]; color=:black)
	
	ax = Makie.Axis(fig[2,3], width=_sz, height=_sz, xlabel=L"\log(P_\mathrm{rbm}(\text{Global}))", ylabel="Class IV. normalized response", yscale=log10, xgridvisible=false, ygridvisible=false)
	Makie.scatter!(ax, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:global, TransitionPaths2024.onehot(response_data_sequences)), [r ≤ 0 ? NaN : r for r = all_response_data.response_data_C4]; color=:black)

	# Correlations barplot
	responses = Dict(
		:I => all_response_data.response_data_C1,
		:II => all_response_data.response_data_C2,
		:IV => all_response_data.response_data_C4
	)

	cat = [i for (i, target) = enumerate((:I, :II, :IV)) for (j, rbm) = enumerate((:global, :I, :II, :IV))]
	grp = [j for (i, target) = enumerate((:I, :II, :IV)) for (j, rbm) = enumerate((:global, :I, :II, :IV))]
	
	# correlations = [
	# 	cor(
	# 		TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(rbm, TransitionPaths2024.onehot(response_data_sequences))[findall(>(0), responses[target])], 
	# 		log.(filter(>(0), responses[target]))
	# 	)
	# 	for target = (:I, :II, :IV) for rbm = (:global, :I, :II, :IV)
	# ]
			
	correlations = [
		cor(TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(rbm, TransitionPaths2024.onehot(response_data_sequences)), responses[target])
		for target = (:I, :II, :IV) for rbm = (:global, :I, :II, :IV)
	]

	colors = [:gray, :cyan, :orange, :green]
		
	ax = Makie.Axis(fig[1,4], width=300, height=_sz, xticks=(1:3, ["I", "II", "IV"]), ylabel="Experimental response corr.", xgridvisible=false, ygridvisible=false)
	Makie.barplot!(ax, cat, correlations; dodge=grp, color=colors[grp], dodge_gap=0.05, width=0.9)
	Makie.hlines!(ax, 0.0; color=:black, linestyle=:dash)
	Makie.ylims!(-0.4, 0.8)

	labels = ["RBM(glob.)", "RBM(I)", "RBM(II)", "RBM(IV)"]
	elements = [Makie.PolyElement(; polycolor = colors[j]) for j = 1:4]
	Makie.Legend(fig[1,5], elements, labels; framed=false)

	
	response_data_v2_sequences = TransitionPaths2024.Exp_20240703_sequences().sequences[
		indexin(TransitionPaths2024.response_data_all_v2_20240708().name, TransitionPaths2024.Exp_20240703_sequences().names)
	]
	rbm_log_likelihoods_v2 = TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:global, TransitionPaths2024.onehot(response_data_v2_sequences))
	
	ax = Makie.Axis(fig[2,4], width=300, height=_sz, xlabel=L"P_{\mathrm{RBM}}(\text{global})", ylabel="Frequency", xgridvisible=false, ygridvisible=false)
	Makie.hist!(ax, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:global, TransitionPaths2024.onehot(TransitionPaths2024.Eugenio_MSA_20230419())); normalization=:pdf, bins=-90:5:10, color=(:green, 0.3), label="MSA")
	Makie.stephist!(ax, rbm_log_likelihoods_v2[TransitionPaths2024.response_data_all_v2_20240708().type .== "Natural"], color=:green, label="Natural", normalization=:pdf, bins=-90:5:10, linewidth=4)
	Makie.stephist!(ax, rbm_log_likelihoods_v2[TransitionPaths2024.response_data_all_v2_20240708().type .== "Designed"], color=:blue, label="Designed", normalization=:pdf, bins=-90:5:10, linewidth=4)
	Makie.stephist!(ax, rbm_log_likelihoods_v2[TransitionPaths2024.response_data_all_v2_20240708().type .== "Shuffled"], color=:red, label="Shuffled", normalization=:pdf, bins=-90:5:10, linewidth=4)
	Makie.xlims!(ax, -90, 0)
	Makie.ylims!(ax, 0, 0.1)
	Makie.axislegend(ax; position=:lt, framevisible=false)

	Makie.resize_to_layout!(fig)
	#Makie.save("/DATA-SSD/cossio/projects/2024/Transition_Paths/TransitionPaths2024.jl/figures/Fig5.pdf", fig)
	fig
end

# ╔═╡ 96a94126-c0c4-4489-8b33-972bcea21c95
md"### I -> II Rep.2"

# ╔═╡ c9e5f0b1-91b5-4d27-9199-c3de59592921
let fig = Makie.Figure()
	ax = Makie.Axis(fig[1,1][1,1], width=250, height=250, xlabel=L"I_{1}", ylabel=L"I_{2}", xgridvisible=false, ygridvisible=false) # h_38 -> h_1, h_36 -> h_2
	Makie.vlines!(ax, -1; color=:gray, linestyle=:dash, linewidth=1)
	Makie.hlines!(ax, -3; color=:gray, linestyle=:dash, linewidth=1)
	for (class, color) = zip((:I, :II, :IV), _class_scatter_colors)
		Makie.scatter!(
			ax,
			-inputs_h_from_v(TransitionPaths2024.Eugenio_RBM_20230419(:global), TransitionPaths2024.onehot(TransitionPaths2024.Eugenio_MSA_20230419()))[38, msa_local_rbm_class_prediction .== class],
			 inputs_h_from_v(TransitionPaths2024.Eugenio_RBM_20230419(:global), TransitionPaths2024.onehot(TransitionPaths2024.Eugenio_MSA_20230419()))[36, msa_local_rbm_class_prediction .== class];
			color=(color, 0.5), markersize=4
		)
	end

	sampled_path = unique(TransitionPaths2024.sampled_path_1to2rep2_20240703())
	path_response_data = TransitionPaths2024.response_data_path_1to2rep2_20240703()
	probed_seqs = TransitionPaths2024.Exp_20240703_sequences().sequences[indexin(path_response_data.ww_names, TransitionPaths2024.Exp_20240703_sequences().names)]
	probed_seqs_idx::Vector{Int} = indexin(probed_seqs, sampled_path)
	
	Makie.scatterlines!(ax, 
		-inputs_h_from_v(TransitionPaths2024.Eugenio_RBM_20230419(:global), TransitionPaths2024.onehot(sampled_path))[38, :],
		 inputs_h_from_v(TransitionPaths2024.Eugenio_RBM_20230419(:global), TransitionPaths2024.onehot(sampled_path))[36, :];
		color=:black
	)

	Makie.scatter!(ax, 
		-inputs_h_from_v(TransitionPaths2024.Eugenio_RBM_20230419(:global), TransitionPaths2024.onehot(probed_seqs))[38, :],
		 inputs_h_from_v(TransitionPaths2024.Eugenio_RBM_20230419(:global), TransitionPaths2024.onehot(probed_seqs))[36, :];
		color=:red, markersize=10
	)
	
	Makie.arrows!(ax,
		-inputs_h_from_v(TransitionPaths2024.Eugenio_RBM_20230419(:global), TransitionPaths2024.onehot(sampled_path))[38, 1:1] .+ 0.35,
		 inputs_h_from_v(TransitionPaths2024.Eugenio_RBM_20230419(:global), TransitionPaths2024.onehot(sampled_path))[36, 1:1] .- 0.35,
		[-0.15], [0.15]; linewidth=2, color=:gray, arrowsize=10
	)
	
	Makie.xlims!(ax, -4.5, 2)
	Makie.ylims!(ax, -7, 3)

	ax = Makie.Axis(fig[1,1][2,1]; width=250, height=250, xgridvisible=false, ygridvisible=false, xlabel="Norm. response I", ylabel="Norm. response II")
	Makie.band!(ax, [-10, thresh_significance], [-10, -10], [thresh_significance, thresh_significance]; color=(:red, 0.3)) # neither specificity is significant region
	Makie.band!(ax, [-10, thresh_significance], [thresh_significance, thresh_significance], [100, 100]; color=(:darkorange, 0.3)) # specific to II
	Makie.band!(ax, [thresh_significance, 100], [-100, -100], [thresh_significance, thresh_significance]; color=(:cyan, 0.3)) # specific to II
	Makie.scatterlines!(ax, path_response_data.C1_response, path_response_data.C2_response; color=:black, linewidth=1)
	Makie.errorbars!(ax, path_response_data.C1_response, path_response_data.C2_response, path_response_data.C1_response_err; color=:gray, direction=:x, linewidth=2, whiskerwidth=5)
	Makie.errorbars!(ax, path_response_data.C1_response, path_response_data.C2_response, path_response_data.C2_response_err; color=:gray, direction=:y, linewidth=2, whiskerwidth=5)
	Makie.scatter!(ax, path_response_data.C1_response, path_response_data.C2_response; color=:red)
	Makie.arrows!(ax, path_response_data.C1_response[1:1] .+ 0.2, path_response_data.C2_response[1:1] .- 0.15, [-0.07], [0.07]; linewidth=2, color=:gray, arrowsize=10)
	Makie.xlims!(ax, -0.5, 10)
	Makie.ylims!(ax, -0.1, 3)

	ax = Makie.Axis(fig[1,2][1,1]; width=300, height=100, xgridvisible=false, ygridvisible=false, ylabel=L"\log P_\mathrm{RBM}(\text{global})")
	Makie.scatterlines!(ax, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:global, TransitionPaths2024.onehot(sampled_path)); color=:black)
	Makie.scatter!(ax, probed_seqs_idx, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:global, TransitionPaths2024.onehot(probed_seqs)); color=:red, markersize=10, label="Designed")
	Makie.hidexdecorations!(ax)

	ax = Makie.Axis(fig[1,2][2,1]; width=300, height=100, xgridvisible=false, ygridvisible=false, ylabel=L"\log P_\mathrm{RBM}(\text{I})")
	Makie.scatterlines!(ax, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:I, TransitionPaths2024.onehot(sampled_path)); color=:black, label="Designed")
	Makie.scatter!(ax, probed_seqs_idx, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:I, TransitionPaths2024.onehot(probed_seqs)); color=:red, markersize=10)
	Makie.hidexdecorations!(ax)
	Makie.axislegend(ax; position=:lb, framevisible=false)

	ax = Makie.Axis(fig[1,2][3,1]; width=300, height=100, xgridvisible=false, ygridvisible=false, ylabel=L"\log P_\mathrm{RBM}(\text{II})")
	Makie.scatterlines!(ax, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:II, TransitionPaths2024.onehot(sampled_path)); color=:black)
	Makie.scatter!(ax, probed_seqs_idx, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:II, TransitionPaths2024.onehot(probed_seqs)); color=:red, markersize=10)
	Makie.hidexdecorations!(ax)

	ax = Makie.Axis(fig[1,2][4,1]; width=300, height=30, xgridvisible=false, ygridvisible=false, yscale=log10, xticks=0:5:50)
	Makie.hidedecorations!(ax)
	Makie.hidespines!(ax)
	
	ax = Makie.Axis(fig[1,2][5,1]; width=300, height=100, xgridvisible=false, ygridvisible=false, ylabel="I/II rel. resp.", yscale=log10, xticks=0:5:50)
	Makie.hspan!(ax, 1, 1e6; color=(:cyan, 0.3))
	Makie.hspan!(ax, 1e-6, 1; color=(:darkorange, 0.3))
	Makie.scatterlines!(ax, probed_seqs_idx, path_response_data.C1_response ./ path_response_data.C2_response; color=:black)
	Makie.scatter!(ax, probed_seqs_idx, path_response_data.C1_response ./ path_response_data.C2_response; color=:red, markersize=10)
	Makie.ylims!(ax, 1e-3, 1e3)

	Makie.resize_to_layout!(fig)
	#Makie.save("/DATA-SSD/cossio/projects/2024/Transition_Paths/TransitionPaths2024.jl/figures/SI_Fig_1.pdf", fig)
	fig
end

# ╔═╡ dee03fac-a312-4e1c-a056-7e2e25a73ce4
md"### I -> IV Direct"

# ╔═╡ f843ff96-94bb-4860-9d6b-0d8292d25159
let fig = Makie.Figure()
	ax = Makie.Axis(fig[1,1][1,1], width=250, height=250, xlabel=L"I_{38}", ylabel=L"I_{36}", xgridvisible=false, ygridvisible=false)
	Makie.vlines!(ax, -1; color=:gray, linestyle=:dash, linewidth=1)
	Makie.hlines!(ax, -3; color=:gray, linestyle=:dash, linewidth=1)
	for (class, color) = zip((:I, :II, :IV), _class_scatter_colors)
		Makie.scatter!(
			ax,
			-inputs_h_from_v(TransitionPaths2024.Eugenio_RBM_20230419(:global), TransitionPaths2024.onehot(TransitionPaths2024.Eugenio_MSA_20230419()))[38, msa_local_rbm_class_prediction .== class],
			 inputs_h_from_v(TransitionPaths2024.Eugenio_RBM_20230419(:global), TransitionPaths2024.onehot(TransitionPaths2024.Eugenio_MSA_20230419()))[36, msa_local_rbm_class_prediction .== class];
			color=(color, 0.5), markersize=4
		)
	end

	sampled_path = unique(TransitionPaths2024.sampled_path_1to4direct_20240703())
	path_response_data = TransitionPaths2024.response_data_path_1to4direct_20240703()
	probed_seqs = TransitionPaths2024.Exp_20240703_sequences().sequences[indexin(path_response_data.ww_names, TransitionPaths2024.Exp_20240703_sequences().names)]
	probed_seqs_idx::Vector{Int} = indexin(probed_seqs, sampled_path)

	Makie.scatterlines!(ax, 
		-inputs_h_from_v(TransitionPaths2024.Eugenio_RBM_20230419(:global), TransitionPaths2024.onehot(sampled_path))[38, :],
		 inputs_h_from_v(TransitionPaths2024.Eugenio_RBM_20230419(:global), TransitionPaths2024.onehot(sampled_path))[36, :];
		color=:black
	)

	Makie.scatter!(ax, 
		-inputs_h_from_v(TransitionPaths2024.Eugenio_RBM_20230419(:global), TransitionPaths2024.onehot(probed_seqs))[38, :],
		 inputs_h_from_v(TransitionPaths2024.Eugenio_RBM_20230419(:global), TransitionPaths2024.onehot(probed_seqs))[36, :];
		color=:red, markersize=10
	)
	
	Makie.arrows!(ax,
		-inputs_h_from_v(TransitionPaths2024.Eugenio_RBM_20230419(:global), TransitionPaths2024.onehot(sampled_path))[38, 1:1] .+ 0.35,
		 inputs_h_from_v(TransitionPaths2024.Eugenio_RBM_20230419(:global), TransitionPaths2024.onehot(sampled_path))[36, 1:1] .- 0.35,
		[-0.15], [0.15]; linewidth=2, color=:gray, arrowsize=10
	)
	
	Makie.xlims!(ax, -4.5, 2)
	Makie.ylims!(ax, -7, 3)

	ax = Makie.Axis(fig[1,1][2,1]; width=250, height=250, xgridvisible=false, ygridvisible=false, xlabel="Norm. response I", ylabel="Norm. response IV")
	Makie.band!(ax, [-10, thresh_significance], [-10, -10], [thresh_significance, thresh_significance]; color=(:red, 0.3)) # neither specificity is significant region
	Makie.band!(ax, [-10, thresh_significance], [thresh_significance, thresh_significance], [100, 100]; color=(:green, 0.3)) # specific to IV
	Makie.band!(ax, [thresh_significance, 100], [-100, -100], [thresh_significance, thresh_significance]; color=(:cyan, 0.3)) # specific to I
	Makie.scatterlines!(ax, path_response_data.C1_response, path_response_data.C4_response; color=:black, linewidth=1)
	Makie.errorbars!(ax, path_response_data.C1_response, path_response_data.C4_response, path_response_data.C1_response_err; color=:gray, direction=:x, linewidth=2, whiskerwidth=5)
	Makie.errorbars!(ax, path_response_data.C1_response, path_response_data.C4_response, path_response_data.C4_response_err; color=:gray, direction=:y, linewidth=2, whiskerwidth=5)
	Makie.scatter!(ax, path_response_data.C1_response, path_response_data.C4_response; color=:red)
	Makie.arrows!(ax, path_response_data.C1_response[1:1] .+ 0.2, path_response_data.C4_response[1:1] .- 0.15, [-0.07], [0.07]; linewidth=2, color=:gray, arrowsize=10)
	Makie.xlims!(ax, -0.5, 12)
	Makie.ylims!(ax, -0.25, 2)

	ax = Makie.Axis(fig[1,2][1,1]; width=300, height=100, xgridvisible=false, ygridvisible=false, ylabel=L"\log P_\mathrm{RBM}(\text{global})")
	Makie.scatterlines!(ax, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:global, TransitionPaths2024.onehot(sampled_path)); color=:black)
	Makie.scatter!(ax, probed_seqs_idx, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:global, TransitionPaths2024.onehot(probed_seqs)); color=:red, markersize=10, label="Designed")
	Makie.hidexdecorations!(ax)

	ax = Makie.Axis(fig[1,2][2,1]; width=300, height=100, xgridvisible=false, ygridvisible=false, ylabel=L"\log P_\mathrm{RBM}(\text{I})")
	Makie.scatterlines!(ax, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:I, TransitionPaths2024.onehot(sampled_path)); color=:black, label="Designed")
	Makie.scatter!(ax, probed_seqs_idx, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:I, TransitionPaths2024.onehot(probed_seqs)); color=:red, markersize=10)
	Makie.hidexdecorations!(ax)
	Makie.axislegend(ax; position=:lb, framevisible=false)

	ax = Makie.Axis(fig[1,2][3,1]; width=300, height=100, xgridvisible=false, ygridvisible=false, ylabel=L"\log P_\mathrm{RBM}(\text{IV})")
	Makie.scatterlines!(ax, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:II, TransitionPaths2024.onehot(sampled_path)); color=:black)
	Makie.scatter!(ax, probed_seqs_idx, TransitionPaths2024.Eugenio_RBM_20230419_loglikelihood(:II, TransitionPaths2024.onehot(probed_seqs)); color=:red, markersize=10)
	Makie.hidexdecorations!(ax)

	rel_response = path_response_data.C1_response ./ path_response_data.C4_response

	ax = Makie.Axis(fig[1,2][4,1]; width=300, height=30, xgridvisible=false, ygridvisible=false, yscale=log10, xticks=0:5:50)
	Makie.hidedecorations!(ax)
	Makie.hidespines!(ax)

	ax = Makie.Axis(fig[1,2][5,1]; width=300, height=100, xgridvisible=false, ygridvisible=false, ylabel="I/IV rel. resp.", yscale=log10, xticks=0:5:50)
	Makie.hspan!(ax, 1, 1e6; color=(:cyan, 0.3))
	Makie.hspan!(ax, 1e-6, 1; color=(:green, 0.3))
	Makie.scatterlines!(ax, [i for (i,x) = zip(probed_seqs_idx, rel_response) if x > 0], [x for x = rel_response if x > 0]; color=:black)
	Makie.scatter!(ax, [i for (i,x) = zip(probed_seqs_idx, rel_response) if x > 0], [x for x = rel_response if x > 0]; color=:red, markersize=10)
	Makie.ylims!(ax, 1e-3, 5e3)

	Makie.resize_to_layout!(fig)
	#Makie.save("/DATA-SSD/cossio/projects/2024/Transition_Paths/TransitionPaths2024.jl/figures/SI_Fig_2.pdf", fig)
	fig
end

# ╔═╡ Cell order:
# ╠═0089bfee-3916-11ef-0452-d739816a3d18
# ╠═70eb32c0-f459-403f-a1d8-8db8a3f2ca06
# ╠═7ff741cb-ea2d-435a-80fa-9f53391a3117
# ╠═d301243a-53c7-4e7d-bb5f-056a9c7a77f8
# ╠═63e09505-9d7c-4db0-8c54-852ee51626be
# ╠═1e3060e6-f40e-44a9-937a-6965a2fa2134
# ╠═d1b1a928-3daf-4c30-8892-b6ced5c2b3a4
# ╠═1b8d68ed-ad51-4c3b-9023-e6372868e7d4
# ╠═474d459d-083e-4639-9e7f-753bde25c217
# ╠═ec73dce9-638a-4f34-8e52-29d03c83c9f6
# ╠═5f7923a5-47c3-42d3-94ce-89fef8c57871
# ╠═7809f3a5-e704-43da-af08-4e47f412c919
# ╠═8febf5c6-92bb-4b91-a022-eb6080c8f2a0
# ╠═4e3e5c24-41df-4fbd-91d9-42ab34f45cbe
# ╠═b5677fef-d1cd-46f5-8f74-7b32fdd4b030
# ╠═6165d60b-bada-4400-b76b-9deb2bd1450e
# ╠═cd1ca226-00f6-4eb4-b10e-bb2a390b236d
# ╠═d93cf7b6-05a0-40d1-8f7b-905f4cd72487
# ╠═383e921b-47d6-43ef-86f4-278ec9e0ce29
# ╠═c3ae9f17-68e0-45d1-a44c-ee4135ea44d7
# ╠═f07f8761-54e5-41a9-837a-268eb73b5e4e
# ╠═e3307dc4-638b-4cf8-af1e-b71088e9cedb
# ╠═ae2e577d-dfcb-4bd4-9e5c-41de0ce892a0
# ╠═d7959c5e-6c90-4995-9c34-b4de1c6765f1
# ╠═ed61d8ed-751f-4a0b-96af-e5168fff09bc
# ╠═401dddb8-7acd-47d6-ac95-4e2fdce29361
# ╠═f93754f9-437c-486e-bf3f-1f4d433a19dc
# ╠═f5d53e91-81ba-45d1-a92c-99e7813cb04f
# ╠═0de36432-13f5-437a-b007-a23b56bed610
# ╠═7a784b46-a25a-4d0e-8053-8286ec8991a8
# ╠═7ccab879-dd94-48fa-84ac-e4b0ce916ae1
# ╠═78e9c9a7-9280-4d2b-b495-4a725b5eee1c
# ╠═d5fb2c5a-e2dc-4f7b-9b19-e96787ddff71
# ╠═73169a6f-559e-46d3-af29-0e35bd487fb5
# ╠═abe10f1a-c2a3-4596-89d4-80252932e3ca
# ╠═f85a0a90-b33a-4533-ac39-6a497a74f685
# ╠═02a01dae-ac90-4e93-9351-716f9ec8b4f3
# ╠═81b9f641-4a13-42cc-8d43-c1b2ae00fc2a
# ╠═c1a63405-1e5a-4a9a-b7ca-63664460017f
# ╠═8de38839-0452-413b-ba6f-f09016cb4e39
# ╠═b55b9ccd-2f0c-4d66-b929-68a17361dee4
# ╠═0ca0ac10-2c96-4cc5-b6fe-ea2f142cba76
# ╠═449467a6-8cde-4653-9402-557466ffadba
# ╠═6f35366a-b6c5-4b05-a236-4b13b0f3b546
# ╠═606f2eea-b7f0-47b2-93cc-7fa482f342fe
# ╠═71f7257d-6998-4833-98c7-5513e5b23b93
# ╠═9a991c6a-8479-4e3d-90f4-c88f7d298c99
# ╠═e5fb8b94-2b1a-4544-88e1-3636bd4798aa
# ╠═0ae8f8b8-225d-4022-8971-4fb187c726d3
# ╠═0ebe69bb-00c4-4381-ab36-770b27ae356e
# ╠═d4216956-091a-4780-94ca-6077ee29a1e6
# ╠═bb6e6dc0-9e56-4c5f-80b6-89637573b4eb
# ╠═10489d92-333c-41e6-a3b4-5f29e3a3689a
# ╠═cf6d18bc-90da-4576-a8a4-1ff412124bdb
# ╠═92da5b0e-75e1-4aaa-911d-f375a36174f4
# ╠═191d062a-6da4-46f2-8dfe-e70e35adcfd5
# ╠═6887c9cb-f9be-43a7-82a4-654fb6bc008f
# ╠═776e0bf2-2fbb-4377-9ceb-08bbab82f4a6
# ╠═19fa19d9-74b5-44e5-8fdc-dd5f8d1bf8c9
# ╠═76b90e10-3325-49f8-95c0-9d48a2cbdbb7
# ╠═6184d501-0f4e-4d22-8eb0-88a07b0efab8
# ╠═d0cc1b76-9705-47ff-8c08-3402531002f2
# ╠═4354e416-143a-45ca-916e-6304accf1bdf
# ╠═9734267d-e44f-4f1d-9c7b-3464ed5a2677
# ╠═bf689107-5744-4535-b455-f0678071b29e
# ╠═3e59d567-6300-4e01-a3c9-7226acd6174d
# ╠═4a19b62a-8644-4ad1-a516-8b04bd84464a
# ╠═ef668901-d457-42fe-bcb9-0dd53bf04e28
# ╠═d2772e25-a731-433a-b10c-a92901172ea1
# ╠═38d65a29-446d-4e04-b876-b039baa00369
# ╠═340db97a-f8e3-4523-8101-d9cf6bfb549d
# ╠═735b7450-6722-46e8-a1b9-3b963d4324e2
# ╠═3a544f67-efe6-47c3-b522-d7f3a2049b58
# ╠═33cae588-2a3d-4755-a856-517ca4593384
# ╠═2b5cdff8-7b1d-4f4a-8b13-bc7bb5a73071
# ╠═7742b870-4495-498c-b430-f83635a03536
# ╠═54cab97b-82b1-4c07-b480-f2477c4a8ce9
# ╠═f70aa4f9-2b77-48f3-95a3-ee9e0a4054d5
# ╠═ba030665-6baf-4fd1-8847-8df4fa81887b
# ╠═d29bbfc0-d0b7-4188-84b5-24e76ea96b4b
# ╠═f38bfcd8-c5fb-4b66-8a36-779c3ee41981
# ╠═4ac68e13-6137-45e0-a118-c36f14794318
# ╠═e3744752-d14d-4385-b39a-06fd0a6b25e6
# ╠═e7ae3d8f-998e-4422-b235-f57a1308ed6d
# ╠═ab5b254b-43f9-474f-8496-8d708e9d2e82
# ╠═3b596356-3652-42c8-9dcc-1e744225e2bb
# ╠═c55755eb-8fdf-4cdb-9159-2f5fd3ca8dc9
# ╠═2e808e48-8167-451b-a175-535d194292a1
# ╠═df92a6e3-f2ea-4bde-a6b5-e620632ecc26
# ╠═6c0cf19b-79e3-4456-aa25-6d3acbc720bb
# ╠═46bf2dc8-1fd8-4a07-8318-ebfd7146f608
# ╠═1b6ea41a-c1fc-4270-84b7-c8bc35546898
# ╠═ee72ebf0-d19a-4d7d-856d-6a34d39c6f64
# ╠═d66e29e4-9304-484d-8bc9-64d3c372420f
# ╠═ddc3a193-e723-48b7-b751-d5d55a1c668d
# ╠═d280337a-4d54-4037-ab05-87eee053be30
# ╠═ca227bca-6e39-439a-a0b1-314905427d2f
# ╠═f8baf3e0-d13a-4e88-8d45-0db5b27ec77a
# ╠═b54c6ad1-6b53-4d67-89b6-0f3e2b928887
# ╠═3e7fb491-4697-4d83-ac43-8ec002ad801b
# ╠═2103df73-a592-4dda-bdac-df585bfdb3a7
# ╠═99d18896-729b-44f8-b594-b140811141f0
# ╠═c8059e51-6122-430a-98d6-128d7f7aba21
# ╠═c20b9b9e-45b9-47a2-8d99-4cecec1452eb
# ╠═f1180d4d-da49-44a4-b8f2-4f5a83a5c425
# ╠═d800b387-dc5f-4e96-839a-72f18d949a55
# ╠═508812bf-8d45-43d2-9e9f-70d6485b5f2d
# ╠═06c406c0-d8c4-4ca1-8357-fc04b5f4b281
# ╠═f2d74f94-0e12-40dc-b834-55af6fb5076f
# ╠═a7f991c7-f088-402c-951c-74bf7481fb8d
# ╠═ec0eb26b-8818-4197-a75a-bd694728bca0
# ╠═d5f3b04c-87e5-4a8b-9ed3-dce071b9d093
# ╠═8837f58a-73b8-4e6e-b0a7-c7b8faba94be
# ╠═3bd5e451-753d-4fe8-a8f3-ac106cf762c8
# ╠═c9794e23-393e-4ea5-9e67-41ab7442470a
# ╠═a7762c4a-3cf9-4f21-982f-dcadaa4c417c
# ╠═4e5cfd73-e7d9-4a5d-8bc5-3fd8bea45a11
# ╠═3376c2a3-f407-4c69-b391-4958b46658ca
# ╠═2a23e69c-0169-4e54-9f95-da5fbf97a657
# ╠═24493642-3e88-4312-ad51-4345ef6aafd8
# ╠═2fecd867-e4d7-4267-b12e-52fc02cea62d
# ╠═ef351be2-07bb-4f2a-9d0e-80518788c890
# ╠═3ed7f78d-e40a-4c8d-b2db-c04451663b0c
# ╠═f29f4108-aa2a-420c-97cc-34a472f0f62d
# ╠═516f1104-1f8d-41a8-a73b-488526413a2b
# ╠═12ec5e4d-8db0-4c65-880a-e0c7b36ce9bc
# ╠═8fc31f9d-de38-4d8a-90ac-53a716d55da9
# ╠═251a2b2d-3a52-4a65-bf94-2228cf599fc0
# ╠═9d150281-a7b0-449a-af71-ee77993502de
# ╠═538c4027-6cf0-47fc-bbe7-d277df95bf2b
# ╠═00ecfd29-ca78-41d1-95e2-d0ba3bd601fc
# ╠═28c5e03b-b8fa-4896-983f-8a2c5f639ba9
# ╠═343cee63-9ebc-4f07-ac33-c8b23e47ff70
# ╠═2aadc194-117e-4100-abdc-755626d6577a
# ╠═971044ac-8df4-4fec-8d4a-5cf20ab2ece1
# ╠═96a94126-c0c4-4489-8b33-972bcea21c95
# ╠═c9e5f0b1-91b5-4d27-9199-c3de59592921
# ╠═dee03fac-a312-4e1c-a056-7e2e25a73ce4
# ╠═f843ff96-94bb-4860-9d6b-0d8292d25159
