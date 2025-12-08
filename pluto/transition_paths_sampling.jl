### A Pluto.jl notebook ###
# v0.20.13

using Markdown
using InteractiveUtils

# ╔═╡ 70a2b42c-81ac-11f0-0a00-1d2d7e034ef0
import Pkg, Revise; Pkg.activate(Base.current_project())

# ╔═╡ e06d0083-bbaf-4e38-9b66-b7606465befe
using RestrictedBoltzmannMachines: RBM, Potts, sample_h_from_v, inputs_v_from_h, sample_from_inputs

# ╔═╡ 1a1f8fa1-b460-47ac-860d-982823a44979
using LinearAlgebra: I

# ╔═╡ bc2d8ac8-818a-421e-a6a9-04361516f2b2
import PlutoUI, Makie, CairoMakie, TransitionPaths2024

# ╔═╡ 748205d9-ec7c-4c6f-a2c6-ad9d746c0ee5
import RestrictedBoltzmannMachines as RBMs

# ╔═╡ 702986ec-d642-4e6e-afcc-fd4edcb98f89
md"# Test"

# ╔═╡ c19fce37-6c23-46d2-be15-5b7d68e87e91
q = 21; # alphabet

# ╔═╡ 7a4e7841-2ccf-4f8c-8e89-0a4bbe05d373
L = 10; # seq. length

# ╔═╡ 946895b8-0cce-4b47-8350-1ceb0ed0fef9
M = 5; # num. hidden units

# ╔═╡ d85ede4b-672e-4da2-a539-5a70b7cb99ff
T = 50; # path length

# ╔═╡ dc1ce7c2-b972-4d1b-8596-e754d58ae1b7
B = 1; # batch size

# ╔═╡ 73fd9a49-4860-4844-a01d-db9811eb8ffb
rbm = RBMs.RBM(RBMs.Potts(; θ=zeros(q, L)), RBMs.Binary(; θ=zeros(M)), zeros(q, L, M));

# ╔═╡ b53d57cf-11b5-4829-8d07-ec7ffa84c271
mutation_costs = 100I(q);

# ╔═╡ 83ef49bd-79ec-4ff4-8571-95b8791cdf1e
v_path = TransitionPaths2024.onehot(rand(Int8(1):Int8(q), L, T, B));

# ╔═╡ aa933aa7-d2ba-4498-b21d-75dd8ccfa935
size(v_path) == (q,L,T,B)

# ╔═╡ 7227416f-4d8c-4961-a851-d84db9df8554
begin
	v_path_1 = copy(v_path)
	nsteps = 1000
	mut_steps = zeros(Int, T - 1, nsteps)
	for step = 1:nsteps
		v_path_1 .= TransitionPaths2024.sample_path(rbm, v_path_1, mutation_costs; nsteps=1)
		mut_steps[:, step] .= dropdims(sum(v_path_1[:,:,2:end,:] .!= v_path_1[:,:,1:end-1,:]; dims=(1,2)); dims=(1,2))
	end
end

# ╔═╡ 68ec1dad-bf01-40ed-a374-54238c25e696
let fig = Makie.Figure()
	ax = Makie.Axis(fig[1,1]; width=400, height=400)
	Makie.lines!(ax, dropdims(sum(mut_steps; dims=1); dims=1) / 2)
	Makie.hlines!(ax, L; color=:red, linestyle=:dash)
	Makie.resize_to_layout!(fig)
	fig
end

# ╔═╡ cb9dd1fe-32f2-40ea-9a25-88a63a51f436
dropdims(sum(mut_steps; dims=1); dims=1) / 2

# ╔═╡ ba7e8f59-2793-455c-882a-4b5b935d73fa
begin
	total_inputs = zero(float(v_path))
	total_inputs[:, :, 1:T-1, :] .+= reshape(mutation_costs * reshape(v_path[:,:,2:T,:], q, :), q, L, T-1, B) # inputs from the future
	total_inputs[:, :, 2:T, :] .+= reshape(mutation_costs' * reshape(v_path[:,:,1:T-1,:], q, :), q, L, T-1, B) # inputs from the past
end;

# ╔═╡ c71967e7-8132-491c-8715-1f88414fc6eb
mutation_costs_inputs_even_from_odd(v_path, mutation_costs) == total_inputs[:, :, 2:2:end, :]

# ╔═╡ 54d0d741-0909-4545-83bd-752d96fa3d27
mutation_costs_inputs_odd_from_even(v_path, mutation_costs) == total_inputs[:, :, 1:2:end, :]

# ╔═╡ 49b37781-b94e-4d1d-a43e-f3fdc3274991
md"# Test on WW domain"

# ╔═╡ ab576ae0-542a-49ce-b79e-c771c75dc7b8
TransitionPaths2024.sampled_path_1to1_20240703()

# ╔═╡ 25abc77b-0e09-4d09-b97a-9ad9b8b428fb
my_path = reshape(TransitionPaths2024.onehot(TransitionPaths2024.sampled_path_1to1_20240703()), 21, 31, 19, 1);

# ╔═╡ 6d1cd5aa-358e-4f89-98f8-73e77a4927e1


# ╔═╡ 9eb9c810-1d13-4cfd-9e4b-4b48d867c152
eugenio_rbm = TransitionPaths2024.Eugenio_RBM_20230419()

# ╔═╡ 8035cc81-d2a4-47ba-9107-b2fddbab0d7e
my_mutation_costs = 20I(q);

# ╔═╡ b8e87f5a-a04a-4edb-9fc0-50762775d66d
begin
	my_path_1 = copy(my_path)
	#my_T = size(my_path, 3)
	my_T = 10size(my_path, 3)
	my_nsteps = 1000
	my_mut_steps = zeros(Int, my_T - 1, my_nsteps)
	for step = 1:my_nsteps
		my_path_1 .= sample_path(eugenio_rbm, my_path_1, my_mutation_costs; nsteps=1)
		my_mut_steps[:, step] .= dropdims(sum(my_path_1[:,:,2:end,:] .!= my_path_1[:,:,1:end-1,:]; dims=(1,2)); dims=(1,2))
	end
end

# ╔═╡ 90691a4e-1ba7-4959-a3fd-8d22ba1bc2f4
my_path_1

# ╔═╡ 45d16df8-b6a0-40a6-9311-19416aa5c64e
size(my_mut_steps)

# ╔═╡ f4326d3e-6f5a-4998-b6a3-0a5c648bca69
dropdims(sum(my_mut_steps; dims=1); dims=1) / 2

# ╔═╡ 7a2bd268-1ded-4abf-8c26-56ecb357a3b6
my_mut_steps[:,end]

# ╔═╡ Cell order:
# ╠═70a2b42c-81ac-11f0-0a00-1d2d7e034ef0
# ╠═bc2d8ac8-818a-421e-a6a9-04361516f2b2
# ╠═748205d9-ec7c-4c6f-a2c6-ad9d746c0ee5
# ╠═e06d0083-bbaf-4e38-9b66-b7606465befe
# ╠═1a1f8fa1-b460-47ac-860d-982823a44979
# ╠═702986ec-d642-4e6e-afcc-fd4edcb98f89
# ╠═c19fce37-6c23-46d2-be15-5b7d68e87e91
# ╠═7a4e7841-2ccf-4f8c-8e89-0a4bbe05d373
# ╠═946895b8-0cce-4b47-8350-1ceb0ed0fef9
# ╠═d85ede4b-672e-4da2-a539-5a70b7cb99ff
# ╠═dc1ce7c2-b972-4d1b-8596-e754d58ae1b7
# ╠═73fd9a49-4860-4844-a01d-db9811eb8ffb
# ╠═b53d57cf-11b5-4829-8d07-ec7ffa84c271
# ╠═83ef49bd-79ec-4ff4-8571-95b8791cdf1e
# ╠═aa933aa7-d2ba-4498-b21d-75dd8ccfa935
# ╠═7227416f-4d8c-4961-a851-d84db9df8554
# ╠═68ec1dad-bf01-40ed-a374-54238c25e696
# ╠═cb9dd1fe-32f2-40ea-9a25-88a63a51f436
# ╠═ba7e8f59-2793-455c-882a-4b5b935d73fa
# ╠═c71967e7-8132-491c-8715-1f88414fc6eb
# ╠═54d0d741-0909-4545-83bd-752d96fa3d27
# ╠═49b37781-b94e-4d1d-a43e-f3fdc3274991
# ╠═ab576ae0-542a-49ce-b79e-c771c75dc7b8
# ╠═25abc77b-0e09-4d09-b97a-9ad9b8b428fb
# ╠═6d1cd5aa-358e-4f89-98f8-73e77a4927e1
# ╠═9eb9c810-1d13-4cfd-9e4b-4b48d867c152
# ╠═8035cc81-d2a4-47ba-9107-b2fddbab0d7e
# ╠═b8e87f5a-a04a-4edb-9fc0-50762775d66d
# ╠═90691a4e-1ba7-4959-a3fd-8d22ba1bc2f4
# ╠═45d16df8-b6a0-40a6-9311-19416aa5c64e
# ╠═f4326d3e-6f5a-4998-b6a3-0a5c648bca69
# ╠═7a2bd268-1ded-4abf-8c26-56ecb357a3b6
