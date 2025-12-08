"""
    reweights(v)

Computes the reweighting vector for the MSA.
"""
function reweights(v::BitArray{3}; thresh::Real=0.1)
    q, L, B = size(v)
    Z = Int8.(reshape(first.(Tuple.(argmax(v; dims=1))), L, B))
    #thresh = DCAUtils.compute_theta(Z)
    W, Beff = DCAUtils.compute_weights(Z, q, thresh; verbose=true)
    return (; wts = W, Beff, thresh)
end
