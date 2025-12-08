"""
    load_rbm_from_jerome_potts_drelu_as_xrelu(;
        fields, theta_plus, theta_minus, gamma_plus, gamma_minus, weights
    )

Load an RBM with Potts visible units and dReLU hidden units, in Jerome's format.
I assume that the input arrays are taken directly from Python.
The dReLU hidden layer is converted to xReLU.
"""
function load_rbm_from_jerome_potts_drelu_as_xrelu(;
    # Potts visible layer parameters (N visible units); dimensions: q * N
    fields::AbstractMatrix,

    # dReLU hidden layer parameters (M hidden units)
    theta_plus::AbstractVector, theta_minus::AbstractVector,
    gamma_plus::AbstractVector, gamma_minus::AbstractVector,

    # weights; dimensions: q * N * M
    weights::AbstractArray{<:Real,3}
)
    potts = Potts(; θ = fields)
    drelu = dReLU(; θp = -theta_plus, θn = theta_minus, γp = gamma_plus, γn = gamma_minus)
    xrelu = xReLU(drelu)
    return RBM(potts, xrelu, weights)
end
