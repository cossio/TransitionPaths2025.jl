"""
    Eugenio_RBM_20230419()

RBM trained by Eugenio.
"""
Eugenio_RBM_20230419() = HDF5.h5open(Eugenio_RBM_20230419_path(), "r") do hdf
    fields = HDF5.read(hdf, "fields")
    theta_plus = HDF5.read(hdf, "theta_plus")
    gamma_plus = HDF5.read(hdf, "gamma_plus")
    theta_minus = HDF5.read(hdf, "theta_minus")
    gamma_minus = HDF5.read(hdf, "gamma_minus")
    weights = HDF5.read(hdf, "weights")
    return load_rbm_from_jerome_potts_drelu_as_xrelu(
        ; fields, theta_plus, theta_minus, gamma_plus, gamma_minus, weights
    )
end

"""
    Eugenio_MSA_20230419()

MSA used by Eugenio (also by Jerome). Older version of Pfam alignment.
"""
function Eugenio_MSA_20230419()
    path = Eugenio_fasta_20230419_path()
    seqs = FASTX.sequence.(FASTX.FASTA.Reader(open(path)))
    return LongAA.(seqs)
end

# from Eugenio (20230419): Wild-type sequences
const Type_I_wt  = aa"LPPGWERRADSL-GRTYYVDHNTRTTTWTRP"
const Type_IV_wt = aa"LPPGWEKRMSRSSGRVYYFNHITNASQWERP"
const Type_II_wt = aa"AKSMWTEHKSPD-GRTYYYNTETKQSTWEKP"

# Hidden units that clusterize sequences by specificity, as in Fig. 3E of J.Tubiana (2019)
# Note that Eugenio's indices are 0-based, hence we add 1 here.
const Cluster_Hidden_Units = [37, 35] .+ 1

"""
    Eugenio_RBM_specific_20230424()

RBM trained by Eugenio, specific to one of Type I, II, or IV sequences.
"""
function Eugenio_RBM_specific_20230424(type::Symbol)
    if type === :I
        path = Eugenio_RBM_typeI_20230424_path()
    elseif type === :II
        path = Eugenio_RBM_typeII_20230424_path()
    elseif type === :IV
        path = Eugenio_RBM_typeIV_20230424_path()
    else
        throw(ArgumentError("type must be one of :I, :II, or :IV"))
    end

    rbm = HDF5.h5open(path, "r") do hdf
        fields = HDF5.read(hdf, "fields")
        theta_plus = HDF5.read(hdf, "theta_plus")
        gamma_plus = HDF5.read(hdf, "gamma_plus")
        theta_minus = HDF5.read(hdf, "theta_minus")
        gamma_minus = HDF5.read(hdf, "gamma_minus")
        weights = HDF5.read(hdf, "weights")
        return load_rbm_from_jerome_potts_drelu_as_xrelu(
            ; fields, theta_plus, theta_minus, gamma_plus, gamma_minus, weights
        )
    end

    return rbm
end

"""
    Literature_sequences_20230424()

Sequences validated in the literature, originally collected by Jerome.
"""
function Literature_sequences_20230424()
    path = Literature_sequences_20230424_path()
    return HDF5.h5open(path, "r") do hdf
        sequences_1 = HDF5.read(hdf, "sequences_1") # type I
        sequences_2 = HDF5.read(hdf, "sequences_2") # type II
        sequences_3 = HDF5.read(hdf, "sequences_3") # type III
        sequences_4 = HDF5.read(hdf, "sequences_4") # type IV
        sequences_unknown = HDF5.read(hdf, "sequences_unknown") # unknown type ...
        # Convert from Potts (0-based indexing) to LongAA
        return (;
            sequences_I = aaseq(Int8.(sequences_1 .+ 1)),
            sequences_II = aaseq(Int8.(sequences_2 .+ 1)),
            sequences_III = aaseq(Int8.(sequences_3 .+ 1)),
            sequences_IV = aaseq(Int8.(sequences_4 .+ 1)),
            sequences_unknown = aaseq(Int8.(sequences_unknown .+ 1))
        )
    end
end

function Eugenio_Probed_Sequences_202306()
    path = Eugenio_Probed_Sequences_202306_path()
    seqs = FASTX.sequence.(FASTX.FASTA.Reader(open(path)))
    return LongAA.(seqs)
end

"""
    Eugenio_RBM_logZ_20230419(which)

Returns the precomputed logZ values for Eugenio's RBMs.
"""
function Eugenio_RBM_logZ_20230419(which::Symbol)
    if which === :I
        return 259.22567749023437500000
    elseif which === :II
        return 240.27294921875000000000
    elseif which === :IV
        return 323.63775634765625000000
    elseif which === :global
        return 285.41845703125000000000
    else
        throw(ArgumentError("argument must be one of :I, :II, :IV, or :global"))
    end
end

function Eugenio_RBM_20230419(which::Symbol)
    if which ∈ (:I, :II, :IV)
        return Eugenio_RBM_specific_20230424(which)
    elseif which === :global
        return Eugenio_RBM_20230419()
    else
        throw(ArgumentError("argument must be one of :I, :II, :IV, or :global"))
    end
end

function Eugenio_RBM_20230419_loglikelihood(which::Symbol, seqs::AbstractArray)
    Eeff = free_energy(Eugenio_RBM_20230419(which), seqs)
    logZ = Eugenio_RBM_logZ_20230419(which)
    return -Eeff .- logZ
end

function Eugenio_RBM_20230419_probed_sequences_free_energies_eval_from_python(which::Symbol)
    if which === :I
        return vec(readdlm(joinpath(Probed_Sequences_Eeff_from_Python_20240703_path(), "Eeff_I.txt")))
    elseif which === :II
        return vec(readdlm(joinpath(Probed_Sequences_Eeff_from_Python_20240703_path(), "Eeff_II.txt")))
    elseif which === :IV
        return vec(readdlm(joinpath(Probed_Sequences_Eeff_from_Python_20240703_path(), "Eeff_IV.txt")))
    elseif which === :global
        return vec(readdlm(joinpath(Probed_Sequences_Eeff_from_Python_20240703_path(), "Eeff_glob.txt")))
    else
        throw(ArgumentError("argument must be one of :I, :II, :IV, or :global"))
    end
end

"""
    Eugenio_MSA_Predicted_Classes_20240715()

Predicted classes (I, II, or IV) of MSA sequences, according to inputs to two hidden units
of the "global" RBM.
"""
function Eugenio_MSA_Predicted_Classes_20240715()
    msa = Eugenio_MSA_20230419()
    predicted_class = Eugenio_Predict_Class_from_Inputs_20240715(onehot(msa))
    return (; msa, predicted_class)
end

function Eugenio_Predict_Class_from_Inputs_20240715(sequences::AbstractArray{<:Real,3})
    rbm = Eugenio_RBM_20230419(:global)
    Ih = inputs_h_from_v(rbm, sequences)

    # Thresholds are taken from Eugenio's PhD thesis (see Fig. 7.6).
    predicted_class = map(-Ih[38, :], Ih[36, :]) do x, y
        if x < -1
            return :I
        elseif x > -1 && y > -3
            return :II
        elseif x > -1 && y < -3
            return :IV
        else
            return :unknown
        end
    end

    return predicted_class
end
