module RecurrentNetworkAnalysis
using LinearAlgebra

"""
    get_rotational_modes(Wrr::AbstractMatrix)

Get the complex conjugate pairs of eigenvectors  that make up the rotational modes of the matrix `Wrr`
"""
function get_rotational_modes(Wrr::AbstractMatrix)
    ee = eigen(Wrr)
    get_rotational_modes(ee.values)
end

function get_rotational_modes(λ::AbstractVector{T}) where T <: Complex
    ij = [(i,j) for i in 1:length(λ), j in 1:length(λ)]
    pp = filter(k->(λ[k[1]]==conj(λ[k[2]])&&(k[1]>k[2])&&(imag(λ[k[1]])!=0)), ij)
    pp
end

"""
    non_normality(W::AbstractMatrix{T}) where T <: Real

Henrici's index of matrix non-normality. A value of 0 indicates a normal matrix, while a value of 1 indicate
a maximally non-normal matrix.

A normal matrix can be described by an orthonormal eigenvalue basis, while for a non-normal matrix, the eigvectors
do no necessarily form an orthonormal basis.
"""
function non_normality(W::AbstractMatrix{T}) where T <: Real
    nw = sum(abs2, W) 
    ee = eigen(W)
    ns = sum(abs2, ee.values)
    sqrt(nw-ns)./sqrt(nw)
end

end # module RecurrentNetworkAnalysis
