module AlgebraicMeshFEM
using AlgebraicMeshes, ClassicalOrthogonalPolynomials, ContinuumArrays
using ContinuumArrays: Basis, LayoutVector
import Base: getindex, axes

"""
    ElementIndex(dim, elnum, basisind)


"""
struct ElementIndex{N,Ind}
    dim::Int
    elindex::Int
    basisind::Ind
end

"""
   AlgebraicMeshAxis

gives the map between a mesh and the axes of the correspond
coefficient vector.
"""
struct AlgebraicMeshAxis{λ, M<:Mesh, Sz<:Tuple} <: AbstractUnitRange{Int}
    mesh::M
    sizes::Sz
end

"""
    findelementindex(meshaxis, k)

converts an integer index `k` to the corresp[onding `ElementIndex`
"""
function findelementindex(meshaxis, k::Int)

end

struct AlgebraicMeshPolynomial{λ, T, M<:Mesh, Sz<:Tuple} <: Basis{T}
    mesh::M
    sizes::Sz
end


axes(P::AlgebraicMeshPolynomial{λ}) = (Inclusion(P.mesh), AlgebraicMeshAxis{λ}(P.mesh, P.sizes))

getindex(P::AlgebraicMeshPolynomial, 𝐱::SVector, k::Int) = P[𝐱, findelementindex(axes(P,2), k)]
function getindex(P::AlgebraicMeshPolynomial{0,T}, 𝐱::SVector, K::ElementIndex) where T
    # TODO: check axes
    K.dims == 3 || return zero(T) # only elements have a basis attached
    el = elements(P.mesh)[K.elindex]
    𝐱 ∈ el || return zero(T) # zero outside element
    legendre(el)[𝐱, K.basisind]
end



struct AlgebraicMeshCoefficients{λ, T, M<:Mesh, D<:Dict} <: LayoutVector{T}
    mesh::M
    data::D
end

end # module AlgebraicMeshFEM
