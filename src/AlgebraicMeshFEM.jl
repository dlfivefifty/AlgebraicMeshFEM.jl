module AlgebraicMeshFEM
using AlgebraicMeshes, ArrayLayouts, ClassicalOrthogonalPolynomials, ContinuumArrays, StaticArrays
using ContinuumArrays: Basis
import Base: getindex, axes
export AlgebraicMeshVector

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
struct AlgebraicMeshAxis{λ, M<:AlgebraicMesh, Sz<:Tuple} <: AbstractUnitRange{Int}
    mesh::M
    sizes::Sz
end

"""
    findelementindex(meshaxis, k)

converts an integer index `k` to the corresp[onding `ElementIndex`
"""
function findelementindex(meshaxis, k::Int)

end

struct AlgebraicMeshPolynomial{λ, T, M<:AlgebraicMesh, Sz<:Tuple} <: Basis{T}
    mesh::M
    sizes::Sz
end


axes(P::AlgebraicMeshPolynomial{λ}) where λ = (Inclusion(P.mesh), AlgebraicMeshAxis{λ}(P.mesh, P.sizes))

getindex(P::AlgebraicMeshPolynomial, 𝐱::SVector, k::Int) = P[𝐱, findelementindex(axes(P,2), k)]
function getindex(P::AlgebraicMeshPolynomial{0,T}, 𝐱::SVector, K::ElementIndex) where T
    # TODO: check axes
    K.dims == 3 || return zero(T) # only elements have a basis attached
    el = elements(P.mesh)[K.elindex]
    𝐱 ∈ el || return zero(T) # zero outside element
    legendre(el)[𝐱, K.basisind]
end



struct AlgebraicMeshArray{T, N, M<:NTuple{N,AlgebraicMesh}, D<:NTuple{N,Tuple}} <: LayoutArray{T,N}
    mesh::M
    data::D # a tuple of vectors with same structure as mesh, containing vectors/matrices.
            # we interlace these to produce a single vector/matrix.

    function AlgebraicMeshArray{T, N, M, D}(meshes, data) where {T, N, M<:NTuple{N,AlgebraicMesh}, D<:Tuple}    
        # check mesh and data sizes match
        @assert length(meshes) == length(data)
        for (mesh,d) in zip(meshes, data)
            @assert length(mesh.complex) == length(data)
            for (m,v) in zip(mesh.complex,d)
                @assert length(m) == length(v)
            end
        end
        new{T,N,M,D}(mesh, data)
    end
end

const AlgebraicMeshVector{T, M<:Tuple{AlgebraicMesh}, D<:Tuple{Tuple}} = AlgebraicMeshArray{T, 1, M, D}
const AlgebraicMeshMatrix{T, M<:NTuple{2,AlgebraicMesh}, D<:NTuple{2,Tuple}} = AlgebraicMeshArray{T, 2, M, D}

end # module AlgebraicMeshFEM
