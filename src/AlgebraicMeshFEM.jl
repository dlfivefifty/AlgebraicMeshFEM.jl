module AlgebraicMeshFEM
using AlgebraicMeshes, ArrayLayouts, ClassicalOrthogonalPolynomials, ContinuumArrays, StaticArrays
using ContinuumArrays: Basis
import Base: getindex, axes, first, last, size
export AlgebraicMeshVector, AlgebraicMeshAxis, ElementIndex

"""
    ElementIndex(dim, elnum, basisind)


"""
struct ElementIndex{Ind}
    dim::Int
    elindex::Int
    basisind::Ind
end

"""
   AlgebraicMeshAxis

gives the map between a mesh and the axes of the correspond
coefficient vector.
"""
struct AlgebraicMeshAxis{M<:AlgebraicMesh, Ax<:Tuple, LENG<:Integer} <: AbstractUnitRange{Int}
    mesh::M
    axes::Ax
    length::LENG
end

AlgebraicMeshAxis(mesh::AlgebraicMesh, ax::Tuple) = AlgebraicMeshAxis(mesh, ax, mapreduce(d -> mapreduce(length,+,d),+,ax))

first(a::AlgebraicMeshAxis) = 1
last(a::AlgebraicMeshAxis) = a.length
Base.unitrange(a::AlgebraicMeshAxis) = 1:length(a)


"""
    findelementindex(meshaxis, k)

converts an integer index `k` to the corresp[onding `ElementIndex`
"""
function findelementindex(meshaxis::AlgebraicMeshAxis, k::Int)
    
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



struct AlgebraicMeshVector{T, Ax<:AlgebraicMeshAxis, D<:Tuple} <: LayoutVector{T}
    axis::Ax
    data::D # a tuple of vectors with same structure as mesh, containing vectors/matrices.
            # we interlace these to produce a single vector/matrix.

    function AlgebraicMeshVector{T, Ax, D}(ax, data) where {T, Ax<:AlgebraicMeshAxis, D<:Tuple}    
        # check mesh and data sizes match
        @assert length(ax.mesh.complex) == length(data)
        for (m,v) in zip(ax.mesh.complex,data)
            @assert length(m) == length(v)
        end
        new{T,Ax,D}(ax, data)
    end
end

AlgebraicMeshVector(ax::AlgebraicMeshAxis, data) = AlgebraicMeshVector{Float64, typeof(ax), typeof(data)}(ax, data)
AlgebraicMeshVector(mesh::AlgebraicMesh, data) = AlgebraicMeshVector(AlgebraicMeshAxis(mesh, map(d -> map(e -> axes(e,1), d), data)), data)

axes(a::AlgebraicMeshVector) = (a.axis,)
size(a::AlgebraicMeshVector) = (length(a.axis),)

getindex(a::AlgebraicMeshVector, k::Int) = a[findelementindex(axes(a,1),k)]

# struct AlgebraicMeshMatrix{T, Ax, D::Tuple} <: LayoutVector{T}
#     axes::M
#     data::D # a tuple of vectors with same structure as mesh, containing vectors/matrices.
#             # we interlace these to produce a single vector/matrix.

#     function AlgebraicMeshArray{T, M, D}(mesh, data) where {T, N, M::AlgebraicMesh, D<:Tuple}    
#         # check mesh and data sizes match
#         @assert length(mesh.complex) == length(data)
#         for (m,v) in zip(mesh.complex,data)
#             @assert length(m) == length(v)
#         end
#         new{T,M,D}(mesh, data)
#     end
# end


end # module AlgebraicMeshFEM
