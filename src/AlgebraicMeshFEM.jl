module AlgebraicMeshFEM
using AlgebraicMeshes, ArrayLayouts, ClassicalOrthogonalPolynomials, ContinuumArrays, StaticArrays, LinearAlgebra
import ContinuumArrays: Basis, AffineMap, AbstractQuasiVector, AbstractAffineQuasiVector, affine_getindex, measure
import Base: getindex, axes, first, last, size, oneto
import ClassicalOrthogonalPolynomials: legendre
export AlgebraicMeshVector, AlgebraicMeshAxis, AlgebraicMeshPolynomial, ElementIndex, findelementindex


###
# Classical OPs on segments
###

measure(ℓ::Inclusion{<:Any,<:LineSegment}) = norm(ℓ.domain.b - ℓ.domain.a)
legendre(ℓ::LineSegment{d,T}) where {d,T} = Legendre{float(T)}()[affine(ℓ,ChebyshevInterval{T}()), :]
AffineMap(domain::AbstractQuasiVector{T}, range::AbstractQuasiVector{V}) where {d,T<:Number,V<:SVector{d}} = AffineMap{promote_type(SVector{d,T},V), typeof(domain),typeof(range)}(domain,range)
AffineMap(domain::AbstractQuasiVector{T}, range::AbstractQuasiVector{V}) where {T<:SVector,V} = AffineMap{promote_type(eltype(T),V), typeof(domain),typeof(range)}(domain,range)
getindex(A::AbstractAffineQuasiVector{<:Any,<:Any,<:Inclusion{<:SVector{d}}}, k::SVector{d}) where d = affine_getindex(A, k)

###
# Element-based arrays
###
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
struct AlgebraicMeshAxis{Ax<:Tuple, LENG<:Integer} <: AbstractUnitRange{Int}
    axes::Ax
    length::LENG
end

AlgebraicMeshAxis(ax::Tuple) = AlgebraicMeshAxis(ax, mapreduce(d -> mapreduce(length,+,d),+,ax))

first(a::AlgebraicMeshAxis) = 1
last(a::AlgebraicMeshAxis) = a.length
Base.unitrange(a::AlgebraicMeshAxis) = 1:length(a)


"""
    findelementindex(meshaxis, k)

converts an integer index `k` to the corresp[onding `ElementIndex`
"""
function findelementindex(meshaxis::AlgebraicMeshAxis{<:NTuple{2,Any}}, k::Int)
    # put all vertices first
    n_v,n_e = map(length,meshaxis.axes)
    k ≤ n_v && return ElementIndex(1, k, 1)
    k -= n_v
    ind,el = divrem(k-1, n_e)
    ElementIndex(2, el+1, ind+1)
end

struct AlgebraicMeshPolynomial{λ, T, M<:AlgebraicMesh, Ax<:Tuple} <: Basis{T}
    mesh::M
    axes::Ax
end

_mesh2axes() = ()
_mesh2axes(::SVector) = oneto(1)
_mesh2axes(::LineSegment) = oneto(∞)
_mesh2axes(a::Vector, b...) = (map(_mesh2axes,a), _mesh2axes(b...)...)


AlgebraicMeshPolynomial{λ}(mesh, axes) where λ = AlgebraicMeshPolynomial{λ,Float64,typeof(mesh),typeof(axes)}(mesh, axes)
AlgebraicMeshPolynomial{λ}(mesh) where λ = AlgebraicMeshPolynomial{λ}(mesh,_mesh2axes(mesh.complex...))



axes(P::AlgebraicMeshPolynomial{λ}) where λ = (Inclusion(P.mesh), AlgebraicMeshAxis(P.axes))

getindex(P::AlgebraicMeshPolynomial, 𝐱::SVector, k::Int) = P[𝐱, findelementindex(axes(P,2), k)]
function getindex(P::AlgebraicMeshPolynomial{0,T}, 𝐱::SVector, K::ElementIndex) where T
    # TODO: check axes
    el = P.mesh.complex[K.dim][K.elindex]
    𝐱 ∈ el || return zero(T) # zero outside element
    convert(T,legendre(el)[𝐱, K.basisind])::T
end

function getindex(P::AlgebraicMeshPolynomial{1,T,<:Any,<:NTuple{2,Any}}, 𝐱::SVector, K::ElementIndex) where T
    # TODO: check axes
    el = P.mesh.complex[K.dim][K.elindex]
    if K.dim == 1 # vertices
        ne = neighborhood(P.mesh, el)
        𝐱 ∈ ne || return zero(T) # zero outside element
        ed = edges(ne)
        j = findfirst(e -> 𝐱 ∈ e, ed)
        vertexmode(ed[j], el, 𝐱)
    else # K.dim == 2 # elements
        𝐱 ∈ el || return zero(T) # zero outside element
        convert(T,bubble(el)[𝐱, K.basisind])::T
    end
end



struct AlgebraicMeshVector{T, Ax<:AlgebraicMeshAxis, D<:Tuple} <: LayoutVector{T}
    axis::Ax
    data::D # a tuple of vectors with same structure as mesh, containing vectors/matrices.
            # we interlace these to produce a single vector/matrix.

    function AlgebraicMeshVector{T, Ax, D}(ax, data) where {T, Ax<:AlgebraicMeshAxis, D<:Tuple}    
        # check mesh and data sizes match
        @assert length(ax.axes) == length(data)
        for (m,v) in zip(ax.axes,data)
            @assert length(m) == length(v)
        end
        new{T,Ax,D}(ax, data)
    end
end

AlgebraicMeshVector(ax::AlgebraicMeshAxis, data) = AlgebraicMeshVector{Float64, typeof(ax), typeof(data)}(ax, data)
AlgebraicMeshVector(mesh::AlgebraicMesh, data) = AlgebraicMeshVector(AlgebraicMeshAxis(map(d -> map(e -> axes(e,1), d), data)), data)

axes(a::AlgebraicMeshVector) = (a.axis,)
size(a::AlgebraicMeshVector) = (length(a.axis),)

getindex(a::AlgebraicMeshVector{T}, K::ElementIndex) where T = convert(T,a.data[K.dim][K.elindex][K.basisind])::T
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
