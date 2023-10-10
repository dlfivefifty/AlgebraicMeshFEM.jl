module AlgebraicMeshFEM
using AlgebraicMeshes, ArrayLayouts, ClassicalOrthogonalPolynomials, ContinuumArrays, StaticArrays, LinearAlgebra, BlockArrays, MultivariateOrthogonalPolynomials
import ContinuumArrays: Basis, AffineMap, AbstractQuasiVector, AbstractAffineQuasiVector, affine_getindex, measure, affinemap_A, affinemap_b, PaddedLayout
import Base: getindex, axes, first, last, size, oneto
import ArrayLayouts: colsupport, MemoryLayout
import ClassicalOrthogonalPolynomials: legendre
import DomainSets: UnitInterval, Rectangle, leftendpoint, rightendpoint
import BlockArrays: blockcolsupport
export AlgebraicMeshVector, AlgebraicMeshAxis, AlgebraicMeshPolynomial, ElementIndex, findelementindex


###
# Classical OPs on segments
###

measure(ℓ::Inclusion{<:Any,<:LineSegment}) = norm(ℓ.domain.b - ℓ.domain.a)
function affinemap_A(m::AffineMap{<:Number,<:Inclusion{<:SVector}})
    domain, range = getfield(m, :domain), getfield(m, :range)
    ba = last(domain)-first(domain)
    c,d = first(range),last(range)
    (d-c)/norm(ba)^2 * ba'
end
function affinemap_b(m::AffineMap{<:Number,<:Inclusion{<:SVector}})
    domain, range = getfield(m, :domain), getfield(m, :range)
    a,b = first(domain),last(domain)
    c,d = first(range),last(range)
    ba = b-a
    (c-d)/norm(ba)^2 * ba'a + c
end
legendre(ℓ::LineSegment{d,T}) where {d,T} = Legendre{float(T)}()[affine(ℓ,ChebyshevInterval{T}()), :]
AffineMap(domain::AbstractQuasiVector{T}, range::AbstractQuasiVector{V}) where {d,T<:Number,V<:SVector{d}} = AffineMap{promote_type(SVector{d,T},V), typeof(domain),typeof(range)}(domain,range)
AffineMap(domain::AbstractQuasiVector{T}, range::AbstractQuasiVector{V}) where {T<:SVector,V} = AffineMap{promote_type(eltype(T),V), typeof(domain),typeof(range)}(domain,range)
getindex(A::AbstractAffineQuasiVector{<:Any,<:Any,<:Inclusion{<:SVector{d}}}, k::SVector{d}) where d = affine_getindex(A, k)

function legendre(r::Rectangle)
    (a,c) = leftendpoint(r)
    (b,d) = rightendpoint(r)
    RectPolynomial(legendre(a..b), legendre(c..d))
end

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
struct AlgebraicMeshAxis{N, Ax<:NTuple{N,Any}, LENG<:Integer} <: AbstractUnitRange{Int}
    axes::Ax
    length::LENG
end

AlgebraicMeshAxis(ax::Tuple) = AlgebraicMeshAxis(ax, mapreduce(d -> isempty(d) ? 0 : mapreduce(length,+,d),+,ax))

first(a::AlgebraicMeshAxis) = 1
last(a::AlgebraicMeshAxis) = a.length
Base.unitrange(a::AlgebraicMeshAxis) = 1:length(a)


"""
    findelementindex(meshaxis, k)

converts an integer index `k` to the corresp[onding `ElementIndex`
"""
function findelementindex(meshaxis::AlgebraicMeshAxis{2}, k::Int)
    # put all vertices first
    n_v,n_e = map(length,meshaxis.axes)
    k ≤ n_v && return ElementIndex(1, k, 1)
    k -= n_v
    ind,el = divrem(k-1, n_e)
    ElementIndex(2, el+1, ind+1)
end

function findelementindex(meshaxis::AlgebraicMeshAxis{3}, k::Int)
    # put all vertices first
    n_v,n_f,n_e = map(length,meshaxis.axes)
    k ≤ n_v && return ElementIndex(1, k, 1)
    k -= n_v
    # block size is n_f + n_e*K
    # so number of inds up to block N is n_f*N + n_e*N*(1+N)/2
    # solving for N we get
    K = (-n_e - 2n_f + isqrt(8*(k-1)n_e + (n_e+2n_f)^2)) ÷ (2n_e) + 1
    # subtract previous blocks
    prev_inds = ((K-1)*K)÷2
    k -= n_f*(K-1) + n_e*prev_inds

    k ≤ n_f && return ElementIndex(2, k, K)
    k -= n_f


    el,ind = divrem(k-1, K)
    # want to do Block(K)[ind+1] but for now we keep type instability
    ElementIndex(3, el+1, prev_inds + ind+1)
end

struct AlgebraicMeshPolynomial{λ, T, M<:AlgebraicMesh, Ax} <: Basis{T}
    mesh::M
    caxis::Ax
end

_mesh2axes() = ()
_mesh2axes(::SVector) = oneto(1)
_mesh2axes(::LineSegment) = oneto(∞)
_mesh2axes(::Rectangle) = blockedrange(oneto(∞))
_mesh2axes(a::Vector, b...) = (map(_mesh2axes,a), _mesh2axes(b...)...)


AlgebraicMeshPolynomial{λ}(mesh, axes) where λ = AlgebraicMeshPolynomial{λ,Float64,typeof(mesh),typeof(axes)}(mesh, axes)
AlgebraicMeshPolynomial{λ}(mesh) where λ = AlgebraicMeshPolynomial{λ}(mesh, AlgebraicMeshAxis(_mesh2axes(mesh.complex...)))



axes(P::AlgebraicMeshPolynomial{λ}) where λ = (Inclusion(P.mesh), P.caxis)

getindex(P::AlgebraicMeshPolynomial, 𝐱::SVector, k::Int) = P[𝐱, findelementindex(axes(P,2), k)]
function getindex(P::AlgebraicMeshPolynomial{0,T}, 𝐱::SVector, K::ElementIndex) where T
    # TODO: check axes
    el = P.mesh.complex[K.dim][K.elindex]
    𝐱 ∈ el || return zero(T) # zero outside element
    convert(T,legendre(el)[𝐱, K.basisind])::T
end

function vertexmode(vertex, edge::LineSegment, 𝐱)
    if vertex == edge.a
        1-affine(edge, UnitInterval())[𝐱]
    else
        @assert vertex == edge.b
        affine(edge, UnitInterval())[𝐱]
    end
end

bubble(edge::LineSegment) = Weighted(Jacobi(1,1))[affine(edge, ChebyshevInterval()), :]

function getindex(P::AlgebraicMeshPolynomial{1,T,<:Any,<:AlgebraicMeshAxis{2}}, 𝐱::SVector, K::ElementIndex) where T
    # TODO: check axes
    el = P.mesh.complex[K.dim][K.elindex]
    if K.dim == 1 # vertices
        ne = neighborhood(P.mesh, el)
        𝐱 ∈ ne || return zero(T) # zero outside element
        ed = edges(ne)
        j = findfirst(e -> 𝐱 ∈ e, ed)
        vertexmode(el, ed[j], 𝐱)
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

# MemoryLayout(::Type{<:AlgebraicMeshVector}) = PaddedLayout{UnknownLayout}()

axes(a::AlgebraicMeshVector) = (a.axis,)
size(a::AlgebraicMeshVector) = (length(a.axis),)

getindex(a::AlgebraicMeshVector{T}, K::ElementIndex) where T = convert(T,a.data[K.dim][K.elindex][K.basisind])::T
getindex(a::AlgebraicMeshVector, k::Int) = a[findelementindex(axes(a,1),k)]

function colsupport(a::AlgebraicMeshVector{<:Any,<:AlgebraicMeshAxis{2}}, j)
    n_v,n_e = map(length,a.axis.axes)
    r = n_v
    for k = 1:n_e
        N = last(colsupport(a.data[2][k]))
        r = max(r, (N-1) * n_e + k + n_v)
    end
    oneto(r)
end

function blockcolsupport(a::AlgebraicMeshVector{<:Any,<:AlgebraicMeshAxis{3}}, j)
    n_v,n_f,n_e = map(length,a.axis.axes)
    r = 1
    for k = 1:n_f
        N = last(colsupport(a.data[2][k]))
        r = max(r,last(colsupport(a.data[2][k])))
    end
    for k = 1:n_e
        r = max(r,Int(last(blockcolsupport(a.data[3][k]))))
    end

    Block.(oneto(r))
end

function colsupport(a::AlgebraicMeshVector{<:Any,<:AlgebraicMeshAxis{3}}, j)
    n_v,n_f,n_e = map(length,a.axis.axes)
    N = Int(last(blockcolsupport(a, j)))
    oneto(n_v + N*n_f + (N*(N+1)÷2)*n_e)
end

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
