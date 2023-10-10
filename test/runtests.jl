using AlgebraicMeshFEM, AlgebraicMeshes, DomainSets, StaticArrays, InfiniteArrays, Test
using ContinuumArrays: affine

@testset "line segments" begin
    ℓ = LineSegment(SVector(0,0), SVector(1,0))
    a = affine(ℓ, ChebyshevInterval())
    a[SVector(0.1,0)]
end

@testset "1D" begin
    mesh = interior(AlgebraicMesh([LineSegment(SVector(0,0), SVector(1,0)),LineSegment(SVector(1,0), SVector(2,0))]))
    P = AlgebraicMeshPolynomial{0}(mesh)
    P[SVector(0.1,0), 2]

    data = ([1], [[2; zeros(∞)], [2;zeros(∞)]])

    c = AlgebraicMeshVector(mesh, data)
    P[SVector(0.1,1),1]
end

mesh = interior(AlgebraicMesh([(0..1) × (1..2) , (0..1) × (0..1), (1..2) × (0..1)]))


using CairoMakie



