using AlgebraicMeshFEM, AlgebraicMeshes, DomainSets, StaticArrays, InfiniteArrays, ClassicalOrthogonalPolynomials, Test
using ContinuumArrays: affine

@testset "line segments" begin
    ℓ = LineSegment(SVector(0,0), SVector(1,0))
    a = affine(ℓ, ChebyshevInterval())
    @test a[SVector(0.1,0)] ≈ -0.8

    ℓ = LineSegment(SVector(1,2), SVector(2,4))
    a = affine(ℓ, ChebyshevInterval())
    @test a[SVector(1.1, 2.2)] ≈ -0.8
end

@testset "1D" begin
    mesh = interior(AlgebraicMesh([LineSegment(SVector(0,0), SVector(1,0)),LineSegment(SVector(1,0), SVector(2,0))]))
    P = AlgebraicMeshPolynomial{0}(mesh)
    @test P[SVector(0.1,0), 2] ≈ 1
    @test P[SVector(1.1,0), 2] ≈ 0
    @test P[SVector(0.1,0), 3] ≈ 0
    @test P[SVector(1.1,0), 3] ≈ 1

    data = ([1], [[2; zeros(∞)], [3;zeros(∞)]]);

    c = AlgebraicMeshVector(mesh, data);
    f = P*c
    @test f[SVector(0.1,0)] ≈ 2
    @test f[SVector(1.1,0)] ≈ 3

    P = AlgebraicMeshPolynomial{1}(mesh)
    @test P[SVector(0.1,0), 1] ≈ 0.1
    @test P[SVector(1.1,0), 1] ≈ 0.9
    @test P[SVector(0.1,0), 2] ≈ 0.1

    f = P * c
    W = Weighted(Jacobi(1,1))
    @test f[SVector(0.1,0)] ≈ 0.1 + 2*W[2*0.1-1,1]
    @test f[SVector(1.1,0)] ≈ 0.9 + 3*W[2*0.1-1,1]
end

@testset "2D" begin
    mesh = interior(AlgebraicMesh([(0..1) × (1..2) , (0..1) × (0..1), (1..2) × (0..1)]))
    P = AlgebraicMeshPolynomial{0}(mesh)    
end


