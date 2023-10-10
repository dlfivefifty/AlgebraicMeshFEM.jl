using AlgebraicMeshFEM, AlgebraicMeshes, DomainSets, StaticArrays, InfiniteArrays, ClassicalOrthogonalPolynomials, Test
using ContinuumArrays: affine
using MultivariateOrthogonalPolynomials: DiagTrav

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
    W = Weighted(Jacobi(1,1))

    @test P[SVector(0.1,0), 1] ≈ 0.1
    @test P[SVector(1.1,0), 1] ≈ 0.9
    @test P[SVector(0.1,0), 2] ≈ W[2*0.1-1,1]
    @test P[SVector(1.1,0), 2] == 0
    @test P[SVector(0.1,0), 3] == 0
    @test P[SVector(1.1,0), 3] ≈ W[2*0.1-1,1]

    f = P * c
    @test f[SVector(0.1,0)] ≈ 0.1 + 2*W[2*0.1-1,1]
    @test f[SVector(1.1,0)] ≈ 0.9 + 3*W[2*0.1-1,1]
end

@testset "2D" begin
    mesh = interior(AlgebraicMesh([(0..1) × (1..2) , (0..1) × (0..1), (1..2) × (0..1)]))
    P = AlgebraicMeshPolynomial{0}(mesh)
    @test P[SVector(0.1,1.2),1:6] == [0,0,1,0,0,0]
    @test P[SVector(0.1,0.2),1:6] == [0,0,0,1,0,0]
    @test P[SVector(1.1,0.2),1:6] == [0,0,0,0,1,0]

    A = zeros(∞,∞); A[1,1] = 4;
    B = zeros(∞,∞); B[1,1] = 5;
    C = zeros(∞,∞); C[1,1] = 6;

    data = (Int[], [[2; zeros(∞)], [3;zeros(∞)]], [DiagTrav(A),DiagTrav(B),DiagTrav(C)]);
    c = AlgebraicMeshVector(mesh, data);
    f = P*c
    @test f[SVector(0.1,1.2)] ≈ 4
    @test f[SVector(0.1,0.2)] ≈ 5
    @test f[SVector(1.1,0.2)] ≈ 6
end


