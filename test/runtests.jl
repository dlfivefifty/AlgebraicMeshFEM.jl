using AlgebraicMeshFEM, AlgebraicMeshes, DomainSets, StaticArrays, InfiniteArrays, Test


@testset "1D" begin
    mesh = interior(AlgebraicMesh([LineSegment(SVector(0,0), SVector(1,0)),LineSegment(SVector(1,0), SVector(2,0))]))
    data = ([1], [[2; zeros(∞)], [2;zeros(∞)]])
    
    AlgebraicMeshAxis(mesh, map(d -> map(e -> axes(e,1), d), data))

    c = AlgebraicMeshVector(mesh, )
end

mesh = interior(AlgebraicMesh([(0..1) × (1..2) , (0..1) × (0..1), (1..2) × (0..1)]))


using CairoMakie



