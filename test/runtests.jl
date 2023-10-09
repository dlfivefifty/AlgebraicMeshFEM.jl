using AlgebraicMeshFEM, AlgebraicMeshes, DomainSets, Test

mesh = interior(AlgebraicMesh([(0..1) × (1..2) , (0..1) × (0..1), (1..2) × (0..1)]))


using CairoMakie

AlgebraicMeshVector(mesh, 

