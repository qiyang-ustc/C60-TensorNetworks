using Test
include(".\\WangLandau.jl")
# for i in 1:1:60
#     for j in 1:1:3
#         @show i,j,node_index(node_neighbour(index_node(i),j))
#     end
# end
@testset "Graph.jl" begin
node=Node(0,1)

next_node=node_neighbour(node,1)
@test next_node==Node(4,1)

next_node=node_neighbour(node,2)
@test next_node==Node(1,1)

next_node=node_neighbour(node,3)
@test next_node==Node(0,2)

node=Node(2,6)

next_node=node_neighbour(node,3)
@test next_node==Node(2,7)

next_node=node_neighbour(node,2)
@test next_node==Node(2,5)

@test node_index(Node(0,12))==12
@test node_index(Node(2,5))==29
for i in 0:4
    for j in 1:12
        @test node_index(Node(i,j))==12*i+j
        @test node_index(index_node(12*i+j))==12*i+j
    end
end
@test Int(rand(UInt)%3)+1 < 4
@test Int(rand(UInt)%3)+1 > 0
end;

@testset "Energy" begin
test_spins = ones(Int,60)
@test cal_energy(test_spins)==90
for i in 1:1:100
    Random.seed!(123456)
    test_spins = rand([-1,1],60)
    for j in 1:1:60
        test_energy = cal_energy(test_spins)
        test_new_energy = test_energy + delta_energy(j,test_spins)
        test_spins[j] = -test_spins[j]
        @test test_new_energy == cal_energy(test_spins)
    end
end
end