using Random
#point:60  bond:90
#5-5-10-5-
const C60_graph=[
[[-1,1],[1,1],[0,2]],
[[0,1],[0,4],[0,3]],
[[0,2],[0,5],[-1,4]],
[[0,2],[+1,3],[0,6]],
[[0,3],[-1,8],[0,6]],
[[0,4],[0,5],[0,7]],
[[0,6],[0,8],[0,9]],
[[0,7],[1,5],[0,10]],
[[0,7],[0,11],[-1,10]],
[[0,8],[0,11],[+1,9]],
[[0,9],[0,10],[0,12]],
[[-1,12],[0,11],[+1,12]]
]

struct Node
    tree::Int
    index::Int
end
function node_index(node::Node)
    return node.tree*12+node.index  #(0-5)+(1-12)
end
function node_neighbour(node::Node,direction::Int)
    tree,index = C60_graph[node.index][direction]
    new_index = index
    new_tree = (5+tree+node.tree) %5
    return Node(new_tree,new_index)
end
function C60(beta::Float64,steps::Int)
    cosh_beta = cosh(beta)
    tanh_beta = tanh(beta)
    bond_state = zeros(Int,90)
    node = Node(0,1) 
    now_index = node_index(node)
    edges = zeros(Bool,60,60)
    lnZ = 90*(log(cosh_beta))
    number_of_measurement = 0.0
    sum_of_tanb = 0.0
    for i in 1:1:steps
        direction = Int(rand(UInt)%3+1)
        next_node = node_neighbour(node,direction)
        next_index= node_index(next_node)
        if edges[max(now_index,next_index),min(now_index,next_index)]
            edges[max(now_index,next_index),min(now_index,next_index)] = false
            node = next_node
        else    
            if rand()<tanh_beta
                edges[max(now_index,next_index),min(now_index,next_index)] = true
                node = next_node
            end
        end
        # @show node
        if node_index(node)==1
            sum_of_tanb += (-1)^(sum(edges)%2)
            number_of_measurement += tanh_beta^(sum(edges))
        end
    end
    lnZ += log(sum_of_tanb/number_of_measurement)
    return lnZ/60+log(2.0)
end
function C60_RandomLoop(beta::Float64,steps::Int)
    cosh_beta = cosh(beta)
    tanh_beta = tanh(beta)
    bond_state = zeros(Int,90)
    node = Node(0,1) 
    now_index = node_index(node)
    edges = zeros(Bool,60,60)
    lnZ = 90*(log(cosh_beta))
    number_of_measurement = 0.0
    sum_of_tanb = 0.0
    for i in 1:1:steps
        direction = Int(rand(UInt)%3+1)
        next_node = node_neighbour(node,direction)
        next_index= node_index(next_node)
        if edges[max(now_index,next_index),min(now_index,next_index)]
            edges[max(now_index,next_index),min(now_index,next_index)] = false
            node = next_node
        else    
            # if rand()<tanh_beta
                edges[max(now_index,next_index),min(now_index,next_index)] = true
                node = next_node
            # end
        end
        # @show node
        if node_index(node)==1
            sum_of_tanb += (-tanh_beta)^(sum(edges))
            number_of_measurement+=1
        end
    end
    # @show lnZ
    # @show number_of_measurement
    lnZ += log(2^32*sum_of_tanb/number_of_measurement)
    # @show lnZ
    if abs(tanh_beta)<1E-10
        return log(2.0)
    end
    return lnZ/60+log(2.0)
end
print(C60(0.05,100))
# print(C60_RandomLoop(1.0,10000000))

