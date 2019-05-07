using Random
include(".\\Spectrum.jl")
#point:60  bond:90
#5-5-10-5-.....
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
    return node.tree*12+node.index  #(0-4)+(1-12)
end
function node_neighbour(node::Node,direction::Int)
    dtree,index = C60_graph[node.index][direction]
    new_index = index
    new_tree = (5+dtree+node.tree) %5
    return Node(new_tree,new_index)
end
function index_node(index::Int)
    return Node(div(index-1,12),(index-1)%12+1)
end

function cal_energy(spins::Array{Int,1})
    energy = 0
    for pick_index in 1:1:60
        now_node = index_node(pick_index)
        energy_delta = 0
        for inbr in 1:1:3
            nbr_node = node_neighbour(now_node,inbr)
            nbr_index = node_index(nbr_node)
            energy_delta += spins[nbr_index]
        end
        energy_delta *= spins[pick_index]
        energy += energy_delta
    end
    return div(energy,2)
end
function delta_energy(pick_index::Int,spins::Array{Int,1})
    now_node = index_node(pick_index)
    energy_delta = 0
    for inbr in 1:1:3
        nbr_node = node_neighbour(now_node,inbr)
        nbr_index = node_index(nbr_node)
        energy_delta += spins[nbr_index]
    end
    energy_delta *= spins[pick_index]
    energy_delta *= -2
    return energy_delta
end

function ground_g(g::Array{Float64,1})
    temp_ground_g=0
    min_g = minimum(g)
    for i in 1:length(g)
        if g[i] > min_g + 1
            temp_ground_g = g[i]
            break
        end
    end
    return temp_ground_g  
end 

function C60_WangLandau(factor::Float64,steps::Int128)
    node = Node(0,1) 
    spins = ones(Int,60)
    spectrum = Spectrum(200)
    set_ground!(spectrum,-100)
    energy = 90   #(As we set all spin to 1 so that the energy should be +90)
    for i in 1:1:steps
        pick_index = Int(rand(UInt)%(60))+1
        if i%1000==0
            normalize_factor(spectrum)
        end
        new_energy = energy + delta_energy(pick_index,spins)
        if rand()<exp(get_factor(spectrum,energy)-get_factor(spectrum,new_energy)) #accept the move 
            spins[pick_index] = - spins[pick_index] # flip spin
            add_factor!(spectrum,new_energy,factor)
            energy = new_energy
        else 
            add_factor!(spectrum,energy,factor)
        end
    end
    energy_level = collect(-99:1:100)
    min_s = minimum(spectrum.factor)
    for i in 1:length(spectrum.factor)
        if spectrum.factor[i] < min_s + 1
            spectrum.factor[i] = spectrum.factor[i] - 100000.0
        end
    end
    g = exp.(spectrum.factor)
    g = log.(g/sum(g)) #normalization
    beta = 1.0
    lnZ = log(2.0)+  log(sum(exp.(g-beta.*energy_level)))/60
    min_g = minimum(g) 
    @show g
    @show ground_g(g)
    print("Degeneracy of Gound State is ",2^(60.0+ground_g(g)/log(2)))

    # for i in collect(1:1:50)
    #     print(i,' ',log(2.0)+  sum(g.*exp.(-beta.*energy_level))/60,"\n")
    # end
    return lnZ
end
#should be tuned carefully
print('\n',"lnZ=",C60_WangLandau(0.000001,Int128(10000000000)))

