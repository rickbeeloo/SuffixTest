

using Random 
Random.seed!(1234)
using BenchmarkTools

const LIBSAIS = "libsais.so.2"
const MASK = Int32(1<<30) 


struct Color
    origin::Vector{Int32}
    max_len::Vector{Int32}
end

function create_suffix_array(in_vector::Vector{Int32}, free_space::Int32)
    out_vector = zeros(Int32, length(in_vector) + free_space)
    d = Dict(Iterators.map(reverse,pairs(sort(unique(in_vector)))))
    in_vector = Int32.(get.(Ref(d), in_vector, 0))
    n = length(in_vector)
    k = maximum(in_vector) +1 
    @ccall LIBSAIS.libsais_int(in_vector::Ptr{Int32}, out_vector::Ptr{Int32}, n::Int32, k::Int32, free_space::Int32)::Int32
    return out_vector
end

flipnode(n::Int32) = n ⊻ MASK
isflipped(n::Int32) = ifelse(n & MASK != 0, true, false)
noflip(n::Int32) = n & (MASK - 1)
convert_node(n::Int32) = n < 0 ? flipnode(abs(n)) : noflip(n) 

function convert_nodes!(in_vector::Vector{Int32})
    for i in eachindex(in_vector)
        in_vector[i] = convert_node(in_vector[i])
    end
end

function concat_with_seperator(vectors::Vector{Vector{Int32}})
    # Probably a cleaner way to do this :) Like using map 
    # to flip the node ids in slices and copying the slices to the output
    total_size = sum(map(length, vectors))
    concat_arr = zeros(Int32, total_size + length(vectors)) 
    vect_id = -1 * length(vectors) # to have it decending for sorting
    for i in eachindex(vectors)
        convert_nodes!(vectors[i])
    end
    # Concat with seperator + store sign in most significant bit
    i = 1
    @inbounds for v in vectors
        for node_id in v 
            concat_arr[i] = node_id
            i +=1
        end 
        concat_arr[i] = vect_id 
        vect_id +=1
        i +=1
    end
    return concat_arr
end

function create_k_suffix_array(vectors::Vector{Vector{Int32}}, free_space::Int32)
    concat_array = concat_with_seperator(vectors)
    suffix_array = create_suffix_array(concat_array, free_space)
    suffix_array .+= 1
    return concat_array, suffix_array
end


function locate_insert_point(sa::Vector{Int32}, concat_arr::Vector{Int32}, ref::AbstractVector{Int32})
    low = 1 
    high = length(sa)
    suffix_arr_len = length(sa)
    @inbounds while low <= high 
        mid = low + ((high - low) >>> 0x01) 
        suffix_start = sa[mid]
        suffix = view(concat_arr, suffix_start:suffix_arr_len)
        if ref < suffix 
            high = mid - 1
        elseif ref > suffix
            low = mid + 1
        else 
            # then they should be equal, exact match
            return mid 
        end
    end
    return low
end

function matches_till(arr1::AbstractVector{Int32}, arr2::AbstractVector{Int32})
    smallest = min(length(arr1), length(arr2))
    @inbounds for i in 1:smallest
        if arr1[i] != arr2[i]
            return i - 1
        end 
    end
    return smallest
end

function get_suffix_match_location(sa::Vector{Int32}, concat_arr::Vector{Int32}, array_to_find::AbstractVector{Int32}, point_in_sa::Int)
    # Get the suffix start position 
    suffix_start = sa[point_in_sa]
    # Get a view of the whole suffix to iterate over
    suffix = view(concat_arr, suffix_start:length(concat_arr))
    # Check for how long the suffix is the same as the array_to_find
    match_size = matches_till(array_to_find, suffix)
    match_size == 0 && return nothing 
    return suffix_start:suffix_start+match_size-1
end


function scan_flanks(sa::Vector{Int32}, insert_point::Int, array_to_find::AbstractVector{Int32}, concat_arr::AbstractVector{Int32})
    # Now we know where to look in the suffix area, check if we can prefix match the suffix left/right from this point
    concat_matching_indexes = []

    # Scan left of insert point
    left_shift = insert_point
    while left_shift - 1 > 0
        left_match_location = get_suffix_match_location(sa, concat_arr, array_to_find, left_shift - 1)
        isnothing(left_match_location) && break
        push!(concat_matching_indexes, left_match_location)
        left_shift -= 1
    end

    # Scan right of insert point (including the insert point index itself)
    right_shift = insert_point
    while right_shift <= length(sa)
        right_match_location = get_suffix_match_location(sa, concat_arr, array_to_find, right_shift)
        isnothing(right_match_location) && break 
        push!(concat_matching_indexes, right_match_location)
        right_shift += 1
    end

    return concat_matching_indexes

end

function find_longest_match(ref::Vector{Int32}, position::Int64, sa::Vector{Int32}, query_concat::Vector{Int32})
    array_to_find = view(ref, position:length(ref))
    insert_point = locate_insert_point(sa, query_concat, array_to_find)
    # Scan the left and right flank of the suffix array for matches
    concat_matching_indexes = scan_flanks(sa, insert_point, array_to_find, query_concat)
    return concat_matching_indexes
end

function slide_over_ref(ref_id::Int32, ref::Vector{Int32}, sa::Vector{Int32}, query_concat::Vector{Int32}, query_colors::Color, size_map::Dict{Int32,Int64})
    @inbounds for i in 1:length(ref)
        # Find all matches with the Qs for this sub-ref region 
        query_match_locations = find_longest_match(ref, i, sa, query_concat)
        if length(query_match_locations) > 0
            # For each location obtain the nt sizes
            for m in query_match_locations
                q_match_region = view(query_concat, m)
                # Get the nt sizes from the size map 
                q_match_size = sum(get.(Ref(size_map), q_match_region, 0))
                # If the sizes are bigger (based on the colors) replace origin and max size
                for i in m
                    if query_colors.max_len[i] < q_match_size
                        query_colors.max_len[i] = q_match_size
                        query_colors.origin[i] = ref_id
                    end
                end
            end
        end
    end
end

function find_longest_matches!(ref_id::Int32, ref::Vector{Int32}, sa::Vector{Int32}, query_concat::Vector{Int32}, query_colors::Color, size_map::Dict{Int32,Int64})
    # Convert the ref nodes to search forward matches
    convert_nodes!(ref)
    slide_over_ref(ref_id, ref, sa, query_concat, query_colors, size_map)

    # Then flip for reverse matches 
    ref = flipnode.(reverse!(ref))
    slide_over_ref(ref_id, ref, sa, query_concat, query_colors, size_map)

end



function start()
    # Test queries
    q1 = Int32[6, 100, 101, -5, -3]
    q2 = Int32[100, 101, 1,2,3,4]
    queries = [q1, q2]

    # The refs 
    ref1 = Int32[100, 101, 2,3,4, 3, 5]
    ref2 = Int32[100, -4, -3, -2, -1, 90]

    # Concat queries and build suffix array
    concat_arr, sa = create_k_suffix_array(queries, Int32(0))
    
    # Make a color vector holding the max length + origin 
    query_colors = Color(zeros(Int32, length(concat_arr)), zeros(Int32, length(concat_arr)))

    # Generate test nucleotide size mapping
    unique_nodes = Set(reduce(vcat, queries)) 
    #size_map = Dict(unique_nodes .=> rand(1:100,length(unique_nodes)))
    # Just easier to check for now if sizes are 1:
    size_map = Dict(unique_nodes .=> ones(Int64, length(unique_nodes))) 
    println(size_map)

    # Find the matches for the Qs with Ref1 and Ref2
    find_longest_matches!(Int32(1), ref1, sa, concat_arr, query_colors, size_map)
    find_longest_matches!(Int32(2), ref2, sa, concat_arr, query_colors, size_map)

    # Some output printing
    println("Suffix array")
    for (i, si) in enumerate(sa)
        println(i, " -> ", view(concat_arr, si:length(concat_arr)))
    end
    
    println("Queries")
    for query in queries
        println("Q: ", query)
    end 
    println()

    println("Concat array")
    println(concat_arr)
    println()
    println("Ref1: ", ref1)
    println("Ref2: ", ref2)
    println() 
    
    # Mask to add back the split identifiers from the concat 
    m = concat_arr .< Int32(0) 
    query_colors.origin[m] = concat_arr[m]
    query_colors.max_len[m] = concat_arr[m]
    println("ori\tlen")
    for (i, (origin, size)) in enumerate(zip(query_colors.origin,query_colors.max_len))
        if concat_arr[i] < 0
            println()
        else
            println(origin, "\t", size)
        end
    end
    #println("Origins: ",  query_colors.origin)
    #println("Lens: ", query_colors.max_len)
end 

start()


#println(matches_till(Int32[10,7,3], Int32[1,2,3,4,5]))
