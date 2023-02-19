const LIBSAIS = "/home/rickb/tools/suffix/libsais/libsais.so.2"
const MASK = Int32(1<<30) 

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
convert_node(n::Int32) = n < 0 ? flipnode(abs(n)) : noflip(n) # not sure what to name this fn

function convert_nodes(in_vector::Vector{Int32})
    return convert_node.(in_vector)
end

function concat_with_seperator(vectors::Vector{Vector{Int32}})
    # Probably a cleaner way to do this :) Like using map 
    # to flip the node ids in slices and copying the slices to the output
    total_size = sum(map(length, vectors))
    concat_arr = zeros(Int32, total_size + length(vectors)) 
    vect_id = -1 * length(vectors) # to have it decending for sorting
    vectors = map(convert_nodes, vectors)
    # # Concat with seperator + store sign in most significant bit
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
    while low <= high 
        mid = low + ((high - low) >>> 0x01) 
        suffix_start = sa[mid]
        suffix = view(concat_arr, suffix_start:suffix_arr_len)
        if ref < suffix 
            high = mid - 1
        elseif ref > suffix
            low = mid + 1
        else 
            return low # This would be an exact match
        end
    end
    return low
end

function matches_till(arr1::AbstractVector{Int32}, arr2::AbstractVector{Int32})
    smallest = min(length(arr1), length(arr2))
    for i in 1:smallest
        if arr1[i] != arr2[i]
            return i - 1
        end 
    end
    return smallest
end

function get_suffix_match_location(sa::Vector{Int32}, concat_arr::Vector{Int32}, array_to_find::AbstractVector{Int32}, point_in_sa::Int)
    # Get the suffix start position 
    suffix_start = sa[point_in_sa]
    # Get a view of the whole suffix to iteratore over 
    suffix = view(concat_arr, suffix_start:length(concat_arr))
    # Check for how long this equals the array_to_find
    match_size = matches_till(array_to_find, suffix)
    if match_size == 0
        return nothing 
    else
        return suffix_start:suffix_start+match_size-1
    end
end


function scan_flanks(sa::Vector{Int32}, insert_point::Int, array_to_find::AbstractVector{Int32}, concat_arr::AbstractVector{Int32})
    # Now we know where to look in the suffix area we check the midpoint, and left and right flanks for
    # the size of the match untill we don't have a match anymore
    concat_size = length(concat_arr)
    println("> Insert point in SA: ", insert_point)

    # Check the middle one 
    mount_match_location = get_suffix_match_location(sa, concat_arr, array_to_find, insert_point)
    if !isnothing(mount_match_location)
        println("\t== MATCH at insert point: ", concat_arr[mount_match_location], " at: ", insert_point, " in SA  and ", mount_match_location, " in concat array")
    end
    
    # Check the left of the insert point
    left_shift = insert_point
    while left_shift - 1 > 0
        left_match_location = get_suffix_match_location(sa, concat_arr, array_to_find, left_shift - 1)
        if isnothing(left_match_location)
            break 
        else
            println("\t<< MATCH LEFT ", concat_arr[left_match_location], " at: ", left_shift - 1, " in SA  and ", left_match_location, " in concat array")
        end
        left_shift -= 1
    end

    # Check the right flank
    right_shift = insert_point
    while right_shift + 1 <= length(sa)
        left_match_location = get_suffix_match_location(sa, concat_arr, array_to_find, right_shift + 1)
        if isnothing(left_match_location)
            break 
        else
            println("\t>> MATCH RIGHT: ", concat_arr[left_match_location], " at: ", right_shift + 1, " in SA  and ", left_match_location, " in concat array")
        end
        right_shift += 1
    end

end

function find_longest_match(ref::Vector{Int32}, position::Int64, sa::Vector{Int32}, lcp::Vector{Int32}, query_concat::Vector{Int32})
    # search locations in suffix array having the ref[position:] as prefix
    array_to_find = view(ref, position:length(ref))
    println("Ref partition: ", array_to_find )
    # Find the area of suffixes matching the first number of this ref part
    insert_point = locate_insert_point(sa, query_concat, array_to_find)
    # Check if we inserted somewhere in the SA
    if insert_point <= length(sa)
        # Scan the flanks for prefix matches with the suffix
        scan_flanks(sa, insert_point, array_to_find, query_concat)
    end
    println()
end

function find_longest_matches(ref::Vector{Int32}, sa::Vector{Int32}, lcp::Vector{Int32}, query_concat::Vector{Int32} )
    for i in 1:length(ref)
        find_longest_match(ref, i, sa, lcp, query_concat)
    end
end


function start()
    ref = Int32[1,2, 100, 80, 4, 3, 8]
    x = [Int32[1,2,3,4], Int32[5,2,3,1,2,8]]
    concat_arr, sa = create_k_suffix_array(x, Int32(0))

    println("Suffix array")
    for (i, si) in enumerate(sa)
        println(i, " -> ", concat_arr[si:length(concat_arr)])
    end
    println()
    println("Query concat array")
    println(concat_arr)
    println()
    println("Ref array")
    println(ref)
    println()

    find_longest_matches(ref, sa, lcp,  concat_arr)
end 

start()
