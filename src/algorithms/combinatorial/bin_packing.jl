using DataStructures
using DataStructures: SortedDict

function bin_packing(items::AbstractVector{Int64}, bin_capacity::Int64)
    n = length(items)
    if n == 0
        return 0, Vector{Vector{Int64}}()
    end
    
    # Create (original_index, size) pairs and sort descending by size
    items_sizes = [(i, items[i]) for i in 1:n]
    sort!(items_sizes, by = x -> -x[2])
    
    # Track bins: each bin is a list of item indices
    bins = Vector{Vector{Int64}}()
    
    # SortedDict: remaining_capacity -> list of bin indices with that capacity
    rem = SortedDict{Int64, Vector{Int64}}()
    
    for (idx, size) in items_sizes
        # Skip items that are too large (or handle as error)
        if size > bin_capacity
            @warn "Item $idx with size $size exceeds bin capacity $bin_capacity"
            continue
        end
        
        if size <= 0
            continue
        end
        
        # Find first bin with remaining capacity >= size (best fit from available)
        token = DataStructures.searchsortedfirst(rem, size)
        
        if token != DataStructures.pastendsemitoken(rem)
            # Found a bin that can fit this item
            capacity, bin_list = DataStructures.deref((rem, token))
            bin_idx = pop!(bin_list)
            if isempty(bin_list)
                delete!(rem, capacity)
            end
            
            push!(bins[bin_idx], idx)
            new_capacity = capacity - size
            if new_capacity > 0
                push!(get!(rem, new_capacity, Int64[]), bin_idx)
            end
        else
            # No existing bin can fit this item - create new bin
            push!(bins, Int64[idx])
            new_bin_idx = length(bins)
            new_capacity = bin_capacity - size
            if new_capacity > 0
                push!(get!(rem, new_capacity, Int64[]), new_bin_idx)
            end
        end
    end
    
    return length(bins), bins
end

# precompile(bin_packing, (Vector{Int64}, Int64))