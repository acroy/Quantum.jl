module t

type Trie
    value
    children::Dict{Any, Trie}
    is_leaf::Bool
    function Trie()
        self = new()
        self.children = Dict{Any, Trie}()
        self.is_leaf = false
        self
    end
end

function add_leaf(t::Trie, val, key)
    node = t
    for i in key
        if !haskey(node.children, i)
            node.children[i] = Trie()
        end
        node = node.children[i]
    end
    node.is_leaf = true
    node.value = val
end

function subtrie(t::Trie, prefix)
    node = t
    for i in prefix
        if !haskey(node.children, i)
            return nothing
        else
            node = node.children[i]
        end
    end
    node
end

get(t::Trie, key) = get(t, key, nothing)
function get(t::Trie, key, notfound)
    node = subtrie(t, key)
    if node != nothing && node.is_leaf
        return node.value
    end
    notfound
end

function keys(t::Trie, prefix, found)
    if t.is_key
        push!(found, prefix)
    end
    for (char,child) in t.children
        keys(child, ,found)
    end
end
keys(t::Trie, prefix::String) = (found=String[]; keys(t, prefix, found); found)
keys(t::Trie) = keys(t, "")

end#module