# bob = {"A"=>("B", "E"),
# 	   "B"=>("A", "C"),
# 	   "C"=>("B", "E", "D"),
# 	   "D"=>("C",),
# 	   "E"=>("A","F","C"),
# 	   "F"=>("E",)}

# function bfs(dict, start, dest)
# 	next = {}
# 	visited = {}
# 	unshift!(next, init)
# 	push!(visited, init)
# 	while !isempty(next)
# 		current = pop!(next)
# 		println("current is now $current")
# 		if current==dest
# 			return current
# 		end
# 		for i in dict[current]
# 			print("is $i in $(repr(visited))?")
# 			if !in(i, visited)
# 				print(" NO!")
# 				push!(visited, i)
# 				unshift!(next, i)
# 			end
# 			println("")
# 		end
# 	end
# 	return None
# end

# function newbfs(dict, start, dest)
# 	next = {}
# 	visited = {}
# 	path = {}
# 	unshift!(next, init)
# 	push!(visited, init)
# 	push!(path, init)
# 	while !isempty(next)
# 		current = pop!(next)
# 		println("current is now $current")
# 		shared = intersect(dict[current], dict[dest])
# 		if !isempty(shared)
# 			return {}
# 		end
# 		for i in dict[current]
# 			print("is $i in $(repr(visited))?")
# 			if !in(i, visited)
# 				print(" NO!")
# 				push!(visited, i)
# 				unshift!(next, i)
# 			end
# 			println("")
# 		end
# 	end
# 	return None
# end

# function generate_conv{S<:String}(start::S, dest::S, conv::Converter)
# 	if start == dest
# 		error("Source and target bases are the same")
# 	end 
# 	if !in(start, collect(keys(conv.neighbors))) || !in(dest, collect(keys(conv.neighbors)))
# 		error("Either the source or target basis is missing from the given Converter")
# 	end
# 	init = conv.neighbors[start][1]
# 	coeffs = spones(conv.edges["$start_$init"])
# 	dest_reached = false
# 	visited = {}
# 	current = start
# 	while !dest_reached
# 		shared = intersect(conv.neighbors[current], conv.neighbors[dest]) 
# 		if !isempty(shared)
# 			coeffs = coeffs * conv.edges["$current_$(shared[1])"] * conv.edges["$(shared[1])_$dest"] 
# 			add_conversion!(start, dest, coeffs, conv)
# 			dest_reached=true
# 		end
# 		push!(visited, current)
# 		current = setdiff(conv.neighbors[current], visited)[1]
# 	end
# end

inds = [1:4]
bob = [[1:5], [6:11], [12:18], [19:24]]
#f(inds, bob)==[1, 7, 14, 22]

function map_inds(inds, bob)
	result = Array(Any, length(bob))
	for i=1:length(bob)
		result[i] = bob[i][inds[i]]
	end
	return result
end


