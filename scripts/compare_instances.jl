using GuillotineModels.InstanceReader

# What do we want to test:
# 0) If two instances have the same number of pieces, L, and W.
# 1) If the two instances have the same plate types.
# 1.1) If they also share profit.
# 1.2) If they also share demand.
# 2) If reversing L and W (l and w) we then get the same as (1).

function b2s(b)
	return b ? "EQUAL" : "different"
end

macro tell_if_equal(a, b)
	v_a, v_b = esc(a), esc(b)
	s_a, s_b = string(a), string(b)
	return quote
		equal = $v_a == $v_b
		println($s_a * " == " * $s_b * ": " * b2s(equal))
		equal
	end
end

# NOTE: assume the pieces are sorted by their fields, in order.
function check_ptype_uniqueness(l, w, p, d)
	@assert isone(length(unique(length.([l, w, p, d]))))
	@assert length(l) > 0
	is_first = true
	llv, lwv, lpv, ldv = first.((l, w, p, d))
	for (lv, wv, pv, dv) in zip(l, w, p, d)
		if is_first
			is_first = false
			continue
		end
		if lv == llv && wv == lwv
			if pv == lpv
				if dv == ldv
					println("WARNING: two or more pieces have the same length, width," *
					" profit, and demand.")
				else
					println("WARNING: two or more pieces have the same length, width," *
					" and profit.")
				end
			else
				println("WARNING: two or more pieces have the same length and width.")
			end
			return
		end
		llv, lwv, lpv, ldv = lv, wv, pv, dv
	end
	return
end

# NOTE:
# We disregard:
# 1. The order of the pieces inside the instance.
# 2. The fact what one considers to be length the other considers to be width.
# We should also:
# 1. If the same piece type appears multiple times the lines should be combined.
function compare(path_1, path_2)
	N1, L1, W1, l1, w1, p1, q1 = read_from_file(path_1)
	N2, L2, W2, l2, w2, p2, q2 = read_from_file(path_2)
	# Sort the pieces in both instances.
	l1, w1, p1, q1 = map(i -> getindex.(sort(collect(zip(l1, w1, p1, q1))), i), 1:4)
	l2, w2, p2, q2 = map(i -> getindex.(sort(collect(zip(l2, w2, p2, q2))), i), 1:4)

	println("Comparing $(basename(path_1)) with $(basename(path_2))")
	check_ptype_uniqueness(l1, w1, p1, q1)
	check_ptype_uniqueness(l2, w2, p2, q2)
	if @tell_if_equal N1 N2
		if @tell_if_equal (L1, W1) (L2, W2)
			if @tell_if_equal (l1, w1) (l2, w2)
				@tell_if_equal p1 p2
				@tell_if_equal q1 q2
			end
		elseif @tell_if_equal (L1, W1) (W2, L2)
			if @tell_if_equal (l1, w1) (w2, l2)
				@tell_if_equal p1 p2
				@tell_if_equal q1 q2
			end
		end
	end
end

function compare_pairwise(folder_path)
	list = readdir(folder_path; join = true, sort = true)
	n = length(list)
	for i = 1:n, j = (i+1):n
		compare(list[i], list[j])
	end
	return
end

#compare(ARGS[1], ARGS[2])
compare_pairwise(ARGS[1])

