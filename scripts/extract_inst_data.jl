import GuillotineModels
const IR = GuillotineModels.InstanceReader

function kp_ub(e, d_, p_, a_, A)
	se = sortperm(e; rev = true)
	d, p, a = getindex.((d_, p_, a_), (se, se, se))
	rem_A = A
	ub = 0.0
	for (i, qt) in pairs(IndexLinear(), d)
		for _ in 1:qt
			if rem_A < a[i]
				ub += floor((rem_A  * p[i]) / a[i])
				@goto out
			end
			rem_A -= a[i]
			ub += p[i]
		end
	end
	@label out
	return ub
end

function extract_inst_data(inst_name)
	N, L, W, l, w, p, d = IR.read_from_file(inst_name)
	n_ = sum(d)
	a = l .* w
	A = L .* W
	csv_lb = sum((d .* a) ./ A)
	e = p ./ (l .* w)
	kp_ub_ = kp_ub(e, d, p, a, A)
	min_l, min_w, min_a = minimum.((l, w, a))
	return inst_name, N, n_, L, W, A, min_l, min_w, min_a, kp_ub_, csv_lb
end

function create_csv_for(inst_names)
	header = "instance;N;n_;L;W;A;min_l;min_w;min_a;kp_ub;csv_lb\n"
	lines = map(inst_names) do inst
		values = extract_inst_data(inst)
		join(map(string, values), ";")
	end
	body = join(lines, "\n")
	csv = header * body

	return csv
end

const DATASET_C = vcat(
	"A" .* string.(1:5),
	"CHL" .* string.([1, 2, 5, 6, 7]),
	"CU" .* ["1", "2"],
	"CW" .* string.(1:3),
	"Hchl" .* ["2", "9"],
	"Hchl" .* string.([3, 4, 6, 7, 8]) .* "s",
	"HH",
	"OF" .* ["1", "2"],
	"STS" .* ["2", "4"],
	"W",
	"2",
	"3"
)

println(create_csv_for("../instances/" .* DATASET_C))

