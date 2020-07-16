# The basics:
# 1) If the type of the extracted data is not defined, then it is String.
# 2) If not provided, the custom parser is `identity` for String, and
#    `parse` for any type `<: Number`.
# 3) If the `parse` throws, the extractor does not try to be smart.
# 4) In the case `extraction` returns NoDataFound, and the default is
#    NoDefault, then throw an error; otherwise return the default.
# 5) If no default is given, then NoDefault is used.

abstract struct AbstractExtractor{T} end
struct GenericExtractor{T} <: AbstractExtractor{T}
	extract
	default :: Union{T, NoDefault}
end
struct RegexExtractor{T} <: AbstractExtractor{T}
	regex :: Regex
	default :: Union{T, NoDefault}
end
function default(e :: GenericExtractor{T}) :: Union{T, NoDefault}
	return e.default
end
function default(e :: RegexExtractor{T}) :: Union{T, NoDefault}
	return e.default
end
function raw_extractor(e, ::Type{T}, default :: Union{T, NoDefault})
	
end
