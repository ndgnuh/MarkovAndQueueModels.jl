"""
# Cách dùng
- `m = MarkovChain(P)`: Tạo xích Markov với ma trận chuyển `P`
- `m = MarkovProcess(G)`: Tạo quá trình Markov với ma trận sinh `G`
- `stationarydist(m)`: Phân phối dừng của `m`, với `m` là xích hoặc quá trình Markov
- `cost(m, ...)`: Chi phí lâu dài của quá trình Markov `m`, `?cost` để xem thêm.
- `imbeddedchain(m)`: Xích Markov nhúng của quá trình Markov `m`.
- `steadydist(m)`: Phân phối ổn định của quá trình Markov `m`.
"""
module MarkovStuff
using LinearAlgebra

export MarkovProcess, MarkovChain, stationarydist, cost,
	imbeddedchain, steadydist


"""
	MarkovChain(G::AbstractMatrix)

Xích Markov với ma trận chuyển P
"""
struct MarkovChain
	P::AbstractMatrix
	function MarkovChain(P::AbstractMatrix)
		if any(1 .≉ round.(sum(P; dims=2); digits=12))
			@error "Tổng hàng khác 1"
		else
			new(P)
		end
	end
end


"""
	MarkovProcess(G::AbstractMatrix)

Quá trình Markov với ma trận sinh G
"""
struct MarkovProcess
	G::AbstractMatrix
	function MarkovProcess(G::AbstractMatrix)
		if any(0 .≉ round.(sum(G; dims=2); digits=12))
			@error "Tổng hàng khác 0"
		else
			new(G)
		end
	end
end

include(joinpath(@__DIR__, "markovchain.jl"))
include(joinpath(@__DIR__, "markovprocess.jl"))

end
