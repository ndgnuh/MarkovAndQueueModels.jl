"""
# Cách dùng
- `m = MarkovChain(P)`: Tạo xích Markov với ma trận chuyển `P`
- `m = MarkovProcess(G)`: Tạo quá trình Markov với ma trận sinh `G`
- `stationarydist(m)`: Phân phối dừng của `m`, với `m` là xích hoặc quá trình Markov
- `cost(m, ...)`: Chi phí lâu dài của quá trình Markov `m`, `?cost` để xem thêm.
- `imbeddedchain(m)`: Xích Markov nhúng của quá trình Markov `m`.
- `transition_matrix(m, t)`: Tạo ma trận chuyển `P(t)` tại thời điểm `t`
"""
module MarkovStuff
using LinearAlgebra

export MarkovProcess, MarkovChain, stationarydist, cost,
	imbeddedchain, steadydist, transition_matrix


"""
	MarkovChain(G::AbstractMatrix)

Xích Markov với ma trận chuyển P
"""
struct MarkovChain
	P::AbstractMatrix
	function MarkovChain(P::AbstractMatrix{Real})
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

	function MarkovProcess(G::AbstractMatrix{Real})
		if any(0 .≉ round.(sum(G; dims=2); digits=12))
			@error "Tổng hàng khác 0"
		else
			new(G)
		end
	end

	function MarkovProcess(G::AbstractMatrix)
		new(G)
	end

	function MarkovProcess(m::MarkovChain, λ::AbstractVector)
		G = copy(m.P)
		G = G .* λ
		G[diagind(G)] .= .-λ
		return new(G)
	end
end

include(joinpath(@__DIR__, "markovchain.jl"))
include(joinpath(@__DIR__, "markovprocess.jl"))

end
