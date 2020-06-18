"""
	stationarydist(m::MarkovProcess)

Phân phối dừng của quá trình Markov
"""
function stationarydist(m::MarkovProcess)
	G = m.G
	A = vcat(transpose(G)[1:end-1, :], ones(1, size(G, 1)))
	B = vcat(zeros(size(G, 1) - 1), 1)
	inv(A) * B
end

"""
	imbeddedchain(m::MarkovProcess)

Xích Markov nhúng của một quá trình Markov
"""
function imbeddedchain(m::MarkovProcess)
	λ = -diag(m.G)
	p = m.G ./ λ - I
	p[diagind(p)] .= 0
	p[diagind(p)] = 1 .- sum(p; dims=2)
	MarkovChain(p) 
end

"""
	transition_matrix(m::MarkovProcess, t)

Tạo ma trận chuyển tại thời điểm t
"""
function transition_matrix(m::MarkovProcess, t)
	D, Q = eigen(m.G)
	Q * diagm(exp.(D * t)) * inv(Q)
end

function Base.show(io::IO, m::MarkovProcess)
	display("Quá trình Markov với ma trận sinh")
	display(m.G)
end
