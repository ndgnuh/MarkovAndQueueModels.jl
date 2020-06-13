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
	steadydist(m::MarkovProcess)

Phân phối ổn định của quá trình Markov
"""
function steadydist(m::MarkovProcess)
	stationarydist(imbeddedchain(m))
end

function Base.show(io::IO, m::MarkovProcess)
	display("Quá trình Markov với ma trận sinh")
	display(m.G)
end
