"""
	stationarydist(c::MarkovChain)

Phân phối dừng của xích markov
"""
function stationarydist(c::MarkovChain)
	P = m.P
	n = size(P, 1)
	A = [(transpose(p) - I)[1:end-1, :]; ones(1, n)]
	B = [zeros(n - 1); 1]
	inv(A) * B
end

function Base.show(io::IO, m::MarkovChain)
	display("Xích Markov với ma trận chuyển")
	display(m.P)
end
