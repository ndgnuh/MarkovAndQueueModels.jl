"""
	MM1K(λ, μ, k)

Tạo mô hình `M/M/1/K`.
"""
struct MM1K{T} <: AbstractMMCK
	λ::Union{T, Real}
	μ::Union{T, Real}
	k::Union{T, Integer}
	ρ::Union{T, Real}
	function MM1K(λ, μ, k)
		T = Union{typeof(μ), typeof(λ)}
		new{T}(λ, μ, k, λ/μ)
	end
end

function pn(m::MM1K, n::Int)
	ρ = m.ρ
	k = m.k
	if ρ == 1
		1 / (k + 1)
	else
		ρ^n * (1 - ρ) / (1 - ρ^(k+1))
	end
end

function L(m::MM1K)
	ρ = m.ρ
	k = m.k
	if ρ == 1
		k / 2
	else
		l = k + 1
		ρ * (1 + k*ρ^l - l*ρ^k) / (1 - ρ) / (1 - ρ^l)
	end
end

function Lq(m::MM1K)
	ρ = m.ρ
	k = m.k
	if ρ == 1
		k * (k - 1) / 2 / (k + 1)
	else
		L(m) - ρ * (1 - ρ^k) / (1 - ρ^(k+1))
	end
end

function W(m::MM1K)
	L(m) / λe(m)
end

function Wq(m::MM1K)
	W(m) - 1/m.μ
end
