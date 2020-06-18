"""
	MMC(μ, λ, c)

Tạo mô hình M/M/c
"""
struct MMC{T} <: AbstractMMCK
	μ::Union{T, Real}
	λ::Union{T, Real}
	c::Union{T, Integer}
	r::Union{T, Real}
	ρ::Union{T, Real}
	function MMC(λ, μ, c)
		r = λ / μ
		ρ = r / c
		T = Union{typeof(μ), typeof(λ)}
		new{T}(μ, λ, c, r, ρ) 
	end
end

function pn(m::MMC, n::Integer)
	μ = m.μ
	λ = m.λ
	c = m.c
	r = m.r
	ρ = m.ρ
	if n == 0
		d1 = (c * r^c) / factorial(c) / (c - r)
		d2 = sum(r^i / factorial(i) for i = 0:c-1)
		1 / (d1 + d2)
	elseif n < c
		pn(m, 0) * r^n / factorial(n)
	elseif n <= m.k
		pn(m, 0) * r^n / factorial(c) / c^(n-c)
	else 0
	end
end

function Lq(m::MMC)
	r = m.r
	ρ = m.ρ
	c = m.c
	pn(m, 0) * r^c * ρ / factorial(m.c) / (1 - ρ)^2
end

function Wq(m::MMC)
	Lq(m) / m.λ
end

function W(m::MMC)
	Wq(m) + 1 / m.μ
end

function L(m::MMC)
	Lq(m) + m.r
end
