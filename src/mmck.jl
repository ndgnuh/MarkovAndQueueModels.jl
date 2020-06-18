"""
	MMCK(μ, λ, c, k)

Tạo mô hình M/M/c/k
"""
struct MMCK{T} <:AbstractMMCK
	λ::T
	μ::T
	c::Union{T, Integer}
	k::Union{T, Integer}
	r::Union{T, Real}
	ρ::Union{T, Real}
	function MMCK(λ, μ, c, k)
		T = Union{typeof(μ), typeof(λ)}
		new{T}(λ, μ, c,	k, λ/μ,	λ/μ/c)
	end
end

function pn(m::MMCK, n::Integer)
	r = m.r
	k = m.k
	c = m.c
	if n == 0
		1 / (
		     sum(r^i / factorial(i) for i = 0:c-1) +
		     sum(r^i / (c^(i-c) * factorial(c)) for i = c:k)
		     )
	elseif n < c
		pn(m, 0) * r^n / factorial(n)
	elseif n <= k
		pn(m, 0) * r^n / (c^(n-c) * factorial(c))
	else 0
	end
end

function L(m::MMCK)
	sum(i * pn(m, i) for i=0:m.k)
end

function Lq(m::MMCK)
	sum((i - m.c) * pn(m, i) for i=m.c:m.k)
end

function W(m::MMCK)
	L(m) / λe(m)
end

function Wq(m::MMCK)
	Lq(m) / λe(m)
end
