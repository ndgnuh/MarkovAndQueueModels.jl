""" 
	MM1(μ, λ)

Tạo mô hình M/M/1
"""
struct MM1{T<:Any} <: AbstractMMCK
	λ::Union{T, Real}
	μ::Union{T, Real}
	ρ::Union{T, Real}
	function MM1(λ, μ)
		new{typeof(λ)}(λ, μ, λ/μ)
	end
	function MM1(λ::Real, μ::Real)
		if λ/μ >= 1
			@error "ρ phải < 1"
		end
		new{typeof(λ)}(λ, μ, λ/μ)
	end
end

pn(m::MM1, n::Int) = (1 - m.ρ) * m.ρ^n
L(m::MM1) = m.ρ / (1 - m.ρ)
Lq(m::MM1) = m.ρ ^ 2 / (1 - m.ρ)
W(m::MM1) = L(m) / m.λ
Wq(m::MM1) = Lq(m) / m.λ
