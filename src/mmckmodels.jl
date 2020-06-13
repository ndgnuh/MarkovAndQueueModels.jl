"""
Module này thực hiện tính toán với hàng đợi

# Cách dùng

- `m = TênMôHình(tham_số...)`: Tạo một mô hình, danh sách mô hình xem ở dưới
- `λe(m)`: Tính cường độ đến hiệu quả ``\\lambda_e``.
- `pn(m)`: Tính xác suất ``P\\left[N(t) = n\\right]``.
- `L(m)`: Lượng khách trung bình trong hệ thống
- `Lq(m)`: Lượng khách trung bình trong hàng đợi
- `W(m)`: Thời gian trung bình trong hệ thống
- `Wq(m)`: Thời gian trung bình trong hàng đợi

## Ví dụ
~~~julia
julia> include("$(@__FILE__)")
Main.MMCKModels

julia> using .MMCKModels

julia> m = MM1(6, 8)
MM1(6, 8, 0.75)

julia> L(m)
3.0

julia> Lq(m)
2.25

julia> W(m)
0.5
~~~

Các loại mô hình hàng đợi:
- Hiện có: `MM1`, `MM1K`
- Dự kiến: `MMC`

Viết `?TênMôHình` để xem cách dùng.
"""
module MMCKModels

export MM1, MM1K, MMC,
	L, Lq, W, Wq, λe, pn

abstract type AbstractMMCK end;

""" 
	MM1(μ, λ)

Tạo mô hình M/M/1
"""
struct MM1 <: AbstractMMCK
	λ; μ; ρ
	function MM1(λ::Real, μ::Real)
		if λ/μ >= 1
			@error "ρ phải < 1"
		end
		new(λ, μ, λ/μ)
	end
	function MM1(λ, μ)
		new(λ, μ, λ/μ)
	end
end

pn(m::MM1, n::Int) = (1 - m.ρ) * m.ρ^n
L(m::MM1) = m.ρ / (1 - m.ρ)
Lq(m::MM1) = m.ρ ^ 2 / (1 - m.ρ)
W(m::MM1) = L(m) / m.λ
Wq(m::MM1) = Lq(m) / m.λ

"""
	MM1K(λ, μ, k)

Tạo mô hình `M/M/1/K`.
"""
struct MM1K <: AbstractMMCK
	λ; μ; k; ρ
	MM1K(λ, μ, k) = new(λ, μ, k, λ/μ)
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

λe(m::MM1K) = m.λ * (1 - pn(m, m.k))
L(m::MM1K) = let ρ = m.ρ, k = m.k
	if ρ == 1
		k / 2
	else
		l = k + 1
		ρ * (1 + k*ρ^l - l*ρ^k) / (1 - ρ) / (1 - ρ^l)
	end
end
Lq(m::MM1K) = let ρ = m.ρ, k = m.k
	if ρ == 1
		k * (k - 1) / 2 / (k + 1)
	else
		L(m) - ρ * (1 - ρ^k) / (1 - ρ^(k+1))
	end
end
W(m::MM1K) = L(m) / λe(m)
Wq(m::MM1K) = W(m) - 1/μ

"""
	MMC(μ, λ, c)

Tạo mô hình M/M/c
"""
struct MMC <: AbstractMMCK
	μ; λ; c; r; ρ
	function MMC(μ, λ, c)
		new(μ, λ, c, λ/μ, λ/μ/c)
	end
end

function pn(m::MMC, n::Integer)
	μ = m.μ
	λ = m.λ
	c = m.c
	r = m.r
	ρ = m.ρ
	if n == 0
		1 / ((c * r^c) / factorial(c) / (c - r) + sum(r^n/factorial(i) for i = 0:c-1))
	elseif n < c
		pn(m, 0) * r^n / factorial(n)
	else
		pn(m, 0) * r^n / factorial(c) / c^(n-c)
	end
end

end
