"""
Module này chứa các hàm tính chi phí. Chạy `?cost` để xem thêm.
"""
module CostFunctions
using LinearAlgebra
using ..MMs
using ..MarkovStuff

"""
	cost(m::AbstractMMCK, daily)

Tính chi phí trung bình trên một đơn vị thời gian của dãy xếp hàng
- `daily` là chi phí tốt cho một khách trên một đơn vị thời gian
"""
function cost(m::AbstractMMCK, daily)
	L(m) * daily
end

"""
	cost(m::AbstractMMCK, daily, outsource)

Tính chi phí trung bình trên một đơn vị thời gian của dãy xếp hàng
trong trường hợp thuê công nhân ngoài hệ thống phục vụ, chỉ áp dụng với
MMCK với `k` hữu hạn.
- `daily` là chi phí tốt cho một khách trên một đơn vị thời gian
- `outsource` là chi phí thuê công nhân ngoài hệ thống phục vụ
"""
function cost(m::AbstractMMCK, daily, outsource)
	cost(m, daily) + pn(m, m.k) * m.λ * outsource
end

"""
	cost(m::MarkovProcess, r::AbstractVector)

Chi phí trung bình dài hạn trên đơn vị thời gian của quá trình Markov `m`, tính theo vector lợi nhuận `r`.
"""
function cost(m::MarkovProcess, r::AbstractVector)
	Π = stationarydist(m)
	sum(r .* Π)
end

"""
	cost(m::MarkovProcess, r::AbstractVector, h::AbstractMatrix)

Chi phí trung bình dài hạn trên đơn vị thời gian của quá trình Markov `m`, tính theo vector lợi nhuận `r` và chi phí nhảy `h`.
"""
function cost(m::MarkovProcess,
	      r::AbstractVector,
	      h::AbstractMatrix)
	G = m.G
	n = size(G, 1)
	Π = stationarydist(m)
	sum(Π[i] * (r[i] + sum(G[i,k] * h[i,k] for k = 1:n)) for i = 1:n)
end

"""
	cost(m::MarkovProcess, r::AbstractVector, β::Real)

Chi phí trung bình dài hạn trên đơn vị thời gian của quá trình Markov `m`, tính theo vector lợi nhuận `r` và có chiết khấu `β`.
"""
function cost(m::MarkovProcess,
	      r::AbstractVector,
	      β::Real)
	inv(β * I - m.G) * r
end

export cost
end
