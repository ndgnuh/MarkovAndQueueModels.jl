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
- Hiện có: `MM1`, `MM1K`, `MMC`, `MMCK`

Viết `?TênMôHình` để xem cách dùng.
"""
module MMs

#export MM1, MM1K, MMC,
#	
export AbstractMMCK, MM1, MM1K, MMC, MMCK,
	L, Lq, W, Wq, λe, pn, cost

abstract type AbstractMMCK end;

include(joinpath(@__DIR__, "mm1.jl"))
include(joinpath(@__DIR__, "mm1k.jl"))
include(joinpath(@__DIR__, "mmc.jl"))
include(joinpath(@__DIR__, "mmck.jl"))

function λe(m::AbstractMMCK)
	m.λ * (1 - pn(m, m.k))
end

end # module
