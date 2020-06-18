using MarkovModelAndMMs.MarkovStuff
using MarkovModelAndMMs.MMs
using MarkovModelAndMMs.CostFunctions
using LinearAlgebra
using Test


let m = MarkovProcess([-1//2 1//4 1//4; 3//4 -1//1 1//4; 1//2 1//6 -2//3]),
	c = MarkovChain([0 1//2 1//2; 3//4 0 1//4; 3//4 1//4 0]),
	λ = [1//2, 1//1, 2//3]
	@test MarkovProcess(c, λ).G == m.G
	@test imbeddedchain(m).P == c.P
end
