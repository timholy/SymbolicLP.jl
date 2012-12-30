require("SymbolicLP")
using SymbolicLP
using GLPK

# Soft constraints
lpb = LPBlock([:x])
addconstraints(lpb, :(soft(x < -3, 1)))
z, x, flag = lpsolve(lpb)
@assert x[1] <= -3.0
lpb = LPBlock([:x])
addconstraints(lpb, :(soft(x > 3, 1)))
z, x, flag = lpsolve(lpb)
@assert x[1] >= 3.0
lpb = LPBlock([:x])
addconstraints(lpb, :(soft(x == -5, 1)))
z, x, flag = lpsolve(lpb)
@assert abs(x[1] + 5) < sqrt(eps())

# A simple case
lpb = LPBlock([:left, :middle, :right], Expr[:(left < middle < right)])
addconstraints(lpb,
               :(soft(middle-left==2*(right-middle), 5)),
               :(left == 0),
               :(right == 100))
lpp, chunks = lpparse(Float64, lpb)
lpd = lpeval(lpp)
z, x, flag = lpsolve(lpd)
x = x[chunks[1]]
@show x
@assert all(abs(x-[0,200/3,100]) .< sqrt(eps()))

# A more complex case: two blocks, with cross-references and a function call
lpb1 = LPBlock([:pigs, :cows, :chickens])
lpb2 = LPBlock([:stalls, :nests])
barndims = [5,10]  # in units of stall dimensions
nshelves = 4
nests_per_shelf = 10
addconstraints(lpb1, :(pigs >= 0), :(cows >= 0), :(chickens >= 0))
addconstraints(lpb1, :(pigs+cows<$lpb2[stalls]), :(chickens<$lpb2[nests]))
addconstraints(lpb2, :(nests <= $nshelves*$nests_per_shelf),
               :(stalls < prod($barndims)),
               :(stalls >= 0),
               :(nests >= 0))
addconstraints(lpb1, :(soft(pigs == 2*cows, 7)))
addobjective(lpb1, :(-80*pigs - 100*cows - 10*chickens))
lpp, chunks = lpparse(Float64, lpb1, lpb2)
lpd = lpeval(lpp)
z, x, flag = lpsolve(lpd, "Int", 1:5) #, "msg_lev", GLPK.MSG_ALL)
x1 = x[chunks[1]]
x2 = x[chunks[2]]
@show x1
@show x2
barndims[2] = 6
lpd = lpeval(lpp)
z, x, flag = lpsolve(lpd, "Int", 1:5)
x1 = x[chunks[1]]
x2 = x[chunks[2]]
@show x1
@show x2
