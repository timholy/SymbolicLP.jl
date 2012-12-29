SymbolicLP.jl
=============

This package supports symbolic linear programming for the [Julia][Julia] programming language. It allows you to define the objective function and any constraints using symbols for the different variables, which in some cases may be easier than building matrices directly.

This package was written to support the needs of a linear programming-based layout manager. Consequently, it supports a number of "advanced" features:

- "Soft" constraints
- A syntax for building sub-blocks of linear programs and links between variables
- Facilities for delayed evaluation that allow solution of related linear programs differing in their constants.


## Installation

Within Julia, use the package manager:
```julia
load("pkg.jl")
Pkg.init()     # if you've never installed a package before
Pkg.add("SymbolicLP")
```

## Simple usage

To use this module, begin with
```julia
require("SymbolicLP")
using SymbolicLP
```

Let's take a very simple (in fact, trivial) example as a first demonstration. Let's say you have 3 variables, which we'll call `left`, `middle`, and `right`. We set up a linear programming problem for these variables in this way:

```julia
lpb = LPBlock([:left, :middle, :right])
```

An `LPBlock` means "linear programming block". The purpose of this notation will become clearer in the next example. We've initialized this block by saying that our program will be in terms of these three variables.

Our next step is to add some constraints. First, `left` should be to the left of `middle`, and `middle` should be to the left of `right`, which we can phrase as

```julia
addconstraints(lpb, :(left < middle < right))
```

To keep things simple, we're going to provide fixed values for `left` and `right`:

```julia
addconstraints(lpb, :(left == 0),
                    :(right == 100))
```

Finally, we'll say that the middle had best be in the middle:

```julia
addconstraints(lpb, :(middle-left==right-middle))
```

Now we solve this linear programming problem,

```julia
julia> z, x, flag = lpsolve(lpb);

julia> x
3-element Float64 Array:
  -0.0
  50.0
 100.0
```

If it's your first time calling `lpsolve`, there will be a substantial delay as Julia JIT-compiles all the linear programming routines. Subsequent calls will be much faster. However, suppose you'll be solving "similar" linear programs many times, perhaps with slightly different constants each time. If speed also matters, then this simple approach may not be ideal: Julia is having to parse the symbolic program into a numeric problem, then solve the numeric program.

It would be nice if there were a convenient way to do most of the symbolic parsing ahead of time, and then call the numeric solver with the new constants.

## Advanced usage

Let's say you're a farmer raising pigs, cows, and chickens. You get a payout of $80 per pig, $100 per cow, and $10 per chicken. However, your community values its historical balance of having approximately two pigs for every cow and imposes a tax if this balance is violated, in the amount of `$8*|pigs-2cows|`. You are also resource-constrained, in that your barn has a certain footprint, and this limits the number of nests you can have for the chickens and the number of stalls you can build for the pigs and cows.

Just for fun, let's build this program in two blocks:
```julia
lpb1 = LPBlock([:pigs, :cows, :chickens])
lpb2 = LPBlock([:stalls, :nests])
```
Now we have to specify relationships:
```julia
barndims = [5,10]  # in units of stall dimensions
nshelves = 4
nests_per_shelf = 10
```
and add constraints, such as the fact that there cannot be a negative number of pigs and cows, that each pig or cow needs a stall, and each chicken needs a nest:

```julia
addconstraints(lpb1, :(pigs >= 0), :(cows >= 0), :(chickens >= 0))
addconstraints(lpb2, :(stalls >= 0), :(nests >= 0))
addconstraints(lpb1, :(pigs+cows<$lpb2[stalls]), :(chickens<$lpb2[nests]))
```
Here you can see the syntax for making a cross-reference, `$lpb2[stalls]` rather than just `stalls`, so that linear program 1 can refer to variables in linear program 2.

We need space for our nests and stalls:
```julia
addconstraints(lpb2, :(nests <= $nshelves*$nests_per_shelf),
               :(stalls < prod($barndims)))
```
The key new feature here is the function call `prod($barndims)`. We could just supply this as a constant, but suppose the farmer hasn't decided on a particular barn. It might be nice to re-run the same linear program with different dimensions. The function call is not evaluated immediately, but instead after the symbolic part is parsed.

Finally, let's incorporate the tax and the payout:
```julia
addconstraints(lpb1, :(soft(pigs == 2*cows, 8)))
addobjective(lpb1, :(-80*pigs - 100*cows - 10*chickens))
```
Note that the tax is implemented as a "soft" constraint. The first argument of `soft` says what you want to be "approximately" true; however, this differs from "hard" constraints in that it is not required to be strictly observed. The second argument specifies how much of a price you pay when it's violated. Note that soft constraints can also take the form of inequalities.

Now, because we want to efficiently run this "same" linear program for different barn sizes, we first perform all the parsing:
```julia
lpp, chunks = lpparse(Float64, lpb1, lpb2)
```

`lpp` is a "parsed" linear program, but it hasn't yet evaluated that function call. Let's do that now:
```julia
lpd = lpeval(lpp)
```

Now we can solve the program. In this case, let's use integers for all variables:
```julia
z, x, flag = lpsolve(lpd, "Int", vcat(chunks...))
```
The `chunks` variable is an array of index vectors specifying how variables in the combined linear program map back into the separate programs `lpb1` and `lpb2`. If you inspect the contents of `chunks`, you'll notice that there have been two new variables added; these are required to implement the soft constraint.

If we want to change the barn dimensions, we just repeat the last few steps:
```julia
barndims[2] = 6
lpd = lpeval(lpp)
z, x, flag = lpsolve(lpd, "Int", vcat(chunks...))
```
Note that you can't just say `barndims = [5,6]`, because `prod($barndims)` was written in terms of a particular object (i.e., memory location), and replacing `barndims` disassociates the named variable from the original memory location. However, adjusting individual elements does not disassociate these two, and hence is one good way to support re-evaluation with alternate constant values.

Another good strategy uses types and accessor functions:
```julia
type Taxes
    tax1::Float64
    tax2::Float64
end
tax1(t::Taxes) = t.tax1
tax2(t::Taxes) = t.tax2

addconstraints(lpb1, :(soft(pigs == 2*cows, tax1($t))))
```
With this strategy, you can adjust `t.tax1 = 6` and re-run the program (starting with `lpeval` rather than `lpparse`) to see whether it's still advantageous to strive to meet the target balance.

A third strategy uses anonymous functions:
```julia
tax = 8
taxfun = () -> tax
addconstraints(lpb1, :(soft(pigs == 2*cows, taxfun())))
```
Adjusting the `tax` variable will now cause the new value to be used by `lpeval`.

[Julia]: http://julialang.org "Julia"
