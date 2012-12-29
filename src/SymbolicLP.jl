# Utilities for linear programming

require("linprog.jl")
require("Debug")

module SymbolicLP
import GLPK
using LinProgGLPK  # for lpsolve
using Debug

import Base.(+), Base.(-), Base.(*), Base.(/), Base.split

export LPBlock, LPData, addconstraints, addobjective, lpparse, lpeval, lpsolve

## Problem specification syntax
# The total LP problem can be split up into blocks. The LPBlock type
#   stores data for a specific block.
# Variables within a block are referenced by a symbol, e.g., :height, :width
# Variables can be referenced across blocks as an expression :($obj[height]),
#   where obj is an LPBlock object
# Constraints can be expressed as :(height > 5) or :(height == 2*width) or
#   :($xlabel1[baseline] == $xlabel2[baseline]). Comparisons can be chained,
#   e.g., :(left < middle < right).
# The objective function to be minimized can be written :(6*width - height)
# There is automatic support for simple "soft" constraints,
#   :(soft(height >= 5, 7))
#   This allows height values less than 5, but any such values are penalized.
#   This is represented in the following way:
#     A new variable, which here will be called :delta
#     New "hard" (strict) constraints: :(height + delta >= 5), :(delta >= 0)
#     New term added to objective function: :(7*delta)
#   Soft constraints can be expressed as equalities,
#   :(soft(height == 2*width, 12))
#   This introduces two new variables, :deltapos and :deltaneg, four new
#     constraints, and an additional sum of two terms to the objective.
#   Soft constraints cannot be chained.
#   There cannot be any cross-references to auxillary variables for
#     soft constraints
# The LP is "compiled" (parsed) from LPBlocks using
#   lpparse(T, block1, block2, ...), merging all mutually-referencing blocks into
#   a single LP. This ultimately produces f, A, b, Aeq, ..., with the possible
#   exception of function calls (see next point). A parsed LP is stored in
#   the LPParsed type.
# Constraint expressions can include function calls, e.g., :(x <= sin(pi/2)).
#   These function calls are not evaluated when the LP problem is compiled; that
#   is done on a separate lpeval(parsedlp) step, resulting in an object of
#   type LPData. This has the actual, f, A, b, Aeq, ... that one feeds to the
#   LP solver.

# TODO? efficient enabling/disabling particular constraints. Do this by adding names for individual constrains in LPBlock.

type LPBlock
    syms::Vector{Symbol}
    constraints::Vector{Expr}
    objective::Vector{Expr}
end
LPBlock(syms::Vector{Symbol}, constraints::Vector{Expr}, obj::Expr) = LPBlock(syms, constraints, [obj])
LPBlock(syms::Vector{Symbol}, constraints::Vector{Expr}) = LPBlock(syms, constraints, Expr[])
LPBlock(syms::Vector{Symbol}) = LPBlock(syms, Expr[], Expr[])

function addconstraints(lpb::LPBlock, exa::Expr...)
    for ex in exa
        push(lpb.constraints, ex)
    end
    lpb
end

function addobjective(lpb::LPBlock, ex::Expr)
    push(lpb.objective, ex)
    lpb
end

type LPData{T}
    f::Vector{T}
    A::AbstractSparseMatrix{T}
    b::Vector{T}
    Aeq::AbstractSparseMatrix{T}
    beq::Vector{T}
    lb::Vector{T}
    ub::Vector{T}
end

type LPParsed{T}
    # Storage for explicit values
    vals::LPData{T}
    # Storage for function calls. Only one of vals.b[i] and funs.b[i] can
    # be non-zero.
    funs::LPData{Function}
    findex::Vector{Int}
    rowineqindex::Vector{Int}
    roweqindex::Vector{Int}
end

function lpeval{T}(::Type{T}, M::SparseMatrixCSC{Function})
    n = length(M.nzval)
    nzval = Array(T, n)
    nzvalfun = M.nzval
    for i = 1:n
        nzval[i] = nzvalfun[i]()
    end
    SparseMatrixCSC(M.m, M.n, copy(M.colptr), copy(M.rowval), nzval)
end
lpeval{T}(::Type{T}, fa::Vector{Function}) = T[convert(T, f()) for f in fa]
function stuff{T}(v::Vector{T}, indx::Vector{Int}, n::Int)
    a = zeros(T, n)
    a[indx] = v
    a
end

lpeval{T}(lpp::LPParsed{T}) = LPData{T}(
    isempty(lpp.findex) ? lpp.vals.f : lpp.vals.f + stuff(lpeval(T, lpp.funs.f), lpp.findex, length(lpp.vals.f)),
    lpp.vals.A + lpeval(T, lpp.funs.A),
    isempty(lpp.rowineqindex) ? lpp.vals.b : lpp.vals.b + stuff(lpeval(T, lpp.funs.b), lpp.rowineqindex, length(lpp.vals.b)),
    lpp.vals.Aeq + lpeval(T, lpp.funs.Aeq),
    isempty(lpp.roweqindex) ? lpp.vals.beq : lpp.vals.beq + stuff(lpeval(T, lpp.funs.beq), lpp.roweqindex, length(lpp.vals.beq)),
    lpp.vals.lb,
    lpp.vals.ub)


function lpsolve(lpd::LPData, pvtuple...)
    pv = {pvtuple...}
    intindx = find(pv[1:2:end] .== "Int")
    useint = !isempty(intindx)
    local param
    if useint
        param = GLPK.IntoptParam()
        colkind = fill(GLPK.CV, length(lpd.f))
        colkind[pv[intindx[end]+1]] = GLPK.IV
        keep = trues(length(pv))
        keep[intindx] = false
        keep[intindx+1] = false
        pv = pv[keep]
    else
        param = GLPK.SimplexParam()
    end
    param["msg_lev"] = GLPK.MSG_ERR
    param["presolve"] = GLPK.ON
    for i = 1:2:length(pv)
        param[pv[i]] = pv[i+1]
    end
    if useint
        return mixintprog(lpd.f, lpd.A, lpd.b, lpd.Aeq, lpd.beq, lpd.lb, lpd.ub, colkind, param)
    else
        return linprog_simplex(lpd.f, lpd.A, lpd.b, lpd.Aeq, lpd.beq, lpd.lb, lpd.ub, param)
    end
end

function lpsolve(lpb::LPBlock, pv...)
    lpp, chunks = lpparse(Float64, lpb)
    lpd = lpeval(lpp)
    z, x, flag = lpsolve(lpd, pv...)
    x = x[chunks[1]]
    z, x, flag
end

# Cannonical forms of comparison operators.
# For LP, <= is the same as <, >= is the same as >.
function getop(s::Symbol)
    s == :(<=) ? :(<) :
    s == :(>=) ? :(>) : s
end

# Parsing the LP problem
function parsesoft(blk::LPBlock)
    out = LPBlock(copy(blk.syms), Expr[], copy(blk.objective))
    for ex in blk.constraints
        if ex.head == :call && ex.args[1] == :soft
            exc::Expr = ex.args[2]
            pval = ex.args[3]
            if exc.head != :comparison || length(exc.args) != 3
                error("soft constraints must be simple comparisons (with <, >, <=, >=, or ==), but this constraint does not qualify: ", exc)
            end
            args = exc.args
            op = getop(args[2])
            if op == :(==)
                deltaneg = gensym()
                deltapos = gensym()
                push(out.syms, deltaneg)
                push(out.syms, deltapos)
                push(out.constraints, :($deltaneg < 0))
                push(out.constraints, :($deltapos > 0))
                push(out.constraints, :($(args[1]) - $deltaneg > $(args[3])))
                push(out.constraints, :($(args[1]) - $deltapos < $(args[3])))
                push(out.objective, :($pval*($deltapos-$deltaneg)))
            else
                delta = gensym()
                push(out.syms, delta)
                push(out.constraints, expr(:comparison, Any[delta, op, 0])) # :($delta $op 0))
                push(out.constraints, :($(args[1]) + $delta $op $(args[3])))
                if op == :(<)
                    push(out.objective, :(-$pval*$delta))
                else
                    push(out.objective, :($pval*$delta))
                end
            end
        else
            # This is not a soft constraint, just copy it
            push(out.constraints, ex)
        end
    end
    out
end

# Type for algebra on numbers and thunks
typealias NumOrThunk Union(Number, Function)

type LinearIndexExpr{T}
    indx::Vector{Int}
    coef::Vector{T}
    rhs::T
end
LinearIndexExpr{T}(::Type{T}, indx::Vector{Int}, coef::Vector, rhs) = LinearIndexExpr{T}(indx, convert(Array{T}, coef), convert(T, rhs))
indx(l::LinearIndexExpr) = l.indx
indx(l::Nothing) = Array(Int, 0)
coef{T}(::Type{T}, l::LinearIndexExpr{T}) = l.coef
coef{T}(::Type{T}, l::Nothing) = Array(T, 0)

@debug begin
function lpparse{T<:FloatingPoint}(::Type{T}, blocks::LPBlock...)
    nblocks = length(blocks)
    # Parse soft constraints
    sblocks = LPBlock[parsesoft(b) for b in blocks]
    # Build dictionaries to resolve references, including cross-block references
    len = Int[length(b.syms) for b in sblocks]
    offset = [0, cumsum(len)]
    blockindex = Dict{LPBlock, Int}(blocks, 1:nblocks)
    symindexes = Array(Dict{Symbol, Int}, nblocks)
    chunks = Array(Range1{Int}, nblocks)
    for i = 1:nblocks
        symindexes[i] = Dict(sblocks[i].syms, offset[i]+1:offset[i+1])
        chunks[i] = offset[i]+1:offset[i]+length(blocks[i].syms)
    end
    n = offset[end]  # total number of variables in LP problem
    
    eqrows = Array(LinearIndexExpr{NumOrThunk}, 0)
    ineqrows = Array(LinearIndexExpr{NumOrThunk}, 0)
    lb = fill(-inf(T), n)
    ub = fill(inf(T), n)
    objective = nothing
    for iblock in 1:nblocks
        block = sblocks[iblock]
        # Parse the constraints
        for ex in block.constraints
#             @bp
            if ex.head == :(=)
                error("Use ==, not =, in comparisons: ", ex)
            end
            if ex.head != :comparison
                error("Expression is not an equality or inequality:\n", ex)
            end
            # Parse the individual expressions
            pargs = Array(Any, length(ex.args))
            for i = 1:length(ex.args)
                if iseven(i)
                    pargs[i] = getop(ex.args[i])  # copy the comparison operator
                else
                    pargs[i] = resolve(NumOrThunk, ex.args[i], iblock, blockindex, symindexes)
                end
            end
            # Create the individual constraints from chained comparisons
            for i = 2:2:length(ex.args)
                op = pargs[i]
                darg = pargs[i-1] - pargs[i+1]
                if op == :(==)
                    # equality
                    push(eqrows, darg)
                else
                    # inequality
                    gflag = op == :(>)
                    if length(darg.indx) == 0
                        println("Warning: expression that uses no variables: ", ex.args[i-1:i+1])
                    elseif length(darg.indx) == 1 && isa(darg.coef[1], Number) && isa(darg.rhs, Number)
                        # upper/lower bound
#                         @bp
                        if darg.coef[1] < 0
                            gflag = !gflag
                        end
                        val = -darg.rhs/darg.coef[1]
                        j = darg.indx[1]
                        if gflag
                            lb[j] = max(val, lb[j])
                        else
                            ub[j] = min(val, ub[j])
                        end
                    else
                        # inequality
                        if gflag
                            darg.coef = _neg(darg.coef)
                            darg.rhs = _neg(darg.rhs)
                        end
                        push(ineqrows, darg)
                    end
                end
            end
        end
        # Parse the objective
        for term in block.objective
            objective += resolve(NumOrThunk, term, iblock, blockindex, symindexes)
        end
    end
    # Split the thunks from the values
    # We maintain thunks as a smaller array, so we don't have so many ()->0 to evaluate
    valineqrows = Array(LinearIndexExpr{T}, length(ineqrows))
    thunkineqrows = Array(LinearIndexExpr{Function}, 0)
    thunkineqindx = Array(Int, 0)  # stores the row#s that actually have thunks
    valeqrows = Array(LinearIndexExpr{T}, length(eqrows))
    thunkeqrows = Array(LinearIndexExpr{Function}, 0)
    thunkeqindx = Array(Int, 0)
    for i = 1:length(ineqrows)
        valineqrows[i], tmp = split(T, ineqrows[i])
        if !(tmp === nothing)
            push(thunkineqrows, tmp)
            push(thunkineqindx, i)
        end
    end
    for i = 1:length(eqrows)
        valeqrows[i], tmp = split(T, eqrows[i])
        if !(tmp === nothing)
            push(thunkeqrows, tmp)
            push(thunkeqindx, i)
        end
    end
    valf, thunkf = split(T, objective)
    # Create the matrix representation
    rows, cols, val, valb = sparsedata(valineqrows)
    valA = sparse(rows, cols, val, length(ineqrows), n)
    rows, cols, val, valbeq = sparsedata(valeqrows)
    valAeq = sparse(rows, cols, val, length(eqrows), n)
    lpval = LPData(tovector(T, valf, n), valA, valb, valAeq, valbeq, lb, ub)
    rows, cols, val, thunkb = sparsedata(thunkineqrows)
    thunkA = sparse(thunkineqindx[rows], cols, val, length(ineqrows), n)
    rows, cols, val, thunkbeq = sparsedata(thunkeqrows)
    thunkAeq = sparse(thunkeqindx[rows], cols, val, length(eqrows), n)
    lpthunk = LPData(coef(Function, thunkf), thunkA, thunkb, thunkAeq, thunkbeq, Function[], Function[])
    LPParsed(lpval, lpthunk, indx(thunkf), thunkineqindx, thunkeqindx), chunks
end
end # debug

# Resolve references and encode linear operations
resolve{T}(::Type{T}, arg::Number, ind::Int, blockindex::Dict{LPBlock, Int}, symindexes::Vector{Dict{Symbol, Int}}) = convert(T, arg)
resolve{T, Tv<:Number}(::Type{T}, arg::Vector{Tv}, ind::Int, blockindex::Dict{LPBlock, Int}, symindexes::Vector{Dict{Symbol, Int}}) = convert(Array{T}, arg)
resolve{T}(::Type{T}, arg::Function, ind::Int, blockindex::Dict{LPBlock, Int}, symindexes::Vector{Dict{Symbol, Int}}) = get(opmap, arg, arg)
function resolve{T}(::Type{T}, arg::Symbol, ind::Int, blockindex::Dict{LPBlock, Int}, symindexes::Vector{Dict{Symbol, Int}})
    i = get(symindexes[ind], arg, -1)
    if i == -1
        return eval(arg)  # try to evaluate it, perhaps it's a constant
    end
    LinearIndexExpr(T, [i], [one(T)], zero(T))
end
function resolve{T}(::Type{T}, ex::Expr, ind::Int, blockindex::Dict{LPBlock, Int}, symindexes::Vector{Dict{Symbol, Int}})
    local ret
    if ex.head == :ref
        i = get(blockindex, ex.args[1], -1)
        if i == -1
            error("Referenced object not resolvable. Is the list of blocks complete?")
        end
        ret = resolve(T, ex.args[2], i, blockindex, symindexes)
    elseif ex.head == :call
        if contains((:(+), :(-), :(*), :(/)), ex.args[1])
            n = length(ex.args)
            args = Array(Any, n)
            for i = 1:n
                args[i] = resolve(T, ex.args[i], ind, blockindex, symindexes)
            end
            ret = eval(expr(ex.head, args))
        else
            ret = () -> eval(ex)
        end
    end
    ret
end

# Types and operations for intermediate storage during parsing of constraints, references, and function calls
# Algebra with numbers and thunks. We don't want these exported, so give them private names.
_add(n1::Number, n2::Number) = n1+n2
_add(f::Function, n::Number) = () -> f()+n
_add(n::Number, f::Function) = () -> n+f()
_add(f1::Function, f2::Function) = () -> f1()+f2()
_sub(n1::Number, n2::Number) = n1-n2
_sub(f::Function, n::Number) = () -> f()-n
_sub(n::Number, f::Function) = () -> n-f()
_sub(f1::Function, f2::Function) = () -> f1()-f2()
_mul(n1::Number, n2::Number) = n1*n2
_mul(f::Function, n::Number) = () -> f()*n
_mul(n::Number, f::Function) = () -> n*f()
_mul(f1::Function, f2::Function) = () -> f1()*f2()
_div(n1::Number, n2::Number) = n1/n2
_div(f::Function, n::Number) = () -> f()/n
_div(n::Number, f::Function) = () -> n/f()
_div(f1::Function, f2::Function) = () -> f1()/f2()
_neg(n::Number) = -n
_neg(f::Function) = () -> -f()
_mul{T}(n::NumOrThunk, v::Vector{T}) = T[_mul(n, item) for item in v]
_neg{T}(v::Vector{T}) = T[_neg(item) for item in v]
_div{T}(v::Vector{T}, n::NumOrThunk) = T[_div(item, n) for item in v]
opmap = Dict(((+),(-),(*),(/)), (_add, _sub, _mul, _div))

(+){T}(l::LinearIndexExpr{T}, n::NumOrThunk) = LinearIndexExpr(T, l.indx, l.coef, _add(l.rhs, n))
(+)(n::NumOrThunk, l::LinearIndexExpr) = l+n
(+)(l::LinearIndexExpr, n::Nothing) = l
(+)(n::Nothing, l::LinearIndexExpr) = l
(*){T}(l::LinearIndexExpr{T}, n::NumOrThunk) = LinearIndexExpr(T, l.indx, _mul(n, l.coef), _mul(n, l.rhs))
(*)(n::NumOrThunk, l::LinearIndexExpr) = l*n
function (+){T}(l1::LinearIndexExpr{T}, l2::LinearIndexExpr{T})
    indx = Array(Int, 0)
    coef = Array(T, 0)
    ind1 = l1.indx
    ind2 = l2.indx
    i1 = 1
    i2 = 1
    ii1 = ind1[i1]
    ii2 = ind2[i2]
    while i1 <= length(ind1) || i2 <= length(ind2)
        if ii1 < ii2
            push(indx, ii1)
            push(coef, l1.coef[i1])
            i1 += 1
            ii1 = i1 > length(ind1) ? typemax(Int) : ind1[i1]
        elseif ii1 > ii2
            push(indx, ii2)
            push(coef, l2.coef[i2])
            i2 += 1
            ii2 = i2 > length(ind2) ? typemax(Int) : ind2[i2]
        else
            push(indx, ii1)
            push(coef, _add(l1.coef[i1], l2.coef[i2]))
            i1 += 1
            i2 += 1
            ii1 = i1 > length(ind1) ? typemax(Int) : ind1[i1]
            ii2 = i2 > length(ind2) ? typemax(Int) : ind2[i2]
        end
    end
    LinearIndexExpr(T, indx, coef, _add(l1.rhs,l2.rhs))
end
(-){T}(l::LinearIndexExpr{T}) = LinearIndexExpr(T, l.indx, _neg(l.coef), _neg(l.rhs))
(-)(l1::LinearIndexExpr, l2::LinearIndexExpr) = l1 + (-l2)
(-){T}(l::LinearIndexExpr{T}, n::NumOrThunk) = LinearIndexExpr(T, l.indx, l.coef, _sub(l.rhs,n))
(-){T}(n::NumOrThunk, l::LinearIndexExpr) = LinearIndexExpr(T, l.indx, _neg(l.coef), _sub(n,l.rhs))
(/){T}(l::LinearIndexExpr{T}, n::NumOrThunk) = LinearIndexExpr(T, l.indx, _div(l.coef, n), _div(l.rhs, n))

# Separate thunks from explicit values
getnum(x::Number) = x
getnum(x::Function) = 0
getthunk(x::Number) = () -> 0
getthunk(x::Function) = x
function split{T}(::Type{T}, l::LinearIndexExpr{NumOrThunk})
    retnum = LinearIndexExpr(T, l.indx, [getnum(c) for c in l.coef], getnum(l.rhs))
    isthunkcoef = [isa(c, Function) for c in l.coef]
    local retthunk
    if any(isthunkcoef) || isa(l.rhs, Function)
        coefthunk = [getthunk(c) for c in l.coef[isthunkcoef]]
        retthunk = LinearIndexExpr(Function, l.indx[isthunkcoef], coefthunk, getthunk(l.rhs))
    else
        retthunk = nothing
    end
    retnum, retthunk
end
# split{T, Tl<:Number}(::Type{T}, l::LinearIndexExpr{Tl}) = l, nothing
split{T}(::Type{T}, l::Nothing) = nothing, nothing

# Convert a list of LinearIndexExpr to format needed to create the sparse matrix
function sparsedata{T}(la::Vector{LinearIndexExpr{T}})
    len = [length(l.indx) for l in la]
    n = sum(len)
    rows = Array(Int, n)
    cols = Array(Int, n)
    val = Array(T, n)
    m = length(la)
    valb = Array(T, m)
    offset = 0
    for i = 1:m
        l = la[i]
        k = length(l.indx)
        r = offset+(1:k)
        rows[r] = i
        cols[r] = l.indx
        val[r] = l.coef
        valb[i] = _neg(l.rhs)
        offset += k
    end
    rows, cols, val, valb
end

function tovector(::Type{Function}, l::LinearIndexExpr{Function}, n::Int)
    f = fill(()->0, n)
    f[l.indx] = l.coef
    f
end
function tovector{T}(::Type{T}, l::LinearIndexExpr{T}, n::Int)
    f = zeros(T, n)
    f[l.indx] = l.coef
    f
end
tovector(::Type{Function}, l::Nothing, n::Int) = fill(()->0, n)
tovector{T}(::Type{T}, l::Nothing, n::Int) = zeros(T, n)

end # module
