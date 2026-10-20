#!/usr/bin/env julia
#
# src/parsecheck.jl -- syntax check for the Julia sources of this repository
#
#     julia --startup-file=no src/parsecheck.jl file1.jl file2.jl ...
#
# Exits nonzero if any of the files does not parse.  This exists because
# `make` never loads core_entropy_makie.jl -- it is the one Julia file no
# target touches -- and a missing closing parenthesis in it sat in the
# repository for several commits.
#
# The subtlety worth keeping in mind: Meta.parseall does not throw on a syntax
# error, it returns the diagnostic inside the AST as an Expr(:error).  A
# try/catch around it therefore reports success for any input, which is how the
# damaged file got through two checks.  The result has to be walked.

function diagnostics(ex)
    found = String[]
    walk(x) = begin
        if x isa Expr
            if x.head === :error || x.head === :incomplete
                push!(found, x.head === :error ? brief(x.args[1]) : "incomplete expression")
            end
            foreach(walk, x.args)
        end
    end
    walk(ex)
    return found
end

"""The parse error's own message, without its copy of the whole source."""
function brief(err)
    msg = err isa Base.Meta.ParseError ? err.msg : sprint(showerror, err)
    lines = split(msg, '\n')
    keep = String[]
    for l in lines
        push!(keep, l)
        length(keep) == 4 && break
    end
    return strip(join(keep, "\n"))
end

ok = true
for path in ARGS
    diags = diagnostics(Meta.parseall(read(path, String)))
    if isempty(diags)
        println("  parses  ", path)
    else
        global ok = false
        println("  FAILED  ", path)
        println(join(["          " * l for l in split(diags[1], '\n')], "\n"))
    end
end
exit(ok ? 0 : 1)
