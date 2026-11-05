# Command-line front end. `bin/simplexact.jl` calls `Simplexact.main(ARGS)`.

const USAGE = """
usage: simplexact <command> [options] FILE

commands:
  solve FILE     solve an LP-format file (branch and bound if it has integer variables)
  check FILE     parse the file and print a summary
  format FILE    print the file in canonical LP format

solve options:
  --relax                 ignore integrality and solve the LP relaxation
  --rule auto|bland|dantzig
  --strategy best_first|depth_first
  --node-limit N
  --decimal               also print decimal approximations
  --duals                 print shadow prices and reduced costs (LPs only)
  --no-verify             skip the independent certificate check
"""

function show_value(io, q::Q, dec::Bool)
    print(io, fmt(q))
    dec && !isone(denominator(q)) && print(io, "  (≈ ", decimal(q), ")")
end

"""
    main(args; io = stdout) -> Int

Run the CLI with `args` and return the process exit code: 0 on success
(whatever the status of the problem), 1 on usage, parse or verification errors.
"""
function main(args::Vector{String}; io::IO = stdout, err::IO = stderr)
    isempty(args) && (print(err, USAGE); return 1)
    cmd = args[1]
    cmd in ("-h", "--help", "help") && (print(io, USAGE); return 0)
    opts = Dict{String,String}()
    flags = Set{String}()
    files = String[]
    i = 2
    while i <= length(args)
        a = args[i]
        if a in ("--rule", "--strategy", "--node-limit")
            i < length(args) || (println(err, "missing value for $a"); return 1)
            opts[a] = args[i+1]
            i += 2
        elseif a in ("--relax", "--decimal", "--duals", "--no-verify")
            push!(flags, a)
            i += 1
        elseif startswith(a, "--")
            println(err, "unknown option $a"); return 1
        else
            push!(files, a)
            i += 1
        end
    end
    length(files) == 1 || (print(err, USAGE); return 1)
    lp = try
        read_lp(files[1])
    catch e
        e isa LPParseError || e isa SystemError || e isa ArgumentError || rethrow()
        print(err, files[1], ": "); showerror(err, e); println(err)
        return 1
    end
    if cmd == "format"
        write_lp(io, lp)
        return 0
    elseif cmd == "check"
        nint = count(lp.integer)
        println(io, files[1], ": ", lp.sense === :max ? "maximize" : "minimize", ", ",
            nvars(lp), " variables (", nint, " integer), ", nrows(lp), " rows")
        return 0
    elseif cmd != "solve"
        println(err, "unknown command $cmd"); print(err, USAGE); return 1
    end

    rule = Symbol(get(opts, "--rule", "auto"))
    strategy = Symbol(get(opts, "--strategy", "best_first"))
    node_limit = tryparse(Int, get(opts, "--node-limit", "100000"))
    node_limit === nothing && (println(err, "--node-limit needs an integer"); return 1)
    dec = "--decimal" in flags
    check = !("--no-verify" in flags)
    mip = any(lp.integer) && !("--relax" in flags)
    try
        if mip
            sol = solve_mip(lp; rule, strategy, node_limit)
            println(io, "status: ", sol.status)
            if sol.objective !== nothing
                print(io, "objective: "); show_value(io, sol.objective, dec); println(io)
            end
            if sol.status === :node_limit && sol.bound !== nothing
                print(io, "bound: "); show_value(io, sol.bound, dec); println(io)
            end
            println(io, "nodes: ", sol.nodes, ", simplex iterations: ", sol.lp_iterations)
            for j in eachindex(sol.x)
                print(io, "  ", lp.names[j], " = "); show_value(io, sol.x[j], dec); println(io)
            end
            if check && sol.status in (:optimal, :infeasible)
                reason = check_mip_certificate(lp, sol)
                reason === nothing || (println(io, "proof: INVALID (", reason, ")"); return 1)
                println(io, "proof: verified (branch-and-bound tree checked exactly)")
            end
        else
            sol = solve(lp; rule)
            println(io, "status: ", sol.status)
            if sol.status === :optimal
                print(io, "objective: "); show_value(io, sol.objective, dec); println(io)
            end
            println(io, "simplex iterations: ", sol.iterations)
            if sol.status === :optimal
                for j in eachindex(sol.x)
                    print(io, "  ", lp.names[j], " = "); show_value(io, sol.x[j], dec)
                    "--duals" in flags && (print(io, "    reduced cost "); show_value(io, sol.reduced_costs[j], dec))
                    println(io)
                end
                if "--duals" in flags
                    println(io, "shadow prices:")
                    for (i, r) in enumerate(lp.rows)
                        print(io, "  ", r.name, ": "); show_value(io, sol.duals[i], dec); println(io)
                    end
                end
            elseif sol.status === :unbounded
                println(io, "improving ray:")
                for j in eachindex(sol.ray)
                    iszero(sol.ray[j]) && continue
                    print(io, "  ", lp.names[j], ": "); show_value(io, sol.ray[j], dec); println(io)
                end
            elseif sol.status === :infeasible && !isempty(sol.farkas)
                println(io, "Farkas multipliers (rows combined into a contradiction):")
                for (i, r) in enumerate(lp.rows)
                    iszero(sol.farkas[i]) && continue
                    print(io, "  ", r.name, ": "); show_value(io, sol.farkas[i], dec); println(io)
                end
            end
            if check
                reason = check_certificate(lp, sol)
                reason === nothing || (println(io, "certificate: INVALID (", reason, ")"); return 1)
                println(io, "certificate: verified")
            end
        end
    catch e
        e isa ArgumentError || rethrow()
        print(err, "error: "); showerror(err, e); println(err)
        return 1
    end
    0
end
