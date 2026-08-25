#!/usr/bin/env julia
# Verify group operations from Windows Julia via Sage.
#
#   julia --project=tools/sage tools/sage/scripts/verify_groups.jl
#
# Prefers a persistent WSL Sage session; falls back to oneshot sage -python.
# Evidence status on success: verified-computationally for the finite cases below.

using Test
using JSON
using Pkg

Pkg.instantiate()
using SageBridge

function main()
    mkpath(joinpath(@__DIR__, "..", "artifacts"))
    println("=== Julia → Sage group-ops verification ===")

    sess = SageBridge.SageSession()  # Windows default: :wsl_file
    println("transport(initial)=$(sess.transport) url=$(sess.base_url)")
    SageBridge.ensure_server!(sess)
    println("transport(active)=$(sess.transport)")
    ver = try
        SageBridge.sage_version(sess)
    catch e
        string("unknown (", e, ")")
    end
    println("sage version: $ver")

    results = Dict{String,Any}(
        "sage_version" => ver,
        "transport" => string(sess.transport),
        "base_url" => sess.base_url,
        "cases" => Dict{String,Any}(),
    )

    # For oneshot mode, each eval carries its own setup code.
    oneshot = sess.transport === :wsl_oneshot

    @testset "Sage group ops via Julia bridge" begin
        @testset "SymmetricGroup(4)" begin
            if oneshot
                ord = Int(SageBridge.eval!(sess, "Integer(SymmetricGroup(4).order())")["result"])
                ab = Bool(SageBridge.eval!(sess, "bool(SymmetricGroup(4).is_abelian())")["result"])
                prod = string(SageBridge.eval!(
                    sess,
                    "str(p)";
                    code="""
S4 = SymmetricGroup(4)
p = S4("(1,2)") * S4("(2,3)")
""",
                )["result"])
                pord = Int(SageBridge.eval!(
                    sess,
                    "Integer(p.order())";
                    code="""
S4 = SymmetricGroup(4)
p = S4("(1,2)") * S4("(2,3)")
""",
                )["result"])
            else
                SageBridge.symmetric_group!(sess, 4; name="S4")
                ord = SageBridge.group_order(sess, "S4")
                ab = SageBridge.group_is_abelian(sess, "S4")
                prod = SageBridge.perm_product_repr(sess, "(1,2)", "(2,3)"; group_expr="S4")
                SageBridge.exec!(sess, "p = S4(\"(1,2)\") * S4(\"(2,3)\")")
                pord = SageBridge.element_order(sess, "p")
            end
            @test ord == 24
            @test ab == false
            @test occursin("1,3,2", prod)
            @test pord == 3
            results["cases"]["SymmetricGroup(4)"] = Dict(
                "order" => ord,
                "is_abelian" => ab,
                "product_(1,2)*(2,3)" => prod,
                "product_order" => pord,
                "ok" => true,
            )
            println("  S4.order() = $ord")
            println("  S4.is_abelian() = $ab")
            println("  (1,2)*(2,3) = $prod ; order = $pord")
        end

        @testset "AlternatingGroup(5)" begin
            ord = Int(SageBridge.eval!(sess, "Integer(AlternatingGroup(5).order())")["result"])
            @test ord == 60
            results["cases"]["AlternatingGroup(5)"] = Dict("order" => ord, "ok" => true)
            println("  A5.order() = $ord")
        end

        @testset "DihedralGroup(4) [order 8]" begin
            # Sage: DihedralGroup(n) = symmetries of the regular n-gon, order 2n.
            ord = Int(SageBridge.eval!(sess, "Integer(DihedralGroup(4).order())")["result"])
            ab = Bool(SageBridge.eval!(sess, "bool(DihedralGroup(4).is_abelian())")["result"])
            @test ord == 8
            @test ab == false
            results["cases"]["DihedralGroup(4)"] = Dict("order" => ord, "is_abelian" => ab, "ok" => true)
            println("  D4.order() = $ord ; abelian = $ab")
        end

        @testset "CyclicPermutationGroup(5)" begin
            ord = Int(SageBridge.eval!(
                sess,
                "Integer(C.order())";
                code="C = CyclicPermutationGroup(5)",
            )["result"])
            gord = Int(SageBridge.eval!(
                sess,
                "Integer(list(C.gens())[0].order())";
                code="C = CyclicPermutationGroup(5)",
            )["result"])
            @test ord == 5
            @test gord == 5
            results["cases"]["C5"] = Dict("order" => ord, "generator_order" => gord, "ok" => true)
            println("  C5.order() = $ord ; gen.order() = $gord")
        end

        @testset "GAP via Sage libgap" begin
            n = Int(SageBridge.eval!(sess, "Integer(libgap.SymmetricGroup(4).Size())")["result"])
            @test n == 24
            results["cases"]["libgap.SymmetricGroup(4).Size"] = Dict("order" => n, "ok" => true)
            println("  libgap SymmetricGroup(4).Size() = $n")
        end

        @testset "oneshot direct API" begin
            # Always exercise the oneshot path once, independent of session mode.
            r = SageBridge.oneshot_eval("Integer(SymmetricGroup(3).order())")
            @test Int(r["result"]) == 6
            results["cases"]["oneshot_S3"] = Dict("order" => 6, "ok" => true)
            println("  oneshot S3.order() = 6")
        end
    end

    out = joinpath(@__DIR__, "..", "artifacts", "verify_groups_result.json")
    open(out, "w") do io
        JSON.print(io, results, 2)
    end
    println("wrote $out")
    println("=== verification complete ===")
    return results
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
