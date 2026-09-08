using Test

const PC = AkaghefAlgebra.Groups.PresentationCertificates

@testset "finite-presentation word-map certificates" begin
    @test_throws ArgumentError PC.Presentation([1, 1], Vector{Vector{Int}}())
    @test PC.free_reduce([1, -1, 2, -2, 3]) == [3]
    @test PC.inverse_word([1, 2, -3]) == [3, -2, -1]
    images = Dict(1 => [2, 3], 2 => [4], 3 => [-1])
    @test PC.substitute_words([1, -2, 3], images) == [2, 3, -4, -1]

    noncommutative = PC.Presentation([1, 2, 3, 4],
        [PC.Relator(7, [1, 3, 2, -4, 3]), PC.Relator(8, [2, 4, -1])])
    step, after = PC.eliminate_once(noncommutative, 7, 2)
    @test step.replacement == [-3, -1, -3, 4]
    @test 2 in noncommutative.generators && !(2 in after.generators)
    certificate = PC.reduce_presentation(noncommutative; choices=[(7, 2)])
    @test PC.verify_certificate(certificate)
    @test PC.verify_reduction(4, noncommutative.relators, certificate)
    tampered = deepcopy(certificate)
    tampered.steps[1].replacement[1] = -tampered.steps[1].replacement[1]
    @test_throws ArgumentError PC.verify_certificate(tampered)
    @test_throws ArgumentError PC.verify_certificate(
        certificate, PC.Presentation(3, [[1, 2]]))

    # A relator involving one generator does not consume unrelated free
    # generators; the automatic reducer leaves the latter represented.
    free = PC.reduce_presentation(3, [[1]])
    @test PC.verify_certificate(free)
    @test free.final_presentation.generators == [1, 2]
    @test isempty(free.final_presentation.relators)

    duplicate = PC.reduce_presentation(
        PC.Presentation([1, 2, 3],
            [PC.Relator(11, [1, 2]), PC.Relator(12, [-2, -1]),
             PC.Relator(13, [2, 1])]); choices=Tuple{Int,Int}[])
    unique = PC.reduce_presentation(3, [[1, 2]]; choices=Tuple{Int,Int}[])
    @test length(duplicate.final_presentation.relators) == 1
    @test PC.compare_reductions(duplicate, unique).matched

    poincare_trace = PC.reduce_presentation(2, [
        [-2, -2, 1, 2, -1, 2, 1],
        [-2, -1, 2, -1, -2, 1, 1]]; choices=Tuple{Int,Int}[])
    poincare_disk = PC.reduce_presentation(2, [
        [-2, -2, 1, 2, -1, 2, 1],
        [-2, -1, 2, -1, -2, 1, 1],
        [-2, -1, -2, 1, 2, 1]]; choices=Tuple{Int,Int}[])
    comparison = PC.compare_reductions(poincare_trace, poincare_disk)
    @test comparison.matched
    @test comparison.reason == :matched_bounded_relator_products
    @test length(comparison.witnesses.left_to_right) == 3
    @test length(comparison.witnesses.right_to_left) == 2
    @test PC.verify_witnesses(poincare_trace.final_presentation,
                              poincare_disk.final_presentation,
                              comparison.witnesses.left_to_right)
    @test PC.verify_witnesses(poincare_disk.final_presentation,
                              poincare_trace.final_presentation,
                              comparison.witnesses.right_to_left)
    witness = comparison.witnesses.right_to_left[1]
    bad_witness = PC.RelatorWitness(witness.source_relator_ids,
                                    witness.signs, witness.rotations,
                                    vcat(witness.product, [1]),
                                    witness.target_relator_id)
    @test_throws ArgumentError PC.verify_witnesses(
        poincare_disk.final_presentation, poincare_trace.final_presentation,
        [bad_witness; comparison.witnesses.right_to_left[2:end]])
    bad_id = PC.RelatorWitness([999], [1], [0], witness.product,
                                witness.target_relator_id)
    @test_throws ArgumentError PC.verify_witnesses(
        poincare_disk.final_presentation, poincare_trace.final_presentation,
        [bad_id; comparison.witnesses.right_to_left[2:end]])
    bad_sign = PC.RelatorWitness(witness.source_relator_ids,
                                  [0 for _ in witness.signs], witness.rotations,
                                  witness.product, witness.target_relator_id)
    @test_throws ArgumentError PC.verify_witnesses(
        poincare_disk.final_presentation, poincare_trace.final_presentation,
        [bad_sign; comparison.witnesses.right_to_left[2:end]])
    bad_rotation = PC.RelatorWitness(witness.source_relator_ids, witness.signs,
                                      [-1 for _ in witness.rotations],
                                      witness.product, witness.target_relator_id)
    @test_throws ArgumentError PC.verify_witnesses(
        poincare_disk.final_presentation, poincare_trace.final_presentation,
        [bad_rotation; comparison.witnesses.right_to_left[2:end]])

    unresolved = PC.compare_reductions(
        poincare_trace, PC.reduce_presentation(3, [[1, 2]]; choices=Tuple{Int,Int}[]))
    @test !unresolved.matched
    @test unresolved.isomorphism === nothing
end
