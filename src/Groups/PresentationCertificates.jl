module Groups
module PresentationCertificates

export Relator, Presentation, EliminationStep, EliminationCertificate,
       PresentationIsomorphism, RelatorWitness, ReductionComparison,
       free_reduce, inverse_word, substitute_words,
       eliminate_once, reduce_presentation, verify_certificate,
       verify_reduction, verify_witnesses, compare_reductions

"""A reduced signed word in a free group is a vector of nonzero integers."""
free_reduce(word::AbstractVector{<:Integer}) = begin
    stack = Int[]
    for letter in word
        letter == 0 && throw(ArgumentError("free-group words cannot contain zero"))
        value = Int(letter)
        if !isempty(stack) && last(stack) == -value
            pop!(stack)
        else
            push!(stack, value)
        end
    end
    stack
end

"""Return the inverse word, with order reversed and signs changed."""
inverse_word(word::AbstractVector{<:Integer}) =
    Int[-Int(letter) for letter in reverse(word)] |> free_reduce

function _word_less(left::Vector{Int}, right::Vector{Int})
    n = min(length(left), length(right))
    for i in 1:n
        left[i] == right[i] && continue
        return left[i] < right[i]
    end
    length(left) < length(right)
end

function _rotations(word::Vector{Int})
    isempty(word) && return Vector{Int}[Int[]]
    [vcat(word[i:end], word[1:i-1]) for i in eachindex(word)]
end

"""Canonical representative of a cyclic relator, also identifying inversion."""
function _cyclic_normalize(word::AbstractVector{<:Integer})
    reduced = free_reduce(word)
    isempty(reduced) && return Int[]
    candidates = Vector{Int}[]
    for rotation in vcat(_rotations(reduced), _rotations(inverse_word(reduced)))
        push!(candidates, free_reduce(rotation))
    end
    result = first(candidates)
    for candidate in Iterators.drop(candidates, 1)
        _word_less(candidate, result) && (result = candidate)
    end
    copy(result)
end

function _word_equal(left::AbstractVector{<:Integer}, right::AbstractVector{<:Integer})
    Int.(left) == Int.(right)
end

function _word_dict_equal(left::Dict{Int,Vector{Int}}, right::Dict{Int,Vector{Int}})
    Set(keys(left)) == Set(keys(right)) &&
        all(_word_equal(left[key], right[key]) for key in keys(left))
end

"""A presentation relator retaining its original stable identifier."""
struct Relator
    id::Int
    word::Vector{Int}
end

function Relator(id::Integer, word::AbstractVector{<:Integer})
    id > 0 || throw(ArgumentError("relator ids must be positive"))
    Relator(Int(id), Int.(word))
end

"""A finite presentation with stable positive generator and relator ids."""
struct Presentation
    generators::Vector{Int}
    relators::Vector{Relator}
    function Presentation(generators::Vector{Int}, relators::Vector{Relator})
        length(Set(generators)) == length(generators) ||
            throw(ArgumentError("generator ids must be unique"))
        all(>(0), generators) || throw(ArgumentError("generator ids must be positive"))
        ids = [relator.id for relator in relators]
        length(Set(ids)) == length(ids) || throw(ArgumentError("relator ids must be unique"))
        all(>(0), ids) || throw(ArgumentError("relator ids must be positive"))
        allowed = Set(generators)
        for relator in relators
            all(letter -> abs(letter) in allowed, relator.word) ||
                throw(ArgumentError("relator uses an inactive generator"))
        end
        new(copy(generators), [Relator(r.id, r.word) for r in relators])
    end
end

function Presentation(generators::AbstractVector{<:Integer}, words::AbstractVector)
    ids = Int.(generators)
    relators = Relator[Relator(i, Int.(word)) for (i, word) in enumerate(words)]
    Presentation(ids, relators)
end

Presentation(generator_count::Integer, words::AbstractVector) =
    generator_count >= 0 ? Presentation(collect(1:Int(generator_count)), words) :
    throw(ArgumentError("generator count must be nonnegative"))

function _normalize_presentation(presentation::Presentation)
    generators = sort(unique(presentation.generators))
    relators = Relator[]
    for relator in presentation.relators
        word = _cyclic_normalize(relator.word)
        isempty(word) && continue
        push!(relators, Relator(relator.id, word))
    end
    sort!(relators, by=relator -> relator.id)
    Presentation(generators, relators)
end

function _presentation_equal(left::Presentation, right::Presentation)
    left.generators == right.generators || return false
    length(left.relators) == length(right.relators) || return false
    all(a.id == b.id && a.word == b.word for (a, b) in zip(left.relators, right.relators))
end

function _canonical_final(presentation::Presentation)
    active = sort(presentation.generators)
    renumber = Dict{Int,Int}(old => i for (i, old) in enumerate(active))
    words = Vector{Vector{Int}}()
    for relator in presentation.relators
        word = Int[sign(letter) * renumber[abs(letter)] for letter in relator.word]
        word = _cyclic_normalize(word)
        # Stable relator ids are retained in every elimination step.  Only
        # this final presentation normalization removes redundant copies of
        # the same cyclic relator (including inverse and cyclic conjugates).
        isempty(word) && continue
        any(existing -> existing == word, words) || push!(words, word)
    end
    sort!(words, lt=_word_less)
    final_relators = Relator[Relator(i, word) for (i, word) in enumerate(words)]
    Presentation(collect(1:length(active)), final_relators), renumber, active
end

function _image(images::AbstractDict{<:Integer,<:AbstractVector}, generator::Int)
    haskey(images, generator) || throw(ArgumentError("missing image for generator $generator"))
    Int.(images[generator])
end

"""Substitute signed generator images into a free-group word."""
function substitute_words(word::AbstractVector{<:Integer},
                          images::AbstractDict{<:Integer,<:AbstractVector})
    result = Int[]
    for letter in word
        letter == 0 && throw(ArgumentError("free-group words cannot contain zero"))
        image = _image(images, abs(Int(letter)))
        append!(result, letter > 0 ? image : inverse_word(image))
    end
    free_reduce(result)
end

function _substitute_one(word::AbstractVector{<:Integer}, generator::Int,
                         replacement::Vector{Int})
    result = Int[]
    for letter in word
        value = Int(letter)
        value == 0 && throw(ArgumentError("free-group words cannot contain zero"))
        if abs(value) == generator
            append!(result, value > 0 ? replacement : inverse_word(replacement))
        else
            push!(result, value)
        end
    end
    free_reduce(result)
end

struct EliminationStep
    selected_relator_id::Int
    eliminated_generator::Int
    epsilon::Int
    A::Vector{Int}
    B::Vector{Int}
    replacement::Vector{Int}
    before::Presentation
    after::Presentation
    images_after::Dict{Int,Vector{Int}}
end

struct EliminationCertificate
    original::Presentation
    steps::Vector{EliminationStep}
    final_current::Presentation
    final_presentation::Presentation
    forward_images::Dict{Int,Vector{Int}}
    backward_images::Dict{Int,Vector{Int}}
    final_renumber::Dict{Int,Int}
end

struct PresentationIsomorphism
    source::Presentation
    target::Presentation
    forward::Dict{Int,Vector{Int}}
    backward::Dict{Int,Vector{Int}}
end

"""A bounded relator-normal-closure witness on a common final generator set.

`source_relator_ids`, `signs`, and `rotations` describe one or two factors;
`product` is the resulting signed free word before cyclic normalization.
The witness is intentionally small: it certifies only products of at most two
cyclic conjugates (possibly inverted), never arbitrary normal-closure search.
"""
struct RelatorWitness
    source_relator_ids::Vector{Int}
    signs::Vector{Int}
    rotations::Vector{Int}
    product::Vector{Int}
    target_relator_id::Int
end

struct ReductionComparison
    matched::Bool
    isomorphism::Union{Nothing,PresentationIsomorphism}
    reason::Symbol
    # Directions are explicit because each final relator must be certified in
    # both normal-closure inclusions.
    witnesses::NamedTuple{(:left_to_right, :right_to_left),
                          Tuple{Vector{RelatorWitness},Vector{RelatorWitness}}}
end

ReductionComparison(matched::Bool,
                    isomorphism::Union{Nothing,PresentationIsomorphism},
                    reason::Symbol) =
    ReductionComparison(matched, isomorphism, reason,
                        (left_to_right=RelatorWitness[],
                         right_to_left=RelatorWitness[]))

function _occurrence(relator::Relator, generator::Int)
    positions = findall(letter -> abs(letter) == generator, relator.word)
    length(positions) == 1 || return nothing
    position = only(positions)
    position, sign(relator.word[position]), relator.word[1:position-1], relator.word[position+1:end]
end

function _candidate_pairs(presentation::Presentation)
    candidates = Tuple{Int,Int}[]
    for relator in sort(presentation.relators, by=r -> r.id)
        for generator in sort(presentation.generators)
            _occurrence(relator, generator) === nothing || push!(candidates, (relator.id, generator))
        end
    end
    candidates
end

function _eliminate(presentation::Presentation, selected_relator_id::Int,
                    generator::Int, images::Dict{Int,Vector{Int}})
    before = _normalize_presentation(presentation)
    generator in before.generators || throw(ArgumentError("selected generator is inactive"))
    index = findfirst(relator -> relator.id == selected_relator_id, before.relators)
    index === nothing && throw(ArgumentError("selected relator is inactive"))
    selected = before.relators[index]
    occurrence = _occurrence(selected, generator)
    occurrence === nothing && throw(ArgumentError("selected generator must occur exactly once"))
    _, epsilon, A, B = occurrence
    replacement = epsilon == 1 ? inverse_word(vcat(B, A)) : free_reduce(vcat(B, A))
    active_generators = [g for g in before.generators if g != generator]
    next_relators = Relator[]
    for relator in before.relators
        relator.id == selected_relator_id && continue
        word = _substitute_one(relator.word, generator, replacement)
        isempty(word) && continue
        push!(next_relators, Relator(relator.id, _cyclic_normalize(word)))
    end
    after = Presentation(active_generators, next_relators)
    next_images = Dict{Int,Vector{Int}}()
    for (original, image) in images
        next_images[original] = _substitute_one(image, generator, replacement)
    end
    step = EliminationStep(selected_relator_id, generator, epsilon, Int.(A), Int.(B),
                           replacement, before, after, next_images)
    step, after, next_images
end

"""Perform one deterministic elementary Tietze elimination."""
function eliminate_once(presentation::Presentation, selected_relator_id::Integer,
                        generator::Integer)
    normalized = _normalize_presentation(presentation)
    images = Dict{Int,Vector{Int}}(g => [g] for g in normalized.generators)
    step, after, _ = _eliminate(normalized, Int(selected_relator_id), Int(generator), images)
    step, after
end

function _parse_choices(choices)
    choices === nothing && return nothing
    Tuple{Int,Int}[(Int(pair[1]), Int(pair[2])) for pair in choices]
end

"""Reduce by deterministic elementary eliminations and retain a replay log."""
function reduce_presentation(presentation::Presentation; choices=nothing)
    original = _normalize_presentation(presentation)
    current = original
    images = Dict{Int,Vector{Int}}(g => [g] for g in original.generators)
    steps = EliminationStep[]
    requested = _parse_choices(choices)
    cursor = 1
    while !isempty(current.generators)
        choice = if requested === nothing
            candidates = _candidate_pairs(current)
            isempty(candidates) ? nothing : first(candidates)
        else
            cursor > length(requested) ? nothing : requested[cursor]
        end
        choice === nothing && break
        step, current, images = _eliminate(current, choice[1], choice[2], images)
        push!(steps, step)
        requested === nothing || (cursor += 1)
    end
    requested === nothing || cursor == length(requested) + 1 ||
        throw(ArgumentError("choice list contains an unreachable elimination"))
    final_presentation, renumber, active = _canonical_final(current)
    forward = Dict{Int,Vector{Int}}()
    for (original_id, image) in images
        forward[original_id] = Int[sign(letter) * renumber[abs(letter)] for letter in image]
    end
    backward = Dict{Int,Vector{Int}}(renumber[old] => [old] for old in active)
    EliminationCertificate(original, steps, current, final_presentation,
                           forward, backward, renumber)
end

reduce_presentation(generator_count::Integer, relators::AbstractVector; choices=nothing) =
    reduce_presentation(Presentation(generator_count, relators); choices=choices)

function _throw_bad(message)
    throw(ArgumentError("invalid elimination certificate: $message"))
end

"""Replay every step and reject any altered step, final presentation, or map."""
function verify_certificate(certificate::EliminationCertificate)
    original = _normalize_presentation(certificate.original)
    _presentation_equal(original, certificate.original) || _throw_bad("original presentation")
    current = original
    images = Dict{Int,Vector{Int}}(g => [g] for g in original.generators)
    for step in certificate.steps
        _presentation_equal(step.before, current) || _throw_bad("step before-state")
        step.eliminated_generator in current.generators || _throw_bad("inactive generator")
        index = findfirst(relator -> relator.id == step.selected_relator_id, current.relators)
        index === nothing && _throw_bad("inactive selected relator")
        selected = current.relators[index]
        occurrence = _occurrence(selected, step.eliminated_generator)
        occurrence === nothing && _throw_bad("selected letter is not unique")
        _, epsilon, A, B = occurrence
        epsilon == step.epsilon || _throw_bad("selected letter sign")
        A == step.A && B == step.B || _throw_bad("selected relator split")
        expected = epsilon == 1 ? inverse_word(vcat(B, A)) : free_reduce(vcat(B, A))
        expected == step.replacement || _throw_bad("replacement word")
        recomputed, next, next_images = _eliminate(current,
            step.selected_relator_id, step.eliminated_generator, images)
        _presentation_equal(recomputed.after, step.after) || _throw_bad("step after-state")
        _word_dict_equal(next_images, step.images_after) || _throw_bad("step image map")
        current = next
        images = next_images
    end
    _presentation_equal(current, certificate.final_current) || _throw_bad("final current presentation")
    final, renumber, active = _canonical_final(current)
    _presentation_equal(final, certificate.final_presentation) || _throw_bad("final presentation")
    renumber == certificate.final_renumber || _throw_bad("final renumbering")
    forward = Dict{Int,Vector{Int}}()
    for (original_id, image) in images
        forward[original_id] = Int[sign(letter) * renumber[abs(letter)] for letter in image]
    end
    backward = Dict{Int,Vector{Int}}(renumber[old] => [old] for old in active)
    _word_dict_equal(forward, certificate.forward_images) || _throw_bad("forward map")
    _word_dict_equal(backward, certificate.backward_images) || _throw_bad("backward map")
    true
end

"""Replay a certificate only after binding it to the caller's presentation."""
function verify_certificate(certificate::EliminationCertificate,
                            expected::Presentation)
    _presentation_equal(_normalize_presentation(expected), certificate.original) ||
        _throw_bad("certificate is bound to a different original presentation")
    verify_certificate(certificate)
end

"""Count/relator adapter for trace-presentations producers."""
function verify_reduction(generator_count::Integer, relators::AbstractVector,
                          certificate::EliminationCertificate)
    verify_certificate(certificate, Presentation(generator_count, relators))
end

function _rotate_zero_based(word::Vector{Int}, rotation::Int)
    isempty(word) && return Int[]
    r = mod(rotation, length(word))
    r == 0 ? copy(word) : vcat(word[r + 1:end], word[1:r])
end

function _factor_options(relator::Relator)
    options = NamedTuple{(:sign, :rotation, :word),
                         Tuple{Int,Int,Vector{Int}}}[]
    for sign_value in (-1, 1)
        base = sign_value == 1 ? copy(relator.word) : inverse_word(relator.word)
        for rotation in 0:length(base)-1
            push!(options, (sign=sign_value, rotation=rotation,
                            word=_rotate_zero_based(base, rotation)))
        end
    end
    options
end

function _bounded_relator_witness(source::Presentation, target::Relator)
    source_relators = sort(source.relators, by=r -> r.id)
    options = Dict(r.id => _factor_options(r) for r in source_relators)

    # First try one factor, then two factors.  This is deliberately a finite
    # bounded predicate, rather than a normal-closure or CAS search.
    for first in source_relators
        for option in options[first.id]
            product = free_reduce(option.word)
            _cyclic_normalize(product) == target.word &&
                return RelatorWitness([first.id], [option.sign],
                                      [option.rotation], product, target.id)
        end
    end
    for first in source_relators, second in source_relators
        for option_first in options[first.id], option_second in options[second.id]
            product = free_reduce(vcat(option_first.word, option_second.word))
            _cyclic_normalize(product) == target.word &&
                return RelatorWitness([first.id, second.id],
                                      [option_first.sign, option_second.sign],
                                      [option_first.rotation, option_second.rotation],
                                      product, target.id)
        end
    end
    nothing
end

function _bounded_relator_witnesses(source::Presentation, target::Presentation)
    source.generators == target.generators || return nothing
    witnesses = RelatorWitness[]
    for target_relator in sort(target.relators, by=r -> r.id)
        witness = _bounded_relator_witness(source, target_relator)
        witness === nothing && return nothing
        push!(witnesses, witness)
    end
    witnesses
end

"""Replay bounded relator witnesses and reject malformed or tampered data."""
function verify_witnesses(source::Presentation, target::Presentation,
                          witnesses::AbstractVector{<:RelatorWitness})
    source = _normalize_presentation(source)
    target = _normalize_presentation(target)
    source.generators == target.generators ||
        throw(ArgumentError("witness presentations have different final generators"))
    length(witnesses) == length(target.relators) ||
        throw(ArgumentError("witness coverage is incomplete"))
    target_by_id = Dict(relator.id => relator for relator in target.relators)
    source_by_id = Dict(relator.id => relator for relator in source.relators)
    seen = Set{Int}()
    for witness in witnesses
        haskey(target_by_id, witness.target_relator_id) ||
            throw(ArgumentError("witness targets an inactive relator"))
        witness.target_relator_id in seen &&
            throw(ArgumentError("witness target is duplicated"))
        push!(seen, witness.target_relator_id)
        length(witness.source_relator_ids) in (1, 2) ||
            throw(ArgumentError("bounded witness must have one or two factors"))
        length(witness.source_relator_ids) == length(witness.signs) ==
            length(witness.rotations) ||
            throw(ArgumentError("witness factor metadata has inconsistent lengths"))
        factors = Vector{Vector{Int}}()
        for i in eachindex(witness.source_relator_ids)
            id = witness.source_relator_ids[i]
            haskey(source_by_id, id) ||
                throw(ArgumentError("witness uses an inactive source relator"))
            sign_value = witness.signs[i]
            sign_value in (-1, 1) || throw(ArgumentError("witness sign is invalid"))
            source_word = source_by_id[id].word
            rotation = witness.rotations[i]
            0 <= rotation < length(source_word) ||
                throw(ArgumentError("witness rotation is invalid"))
            base = sign_value == 1 ? source_word : inverse_word(source_word)
            push!(factors, _rotate_zero_based(base, rotation))
        end
        product = Int[]
        for factor in factors
            append!(product, factor)
        end
        product = free_reduce(product)
        product == witness.product || throw(ArgumentError("witness product"))
        target_by_id[witness.target_relator_id].word == _cyclic_normalize(product) ||
            throw(ArgumentError("witness target word"))
    end
    length(seen) == length(target.relators) ||
        throw(ArgumentError("witness coverage is incomplete"))
    true
end

function _identity_on_common_final_generators(left::EliminationCertificate,
                                               right::EliminationCertificate)
    forward = Dict{Int,Vector{Int}}()
    for (generator, image) in left.forward_images
        forward[generator] = substitute_words(image, right.backward_images)
    end
    backward = Dict{Int,Vector{Int}}()
    for (generator, image) in right.forward_images
        backward[generator] = substitute_words(image, left.backward_images)
    end
    PresentationIsomorphism(left.original, right.original, forward, backward)
end

"""Compare reductions; failure is unresolved, never a non-isomorphism claim.

If canonical finals differ, a bounded fallback is attempted only when their
final generator sets agree: every relator on each side must be a cyclic
normalization of a product of one or two signed cyclic conjugates from the
other side.  The returned witnesses replay that exact finite search.
"""
function compare_reductions(left::EliminationCertificate, right::EliminationCertificate)
    verify_certificate(left)
    verify_certificate(right)
    if _presentation_equal(left.final_presentation, right.final_presentation)
        iso = _identity_on_common_final_generators(left, right)
        return ReductionComparison(true, iso, :matched_canonical_output)
    end
    left_to_right = _bounded_relator_witnesses(left.final_presentation,
                                                right.final_presentation)
    left_to_right === nothing &&
        return ReductionComparison(false, nothing, :canonical_output_mismatch)
    right_to_left = _bounded_relator_witnesses(right.final_presentation,
                                                left.final_presentation)
    right_to_left === nothing &&
        return ReductionComparison(false, nothing, :canonical_output_mismatch)
    verify_witnesses(left.final_presentation, right.final_presentation, left_to_right)
    verify_witnesses(right.final_presentation, left.final_presentation, right_to_left)
    iso = _identity_on_common_final_generators(left, right)
    ReductionComparison(true, iso, :matched_bounded_relator_products,
                        (left_to_right=left_to_right,
                         right_to_left=right_to_left))
end

end
end
