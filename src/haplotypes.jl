"""
    mutable struct Haplotype
A struct to hold haplotype data in a compact BitMatrix form.
"""
mutable struct Haplotype
    nlc::Int   # number of loci (rows logically = nlc, stored rows = 64 * cld(nlc, 64))
    nhp::Int   # number of haplotypes (columns in gt)
    gt::BitMatrix
    function Haplotype(nlc::Int, nhp::Int, gt::BitMatrix)
        nlc > 0 || error("nlc must be positive")
        (nhp > 0 && iseven(nhp)) || error("nhp must be positive and even")
        size(gt, 1) == 64 * cld(nlc, 64) || error("Number of rows of gt does not match nlc")
        size(gt, 2) == nhp || error("Number of columns of gt does not match nhp")
        # zero out padded rows beyond nlc
        _mask_last_word_columns!(gt, nlc)
        new(nlc, nhp, gt)
    end
end

"""
    Haplotype(nlc::Int, nhp::Int)
Create a Haplotype struct with `nlc` loci and `nhp` haplotypes. The
corresponding BitMatrix is initialized to all zeros.
"""
function Haplotype(nlc::Int, nhp::Int)
    nlc > 0 || error("nlc must be positive")
    (nhp > 0 && iseven(nhp)) || error("nhp must be positive and even")
    gt = falses(64 * cld(nlc, 64), nhp) # as a chunk in BitMatrix is 64-bit UInt64
    return Haplotype(nlc, nhp, gt)
end
