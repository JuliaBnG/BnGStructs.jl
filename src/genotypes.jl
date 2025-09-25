"""
    mutable struct Genotype
A struct to store SNP genotypes in a BitMatrix, with individuals as rows and
alleles as columns. The 
"""
mutable struct Genotype
    nid::Int   # number of individuals
    nas::Int   # number of alleles
    gt::BitMatrix
    function Genotype(nid::Int, nas::Int, gt::BitMatrix)
        nid > 0 || error("nid must be positive")
        (nas > 0 && iseven(nas)) || error("nas must be positive and even")
        size(gt, 1) == 64 * cld(nid, 64) || error("Number of rows of gt does not match nid")
        size(gt, 2) == nas || error("Number of columns of gt does not match nas")
        # zero out padded rows beyond nid
        _mask_last_word_columns!(gt, nid)
        new(nid, nas, gt)
    end
end

"""
    Genotype(nid::Int, nas::Int)
Create a Genotype struct with `nid` individuals and `nas` alleles. The
corresponding BitMatrix is initialized to all zeros.
"""
function Genotype(nid::Int, nas::Int)
    (nas > 0 && iseven(nas)) || error("nas must be positive even")
    nid > 0 || error("nid must be positive")
    gt = falses(64 * cld(nid, 64), nas) # as a chunk in BitMatrix is 64-bit UInt64
    return Genotype(nid, nas, gt)
end
