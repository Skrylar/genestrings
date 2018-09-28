use constants::*;

/// Calculates the number of pieces required to hold this many bits in a genestring.
pub fn part_count_for_bits(bits: u64) -> u64 {
    if bits == 0 {
        1
    } else if bits % PIECE_SIZE_IN_BITS == 0 {
        bits / PIECE_SIZE_IN_BITS
    } else {
        (bits / PIECE_SIZE_IN_BITS) + 1
    }
}
