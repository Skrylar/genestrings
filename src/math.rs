use constants::*;

pub fn part_count_for_bits(bits: u64) -> u64 {
    (bits / PIECE_SIZE_IN_BITS) + 1
}
