#[macro_use] extern crate proptest;

pub mod constants;
pub mod math;

use constants::*;
use math::*;

const PANIC_OUT_OF_BOUNDS: &str = "Requested bits are out of bounds!";

#[derive(Default, Debug)]
pub struct Genestring {
    pieces: Vec<u64>,
}

impl Genestring {
    /// Creates a gene string capable of holding at least `count` bits.
    pub fn with_bits(count: u64) -> Genestring {
        let mut result = Genestring{
            pieces: Vec::with_capacity(count as usize),
        };
        result.pieces.resize(part_count_for_bits(count) as usize, 0);
        result
    }

    // Returns the number of bits in the gene string.
    pub fn bit_len(&self) -> usize {
        self.pieces.len() * PIECE_SIZE_IN_BITS as usize
    }

    // Returns the number of bytes in the gene string.
    pub fn byte_len(&self) -> usize {
        self.pieces.len() * PIECE_SIZE_IN_BYTES as usize
    }

    // Returns the number of integer parts of the gene string.
    pub fn len(&self) -> usize {
        self.pieces.len()
    }

    // Retrieves `bits` number of bits from the string, starting at a given `offset`. Panics if
    // `bits` is larger than 64 or would otherwise go outside the bounds of the string.
    pub fn get(&self, offset: u64, bits: u64) -> u64 {
        if bits == 0 {
            return 0
        }
        
        // safety dance
        if bits > 64 {
            panic!("Can only obtain 64 bits at a time!");
        }

        if bits + offset > self.bit_len() as u64 {
            panic!(PANIC_OUT_OF_BOUNDS);
        }

        // safety dance complete, now figure out which pieces have our bits
        let first_half_idx  = (offset / PIECE_SIZE_IN_BITS) as usize;
        let second_half_idx = ((bits + offset) / PIECE_SIZE_IN_BITS) as usize;

        let mut result: u64 = 0;
        let mut result_b: u64 = 0;

        if first_half_idx != second_half_idx {
            let first_half  = self.pieces[first_half_idx];
            let second_half = self.pieces[second_half_idx];

            let start = offset % PIECE_SIZE_IN_BITS;
            let stop  = PIECE_SIZE_IN_BITS;

            let piece = first_half;
            for i in start..stop {
                let mask = 1 << i;
                result += piece & mask;
            }

            result >>= start;

            let stop  = bits - (stop - start);
            let start = 0;

            let piece = second_half;
            for i in start..stop {
                let mask = 1 << i;
                result_b += piece & mask;
            }

            result_b <<= bits - stop;

            result |= result_b;
        } else {
            let first_half  = self.pieces[first_half_idx];

            let start = offset % PIECE_SIZE_IN_BITS;
            let stop  = start + bits;

            let piece = first_half;
            for i in start..stop {
                let mask = 1 << i;
                result += piece & mask;
            }

            result >>= start;
        }

        result
    }

    // Fills each piece of the genestring from a supplied fill function.
    // The assumed usage of this function is for inserting random values for new DNA.
    pub fn fill<F>(&mut self, mut filler: F)
        where F: FnMut() -> u64
    {
        for i in self.pieces.iter_mut() {
            *i = filler();
        }
    }

    // Assigns bits at a given offset through offset+bits to the given value.
    // The assumed usage of this function is to implement mutation.
    pub fn set(&mut self, offset: u64, bits: u64, value: u64) {
        if bits == 0 {
            return
        }

        // safety dance
        if bits > 64 {
            panic!("Can only write 64 bits at a time!");
        }

        if bits + offset > self.bit_len() as u64 {
            panic!(PANIC_OUT_OF_BOUNDS);
        }

        let first_half_idx  = part_for_bit(offset) as usize;
        let second_half_idx = part_for_bit(offset + bits) as usize;

        let mut source_mask = 0;

        let offset_modulo = offset % PIECE_SIZE_IN_BITS;

        if first_half_idx == second_half_idx {
            // in this path, we are just writing to bits inside the same integer
            for i in 0..bits {
                source_mask += 1 << i;
            }

            let destination_mask = source_mask << offset_modulo;

            self.pieces[first_half_idx] = (self.pieces[first_half_idx] & !destination_mask) + ((value as u64 & source_mask) << offset_modulo);
        } else {
            let first_half  = self.pieces[first_half_idx];
            let second_half = self.pieces[second_half_idx];

            let start = offset % PIECE_SIZE_IN_BITS;
            let stop  = PIECE_SIZE_IN_BITS;

            let piece = first_half;
            let mut piece_mask = 0;
            let mut value_mask = 0;

            for i in start..stop {
                piece_mask += 1 << i;
                value_mask = (value_mask << 1) + 1;
            }

            self.pieces[first_half_idx] = (piece & !piece_mask) | ((value as u64 & value_mask) << start);

            let stop  = bits - (stop - start);
            let start = 0;

            let piece = second_half;
            piece_mask = 0;
            value_mask = 0;
            for i in start..stop {
                piece_mask += 1 << i;
                value_mask = (value_mask << 1) + 1;
            }
            value_mask <<= bits - stop;

            self.pieces[second_half_idx] = (piece & !piece_mask) | ((value as u64 & value_mask) >> (bits - stop));
        }
    }

    // Copies bits starting from a given offset, up to offset+bits, from a donor genestring to this one.
    // The assumed usage of this function is to implement crossover between generations.
    pub fn transplant(&mut self, donor: &Genestring, offset: u64, bits: u64) {
        let end = bits + offset;

        if end > self.bit_len() as u64 || end > donor.bit_len() as u64 {
            panic!(PANIC_OUT_OF_BOUNDS);
        }

        let first_half_idx  = part_for_bit(offset) as usize;
        let second_half_idx = part_for_bit(offset + bits) as usize;

        if first_half_idx == second_half_idx {
            self.set(offset, bits, donor.get(offset, bits));
        } else {
            // middle sections can just be copied directly
            for i in (first_half_idx+1)..(second_half_idx-1) {
                self.pieces[i as usize] = donor.pieces[i as usize];
            }

            // copy first section
            let bits = PIECE_SIZE_IN_BITS - (offset % PIECE_SIZE_IN_BITS);
            self.set(offset, bits, donor.get(offset, bits));

            // copy end section
            let offset = second_half_idx * PIECE_SIZE_IN_BITS;
            let bits = offset - ((bits + offset) % PIECE_SIZE_IN_BITS);
            self.set(offset, bits, donor.get(offset, bits));
        }
    }
}

#[cfg(test)]
mod tests {
    use ::*;

    #[test]
    fn assume_intdiv_rounds_down() {
        assert_eq!(4/5, 0);
        assert_eq!(7/5, 1);
    }

    #[test]
    fn calculating_bi_offsets() {
        // just making sure the way we do bit offsets is correct
        let offset = 50;
        let bits = 32;
        let mut total = 0;

        let start = offset % PIECE_SIZE_IN_BITS;
        let stop  = PIECE_SIZE_IN_BITS;

        for _ in start..stop {
            total += 1;
        }

        let stop  = bits - (stop - start);
        let start = 0;

        for _ in start..stop {
            total += 1;
        }

        assert_eq!(total, bits);
    }

    // These two tests are very basic idiot tests, but are no means exhaustive.

    #[test]
    fn get_set_same_chunk() {
        assert_eq!(PIECE_SIZE_IN_BITS, 64);
        let mut gs = Genestring::with_bits(32);

        eprintln!("{:?}", gs);
        gs.set(8, 12, 842);
        eprintln!("{:?}", gs);
        assert_eq!(gs.get(8, 12), 842);
    }

    #[test]
    fn get_set_different_chunk() {
        assert_eq!(PIECE_SIZE_IN_BITS, 64);
        let mut gs = Genestring::with_bits(128);

        eprintln!("{:?}", gs);
        gs.set(60, 8, 0xFF);
        eprintln!("{:?}", gs);
        assert_eq!(gs.pieces[0], 0xF000000000000000);
        assert_eq!(gs.pieces[1], 0x000000000000000F);
        assert_eq!(gs.get(60, 8), 0xFF);
    }

    #[test]
    fn string_size_minimum() {
        // just making sure this bit of math works as we expect it to
        assert_eq!(PIECE_SIZE_IN_BITS, 64);
        assert_eq!(part_count_for_bits(0), 1);
    }

    // proptest does some more intensive checks to ensure things like split numbers always work
    // or we don't trample non-overlapping numbers doing arithmetic.

    proptest! {
        #[test]
        fn string_size_blocks(blocks in 1..10) {
            // just making sure this bit of math works as we expect it to
            assert_eq!(PIECE_SIZE_IN_BITS, 64);
            assert_eq!(part_count_for_bits(blocks as u64 * PIECE_SIZE_IN_BITS), blocks as u64);
        }

        #[test]
        fn string_size_subblocks(blocks in 1..10, subblock in 1..32) {
            // just making sure this bit of math works as we expect it to
            assert_eq!(PIECE_SIZE_IN_BITS, 64);
            assert_eq!(part_count_for_bits((blocks as u64 * PIECE_SIZE_IN_BITS) + subblock as u64), blocks as u64 + 1);
        }

        #[test]
        fn get_set_single(start in 0..192, len in 1..64, value: u64) {
            // we're going to get and set values at various offsets and make sure we always get
            // back the thing we wanted to start with

            assert_eq!(PIECE_SIZE_IN_BITS, 64);
            let mut gs = Genestring::with_bits(256);

            let mut band: u64 = 0;
            for i in 0..len {
                band += 1 << i;
            }

            let banded_value = value as u64 & band;
            gs.set(start as u64, len as u64, banded_value);
            prop_assert_eq!(gs.get(start as u64, len as u64), banded_value);
        }

        #[test]
        fn get_set_1piece_duo(a in 0..16, b in 32..64, value_a: u16, value_b: u16) {
            // We are going to store two values within the same piece, guaranteed not to overlap,
            // and ensure they do not trample one another.

            assert_eq!(PIECE_SIZE_IN_BITS, 64);
            let mut gs = Genestring::with_bits(64);

            gs.set(a as u64, 16, value_a as u64);
            gs.set(b as u64, 16, value_b as u64);

            prop_assert_eq!(gs.get(a as u64, 16), value_a as u64);
            prop_assert_eq!(gs.get(b as u64, 16), value_b as u64);
        }

        #[test]
        fn get_set_multibinning(a in 0..16, b in 32..100, value_a: u16, value_b: u16) {
            // We have one value which is always in the first piece, and a second value which
            // can span any non-overlap location in either piece. Ensures our single and double
            // piece logics don't conflict.

            assert_eq!(PIECE_SIZE_IN_BITS, 64);
            let mut gs = Genestring::with_bits(128);

            gs.set(a as u64, 16, value_a as u64);
            gs.set(b as u64, 16, value_b as u64);

            prop_assert_eq!(gs.get(a as u64, 16), value_a as u64);
            prop_assert_eq!(gs.get(b as u64, 16), value_b as u64);
        }

        #[test]
        fn transplanting_small_ranges(a in 0..16, b in 32..100, value_a: u16, value_b: u16) {
            assert_eq!(PIECE_SIZE_IN_BITS, 64);
            let mut gs  = Genestring::with_bits(128);
            let mut gs2 = Genestring::with_bits(128);

            gs.set(a as u64, 16, value_a as u64);
            gs.set(b as u64, 16, value_b as u64);

            gs2.transplant(&gs, a as u64, 16);
            gs2.transplant(&gs, b as u64, 16);

            prop_assert_eq!(gs2.get(a as u64, 16), value_a as u64);
            prop_assert_eq!(gs2.get(b as u64, 16), value_b as u64);
        }
    }
}
