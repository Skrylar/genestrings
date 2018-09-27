#[macro_use] extern crate proptest;

pub const PIECE_SIZE_IN_BYTES: u64 = 8;
pub const PIECE_SIZE_IN_BITS: u64  = 64;

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
        result.pieces.resize(((count / PIECE_SIZE_IN_BITS) + 1) as usize, 0);
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
        // safety dance
        if bits > 64 {
            panic!("Can only write 64 bits at a time!");
        }

        if bits + offset > self.bit_len() as u64 {
            panic!(PANIC_OUT_OF_BOUNDS);
        }

        let first_half_idx  = (offset / PIECE_SIZE_IN_BITS) as usize;
        let second_half_idx = ((bits + offset) / PIECE_SIZE_IN_BITS) as usize;

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
    pub fn transplant(&mut self, donor: &mut Genestring, offset: u64, bits: u64) {
        let end = bits + offset;

        if end > self.bit_len() as u64 || end > donor.bit_len() as u64 {
            panic!(PANIC_OUT_OF_BOUNDS);
        }

        let first_half_idx  = offset / PIECE_SIZE_IN_BITS;
        let second_half_idx = (bits + offset) / PIECE_SIZE_IN_BITS;

        unimplemented!();
    }
}

#[cfg(test)]
mod tests {
    use ::*;
    #[test]
    fn calculating_string_size() {
        // just making sure this bit of math works as we expect it to
        assert_eq!(((7 / PIECE_SIZE_IN_BITS) + 1), 1);
        assert_eq!(((70 / PIECE_SIZE_IN_BITS) + 1), 2);
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

    proptest! {
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
    }
}
