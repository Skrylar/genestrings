
pub const PIECE_SIZE_IN_BYTES: usize = 8;
pub const PIECE_SIZE_IN_BITS: usize  = 64;

#[derive(Default)]
pub struct Genestring {
    pieces: Vec<u64>,
}

// XXX i suppose our 'as usize' stuff might not be platform safe; 32-bit platforms might have issues?

impl Genestring {
    /// Creates a gene string capable of holding at least `count` bits.
    pub fn with_bits(count: u64) -> Genestring {
        let mut result = Genestring{
            pieces: Vec::with_capacity(count as usize),
        };
        result.pieces.resize(((count as usize / PIECE_SIZE_IN_BITS) + 1) as usize, 0);
        result
    }

    // Returns the number of bits in the gene string.
    pub fn bit_len(&self) -> usize {
        self.pieces.len() * PIECE_SIZE_IN_BITS
    }

    // Returns the number of bytes in the gene string.
    pub fn byte_len(&self) -> usize {
        self.pieces.len() * PIECE_SIZE_IN_BYTES
    }

    // Returns the number of integer parts of the gene string.
    pub fn len(&self) -> usize {
        self.pieces.len()
    }

    // Retrieves `bits` number of bits from the string, starting at a given `offset`. Panics if
    // `bits` is larger than 64 or would otherwise go outside the bounds of the string.
    pub fn get(&self, offset: usize, bits: usize) -> u64 {
        // safety dance
        if bits > 64 {
            panic!("Can only obtain 64 bits at a time!");
        }

        if bits + offset > self.bit_len() {
            panic!("Requested bits are out of bounds!");
        }

        // safety dance complete, now figure out which pieces have our bits
        let first_half_idx  = offset / PIECE_SIZE_IN_BYTES;
        let second_half_idx = (bits + offset) / PIECE_SIZE_IN_BYTES;

        let mut result: u64 = 0;

        if first_half_idx != second_half_idx {
            let first_half  = self.pieces[first_half_idx];
            let second_half = self.pieces[second_half_idx];

            let start = offset % PIECE_SIZE_IN_BITS;
            let stop  = PIECE_SIZE_IN_BITS;

            let piece = first_half;
            for i in start..stop {
                let mask = 1 << i;
                result <<= 1;
                result += (piece & mask) >> i;
            }

            let stop  = bits - (stop - start);;
            let start = 0;

            let piece = second_half;
            for i in start..stop {
                let mask = 1 << i;
                result <<= 1;
                result += (piece & mask) >> i;
            }
        } else {
            // in this path, all requested bits are within a single piece
            let piece = self.pieces[first_half_idx];

            for i in offset..(offset+bits) {
                let mask = 1 << i;
                result <<= 1;
                result += (piece & mask) >> i;
            }
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
    pub fn set(&mut self, offset: usize, bits: usize, value: usize) {
        unimplemented!();
    }

    // Copies bits starting from a given offset, up to offset+bits, from a donor genestring to this one.
    // The assumed usage of this function is to implement crossover between generations.
    pub fn transplant(&mut self, donor: &mut Genestring, offset: usize, bits: usize) {
        unimplemented!();
    }
}

#[cfg(test)]
mod tests {
    use ::*;
    #[test]
    fn calculating_string_size() {
        // just making sure this bit of math works as we expect it to
        assert_eq!(((7 as usize / PIECE_SIZE_IN_BITS) + 1), 1);
        assert_eq!(((70 as usize / PIECE_SIZE_IN_BITS) + 1), 2);
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
}
