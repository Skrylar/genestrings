Doing genetic programming work in Rust? Need a more compact representation of your genome than vectors of integers? Enter the Genestring.

# Genestring

Genestrings are a container that allow you to define a set number of bits, then slice in to those bits like they were plain old integers all along. This means your moderator genes can truly be one bit, and the level of precision needed for each field can be adjusted to only what is absolutely necessary.

When your genetic fuzzy trees or cartesian genetic programming starts using up too much memory, consider migrating to genestrings for a more compact representation of bits.

## Phenotypes

Genestrings know only of bits and ranges of bits, they do not know anything of your underlying genotypes or phenotypes. You will have to mantain those abstractions above a gene string, should you require things like phenotype-aware mutation functions.

# Tested? Tested.

Genestrings are thoroughly tested with the `proptest` framework, using nightly Rust compilers.

# Licensed for anything
Genestrings are available under the `BSD-3 License`, provided under the `LICENSE.bsd` file in this repository.
