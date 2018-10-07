# Genestring
Genestrings are like bit vectors, but provide access to anywhere from one to sixty-four bits at a time. They are designed for use with genetic programming projects that require a more compact representation of genomes than vectors of integers or vectors of floats.

Encoding and decoding logic must be built on top of a gene string, since all that is provided is a read/write interface and an interface for performing crossover on two strings. This should be sufficient for low level genetic programming work for AI code.
# Motivation
Some genomes can be more efficiently represented by using less than an entire word or floating point value per chromosome. A gene string allows you to fine tune the number of bits particular chromosomes have in resolution, potentially allowing you to fit more in to RAM or make more efficient use of crossover functions. This is similar to work in artificial neural networks, where lower resolution floating point numbers (eg. flexnets) have been used successfully.

Not every genome benefits from using a gene string. Bit shuffling is involved to decode values, so you may need to benchmark to see what your individual CPU/RAM tradeoffs are for a given project.

As an unintended but completely workable alternative, you can use gene strings to prepare network packets of tightly packed numbers.
# Tests
Testing is performed with the `proptest` crate on nightly builds. They should be quite thorough but defect reports are welcome.
# License
Genestrings are available under the `BSD-3 License`, provided under the `LICENSE.bsd` file in this repository.
