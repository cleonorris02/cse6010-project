# CSE 6010 Group 2 Project

This project explores a novel approach to genomic data security by combining DNA sequence representation with modern cryptographic methods. We build upon the methodology introduced in Dokyun Na’s paper “DNA steganography: hiding undetectable secret messages within the single nucleotide polymorphisms of a genome and detecting mutation-induced errors”, published in Microbial Cell Factories, by adapting the framework to incorporate the XChaCha20 encryption algorithm in place of the simple substitution cipher originally used. XChaCha20 offers significant advantages in speed, scalability, and resilience against timing attacks, making it a compelling candidate for protecting sensitive biological information. Using a publicly available DNA sequence dataset from Kaggle, our system encodes encrypted data into genomic sequences and evaluates whether the efficiency and robustness of XChaCha20 can be effectively applied in this context. The goal of this project is to assess the feasibility of stronger encryption techniques within DNA-based data hiding, thereby extending the state of the art in secure storage and transmission of biological information. 

## Error Detection/Correction Module

A C implementation of the block sum parity scheme described in Na et al. (2020) is available in `c/`. The module constructs a two-dimensional DNA block augmented with parity nucleotides for each row and column, enabling detection and correction of single-nucleotide mutations in either the data or parity regions.

### Building and running the demo

```
cd c
make
./error_detection
```

The demonstration program encrypts three DNA words into a parity-protected block, introduces sample mutations (data, row parity, and column parity), and invokes the recovery routine to restore the correct nucleotides.
