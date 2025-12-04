# XChaCha20 Verification Tests

This directory contains tests to verify that the XChaCha20 encryption algorithm is working correctly.

## Building and Running Tests

```bash
make test
```

## What the Tests Verify

The test suite validates the following aspects of the XChaCha20 implementation:

### 1. Encryption/Decryption Round-trip
Verifies that encrypting data and then decrypting it with the same key and nonce returns the original plaintext. This is the fundamental property of symmetric encryption.

### 2. Binary-to-DNA Encoding
Tests the encoding scheme that converts binary ciphertext to DNA nucleotide sequences:
- `00` → `A`
- `01` → `C`  
- `10` → `G`
- `11` → `T`

Each byte (8 bits) becomes 4 nucleotides.

### 3. Deterministic Encryption
Confirms that using the same key and nonce always produces the same ciphertext. This is essential for reproducibility and debugging.

### 4. Full Pipeline Test
Tests the complete workflow:
1. Encrypt DNA sequence with XChaCha20
2. Encode ciphertext as DNA nucleotides
3. Decode DNA back to binary
4. Decrypt to recover original sequence

### 5. Nonce Uniqueness
Verifies that different nonces produce different ciphertext for the same plaintext, which is critical for semantic security.

## Expected Output

A successful test run should show:
```
*** ALL TESTS PASSED ***

This verifies that:
1. XChaCha20 encryption/decryption works correctly (round-trip)
2. Binary-to-DNA encoding is reversible and correct
3. The encryption is deterministic with fixed key/nonce
4. The full pipeline (encrypt -> DNA encode -> DNA decode -> decrypt) works
5. Different nonces produce different ciphertext (semantic security)
```

## Dependencies

- libsodium (for XChaCha20 implementation)
- Standard C library

Install libsodium:
```bash
# Ubuntu/Debian
sudo apt-get install libsodium-dev

# macOS
brew install libsodium
```
