# Hotspot Encryption Utility

This directory contains a standalone C program that parses hotspot sections from a `dog.txt` file and encrypts the hotspot data using the XChaCha20 stream cipher from [libsodium](https://libsodium.gitbook.io/doc/advanced/xchacha20). The workload is distributed across seven OpenMP threads so that multiple hotspot sections can be encrypted in parallel.

## Features
- Parses hotspot sections in the format:
  ```
  Hotspot Positions: 123,456,789
  Reference: ACTG...
  Alternate: ... (ignored)
  ```
- Generates an individual output file for every hotspot section. The file contains a 24-byte XChaCha20 nonce followed by the ciphertext that corresponds to the textual description of the hotspots.
- Uses libsodium's `crypto_stream_xchacha20_xor` for encryption.
- Encrypts sections concurrently using 7 OpenMP threads.

## Building

The program depends on both libsodium and OpenMP. Assuming the libraries are available on your system, build the utility with:

```sh
cd encryption
make
```

This produces an executable named `hotspot_encrypt` inside `encryption/bin`.

If libsodium is installed in a non-standard location, you may need to adjust the `LIBSODIUM_CFLAGS` and `LIBSODIUM_LDFLAGS` variables in the Makefile.

## Usage

Generate or obtain a 32-byte (256-bit) secret key and store it in hexadecimal form in a file. For example:

```sh
# Generate a random key using libsodium's `libsodium-utils`
openssl rand -hex 32 > secret.key
```

Run the program by specifying the input file, the key file, and the directory where you would like the encrypted blocks to be written:

```sh
./bin/hotspot_encrypt dog.txt secret.key encrypted_blocks
```

Each hotspot section produces an output file named `hotspot_<index>.bin` (indices start at 0). The file layout is:

- Bytes 0-23: XChaCha20 nonce
- Remaining bytes: Ciphertext of the hotspot description (`Hotspot Positions:` line + `Reference:` line).

The program will also emit a companion file `hotspot_<index>.meta` that repeats the plaintext metadata in human-readable form to help track which block corresponds to which ciphertext.

## Notes
- The OpenMP loop uses `num_threads(7)` explicitly to satisfy the request for seven threads. If fewer hardware threads are available, OpenMP may still schedule work appropriately.
- Any parsing anomalies (missing lines, malformed numbers, etc.) abort the program with an error message.
- The alternate sequence line is parsed but not encrypted. Adjust `build_plaintext_block` in `src/main.c` if you would like to include it in the ciphertext as well.

