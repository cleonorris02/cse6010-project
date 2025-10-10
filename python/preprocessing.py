import random
import os


# 1. Load and clean DNA sequence


def load_dna(file_path):
    with open(file_path, "r") as f:
        lines = f.readlines()
    # Remove FASTA headers (lines starting with '>') and newlines
    seq = "".join(line.strip() for line in lines if not line.startswith(">"))
    seq = seq.upper().replace("N", "")  # remove any ambiguous bases
    return seq


# 2. Simulate SNPs (single mutations)


def simulate_snps(sequence, snp_rate=0.01):
    """Introduce SNPs randomly for testing."""
    bases = ["A", "T", "G", "C"]
    snps = []
    seq_list = list(sequence)

    for i in range(len(seq_list)):
        if random.random() < snp_rate:
            original = seq_list[i]
            alt = random.choice([b for b in bases if b != original])
            snps.append((i, original, alt))
    return snps


# 3. Identify hotspot regions


def find_hotspots(snps, window_size=1000, threshold=35):
    """Find windows of the sequence with dense SNP activity."""
    hotspots = []
    for i in range(len(snps)):
        start = snps[i][0]
        window = [s for s in snps if start <= s[0] <= start + window_size]
        if len(window) >= threshold:
            hotspots.append(window)
    return hotspots


# 4. Main function


def main():
    dna_file = "data/dog.txt"
    output_folder = "processed_data"

    # Create processed_data folder if it doesn't exist
    os.makedirs(output_folder, exist_ok=True)

    # Load DNA
    dna_seq = load_dna(dna_file)
    print(f"Loaded DNA sequence length: {len(dna_seq)} bases.")

    # Simulate SNPs
    snps = simulate_snps(dna_seq, snp_rate=0.02)  # 2% mutation rate for demo
    print(f"Simulated {len(snps)} SNPs.")

    # Find hotspots
    hotspots = find_hotspots(snps)
    print(f"Found {len(hotspots)} hotspot regions.")

    # Prepare output file path
    input_filename = os.path.basename(dna_file)
    output_file = os.path.join(
        output_folder, f"processed_hotspots_{input_filename}.txt"
    )

    # Write output
    with open(output_file, "w") as f:
        f.write(f"Input File: {input_filename}\n\n")
        for region in hotspots:
            positions = [str(pos) for pos, _, _ in region]
            refs = "".join([ref for _, ref, _ in region])
            alts = "".join([alt for _, _, alt in region])
            f.write("Hotspot Positions: " + ",".join(positions) + "\n")
            f.write("Reference: " + refs + "\n")
            f.write("Alternate: " + alts + "\n\n")

    print(f"✅ Preprocessing complete. Results saved to {output_file}")


if __name__ == "__main__":
    main()
