# Outline and logic generated with ChatGPT (OpenAI), Oct 2025.         
#  Reviewed and modified by Viru Repalle.         


import requests
import gzip
import csv
from collections import defaultdict

url = "https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.40.gz"
hotspot_file = "snp_hotspots_strings.tsv"


def get_snps(url, limit=20000):
    """Stream and collect SNPs with base pair info (chrom, pos, REF, ALT)."""
    print(f"Streaming {url} ...")
    snps = []
    with requests.get(url, stream=True) as r:
        r.raise_for_status()
        with gzip.open(r.raw, "rt") as f:
            for line in f:
                if line.startswith("#"):
                    continue
                parts = line.strip().split("\t")
                if len(parts) < 5:
                    continue
                chrom, pos, ref, alt = parts[0], int(parts[1]), parts[3], parts[4]
                snps.append((chrom, pos, ref, alt))
                if len(snps) >= limit:
                    break
    print(f"✅ Collected {len(snps)} SNPs with base info.")
    return snps


def find_hotspots(snps, window_size=10000, threshold=15):
    """Find genomic windows with high SNP density and produce single DNA strings."""
    print("Analyzing SNP density ...")
    windows = defaultdict(list)
    for chrom, pos, ref, alt in snps:
        window = (chrom, pos // window_size)
        windows[window].append((ref, alt))

    hotspots = []
    for (chrom, start), base_list in windows.items():
        if len(base_list) >= threshold:
            # Produce a single continuous DNA string: REF + ALT concatenated
            base_string = "".join(
                [ref + "".join(alt.split(",")) for ref, alt in base_list]
            )
            hotspots.append(
                (
                    chrom,
                    start * window_size,
                    (start + 1) * window_size,
                    len(base_list),
                    base_string,
                )
            )

    hotspots.sort(key=lambda x: (x[0], x[1]))
    print(f"🔥 Found {len(hotspots)} SNP hotspots.")
    return hotspots


def save_hotspots_to_file(hotspots, output_file):
    """Write SNP hotspots as single DNA strings to TSV."""
    with open(output_file, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["Chromosome", "Start", "End", "SNP_Count", "DNA_String"])
        writer.writerows(hotspots)
    print(f"📂 Saved hotspots as DNA strings to {output_file}")


# --- Run the pipeline ---
snps = get_snps(url, limit=20000)  # adjust limit for larger sample
hotspots = find_hotspots(snps, window_size=10000, threshold=15)
save_hotspots_to_file(hotspots, hotspot_file)

# Optional: print first 3 hotspots for verification
for chrom, start, end, count, dna_str in hotspots[:3]:
    print(f"{chrom}:{start}-{end} → {count} SNPs, DNA string: {dna_str[:50]}...")
