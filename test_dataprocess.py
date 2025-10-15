"""
Simple test script to verify that dataprocess.py runs without errors
"""

import os
from dataprocess import (
    parse_fasta, parse_sequence_file, filter_by_protein_class,
    find_snps_with_max_variation, extract_unique_snp_contexts,
    identify_editing_hotspots
)

def test_dataprocess():
    print("Testing dataprocess.py functions...")
    
    # Create a simple test dataset
    sequences = [
        "ATGCATGCATGC",
        "ACGCATGCATGC",
        "ATGCATACATGC",
        "ATGCATGCGTGC",
        "ATGCATGCATTC"
    ]
    
    # Test find_snps_with_max_variation
    print("Testing find_snps_with_max_variation...")
    snp_positions = find_snps_with_max_variation(
        sequences, 
        min_variants=2, 
        min_coverage=3, 
        add_test_snps=False
    )
    print(f"Found {len(snp_positions)} SNP positions: {snp_positions}")
    
    # Test extract_unique_snp_contexts
    print("\nTesting extract_unique_snp_contexts...")
    unique_snps, contexts = extract_unique_snp_contexts(
        snp_positions, 
        sequences, 
        context_size=3
    )
    print(f"Found {len(unique_snps)} unique SNPs: {unique_snps}")
    
    # Test identify_editing_hotspots
    print("\nTesting identify_editing_hotspots...")
    hotspots = identify_editing_hotspots(
        unique_snps, 
        window_size=5, 
        min_snps=2
    )
    print(f"Found {len(hotspots)} hotspots")
    
    for i, hotspot in enumerate(hotspots):
        print(f"Hotspot {i+1}: positions {hotspot['start']}-{hotspot['end']} with {hotspot['snp_count']} SNPs")
    
    print("\nAll tests completed successfully!")

def test_with_actual_files():
    """Test using actual human.txt and chimpanzee.txt files from the dataset folder"""
    import os
    import pandas as pd
    
    print("\nTesting with actual human.txt and chimpanzee.txt files...")
    
    # Define file paths to the actual dataset files
    base_dir = os.path.dirname(os.path.abspath(__file__))
    dataset_dir = os.path.join(base_dir, "dna-sequence-dataset")
    human_file = os.path.join(dataset_dir, "human.txt")
    chimp_file = os.path.join(dataset_dir, "chimpanzee.txt")
    
    # Verify files exist
    if not os.path.exists(human_file):
        print(f"Error: {human_file} not found. Skipping test.")
        return
        
    if not os.path.exists(chimp_file):
        print(f"Error: {chimp_file} not found. Skipping test.")
        return
        
    print(f"Found sequence files in {dataset_dir}")
    
    # Import the functions we need
    from dataprocess import (
        parse_sequence_file, find_snps_with_max_variation, 
        extract_unique_snp_contexts, identify_editing_hotspots
    )
    
    try:
        # Parse sequence files
        print("\nParsing sequence files...")
        human_df = parse_sequence_file(human_file)
        chimp_df = parse_sequence_file(chimp_file)
        
        # Extract sequences from dataframes
        human_sequences = []
        chimp_sequences = []
        
        if 'sequence' in human_df.columns:
            # Format where each row is a complete sequence
            human_sequences = human_df['sequence'].tolist()
            print(f"Loaded {len(human_sequences)} human sequences")
        else:
            # Format where file contains nucleotide positions
            print("Human file contains nucleotide position data")
            
        if 'sequence' in chimp_df.columns:
            # Format where each row is a complete sequence
            chimp_sequences = chimp_df['sequence'].tolist()
            print(f"Loaded {len(chimp_sequences)} chimp sequences")
        else:
            # Format where file contains nucleotide positions
            print("Chimpanzee file contains nucleotide position data")
        
        # Combine all sequences for analysis
        all_sequences = human_sequences + chimp_sequences
        
        if not all_sequences:
            print("No complete sequences found. Using dataframes directly for analysis.")
            all_sequences = [human_df, chimp_df]
        
        # Find SNPs with maximum variation
        print("\nFinding SNPs with maximum variation...")
        snp_positions = find_snps_with_max_variation(all_sequences, min_variants=2, min_coverage=3)
        print(f"Found {len(snp_positions)} SNP positions: {snp_positions}")
        
        # Extract unique SNPs
        print("\nExtracting unique SNP contexts...")
        unique_snps, contexts = extract_unique_snp_contexts(snp_positions, all_sequences, context_size=5)
        print(f"Found {len(unique_snps)} unique SNPs: {unique_snps}")
        
        # Identify hotspots
        print("\nIdentifying editing hotspots...")
        hotspots = identify_editing_hotspots(unique_snps, window_size=10, min_snps=2)
        print(f"Found {len(hotspots)} hotspots")
        
        for i, hotspot in enumerate(hotspots):
            print(f"Hotspot {i+1}: positions {hotspot['start']}-{hotspot['end']} with {hotspot['snp_count']} SNPs")
        
        # Additional analysis for DNA steganography potential
        print("\nAnalyzing hotspots for DNA steganography potential...")
        for i, hotspot in enumerate(hotspots[:3]):  # Analyze top 3 hotspots
            snp_count = hotspot['snp_count']
            hotspot_len = hotspot['end'] - hotspot['start']
            density = snp_count / hotspot_len if hotspot_len > 0 else 0
            
            info_capacity = snp_count * 2  # 2 bits per SNP position
            
            print(f"Hotspot {i+1}: positions {hotspot['start']}-{hotspot['end']}")
            print(f"  - Information capacity: {info_capacity} bits ({info_capacity/8:.1f} bytes)")
            print(f"  - SNP density: {density:.3f} SNPs/bp")
            print(f"  - Steganography potential: {'High' if density > 0.3 else 'Medium' if density > 0.1 else 'Low'}")
        
        print("\nSuccessfully analyzed human and chimpanzee datasets!")
        
    except Exception as e:
        print(f"Error testing with real files: {str(e)}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    test_dataprocess()
    test_with_actual_files()
