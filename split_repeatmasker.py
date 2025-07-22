#!/usr/bin/env python3
"""
Split RepeatMasker annotations into L1, SVA, ALU, and HERV BED files
"""

import sys
import os

def categorize_repeat(repeat_class):
    """
    Categorize repeat based on RepeatMasker classification
    """
    original_class = repeat_class
    repeat_class = repeat_class.upper()
    
    # SVA elements
    if repeat_class == 'RETROPOSON/SVA':
        return 'SVA'
    
    # ALU elements
    if repeat_class == 'SINE/ALU':
        return 'ALU'

    # MIR elements
    if repeat_class == 'SINE/MIR':
        return 'MIR'
    
    # HERV elements (LTR retrotransposons)
    if repeat_class in ['LTR/ERV1', 'LTR/ERV1?', 'LTR/ERVK', 'LTR/ERVL', 'LTR/ERVL-MaLR', 'LTR/ERVL?']:
        return 'HERV'
    
    # L1 elements (LINE-1)
    if repeat_class == 'LINE/L1':
        return 'L1'

    # L2 elements (LINE-2)
    if repeat_class == 'LINE/L2':
        return 'L2'
    
    # Return None for other repeat types
    return None

def split_repeatmasker_bed(input_file, output_prefix):
    """
    Split RepeatMasker BED file into separate files for L1, SVA, ALU, and HERV
    """
    
    # Open output files
    output_files = {}
    for category in ['L1', 'L2', 'SVA', 'MIR', 'ALU', 'HERV']:
        filename = f"{output_prefix}_{category.lower()}.bed"
        output_files[category] = open(filename, 'w')
        print(f"Created output file: {filename}")
    
    # Counters for statistics
    counters = {'L1': 0, 'L2': 0, 'SVA': 0, 'MIR': 0, 'ALU': 0, 'HERV': 0, 'OTHER': 0}
    
    try:
        with open(input_file, 'r') as infile:
            for line_num, line in enumerate(infile, 1):
                line = line.strip()
                
                # Skip empty lines and comments
                if not line or line.startswith('#'):
                    continue
                
                # Split the line
                fields = line.split('\t')
                
                # Expect RepeatMasker format with at least 13 columns
                if len(fields) < 13:
                    print(f"Warning: Line {line_num} has fewer than 13 columns, skipping: {line}")
                    continue
                
                # RepeatMasker format: columns 12 and 13 contain class and family
                repeat_class = fields[11]  # e.g., "LINE", "SINE", "LTR"
                repeat_family = fields[12]  # e.g., "L1", "Alu", "ERVL-MaLR"
                repeat_full = f"{repeat_class}/{repeat_family}"  # e.g., "LINE/L1"
                category = categorize_repeat(repeat_full)
                
                if category:
                    # Create BED format: chr, start, end, repeat_info
                    chr_name = fields[5]  # chromosome
                    start = fields[6]     # start position
                    end = fields[7]       # end position
                    bed_line = f"{chr_name}\t{start}\t{end}\t{repeat_full}"
                    output_files[category].write(bed_line + '\n')
                    counters[category] += 1
                else:
                    counters['OTHER'] += 1
                    if counters['OTHER'] <= 10:  # Show first 10 examples of uncategorized entries
                        print(f"Uncategorized repeat: {repeat_full}")
    
    finally:
        # Close all output files
        for file_handle in output_files.values():
            file_handle.close()
    
    # Print statistics
    print("\nSplit statistics:")
    for category, count in counters.items():
        print(f"{category}: {count:,} entries")
    
    total_categorized = sum(counters[cat] for cat in ['L1', 'L2', 'SVA', 'MIR', 'ALU', 'HERV'])
    total_entries = total_categorized + counters['OTHER']
    print(f"Total categorized: {total_categorized:,} / {total_entries:,} ({100*total_categorized/total_entries:.1f}%)")

def main():
    if len(sys.argv) != 3:
        print("Usage: python split_repeatmasker.py <input_bed_file> <output_prefix>")
        print("Example: python split_repeatmasker.py repeatmasker.bed hg38_repeats")
        print("This will create: hg38_repeats_l1.bed, hg38_repeats_sva.bed, etc.")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_prefix = sys.argv[2]
    
    if not os.path.exists(input_file):
        print(f"Error: Input file '{input_file}' not found")
        sys.exit(1)
    
    print(f"Splitting RepeatMasker file: {input_file}")
    print(f"Output prefix: {output_prefix}")
    
    split_repeatmasker_bed(input_file, output_prefix)
    print("Done!")

if __name__ == "__main__":
    main()