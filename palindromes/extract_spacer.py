import sys
from pybedtools import BedTool

#using bed file with palindrome coordinates from palindrover, identify the coordinates of spacers

if len(sys.argv) != 2:
    print("Usage: python script_name.py input_bed_file.bed")
    sys.exit(1)

bed_file_path = sys.argv[1]
bed_entries = BedTool(bed_file_path)

bed_entry_list = list(bed_entries)

output_bed_file = open(bed_file_path+'spacers.bed', 'w')  # Open the output BED file for writing

for i in range(0, len(bed_entry_list) - 1, 2):
    entry1 = bed_entry_list[i]
    entry2 = bed_entry_list[i + 1]
    distance = abs(entry1.end - entry2.start)
    print(f"{entry1.name}\t{distance}")
    output_bed_file.write(str(entry1.chrom) + '\t' + str(entry1.end) + '\t' + str(entry2.start) + '\t' + str(entry1.name) + '\n')
