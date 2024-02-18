#if two consecutive regions have the same annotation, then they should be merged into a single region

def merge_bed(file_path):
    merged_rows = []
    
    with open(file_path, 'r') as file:
        previous_row = None
        for line in file:
            row = line.strip().split('\t')
            #print(row)
            if previous_row is None:
                previous_row = row
                continue
            
            #both species and colors need to match for the coordinates to be merged
            if ((row[0] == previous_row[0]) and (row[3] == previous_row[3])):
                # Merge coordinates
                #print("found match")
                merged_row = [previous_row[0], previous_row[1], row[2], row[3]]
                previous_row = merged_row
            else:
                output=str('\t'.join(previous_row))
                outputfile.write(output + "\n")
                previous_row = row
        output=str('\t'.join(previous_row)) #the last entry
        outputfile.write(output + "\n") #the last entry


outputfile=open('circos.all.sequence.classes.merged.bed', 'w')
merge_bed('circos.all.sequence.classes.imputed.bed')
outputfile.close()
