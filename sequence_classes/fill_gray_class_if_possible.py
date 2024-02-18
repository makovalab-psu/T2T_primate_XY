#if an unannotated region neighbors two regions with the same annotation, then it should also have that annotation

def check_previous_next(file_path):
    lines = []
    with open(file_path, 'r') as file:
        lines = file.readlines()

    outputfile.write(lines[0])
    for i in range(1, len(lines) - 1):
        #extract the fourth column
        #print (i)
        current_line = lines[i].strip().split('\t')[3]
        previous_line = lines[i-1].strip().split('\t')[3]
        next_line = lines[i+1].strip().split('\t')[3]

        # (#d9d8d8==Gray==Other)
        if (previous_line == next_line) and (current_line=="fill_color=#d9d8d8"):
            #print(lines[i-1])
            #print(lines[i])
            #print(lines[i+1])
            new_color=previous_line #the new color
            start=lines[i].strip().split('\t')[0:3]
            new_entry=str('\t'.join(start) + "\t" + new_color)
            #print(new_entry)
            #print the current line that was previously gray
            outputfile.write(new_entry + "\n")
            #print("\n")
        else:
            outputfile.write(lines[i])
    outputfile.write(lines[len(lines)-1]) #write the last line

outputfile=open('circos.all.sequence.classes.imputed.bed', 'w')
check_previous_next('circos.all.sequence.classes.bed')
outputfile.close()
