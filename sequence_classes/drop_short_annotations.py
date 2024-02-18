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

            previous_distance = int(previous_row[2]) - int(previous_row[1])
            distance = int(row[2]) - int(row[1])
            previous_color = previous_row[3]
            color = row[3]
            
            #both species and colors need to match for the coordinates to be merged
            if ((row[0] == previous_row[0]) and ((previous_distance<=100) or (distance<=100))):
                # Merge coordinates
                #print("found match, should be merged")
                #print(previous_row)
                #print(row)
                final_color="fill_color=#d9d8d8"

                # (#d9d8d8==Gray==Other)
                if ((previous_color!="fill_color=#d9d8d8") and (color!="fill_color=#d9d8d8")):
                    #both of them are color other than gray, attach to the longer interval 
                    if (previous_distance>distance):
                        #previous distance is bigger, merge distance to previous distance
                        final_color=previous_color
                    else:
                        #distance is bigger, merge previous distance to distance
                        final_color=color

                else: 
                    #only one of them is color other than gray, choose the one that is not gray (#d9d8d8==Gray==Other)
                    if (previous_color=="fill_color=#d9d8d8"):
                        final_color=color
                    else:
                        final_color=previous_color


                merged_row = [previous_row[0], previous_row[1], row[2], final_color]
                #print(merged_row)
                #print("====")
                previous_row = merged_row
            else:
                output=str('\t'.join(previous_row))
                outputfile.write(output + "\n")
                previous_row = row
        output=str('\t'.join(previous_row)) #the last entry
        outputfile.write(output + "\n") #the last entry


outputfile=open('circos.all.sequence.classes.final.bed', 'w')
merge_bed('circos.all.sequence.classes.merged.bed')
outputfile.close()
