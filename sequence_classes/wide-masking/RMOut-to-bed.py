#!/usr/bin/env python3

#Convert the RepeatMasker output to a zero-based UCSC track
#script provided by Gabby Hartley

import sys

with open(sys.argv[1],'r') as rpmsk:
    with open(sys.argv[2],'w') as bed:
        for line in rpmsk:
            line = line.replace("/","\t")
            spl=line.split()
            zeroStart = str(int(spl[5])-1)
            orient = spl[8]
            sequence = spl[4]
            end = spl[6]
            family = spl[9]
            divergence = spl[1]
            score = spl[0]
            if orient == 'C':
                orient == '-'
            if len(spl) == 15:
                if orient == 'C':
                    bed.write(sequence+'\t'+zeroStart+'\t'+end+'\t'+family+'\t'+score+'\t'+'-'+'\t'+spl[10]+'\t'+'undefined'+'\t'+divergence+'\t'+spl[14]+'\n')
                else:
                    bed.write(sequence+'\t'+zeroStart+'\t'+end+'\t'+family+'\t'+score+'\t'+'+'+'\t'+spl[10]+'\t'+'undefined'+'\t'+divergence+'\t'+spl[14]+'\n')
            if len(spl) == 16:
                if orient == 'C':
                    bed.write(sequence+'\t'+zeroStart+'\t'+end+'\t'+family+'\t'+score+'\t'+'-'+'\t'+spl[10]+'\t'+spl[11]+'\t'+divergence+'\t'+spl[15]+'\n')
                else:
                    bed.write(sequence+'\t'+zeroStart+'\t'+end+'\t'+family+'\t'+score+'\t'+'+'+'\t'+spl[10]+'\t'+spl[11]+'\t'+divergence+'\t'+spl[15]+'\n')
            if len(spl) == 17:
                if orient == 'C':
                    bed.write(sequence+'\t'+zeroStart+'\t'+end+'\t'+family+'\\'+spl[10]+'\t'+score+'\t'+'-'+'\t'+spl[11]+'\t'+spl[12]+'\t'+divergence+'\t'+spl[16]+'\n')
                else:
                    bed.write(sequence+'\t'+zeroStart+'\t'+end+'\t'+family+'\\'+spl[10]+'\t'+score+'\t'+'+'+'\t'+spl[11]+'\t'+spl[12]+'\t'+divergence+'\t'+spl[16]+'\n')
