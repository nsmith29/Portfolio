#!/usr/bin/env python

inp = "A.inp"

rewrite_ = False

counter = -1
lines = []

with open(inp,'r') as file_:
    for line in file_:
        counter += 1
        lines.append(line)
        
        if "SCF_GUESS  ATOMIC" in line:
            line_number = counter
            rewrite_ = True
            final_line = "       SCF_GUESS  RESTART \n"
            print("line",line_number,"was :",lines[line_number])
            lines[line_number] = final_line
            print("is now",lines[line_number])
            
    if rewrite_ is True:
        with open(inp, 'w') as restrt:
            for line_ in lines:
                restrt.write(line_)  