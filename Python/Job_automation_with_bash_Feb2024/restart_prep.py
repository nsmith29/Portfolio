#!/usr/bin/env python

restart = "A.restart"

rewrite_ = False

counter = -1
lines = []

with open(restart,'r') as file_:
    for line in file_:
        counter += 1
        lines.append(line)

        if "WFN_RESTART_FILE_NAME" in line or "COORD_FILE_FORMAT" in line or "COORD_FILE_NAME" in line:
            lin_num = counter
            counter -=1
            removed_ = lines.pop(lin_num)
            print("line removed:",removed_)
        
        elif "SCF_GUESS  ATOMIC" in line:
            line_number = counter
            rewrite_ = True
            final_line = "       SCF_GUESS  RESTART \n"
            print("line",line_number,"was :",lines[line_number])
            lines[line_number] = final_line
            print("is now",lines[line_number])
            
    if rewrite_ is True:
        with open(restart, 'w') as restrt:
            for line_ in lines:
                restrt.write(line_)       
            
            