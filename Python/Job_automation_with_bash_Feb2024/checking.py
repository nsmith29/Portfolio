#!/usr/bin/env python

intermediate  = "A.inter.restart"
xyz_ = "A-pos-1.xyz"
xyz__ = "A-pos-A_test.xyz"

MAX_at_500 = False
STEP_VAL_F = False
rewrite_ = False
lines = []

counter = -1

with open(intermediate,'r') as file_:
    for line in file_:
        counter += 1
        lines.append(line)
            
        if "MAX_ITER " in line:
            line_number = counter
            first_part = str(line.split()[0])
            var_iter = str(line.split()[1])
            if var_iter == "500":
                MAX_at_500 = True
                print("Max iterations set to 500")
            else:
                print("more than 500 iterations set")

        elif MAX_at_500 is True and "STEP_START_VAL" in line and STEP_VAL_F is False:
            STEP_VAL_F = None
            var = str(line.split()[1])
            if int(var) >= 400:
                print("Step value passesd 400 with value of",var)
                rewrite_ = True
                final_line = "     " + str(first_part) + "  800" + "\n"
                
                print("line",line_number,"was :",lines[line_number])
                lines[line_number] = final_line
                print("is now",lines[line_number])
                
                
                with open(xyz_, 'r') as xyz1:
                    xyz_lines = xyz1.readlines()
                    with open(xyz__,'w') as xyz2:
                        xyz2.writelines(xyz_lines)
            else:
                print("Step value below 400")                
                
    if rewrite_ is True:
        with open(intermediate, 'w') as restrt:
            for line_ in lines:
                restrt.write(line_)