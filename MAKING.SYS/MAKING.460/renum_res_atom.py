inp1 = open("mod46_renum.pdb", "r")

atom = []
atom_num = []
atom_name = []
chain_1 = []
chain_2 = []
res_num = []
x_cord = []
y_cord = []
z_cord = []
occu = []
b_fact = []
atom_type = []

for lines in inp1:
    if lines.startswith("ATOM"):
        arr = lines.split()
        atom.append(str(arr[0]))
        atom_num.append(int(arr[1]))
        atom_name.append(str(arr[2]))
        chain_1.append(str(arr[3]))
        chain_2.append(str(arr[4]))
        res_num.append(int(arr[5]))
        x_cord.append(float(arr[6]))
        y_cord.append(float(arr[7]))
        z_cord.append(float(arr[8]))
        occu.append(float(arr[9]))
        b_fact.append(float(arr[10]))
        #atom_type.append(str(arr[11]))
        
new_at_num = []
new_res_num = []
out2 = open("mod46_renum2.pdb", "w")

for j in range(len(atom_num)):
    ab = atom_num[j] - 19
    bc = res_num[j] - 2
    new_at_num.append(ab)
    new_res_num.append(bc)
    ss = '{}{:>7s}{}{:<3s}{}{:>5s}{:>12.3f}{:>8.3f}{:>8.2f}{}{}{}'.format('ATOM', str(new_at_num[j]), ' ', str(atom_name[j]),'   C D', str(new_res_num[j]), x_cord[j], y_cord[j], z_cord[j], '  1.00', '  0.00', "\n")
    out2.writelines(ss)
    
inp1.close()
out2.close()