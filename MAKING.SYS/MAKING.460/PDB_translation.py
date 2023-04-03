input1 = input("Insert the filename:")
inp1 = open("dna.pdb", "r")

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
        atom_type.append(str(arr[11]))

for i in range(len(atom_num)):
    if atom_num[i] == 6039:
        a = x_cord[i]
        b = y_cord[i]
        c = z_cord[i]
        #print(a, b, c)
    if atom_num[i] == 7711:
        m = x_cord[i]
        n = y_cord[i]
        l = z_cord[i]
    aa = a - m
    bb = b - n
    cc = c - l
    
x_mod = []
y_mod = []
z_mod = []

out1 = open("mod_dna2.pdb", "w")

for j in range(len(at_num)):
    p = x_cord[j] + aa
    q = y_cord[j] + bb
    r = z_cord[j] + cc
    x_mod.append(p)
    y_mod.append(q)
    z_mod.append(r)
    
    ss = '{}{:>7s}{}{:<3s}{}{:>4s}{:>12.3f}{:>8.3f}{:>8.2f}{}{}{}'.format('ATOM', str(atom_num[j]), '  ', str(atom_name[j]),'   C D', str(res_num[j]), x_mod[j], y_mod[j], z_mod[j], '  1.00', '  0.00', "\n")
    out1.writelines(ss)
    
inp1.close()
out1.close()