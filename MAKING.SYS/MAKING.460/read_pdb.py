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