inp1 = open("snap_486.pdb", "r")
out1 = open("snap_486_formatted.dat", "w")
atom = []
atom_name = []
res_name = []
res_num = []
x = []
y = []
z = []
#atom_type = []
for line in inp1:
    a = line[4:11]
    b = line[11:15]
    c = line[15:20]
    d = line[22:26]
    e = line[26:38]
    f = line[38:46]
    g = line[46:54]
    #h = line[67:-1]
    #print(h)
    atom.append(int(a))
    atom_name.append(b.strip())
    res_name.append(c.strip())
    res_num.append(int(d))
    x.append(float(e))
    y.append(float(f))
    z.append(float(g))
    #atom_type.append(h)
    
m = 1.0
n = 0
atom_type = "C"
    
for i in range(len(atom)):
    if atom[i] <= 1339:
        #ss = '{:>5d}{:>4d}{:>3s}{:>4s}{:>8.3f}{:>8.3f}{:>8.3f}{:>8.3f}{}'.format(atom[i], res_num[i], atom_name[i], res_name[i], x[i], y[i], z[i], m, "\n")
        ss = '{:>5d}{:>4d}{:>3s}{:>4s}{:>8.3f}{:>8.3f}{:>8.3f}{:>8.3f}{}'.format(atom[i], res_num[i], atom_name[i], res_name[i], x[i], y[i], z[i], m, "\n")
        out1.writelines(ss)
    if atom[i] >= 1340:
        ss1 = '{:>5d}{:>4d}{:>2s}{:>5s}{:>8.3f}{:>8.3f}{:>8.3f}{:>8.3f}{:>3d}{}'.format(atom[i], atom[i], atom_name[i], res_name[i], x[i], y[i], z[i], m, n, "\n")
        out1.writelines(ss1)
        
out1.close()
inp1.close()
        