inp1 = open("test1.dat")

di = []
atom1 = []
atom2 = []
atom3 = []
atom4 = []
ang = []

out1 = open("change_pers.dat", "w")

for line in inp1:
    arr = line.split()
    di.append(int(arr[0]))
    atom1.append(int(arr[1]))
    atom2.append(int(arr[2]))
    atom3.append(int(arr[3]))
    atom4.append(int(arr[4]))
    ang.append(float(arr[5]))
    
m = 0.0    
for i in range(len(di)):    
    ss = '{:>5d}{:>5d}{:>5d}{:>5d}{:>5d}{:>8.3f}{:>8.3f}{:>8.3f}{:8.3f}{}'.format(di[i], atom1[i], atom2[i], atom3[i], atom4[i], ang[i], m, m, m, "\n")
    out1.writelines(ss)
inp1.close()
out1.close()