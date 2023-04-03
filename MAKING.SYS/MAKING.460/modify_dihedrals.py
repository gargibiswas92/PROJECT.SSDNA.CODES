inp1 = open("modrep_rpa.dat", "r")
out1 = open("mod_dihed.dat", "w")
count = 0
di_ind = []
at_ind1 = []
at_ind2 = []
at_ind3 = []
at_ind4 = []
ang = []
strength = []
rand1 = []
rand2 = []
for line in inp1:
    arr = []
    count = count + 1
    if count >= 39199 and count <= 40805:
        arr = line.split()
        di_ind.append(int(arr[0]))
        at_ind1.append(int(arr[1]))
        at_ind2.append(int(arr[2]))
        at_ind3.append(int(arr[3]))
        at_ind4.append(int(arr[4]))
        ang.append(float(arr[5]))
        strength.append(float(arr[6]))
        rand1.append(float(arr[7]))
        rand2.append(float(arr[8]))
        
change_di = []
for j in range(542, 560):
    change_di.append(j)
for j in range(778, 798):
    change_di.append(j)
for j in range(978, 999):
    change_di.append(j)
for j in range(1021, 1038):
    change_di.append(j)    
cn = 0    
st = 0.0
for k in range(len(at_ind1)):
    flag = 0
    for i in range(len(change_di)):
        if (at_ind2[k] == change_di[i]) or (at_ind3[k] == change_di[i]):
            flag = 1
            cn = cn + 1
    if flag == 0:
        ssc = '{:>5d}{:>5d}{:>5d}{:>5d}{:>5d}{:>8.3f}{:>8.3f}{:>8.3f}{:>8.3f}{}'.format(di_ind[k], at_ind1[k], at_ind2[k], at_ind3[k], at_ind4[k], ang[k], strength[k], rand1[k], rand2[k], "\n")
        out1.writelines(ssc)
    if flag == 1:
        ssd = '{:>5d}{:>5d}{:>5d}{:>5d}{:>5d}{:>8.3f}{:>8.3f}{:>8.3f}{:>8.3f}{}'.format(di_ind[k], at_ind1[k], at_ind2[k], at_ind3[k], at_ind4[k], ang[k], st, rand1[k], rand2[k], "\n")
        out1.writelines(ssd)
#        print(di_ind[k], at_ind1[k], at_ind2[k], at_ind3[k], at_ind4[k])
print(cn)
inp1.close()
out1.close()