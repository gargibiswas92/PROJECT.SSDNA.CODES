inp1 = open("dynamic_RPA_ssdna460.dat", "r")
out2 = open("contact.dat", "w")
out3 = open("abc2.dat", "w")
count = 0
ind = []
cont1 = []
cont2 = []
dis = []
strength = []
cn = 0
for line in inp1:
    count = count + 1
    if count >= 40807 and count <= 44627:
        #print(line)
        a, b, c, d, e = line.split()
        ind.append(int(a))
        cont1.append(int(b))
        cont2.append(int(c))
        dis.append(float(d))
        strength.append(float(e))
        cn = cn + 1
resi = []
resi2 = []
for p in range(538, 566, 1):
    resi.append(p)
for p in range(644, 672, 1):
    resi.append(p)
for p in range(769, 804, 1):
    resi.append(p)
for p in range(886, 895, 1):
    resi.append(p)
#print(resi)

for q in range(357, 365, 1):
    resi.append(q)
for q in range(317, 327, 1):
    resi.append(q)
for q in range(404, 426, 1):
    resi.append(q)
for q in range(978, 999, 1):
    resi.append(q)
for q in range(1072, 1139, 1):
    resi.append(q)
for q in range(1271, 1283, 1):
    resi.append(q)
#print(resi2)
m = 0.0
n = 0.1
cn = 0
change_cont1 = []
change_cont2 = []

for item in resi:
    for item2 in resi:
        if item != item2:
            #print(item, item2)
            for i in range(len(cont1)):
                if (cont1[i] == item and cont2[i] == item2) and (strength[i] == 1.00):
                    change_cont1.append(cont1[i])
                    change_cont2.append(cont2[i])
                    ssp = '{:>8d}{:>8d}{}'.format(item, item2, "\n")
                    out3.writelines(ssp)
                    cn = cn + 1
                    

for j in range(len(cont1)):
    flag = 0
    for k in range(len(change_cont1)):
        if (cont1[j] == change_cont1[k] and cont2[j] == change_cont2[k]) and (strength[j] == 1.00):
            flag = 1
            break
    if flag == 0:
        ssd = '{:>5d}{:>5d}{:>5d}{:>10.3f}{:>9.6f}{}'.format(ind[j], cont1[j], cont2[j], dis[j], strength[j], "\n")
        out2.writelines(ssd)
    if flag == 1:
        ssp = '{:>5d}{:>5d}{:>5d}{:>10.3f}{:>9.6f}{}'.format(ind[j], cont1[j], cont2[j], dis[j], m, "\n")
        out2.writelines(ssp)
print(cn)

out2.close()        
inp1.close()
out3.close()