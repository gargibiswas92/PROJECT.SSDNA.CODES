import sys
import math
inp1 = open("test2.pdb", "r")
res_num = []
cord_x = []
cord_y = []
cord_z = []
res_name = []
for line in inp1:
    a = line[7:11]
    a1 = line[14:15]
    b = line[26:38]
    c = line[38:46]
    d = line[46:54]
    res_num.append(int(a))
    res_name.append(a1)
    cord_x.append(float(b))
    cord_y.append(float(c))
    cord_z.append(float(d))
#    print(a, b, c, d)
    
for i in range(len(res_num)):
    if res_name[i] == 'P' and res_num[i] != 1519:
            x1 = cord_x[i]
            y1 = cord_y[i]
            z1 = cord_z[i]
    
            x2 = cord_x[i+1]
            y2 = cord_y[i+1]
            z2 = cord_z[i+1]
            
#print(x1, x2, x11, x22)
#print(y1, y2, y11, y22)
#print(z1, z2, z11, z22)
            r1x = x2 - x1
            r1y = y2 - y1
            r1z = z2 - z1
            r1_m = math.sqrt(r1x**2 + r1y**2 + r1z**2)
            print(r1_m)
#               ssp = [res_num[i], " ", res_num[j], " ", r1_m, "\n"]
#               out1.writelines(ssp)
#print(contact, dist)

#out1 = open("cont.txt", "w")
#for key in contact:
#    for i in range(len(res_prot)):
#        ss ='{:>4s} {:>4s} {:>8.6f} {}'.format(str(key), str(contact[key][i]), dist[key][i], "\n")
#        out1.writelines(ss) 
#ss = ["r1x_f =", str(r1x_f), "  r1y_f =", str(r1y_f), "  r1z_f =", str(r1z_f), "\n"]
#ss2 = ["r2x_f =", str(r2x_f), "  r2y_f =", str(r2y_f), "  r2z_f =", str(r2z_f), "\n"]
#out1 = open("cont.txt", "w")
#out1.writelines(ss2)
#out1.close()
inp1.close()
