import sys
import math
inp1 = open("final.pdb", "r")
out1 = open("write_contacts.txt", "w")
res_num = []
cord_x = []
cord_y = []
cord_z = []
res_prot = []
res_dis = []
contact = {}
dist = {}
for line in inp1:
    a = line[0:5]
    b = line[16:24]
    c = line[24:32]
    d = line[32:40]
    res_num.append(int(a))
    cord_x.append(float(b))
    cord_y.append(float(c))
    cord_z.append(float(d))
#    print(a, b, c, d)
    
for i in range(len(res_num)):
    for j in range(len(res_num)):
        flag = 0
        if (res_num[i] == 1397 or res_num[i] == 1412 or res_num[i] == 1457 or res_num[i] == 1430) and res_num[j] <= 1339:
        #if res_num[i] == 1397 and res_num[j] <= 1339:
            #print(res_num[i])
            x1 = cord_x[i]
            y1 = cord_y[i]
            z1 = cord_z[i]
    
            x2 = cord_x[j]
            y2 = cord_y[j]
            z2 = cord_z[j]
            
#print(x1, x2, x11, x22)
#print(y1, y2, y11, y22)
#print(z1, z2, z11, z22)
            r1x = x2 - x1
            r1y = y2 - y1
            r1z = z2 - z1
            r1_m = math.sqrt(r1x**2 + r1y**2 + r1z**2)
#            print(r1_m)
            if r1_m <= 9.00:
               flag = 1
               res_prot.append(res_num[j])
               res_dis.append(r1_m)
               contact[res_num[i]] = res_prot
               dist[res_num[i]] = res_dis
#               ssp = [res_num[i], " ", res_num[j], " ", r1_m, "\n"]
               ssp ='{:>4s} {:>4s} {:>8.6f} {}'.format(str(res_num[i]), str(res_num[j]), r1_m, "\n")

               out1.writelines(ssp)
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
out1.close()
inp1.close()
