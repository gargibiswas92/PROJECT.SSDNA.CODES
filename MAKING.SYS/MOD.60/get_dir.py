import sys
import math
inp1 = open("final.pdb", "r")
res_num = []
cord_x = []
cord_y = []
cord_z = []
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
        if res_num[i] == 1060 and res_num[j] == 321:
            x1 = cord_x[i]
            y1 = cord_y[i]
            z1 = cord_z[i]
    
            x2 = cord_x[j]
            y2 = cord_y[j]
            z2 = cord_z[j]
            
        if res_num[i] == 782 and res_num[j] == 552:
            x11 = cord_x[i]
            y11 = cord_y[i]
            z11 = cord_z[i]
            
            x22 = cord_x[j]
            y22 = cord_y[j]
            z22 = cord_z[j]
            
#print(x1, x2, x11, x22)
#print(y1, y2, y11, y22)
#print(z1, z2, z11, z22)
r1x = x2 - x1
r1y = y2 - y1
r1z = z2 - z1
r1_m = math.sqrt(r1x**2 + r1y**2 + r1z**2)
r1x_f = r1x/r1_m
r1y_f = r1y/r1_m
r1z_f = r1z/r1_m

r2x = x22 - x11
r2y = y22 - y11
r2z = z22 - z11
r2_m = math.sqrt(r2x**2 + r2y**2 + r2z**2)
r2x_f = r2x/r2_m
r2y_f = r2y/r2_m
r2z_f = r2z/r2_m
        
print(r1x_f, r1y_f, r1z_f) 
print(r2x_f, r2y_f, r2z_f)

ss = ["r1x_f =", str(r1x_f), "  r1y_f =", str(r1y_f), "  r1z_f =", str(r1z_f), "\n"]
ss2 = ["r2x_f =", str(r2x_f), "  r2y_f =", str(r2y_f), "  r2z_f =", str(r2z_f), "\n"]
out1 = open("get_dir.txt", "w")
out1.writelines(ss)
out1.writelines(ss2)
out1.close()
inp1.close()
