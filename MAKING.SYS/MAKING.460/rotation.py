import math

inp1 = open("snap18.pdb", "r")
atom_num = []
atom_name = []
res_name = []
res_num = []
x = []
y = []
z = []
count = 0

for line in inp1:
    a = line[4:11]
    b = line[11:15]
    c = line[15:20]
    d = line[22:26]
    e = line[26:38]
    f = line[38:46]
    g = line[46:54]
    
    atom_num.append(int(a))
    atom_name.append(b.strip())
    res_name.append(c.strip())
    res_num.append(int(d))
    x.append(float(e))
    y.append(float(f))
    z.append(float(g))
    count = count + 1
    
def dir_cosine(x1, y1, z1, x2, y2, z2):
    dx = x2 - x1
    dy = y2 - y1
    dz = z2 - z1
    rm = math.sqrt(dx**2 + dy**2 + dz**2)
    return dx/rm, dy/rm, dz/rm
    
for i in range(count):
    if atom_num[i] == 1340:
        ind1 = i
        x1 = x[i]
        y1 = y[i]
        z1 = z[i]
        print(x1, y1, z1)
    if atom_num[i] == 2042:
        ind2 = i
        x2 = x[i]
        y2 = y[i]
        z2 = z[i]
        print(x2, y2, z2)
    alpha, beta, gamma = dir_cosine(x1, y1, z1, x2, y2, z2)
print(ind1, ind2, alpha, beta, gamma)

def translate(x, y, z, a, b, c):
    xn = []
    yn = []
    zn = []
    for i in range(len(x)):
        xx = x[i] + a
        yy = y[i] + b
        zz = z[i] + c
        xn.append(xx)
        yn.append(yy)
        zn.append(zz)
    return xn, yn, zn
xt, yt, zt = translate(x, y, z, 14.216, 0, 0)
        
maxx = max(xt)
minx = min(xt)
print(maxx, minx)
maxy = max(yt)
miny = min(yt)
print(maxy, miny)
maxz = max(zt)
minz = min(zt)
print(maxz, minz)

m = 57.295
ang1 = math.acos(alpha)
ang2 = math.acos(beta)
ang3 = math.acos(gamma)

print(math.acos(alpha), math.acos(beta), math.acos(gamma))
print(math.acos(alpha)*m, math.acos(beta)*m, math.acos(gamma)*m)

def x_rotation(ang, x, y, z):
    x_n1 = []
    y_n1 = []
    z_n1 = []
    for i in range(len(x)):
        a1 = x[i]
        b1 = (math.cos(ang)*y[i]) - (z[i]*math.sin(ang))
        c1 = (math.sin(ang)*y[i]) + (z[i]*math.cos(ang))
        
        x_n1.append(a1)
        y_n1.append(b1)
        z_n1.append(c1)
        
    return x_n1, y_n1, z_n1

def y_rotation(ang, x, y, z):
    x_n2 = []
    y_n2 = []
    z_n2 = []
    for i in range(len(x)):
        a1 = (x[i]*math.cos(ang)) + (z[i]*math.sin(ang))
        b1 = y[i]
        c1 = (z[i]*math.cos(ang)) - (x[i]*math.sin(ang))
        
        x_n2.append(a1)
        y_n2.append(b1)
        z_n2.append(c1)
        
    return x_n2, y_n2, z_n2

def z_rotation(ang, x, y, z):
    x_n3 = []
    y_n3 = []
    z_n3 = []
    for i in range(len(x)):
        a1 = (x[i]*math.cos(ang)) - (y[i]*math.sin(ang))
        b1 = (x[i]*math.sin(ang)) + (y[i]*math.cos(ang))
        c1 = z[i]
        
        x_n3.append(a1)
        y_n3.append(b1)
        z_n3.append(c1)
        
    return x_n3, y_n3, z_n3

#rotx1, roty1, rotz1 = x_rotation(ang1, xt, yt, zt)
#rotx2, roty2, rotz2 = y_rotation(ang2, rotx1, roty1, rotz1)
#rotx3, roty3, rotz3 = z_rotation(ang3, rotx2, roty2, rotz2)

rotx1, roty1, rotz1 = y_rotation(ang2, xt, yt, zt)
rotx2, roty2, rotz2 = z_rotation(ang3, rotx1, roty1, rotz1)
rotx3, roty3, rotz3 = x_rotation(ang1, rotx2, roty2, rotz2)

#for ii in range(len(rotx)):
#    print(rotx[ii])

for i in range(count):
    if atom_num[i] == 1340:
        ind1 = i
        x1 = rotx3[i]
        y1 = roty3[i]
        z1 = rotz3[i]
        print(x1, y1, z1)
    if atom_num[i] == 2042:
        ind2 = i
        x2 = rotx3[i]
        y2 = roty3[i]
        z2 = rotz3[i]
        print(x2, y2, z2)
    alpha1, beta1, gamma1 = dir_cosine(x1, y1, z1, x2, y2, z2)
print(ind1, ind2, alpha1, beta1, gamma1)

m = 180/(math.pi)
ang11 = math.acos(alpha1)
ang22 = math.acos(beta1)
ang33 = math.acos(gamma1)
print(ang11, ang22, ang33)
print(m*ang11, m*ang22, m*ang33)

m = 1.0
n = 0
atom_type = "C"
out1 = open("rot_mod_18.dat", "w")   
for i in range(len(atom_num)):
    if atom_num[i] <= 1339:
        #ss = '{:>5d}{:>4d}{:>3s}{:>4s}{:>8.3f}{:>8.3f}{:>8.3f}{:>8.3f}{}'.format(atom[i], res_num[i], atom_name[i], res_name[i], x[i], y[i], z[i], m, "\n")
        ss = '{:>5d}{:>4d}{:>3s}{:>4s}{:>8.3f}{:>8.3f}{:>8.3f}{:>8.3f}{}'.format(atom_num[i], res_num[i], atom_name[i], res_name[i], rotx3[i], roty3[i], rotz3[i], m, "\n")
        out1.writelines(ss)
    if atom_num[i] >= 1340:
        ss1 = '{:>5d}{:>4d}{:>2s}{:>5s}{:>8.3f}{:>8.3f}{:>8.3f}{:>8.3f}{:>3d}{}'.format(atom_num[i], atom_num[i], atom_name[i], res_name[i], rotx3[i], roty3[i], rotz3[i], m, n, "\n")
        out1.writelines(ss1)
out1.close()

maxx = max(rotx3)
minx = min(rotx3)
print(maxx, minx)
maxy = max(roty3)
miny = min(roty3)
print(maxy, miny)
maxz = max(rotz3)
minz = min(rotz3)
print(maxz, minz)