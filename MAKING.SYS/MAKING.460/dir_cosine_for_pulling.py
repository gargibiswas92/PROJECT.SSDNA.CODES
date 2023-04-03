import math
inp1 = open("cg.pdb", "r")

bead = []
x_cord = []
y_cord = []
z_cord = []

for lines in inp1:
    a = lines[4:11]
    b = lines[26:38]
    c = lines[38:47]
    d = lines[47:55]
    
    bead.append(int(a))
    x_cord.append(float(b))
    y_cord.append(float(c))
    z_cord.append(float(d))
    
def get_coordinates(ind):
    for k in range(len(bead)):
        if bead[k] == ind:
            l = x_cord[k]
            m = y_cord[k]
            n = z_cord[k]
    return l, m, n
def dir_cosine(x1, y1, z1, x2, y2, z2):
    r1 = (x1 - x2)
    r2 = (y1 - y2)
    r3 = (z1 - z2)
    r_m = math.sqrt(r1**2 + r2**2 + r3**2)
    return (r1/r_m), (r2/r_m), (r3/r_m)

def distance(x1, y1, z1, x2, y2, z2):
    r1 = (x1 - x2)
    r2 = (y1 - y2)
    r3 = (z1 - z2)
    r_m = math.sqrt(r1**2 + r2**2 + r3**2)
    return (r_m)

for i in range(len(bead)):
    if bead[i] == 852:
        ind1 = bead[i]
        x1, y1, z1 = get_coordinates(ind1)
    if bead[i] == 1220:
        ind2 = bead[i]
        x2, y2, z2 = get_coordinates(ind2)
        ll, mm, nn = dir_cosine(x1, y1, z1, x2, y2, z2)
        dd = distance(x1, y1, z1, x2, y2, z2)
#print(ll, mm, nn)

def get_points(x1, y1, z1, l, m, n, dist):
    x2 = x1 + dist*l
    y2 = y1 + dist*m
    z2 = z1 + dist*n
    return x2, y2, z2

x3, y3, z3 = get_points(x1, y1, z1, ll, mm, nn, 200)
#print(x3, y3, z3)#force to be applied on residue 928## 2042

x4, y4, z4 = get_points(x2, y2, z2, -ll, -mm, -nn, 200)
#print(x3, y3, z3)#force to be applied on residue 928## 1340

for i in range(len(bead)):
    if bead[i] == 1340:
        ind3 = bead[i]
        #print(ind3)
        x5, y5, z5 = get_coordinates(ind3)
        #print(x5, y5, z5)
    if bead[i] == 2717:
        ind4 = bead[i]
        x6, y6, z6 = get_coordinates(ind4)
        
l1, m1, n1 = dir_cosine(x4, y4, z4, x5, y5, z5)
l2, m2, n2 = dir_cosine(x3, y3, z3, x6, y6, z6)
print(l1, m1, n1, ind3)
print(l2, m2, n2, ind4)