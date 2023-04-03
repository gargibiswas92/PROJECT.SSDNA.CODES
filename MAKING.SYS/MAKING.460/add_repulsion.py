inp1 = open("bead_inf.dat", "r")
inp2 = open("abc3.dat", "r")
out1 = open("added_repulsion.dat", "w")

index1 = []
index2 = []
bead = []
residue = []

for line in inp1:
    a = line[0:5]
    b = line[5:10]
    c = line[10:13]
    d = line[13:17]
    index1.append(int(a))
    index2.append(int(b))
    bead.append(c.rstrip())
    residue.append(d.rstrip())
    
new_rep1 = []
new_rep2 = []
for line2 in inp2:
    a2, b2 = line2.split()
    new_rep1.append(int(a2))
    new_rep2.append(int(b2))
    
d1=4.410
d2=5.664
d3=7.076
ep=1.000
def bead_identity(ind):
    for m in range(len(residue)):
        if index1[m] == ind:
            be = bead[m]
            return be
    

        m = new_rep1[j]
        n = new_rep2[j]
rep = 3681781
for j in range(len(new_rep1)):
        repb1 = bead_identity(new_rep1[j])
        repb2 = bead_identity(new_rep2[j])
        if repb1 == ' CA' and repb2 == ' CA':
            rep = rep + 1
            ss = '{:>8d}{:>5d}{:>5d}{:>10.3f}{:>10.3f}{}'.format(rep, new_rep1[j], new_rep2[j], d3, ep, "\n")
            out1.writelines(ss)
            
        if repb1 == ' CA' and repb2 == ' CB':
            rep = rep + 1
            ss = '{:>8d}{:>5d}{:>5d}{:>10.3f}{:>10.3f}{}'.format(rep, new_rep1[j], new_rep2[j], d2, ep, "\n")
            out1.writelines(ss)
            
        if repb1 == ' CB' and repb2 == ' CA':
            rep = rep + 1
            ss = '{:>8d}{:>5d}{:>5d}{:>10.3f}{:>10.3f}{}'.format(rep, new_rep1[j], new_rep2[j], d2, ep, "\n")
            out1.writelines(ss)
        
        if repb1 == ' CB' and repb2 == ' CB':
            rep = rep + 1
            ss = '{:>8d}{:>5d}{:>5d}{:>10.3f}{:>10.3f}{}'.format(rep, new_rep1[j], new_rep2[j], d1, ep, "\n")
            out1.writelines(ss)
            
out1.close()
inp1.close()
inp2.close()
    
        