import os
import matplotlib.pyplot as plt
diff_com = []
for file in os.listdir():
    if file.startswith('MSD_NT') and file.endswith('.txt'):
        tt = []
        msd = []
        inp1= open(file, 'r')
        for line in inp1:
            arr = line.split()
            tt.append(int(arr[0]))
            msd.append(float(arr[1]))
        ym = msd[800]-msd[300]
        xm = tt[800]-tt[300]
        dm = ym/xm
        diff_com.append(dm)
        plt.plot(tt, msd)

plt.xlabel('Time steps',fontsize=20)
plt.ylabel('MSD (nt$^2$)',fontsize=20)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.title('Distance = 812 $\AA$', fontsize=20)
 
plt.savefig('dis_812_com.png', bbox_inches='tight', dpi=300)

out1 = open('diff_com_812.txt', 'w')

for j in range(len(diff_com)):
    ss = '{:>8.3f}{}'.format(diff_com[j], '\n')
    out1.writelines(ss)

from statistics import mean
from statistics import stdev
me = mean(diff_com)
st = stdev(diff_com)

sp = '{}{:>8.3f}{:>8.3f}'.format('\n',me, st)
out1.writelines(sp)
out1.close()