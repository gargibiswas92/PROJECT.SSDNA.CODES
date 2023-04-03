inp1 = open('all_int_trans.txt','r')
jump = []
for line in inp1:
    arr = line.split()
    a = float(arr[1])
    if a >= 10.0:
        jump.append(a)
    else:
        continue
import matplotlib.pyplot as plt
import numpy as np

data = np.array(jump)

counts, bins, _ = plt.hist(jump, bins=50, range=(0,200))

constant = 9.0
normalized_counts = counts / constant

import matplotlib.pyplot as plt2

plt2.bar(bins[:-1], normalized_counts, width=(bins[1]-bins[0]), color='orangered')
plt2.xlabel('ΔNT between 1000 consecutive MD steps', fontsize=16)
plt2.ylabel('Frequency', fontsize=16)
plt2.title('Distance= 994 $\AA$, k$_D$=0.7', fontsize=20)
plt2.ylim(0,120)
plt2.savefig('int_trans_994_7.png',dpi=300, bbox_inches='tight')