import numpy as np

#discr
x  = np.zeros(5)
x[0] = 32
for i in range(1,5):
    x[i] = x[i-1] * 2.0

#values for quad conv (q = 2)
q = 2.0

snd_order = np.zeros(5)
snd_order[0] = 0.01
for i in range(1,5):
    snd_order[i] = snd_order[i-1]/np.exp(q*np.log(2.0))

#print
for  i in range(5):
    print(x[i], snd_order[i])

