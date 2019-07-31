from pylab import *

w0 = -1.
wa = 0.0
zmax = 2.5

z = linspace(0,zmax,20)
w = zeros(z.shape)

for i in range(len(z)):
	w[i] = w0 + wa*z[i]/(1.+z[i])

print 'DDE_z = ', 
for i in range(len(z)):
	print '%5.3f, '%(z[i]),


print ''

print 'DDE_w = ', 
for i in range(len(w)):
	print '%5.3f, '%(w[i]),