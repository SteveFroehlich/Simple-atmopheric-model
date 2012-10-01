#!/usr/bin/env python
"""
Example: simple line plot.
Show how to make and save a simple line plot with labels, title and grid
"""
import pylab

data = []

#get data from file
FILES = ["out.txt"]
for file in FILES:
	infile = open(file, "r")

for eachline in infile:
	a = eachline
	data.append(a)

s = data
#t = numpy.arange(0.0, 1.0+0.01, 0.01)
#s = numpy.cos(numpy.pi*t)
#s = [46.0,40.43,35.94,32.29,29.32,26.89,24.9,23.26,21.91,20.8,19.88,19.12,18.49,17.97,17.55]
pylab.plot(s)
 
pylab.xlabel('time (hrs)')
pylab.ylabel('temperature (deg F)')
pylab.title('Overnight temperature')
pylab.grid(True)
pylab.savefig('simple_plot')
   
pylab.show()
