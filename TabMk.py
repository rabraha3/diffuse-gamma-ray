

f = open('CldNum.txt','r')
g = open('Clds_Distance.txt','r')
h = open('Clds_Emissivity.txt','r')


em = h.readlines()
dis= g.readlines()
cor= f.readlines()


"""
Ok, so to do this, I need -- for each iteration,
coordinates, distance, CO (use coordinates to get CO template),
(emissivity +/- uncertainty)*scaling

"""
