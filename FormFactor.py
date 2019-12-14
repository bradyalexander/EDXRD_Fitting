import numpy as np
import scipy as sc
import math


a = 4.1568
b = 4.1568
c = 4.1568
alpha = 90
beta = 90
gamma = 90

cx = math.cos(beta*math.pi/180)
cy = (math.cos(alpha*math.pi/180) - math.cos(gamma*math.pi/180)*math.cos(beta*math.pi/180)) / (math.sin(gamma*math.pi/180))
cz = math.sqrt(1 - cx**2 - cy**2)
aMat = np.asarray([[a, 0, 0],
			[b*math.cos(gamma*math.pi/180), b*math.sin(gamma*math.pi/180), 0],
			[c*cx, c*cy, c*cz]])
sMat = np.asarray([np.cross(aMat[1,:],aMat[2,:])/np.cross(aMat[1,:],aMat[2,:]).dot(aMat[0,:].T),
					np.cross(aMat[0,:],aMat[2,:])/np.cross(aMat[0,:],aMat[2,:]).dot(aMat[1,:].T),
					np.cross(aMat[0,:],aMat[1,:])/np.cross(aMat[0,:],aMat[1,:]).dot(aMat[2,:].T)])

atoms = [ ('La', 0, 0, 0, 1, 0.003),
		  ('B', 0, 0, 0, 1, 0.003)]
asymAtoms = [ ('La', 0, 0, 0, 1, 0.003),
			  ('B', 0.1975, 0.5000, 0.5000, 1, 0.003),
			  ('B', 0.5000, 0.1975, 0.5000, 1, 0.003),
			  ('B', 0.5000, 0.5000, 0.1975, 1, 0.003),
			  ('B', 0.8025, 0.5000, 0.5000, 1, 0.003),
			  ('B', 0.5000, 0.8025, 0.5000, 1, 0.003),
			  ('B', 0.5000, 0.5000, 0.5000, 1, 0.003)]

# TODO: calculate fp and fpp dynamically
# 'El' : [[fa], [fb], [fc], [fp, fpp]]
FTable = {'La':[[20.578, 19.599, 11.3727, 3.28719],[2.94817, 0.244475, 18.7726, 133.124],[2.14678],[-0.436115868526585, 2.16014355564795]],
          'B':[[2.0545, 1.3326, 1.0979, 0.7068],[23.2185, 1.021, 60.3498, 0.1403],[0.1932],[-0.00086710537076686, 0.0000530857369224406]]}

wavelength = 0.2388

hkl = (1,0,0)

def calcFormFactor(hkl, FTable, wavelength, asymAtoms, sMat): # Figure out what variables I need here!
	# calculate the inverse
	HKL = np.asarray(hkl)
	pVec = HKL.dot(sMat)
	inv_D = np.linalg.norm(pVec)
	St_L2 = (inv_D/2)**2
	FTList = {}
	for at, ft in FTable.iteritems():
		f0 = ft[2][0]
		for s, r in zip(ft[0], ft[1]):
			f0 = f0 + s * math.exp(-r*St_L2)
		FTList[at] = f0
 	A0_Sum = 0
	App_Sum = 0
	B0_Sum = 0
	Bpp_Sum = 0
	for at in asymAtoms:
		f0 = FTList[at[0]]
		Tpp = math.exp(-8*math.pi*math.pi*St_L2*at[5])
		HR = hkl[0]*at[1] + hkl[1]*at[2] + hkl[2]*at[3]
		A0_Sum += Tpp*at[4]*math.cos(2*math.pi*HR)*(f0+FTable[at[0]][3][0])
		App_Sum += Tpp*at[4]*math.cos(2*math.pi*HR)*(FTable[at[0]][3][1])
		B0_Sum += Tpp*at[4]*math.sin(2*math.pi*HR)*(f0+FTable[at[0]][3][0])
		Bpp_Sum += Tpp*at[4]*math.sin(2*math.pi*HR)*(FTable[at[0]][3][1])
	return (A0_Sum**2)+(App_Sum**2)+(B0_Sum**2)+(Bpp_Sum**2)

calcFormFactor(hkl, FTable, asymAtoms, sMat)
