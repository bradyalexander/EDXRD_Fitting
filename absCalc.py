import numpy as np
import math

def atomicAbsArray(energyList, atom):
	# Read Data from file
	absData = np.genfromtxt('atomicAbsorption/'+atom+'.abs', delimiter=' ')
	# Calculate Linear approximation
	start = absData[0,0]
	end = absData[-1,0]
	numP = len(absData)-1
	slope = (numP - 1 ) / (end-start)
	intercept = - slope * start
	# Calculate the absorption
	absList = []
	for energy in energyList:
		logE = math.log(energy)
		# Guess the index of the nearest data points
		pos = int(slope * logE + intercept)
		# Move until the index is between the next lowest and the next highest
		if pos < 0:
			pos = 0
		elif pos > numP-1:
			pos = numP-1
		elif logE < absData[pos,0]:
			pos -= 1
			while pos > 0 and logE < absData[pos,0]:
				pos -= 1
		elif logE > absData[pos+1,0]:
			pos += 1
			while pos < numP-1 and logE > absData[pos+1,0]:
				pos += 1
		# Calculate the cubic equation for the adjacent data points
		x1 = absData[pos,0]
		x2 = absData[pos+1,0]
		invMe = np.asarray([[1, x1, x1**2,   x1**3],
							[0,  1,  2*x1, 3*x1**2],
							[1, x2, x2**2,   x2**3],
							[0,  1,  2*x2, 3*x2**2]])
		yVals = np.asarray( [absData[pos,3],absData[pos,4],
							absData[pos+1,1],absData[pos+1,2]])
		coeff = np.linalg.solve(invMe, yVals)
		# Add the absorption coefficients to the list
		absList += [math.exp(coeff[0] + coeff[1]*logE + coeff[2]*(logE**2) + coeff[3]*(logE**3))]
	return absList
