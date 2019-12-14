import numpy as np
import scipy as scp

# TODO: try two extreme wavelengths
#		get an idea of how much it changes

# turn my axes into a useful format
# apply symmetry operations to atoms
# calculate a list of all interesting peaks

#phase   -> information about the material structure
#		-> the information we want
#histogram -> information about the data/experiment
#		  -> information about the energy spectra used
def patternCalculator(phase, histogram):
	X = histogram['Data'][0, :]
	peakFF = []
	peakInt = []
	peakWidth = []
	for peak in phase['Peaks']:
		# calculate center x for peak
		wave = x(peak)
		peakFF += [calcFormFactor(peak['hkl'], FTable, wave,  phase['asymAtoms'], phase['sMat'])]
		peakInt += [peakFF[-1] * histogram['Scale'] * phase['Scale']
		peakWidth += [calcPVparams(f_G, f_L)] # Something
	for x in X:
		# x is in wavelength!
		# calculate scattered intensity
		# apply broadening to each peak
		for each nearby peak
			int_p += peakBroad(peakWidth[n], x-x0)*peakFF[n]
		int_p = back(x) + int_p
		# apply intensity correction
		int_p = intCorrection * int_p

def background(coeff, x):
	accum = 0
	for i in range(len(coeff)):
		accum += coeff[i] * x**i
	return accum

def intCorrection(incE, absorption, x):
	# calculate incident energy
	corr = fastLPV(incE, x)
	# calculate absorption
	for material in absorption:
		mMu = 1
		for element in material['ElList']:
			# TODO: write a function that gives mass absorption coefficients for various elements
			mMu += elementMu(element['name'], x)*element['fract']*element['mass']/material['volume']
		corr *= exp(- mMu * material['dist']
	return corr

# Calculates the relative intensity of a log pseudovoigt.
# This function expects input from the function calcPVparams
# This function saves time by not recalculating the same parameters over and over
def fastLPV (params, x):
	return params[0]/(x*(((math.log(x))**2)+params[1])) + (params[2]/x)*math.exp(-((math.log(x))**2)/(params[3]))

# Calculates the relative intensity of a pseudovoigt.
# This function expects input from the function calcPVparams
# This function saves time by not recalculating the same parameters over and over
def fastPV (params, x):
	return params[0]/((x**2)+params[1]) + params[2]*math.exp(-(x**2)/(params[3]))

# Calculates the 5 constants used to repeatedly calculate a pseudovoigt or log pseudovoigt quickly
# Params must be a list of the following pseudovoigt or log pseudovoigt parameters:
# f_G = gaussian peak width
# f_L = lorentzian peak width
# mu = the peak center
def calcPVparams (f_G, f_L):
	f = ((f_G**5) + 2.69269*(f_G**4)*(f_L) + 4.47163*(f_G**2)*(f_L**3) + 0.07842*(f_G)*(f_L**4) + (f_L**5))
	eta = 1.36603*(f_L/f) - 0.47719*((f_L/f)**2) + 0.11116*((f_L/f)**3)
	# This format is mu, first L coeff, second L coeff, first G coeff, second G coeff
	return [ eta*f_L/math.pi , f_L**2 , ((1-eta)/(math.sqrt(2*math.pi)*f_G)) , 2*(f_G**2) ]
