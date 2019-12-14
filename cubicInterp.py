import numpy as np
import math
import sys

def calcSlopes(logData):
	b = np.zeros((len(logData),1))
	b[0,0]  = 3 * (logData[1,1]-logData[0,1]) / ((logData[1,0]-logData[0,0])**2)
	b[-1,0] = 3 * (logData[-1,1]-logData[-2,1]) / ((logData[-1,0]-logData[-2,0])**2)
	for n in range(1,len(b)-1):
		b[n,0] = 3 * (logData[n,1]-logData[n-1,1]) / ((logData[n,0]-logData[n-1,0])**2) + \
				 3 * (logData[n+1,1]-logData[n,1]) / ((logData[n+1,0]-logData[n,0])**2)

	a = np.zeros((len(logData), len(logData)))
	# First and last
	a[0,0]  = 2 / (logData[1,0]-logData[0,0])
	a[-1,-1] = 2 / (logData[-1,0]-logData[-2,0])
	# Starting diagonal
	a[0,1]  = 1 / (logData[1,0]-logData[0,0])
	a[1,0]  = 1 / (logData[1,0]-logData[0,0])
	# Main Diagonal Loop
	for n in range(1,len(a)-1):
		a[n,n] = 2 / (logData[n,0]-logData[n-1,0]) + \
				 2 / (logData[n+1,0]-logData[n,0])
		a[n,n+1] = 1 / (logData[n+1,0]-logData[n,0])
		a[n+1,n] = 1 / (logData[n+1,0]-logData[n,0])
	return np.linalg.solve(a, b)

# Read data from file
my_data = np.genfromtxt('lanthanumTest.txt', delimiter=' ')

# Convert data into a usable form
dataSetList = []
nextDataSet = np.asarray([[math.log(my_data[0,0]), math.log(my_data[0,1])]])
for n in range(1,len(my_data)):
	nextRow = np.asarray([[math.log(my_data[n,0]), math.log(my_data[n,1])]])
	# If there is a matching X value, start a new data set
	if nextRow[0,0] == nextDataSet[-1,0]:
		dataSetList += [nextDataSet]
		nextDataSet = nextRow
	# Otherwise, add to the end of this data set
	else:
		nextDataSet = np.append(nextDataSet, nextRow, axis=0)
dataSetList += [nextDataSet]

resultList = []
for dataSet in dataSetList:
	returnArray = calcSlopes(dataSet)
	resultSet = []
	for n in range(0, len(dataSet)):
		resultSet += [np.append(dataSet[n], returnArray[n])]
	resultList += [resultSet]

fSet = resultList[0]
for n in range(0,len(fSet)-1):
	row = fSet[n]
	#print str(row[0])+" "+str(row[1])+" "+str(row[2])+" "+str(row[1])+" "+str(row[2])
	print "{0:.5f} {1:.5f} {2:.5f} {3:.5f} {4:.5f}".format(row[0], row[1], row[2], row[1], row[2])
row = fSet[-1]
#sys.stdout.write(str(row[0])+" "+str(row[1])+" "+str(row[2])+" ")
sys.stdout.write("{0:.5f} {1:.5f} {2:.5f} ".format(row[0], row[1], row[2]))
for setN in range(1,len(resultList)):
	resultSet = resultList[setN]
	row = resultSet[0]
	#print str(row[1])+" "+str(row[2])
	print "{0:.5f} {1:.5f}".format(row[1], row[2])
	for n in range(1,len(resultSet)-1):
		row = resultSet[n]
		#print str(row[0])+" "+str(row[1])+" "+str(row[2])+" "+str(row[1])+" "+str(row[2])
		print "{0:.5f} {1:.5f} {2:.5f} {3:.5f} {4:.5f}".format(row[0], row[1], row[2], row[1], row[2])
	row = resultSet[-1]
	#sys.stdout.write(str(row[0])+" "+str(row[1])+" "+str(row[2])+" ")
	sys.stdout.write("{0:.5f} {1:.5f} {2:.5f} ".format(row[0], row[1], row[2]))

lastSet = resultList[-1]
row = lastSet[-1]
#print str(row[1])+" "+str(row[2])
print "{0:.5f} {1:.5f}".format(row[1], row[2])
