import time
import string
def readData(pos,director,fileName):
	file=open(fileName,'r')
	while 1:
		line=file.readline() 
		if not line:
			break
		v=string.split(line) 
		pos.append([eval(v[0]),eval(v[1]),eval(v[2])])
		director.append([eval(v[3]),eval(v[4]),eval(v[5])])
	file.close()
if __name__=='__main__': 
	file_name="c:/mikan4/src/common/resMobTensor/trajectory.dat"
	pos,director=[],[]
	readData(pos,director,file_name)
	polyline(pos,1)
	n=len(pos)
	for i in pos:
		point(i,1)
