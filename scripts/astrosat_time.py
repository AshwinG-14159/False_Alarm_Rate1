import sys
from astropy.time import Time, TimeDelta
import argparse


def isot_to_astrosat(time):
	time0 = '2010-01-01T00:00:00'
	t0 = Time(time0,format='isot',scale='utc')
	t = Time(time,format='isot',scale='utc')
	dt = t - t0
	return dt.sec

def astrosat_to_isot(time):
	time0 = '2010-01-01T00:00:00'
	t0 = Time(time0,format='isot',scale='utc')
	time = float(time)
	dt = TimeDelta(time, format='sec')
	t = t0 + dt
	return t.isot

def toDecimal(time):
	time = str(time)
	date = time[0:10]
	t = time[11:len(time)]
	h = t[0:2]
	m = t[3:5]
	dec = (float(h)+(float(m)/60.0))/24.0
	dec = str(round(dec,3))
	decdate = date+dec[1:]
	return decdate

def main():
	parser=argparse.ArgumentParser()
	parser.add_argument("type", type=str, choices=('sec', 'isot'))
	parser.add_argument("time",type=str)
	args=parser.parse_args()
	global time0 
	time0 = '2010-01-01T00:00:00'
	global t0 
	t0 = Time(time0,format='isot',scale='utc')
	#print args.time

	if(args.type=='isot'):
		astrosat_time = isot_to_astrosat(args.time)
		#print toDecimal(args.time)
		print(astrosat_time)
	elif(args.type=='sec'):
		isot_time = astrosat_to_isot(args.time)
		print(isot_time)
	
if __name__ == '__main__':
    main()
