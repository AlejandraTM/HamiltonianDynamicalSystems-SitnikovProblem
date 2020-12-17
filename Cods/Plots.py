import numpy as np
import math as mt
from scipy import constants as cts
import matplotlib.pyplot as plt

Ya=False
r=0
k=0
while(Ya==False):
	print (k,r)
	option=int(input(" 1. Position vs Time Plot \n 2. Velocity vs Time PLot \n 3. Phase Space Plot \n 4. Poincare Map Plot \n 5. Exit \n"))
	if(option==5):
		break;
	if(option!=4):
		if (r==0):
			#Read File Method Solution
			f=open("exc(0,9)init(1,1).txt","r")
			l=[]
			t=[]
			v=[]
			p=[]
			for x in f.readlines():
			    m=x.split(";")
			    m=[y.rstrip() for y in m]
			    m=[y.replace("\\"," ") for y in m]
			    for z in m:
				try:
				    l.append(float(z))
				except:
				    pass
			f.close()

			for i in range(0,len(l),1):
				if(i%3==0):
					t.append(l[i])
				elif(i%3==1):
					p.append(l[i])
				elif(i%3==2):
					v.append(l[i])

			print("Reading Vector size: ",len(l))
			r=1
		if(option==1):
			plt.plot(t,p,'-')
			plt.xlabel('Time')
			plt.xticks([8*np.pi,16*np.pi,32*np.pi,64*np.pi,128*np.pi,256*np.pi,512*np.pi,1024*np.pi,2048*np.pi],
			  [r'$8\pi$', r'$16\pi$',r'$32\pi$',r'$64\pi$',r'$128\pi$',r'$256\pi$',r'$512\pi$',r'$1024\pi$',r'$2048\pi$'])
			plt.xlim(128*np.pi,2048*np.pi)
			plt.ylim(-20,100)
			plt.grid(True)
			plt.ylabel('Distance from initial position ($z_0=1$)')
			plt.title("Position Evolution for $\epsilon = 0.08$.")
			plt.savefig('Pe(0,08)int(1,1)')
			plt.show()
		elif(option==2):
			plt.plot(t,v,'-')
			plt.xlabel('Time')
			plt.xticks([8*np.pi,16*np.pi,32*np.pi,64*np.pi,128*np.pi,256*np.pi,512*np.pi,1024*np.pi,2048*np.pi],
			  [r'$8\pi$', r'$16\pi$',r'$32\pi$',r'$64\pi$',r'$128\pi$',r'$256\pi$',r'$512\pi$',r'$1024\pi$',r'$2048\pi$'])
			plt.xlim(0*np.pi,64*np.pi)
			#plt.ylim(-1000,1000)
			plt.grid(True)
			plt.ylabel('Velocity')
			#plt.title("Velocity vs Time")
			plt.savefig('Ve(0,9)int(1,1)')
			plt.show()
		elif(option==3):
			plt.plot(p,v,'.',markersize=1.0)
			plt.xlabel('Position')
			plt.ylabel('Velocity')
			plt.xlim(-55,55)
			#plt.ylim(-0.5,0.5)
			plt.grid(True)
			plt.title("Phase Space of $\epsilon =0.9$ ")
			plt.savefig('PS(0,9)int(1,1)')
			plt.show()
	elif(option==4):
		if (k==0):
			#Read File Poincare Map
			f=open("PMexc(0,9)init(1,1).txt","r")
			l1=[]
			t1=[]
			v1=[]
			p1=[]
			for x in f.readlines():
			    n=x.split(";")
			    n=[y.rstrip() for y in n]
			    n=[y.replace("\\"," ") for y in n]
			    for z in n:
				try:
				    l1.append(float(z))
				except:
				    pass
			f.close()

			for i in range(0,len(l1),1):
				if(i%3==0):
					t1.append(l1[i])
				elif(i%3==1):
					p1.append(l1[i])
				elif(i%3==2):
					v1.append(l1[i])

			print("Reading Vector size: ",len(l1))
			k=1
		plt.plot(p1,v1,'.', markersize=1.5)
		plt.xlabel('Position')
		plt.ylabel('Velocity')
		#plt.xlim(-.06,.06)
		plt.xlim(-.03,.03)
		plt.ylim(-2.5,2.5)
		plt.grid(True)
		plt.title("Poincare Map of $\epsilon =0.9$ ")
		plt.savefig('MP(0,9)slope(1)4')
		plt.show()

