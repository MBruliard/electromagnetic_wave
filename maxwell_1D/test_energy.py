__author__ = "BRULIARD Margaux"
__partners__ = "NMPP"
__date__ = "2018-08-29"


import numpy as np
from space.spaces import creationSpaces

from util.energy import *
from mat1d.mass import *
from scheme.eulerexpl import *
from scheme.eulerimpl import *
from scheme.rk2 import *
from scheme.rk4 import *
from proj.projectors import *


if __name__ == "__main__":

	# data
	ne = 10
	grid =  np.linspace(0., 1., ne)
	degree = 3

	
	#epsilon0 = 8.854e-2
	#mu0 = 12.566370e-7
	epsilon0 = 1
	mu0 = 1
	omega = 1.
	
	leftdirichlet = 0.
	rightdirichlet = 0.
	
	Eanalytic = lambda x, t: np.sin(2*np.pi*omega*x)*np.cos(omega*t)
	Banalytic = lambda x, t: np.cos(2*np.pi*omega*x)*np.sin(omega*t)
	
	# creating spaces
	V, W = creationSpaces(degree, grid)
	
	mass = MassMatrix(V).tocsc()
	massD = MassMatrix(W, coef=True).tocsc()
	
	
	# creating startE and startB
	startE = FirstProjector(V, Eanalytic)
	startE = startE.tondarray()
	startB = FirstProjector(W, Banalytic, coef=True)
	startB = startB.tondarray()
	
	exple = startE
	explb = startB
	imple = startE
	implb = startB
	rktwoe = startE
	rktwob = startB
	rkfoure = startE
	rkfourb = startB
	
	deltaT = 1e-6
	maxiter = 30
	finalT = maxiter* deltaT
	
	t= 0.
	Eexpl =[computeEnergy(exple, explb, mass, massD, epsilon0, mu0)]
	Eimpl= [computeEnergy(exple, explb, mass, massD, epsilon0, mu0)]
	Erk2 = [computeEnergy(exple, explb, mass, massD, epsilon0, mu0)]
	Erk4 = [computeEnergy(exple, explb, mass, massD, epsilon0, mu0)]
	Energyex =[computeEnergy(exple, explb, mass, massD, epsilon0, mu0)]
	
	for i in range(0, maxiter):
		
		# computing 
		
		# euler explicit
		expl = EulerExplMaxwell1D(V, W, exple, explb,  deltaT, t + deltaT, epsilon0, mu0,  t=t, dirichlet=True)
		exple = expl.efield()
		explb = expl.bfield()
			
			
		# euler implicit
		impl = EulerImplMaxwell1D(V, W, imple, implb,  deltaT, t + deltaT, epsilon0, mu0,  t=t, dirichlet=True)
		imple = impl.efield()
		implb = impl.bfield()
		
		
		# rk2
		rktwo = SchemeRKtwo1D(V, W, rktwoe, rktwob,  deltaT, t + deltaT, epsilon0, mu0,  t=t, dirichlet=True)
		rktwoe = rktwo.efield()
		rktwob = rktwo.bfield()
		
		# rk4
		rkfour = SchemeRKfour1D(V, W, rkfoure, rkfourb,  deltaT, t + deltaT, epsilon0, mu0,  t=t, dirichlet=True)
		rkfoure = rkfour.efield()
		rkfourb = rkfour.bfield()
		
		# analytic
		efieldex = FirstProjector(V, Eanalytic, t=t)
		efieldex = efieldex.tondarray()
		bfieldex = FirstProjector(W, Banalytic, coef=True, t=t)
		bfieldex = bfieldex.tondarray()
		
		
		# computing energy
		Eexpl.append( computeEnergy(exple, explb, mass, massD, epsilon0, mu0))
		Eimpl.append( computeEnergy(imple, implb, mass, massD, epsilon0, mu0))
		Erk2.append( computeEnergy(rktwoe, rktwob, mass, massD, epsilon0, mu0))		
		Erk4.append( computeEnergy(rkfoure, rkfourb, mass, massD, epsilon0, mu0))
		Energyex.append(computeEnergy(efieldex, bfieldex, mass, massD, epsilon0, mu0))
				
		#update
		t = t +deltaT
	# ...
		
	xt = np.linspace(0., finalT, maxiter+1)
	# ploting ending 
	plt.plot(xt, Eexpl, label='euler explicit', color='y')
	plt.plot(xt, Energyex, label='analytic', color='r')
	plt.legend()
	plt.xlabel('t')
	plt.title('Evolution of the system\ndt='+str(deltaT)+' t='+str(finalT))
	plt.savefig('Images/energy_euler_expl.png')
	plt.clf()
	
	plt.plot(xt, Eimpl, label='euler implicit', color='g')
	plt.plot(xt, Energyex, label='analytic', color='r')
	plt.legend()
	plt.xlabel('t')
	plt.title('Evolution of the system\ndt='+str(deltaT)+' t='+str(finalT))
	plt.savefig('Images/energy_euler_impl.png')
	plt.clf()
	
	plt.plot(xt, Erk2, label='RK2', color='c')
	plt.plot(xt, Energyex, label='analytic', color='r')
	plt.legend()
	plt.xlabel('t')
	plt.title('Evolution of the system\ndt='+str(deltaT)+' t='+str(finalT))
	plt.savefig('Images/energy_rk2.png')
	plt.clf()
	
	plt.plot(xt, Erk4, label='RK4', color='b')
	plt.plot(xt, Energyex, label='analytic', color='r')
	plt.legend()
	plt.xlabel('t')
	plt.title('Evolution of the system\ndt='+str(deltaT)+' t='+str(finalT))
	plt.savefig('Images/energy_rk4.png')
	plt.clf()
	
	
#	plt.plot(xt[1:-1], Eexpl[1:-1], label='euler explicit', color='y')
#	plt.plot(xt[1:-1], Eimpl[1:-1], label='euler implicit', color='g')
#	plt.plot(xt[1:-1], Erk2[1:-1], label='RK2', color='c')
	plt.plot(xt[1:-1], Erk4[1:-1], label='RK4', color='b')
#	plt.plot(xt, Energyex, label='analytic', color='r')
	plt.legend()
	plt.xlabel('t')
	plt.title('Evolution of the system\ndt='+str(deltaT)+' t='+str(finalT))
	plt.savefig('Images/comp_energy.png')
	plt.clf()

