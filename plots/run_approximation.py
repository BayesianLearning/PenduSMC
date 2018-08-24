if __name__ == '__main__' and __package__ is None:
    from os import sys, path
    sys.path.append( path.dirname( path.dirname( path.abspath(__file__) ) ) )

from optparse import OptionParser
import numpy as np
import scipy.integrate as integrate
import scipy.stats as stats
from math import sin, pi, sqrt, exp, log
from particle_approximation import ParticleApproximation

# Parse arguments
parser = OptionParser()
parser.add_option("-m", "--method", dest="method", action="store",
                  choices=["SMC", "SIS"], help="method to use (SMC or SIS)")
parser.add_option("-M", "--particles", dest="M", type="int",
                  help="number of particles to use")
parser.add_option("-f", "--filename", dest="filename",
                  help="output file name")

(options, args) = parser.parse_args()

# Setup the inverse problem
mesh = np.array([0, 1.51, 4.06, 7.06, 9.90, 12.66, 15.40, 15.58, 18.56, 21.38, 24.36])
theta_0 = [5*pi/180, 0]

def pendulum_rhs(theta, t, g):
    d_theta = theta[1]
    d2_theta = -g/7.4 * sin(theta[0])
    return [d_theta, d2_theta]

prior = stats.truncnorm(-10, 10, loc=10, scale=1)
error_model = stats.norm(loc=0,scale=0.05)

def potential(g, n):
    if n == 0:
        return prior.logpdf(g)
    sol = integrate.odeint(pendulum_rhs, theta_0, mesh[0:n+1], args=(g,))
    return error_model.logpdf(sol[1:,0]).sum()

vpotential = np.vectorize(potential)

def gaussian_proposal(x):
    proposals = stats.norm(loc=0, scale=0.25).rvs(size=x.size)
    return np.add(x, proposals, out=proposals)

# Run the approximation
approximation = ParticleApproximation(options.M, prior)

if options.method == "SMC":
    for n in range(mesh.size+1):
        importance_potential = lambda x: vpotential(x, n-1)
        target_potential = lambda x: vpotential(x, n)

        approximation.smc_update(importance_potential, target_potential, gaussian_proposal, correction_steps=5, ess_ratio=1.3)
        
        approximation.save(options.filename + '_' + str(n))
else:
    importance_potential = lambda x: prior.logpdf(x)
    target_potential = lambda x: vpotential(x, mesh.size)

    approximation.reweight(importance_potential, target_potential)
    approximation.save(options.filename)    
