import numpy as np
import matplotlib.pyplot as plt

class ParticleApproximation:
    def __init__(self, num_particles, prior, init=True):
        self.num_particles = num_particles
        self.prior = prior
        
        # create an initial Monte Carlo approximation of the prior
        if init:
            self.particles = prior.rvs(size=num_particles)
            self.weights = np.full(num_particles, 1.0/num_particles)
    
    @staticmethod
    def load(filename, prior):
        """
        load(filename, prior)
        
        Loads a previously saved particle approximation (see ParticleApproximation.save).
        
        Parameters
        ----------
        filename : string
            Name of the file in which the particle approximation was previously saved.
        prior
            Prior distribution used for this particle approximation.
        Returns
        -------
        load(filename, prior) : ParticleApproximation
            The particle approximation that was previously saved into the file.
        """
        data = np.load(filename)
        approximation = ParticleApproximation(data['particles'].size, prior, init=False)
        approximation.particles = data['particles']
        approximation.weights = data['weights']
        return approximation
    
    def save(self, filename):
        """
        save(filename)
        
        Saves the particle and weights of the particle approximation in a file for later use.
        
        Parameters
        ----------
        filename : string
            Name to use for the file in which to save to particle approximation.
        """
        np.savez(filename, particles=self.particles, weights=self.weights)
        
    def hist(self, **kwargs):
        """
        hist(**kwargs)
                
        Plots a histrogram of the approximated posterior distribution.
        
        Parameters
        ----------
        kwargs
            Additional keyword arguments to pass to Matplotlib.hist when plotting the approximation.
        """
        plt.hist(self.particles, weights=self.weights, **kwargs)
             
    def integrate(self, f):
        """
        integrate(f)
        
        Estimates the integral of the function with respect to the particle approximation.
        
        Parameters
        ----------
        f : (float) -> float
            Function to integrate.
        
        Returns
        -------
        integrate(f) : float
            Integral of this function with respect to the particle approximation.
        """
        fv = np.vectorize(f)
        return np.dot(fv(self.particles), self.weights)
    
    def sample(self, size):
        """
        sample(int) : np.ndarray
        
        Generate a pseudo random sample from the approximated distribution.
        
        Parameters
        ----------
        size : int
            Size of the generated sample.
        
        Returns
        -------
        sample(size) : np.ndarray
            A sample generated from the approximated distribution of size `size`.
        """
        return np.random.choice(self.particles, size=size, p=self.weights, replace=True)
    
    def resample(self):
        """
        resample()
        
        Resamples the set of particles from the current approximation of the target distribution. This
        will also reset the weighting of the particles to an uniform weighting.
        """
        self.particles = np.random.choice(self.particles, size=self.num_particles, p=self.weights, replace=True)
        self.weights = np.full(self.num_particles, 1/self.num_particles)

   
    def effective_sample_size(self):
        """
        effective_sample_size() : float
        
        Computes the effective sample size to approximate the effective number of samples of this
        approximation based of the deviation of the approximate weights variance.
        
        Returns
        -------
        effective_sample_size : float
            The effective sample size of the current particle approximation
        """
        return 1 / np.dot(self.weights, self.weights)

    def reweight(self, importance_potential, target_potential):
        """
        reweight(importance_potential, target_potential)

        Given the potential functions of a importance and target distribution, and assuming that the
        current ParticleApproximation approximates the importance distribution, this function updates
        the weights of the particles to approximate the target distribution.

        Parameters
        ----------
        importance_potential : (np.ndarray) -> np.ndarray
            Potential function of the importance distribution.            
        target_potential : (np.ndarray) -> np.ndarray
            Potential function of the target distribution.            
        """
        # Update is done in log-scale
        self.weights = np.log(self.weights)
        
        # Compute log-importance-weights and update current weights
        importance_weights = target_potential(self.particles) - importance_potential(self.particles)
        self.weights += importance_weights
        
        # Return to linear-scale to normalize weights
        self.weights = np.exp(self.weights)
        self.weights /= self.weights.sum()

    def mh_correction(self, target_potential, proposal_kernel, n_steps):
        """
        mh_correction(target_potential, proposal_kernel, n_steps) : float
        
        Improve the particle approximation by sampling several steps of the Metropolis-Hastings
        Markov Chain constructed using a given proposal kernel and the with the given target
        distribution as limiting distribution.
        
        Parameters
        ----------
        target_potential: (np.ndarray) -> np.ndarray
            Potential function of the target distribution.
        proposal_kernel: (np.ndarray) -> np.ndarray
            Function to sample proposals for the Metropolis-Hastings algorithm, conditioned
            on the current particles, passed as a parameter.
        n_steps : int, optional
            Number of steps of sampling from the Markov Chain.
            
        Returns
        -------
        mh_correction : float
            The average acceptance rate over the `n_steps`.
        """
        total_accepted = 0
        
        # Sample from the proposal kernel, conditioned on currect particles
        proposals = proposal_kernel(self.particles)
        
        proposal_potentials = target_potential(proposals)
        current_potentials = target_potential(self.particles)
        
        for i in range(n_steps):
            # Compute the log acceptance ratio
            potential_ratio = proposal_potentials - current_potentials
            prior_ratio = self.prior.logpdf(proposals) - self.prior.logpdf(self.particles)
            acceptance_ratio = np.exp(potential_ratio + prior_ratio)

            # Randomly accept the transitions based on the log acceptance ratio
            accepted = np.random.uniform(size=self.num_particles) < acceptance_ratio
            self.particles[accepted] = proposals[accepted]
            
            total_accepted += np.sum(accepted)
            
            # Recompute necessary potentials for next step
            if i < n_steps - 1:
                # Update current potentials
                current_potentials[accepted] = proposal_potentials[accepted]
                
                # Sample new proposals
                proposals = proposal_kernel(self.particles)
                proposal_potentials = target_potential(proposals)

        return total_accepted / (n_steps * self.num_particles)
                
    def smc_update(self, importance_potential, target_potential, proposal_kernel, correction_steps, ess_ratio):
        """
        smc_update(importance_potential, target_potential, proposal_kernel, correction_steps, ess_ratio)
        
        Updates the current particle approximation using the SMC algorithm with Metropolis-Hastings MCMC correction
        steps.
        
        Parameters
        ----------
        importance_potential : (np.ndarray) -> np.ndarray
            Potential function of the importance distribution.       
        target_potential: (np.ndarray) -> np.ndarray
            Potential function of the target distribution.
        proposal_kernel: (np.ndarray) -> np.ndarray
            Function to sample proposals for the Metropolis-Hastings algorithm, conditioned
            on the current particles, passed as a parameter.
        correction_steps : int, optional
            Number of steps of sampling from the Markov Chain.
        ess_ratio: int, optional
            Specifies the lower bound on the effective sample size. A resampling step
            is then performed if `ESS < self.num_particles/ess_ratio`.
            
        Returns
        -------
        smc_update : (float, float)
            (Acceptance ratio of the correction steps, ESS after reweighting)
        """
        self.reweight(importance_potential, target_potential)
        acceptance_ratio = self.mh_correction(target_potential, proposal_kernel, n_steps=correction_steps)

        ess = self.effective_sample_size()
        if ess < self.num_particles/ess_ratio:
            self.resample()
            
        return (acceptance_ratio, ess)
