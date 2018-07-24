# PenduSMC

This folder contains all of the code necessary to reproduce the numerical experiments presented in 

Matthieu Bulté, Jonas Latz, Elisabeth Ullmann (2018)

__A practical example of non-linear Bayesian filtering of model parameters__

Submitted, ([arXiv preprint](https://arxiv.org/abs/1807.08713)).
> __Abstract:__ In this tutorial we consider the non-linear Bayesian filtering of static parameters in a time-dependent model. We outline the theoretical background and discuss appropriate solvers. We focus on particle-based filters and present Sequential Importance Sampling (SIS) and Sequential Monte Carlo (SMC). Throughout the paper we illustrate the concepts and techniques with a practical example using real-world data. The task is to estimate the gravitational acceleration of the Earth g by using observations collected from a simple pendulum. Importantly, the particle filters enable the adaptive updating of the estimate for g as new observations become available. For tutorial purposes we provide the data set and a Python implementation of the particle filters. 

To access the experiments presented in this repository, you can either:
1. Install locally the source code and required packages. How to do this is explained in the `1. Installation` and `2. Running the code` sections of this document.
2. Explore, modify and run the source code in the cloud by clicking the following badge and start a Binder session [![Binder](https://mybinder.org/badge.svg)](https://mybinder.org/v2/gh/BayesianLearning/PenduSMC/master).

## 1. Installation

### 1.1 Source code

To obtain the source code, you can clone this repository using `git` by typing the following command in your terminal
```
git clone git@github.com:BayesianLearning/PenduSMC.git
```

If you do not want to use git, GitHub provides the option to download a zip file of the repository. A link to download the zip file can be found on the home page of the repository or simply click the following link: [Dowload ZIP](https://github.com/BayesianLearning/PenduSMC/archive/master.zip)

### 1.2 Dependencies
Dependencies are managed using the `anaconda` package manager. Anaconda can be installed for any major operating system from the [Anaconda website](https://anaconda.org/).

Anaconda allows creating isolated python environment in which packages can be installed without affecting your global python installation. We provide the configuration for such an environment containing all the dependencies required for running the experiments. You can install this environment by typing the following command in your terminal
```
conda create --name pendusmc --file environment.yml
```

You can then activate the environment by typing the following command in your terminal.
```
source activate pendusmc
```
This will activate the environment for the current terminal session and won't affect any other running or future terminal sessions.

If you do not want to use Anaconda to install the project's dependencies, you can find these dependencies and their required version in the `environment.yml` file and install them on your own.

## 2. Running the code
	
This repository contains several runnable resources:
1. the `PenduSMC.ipynb` notebook contains an interactive version of the pendulum case study presented in the paper,
2. the `plots/` directory contains all the code necessary to generate the figures presented in the paper. This includes the code necessary to run the benchmarks and a notebook for creating the figures themselves. 

In both cases, you will need to activate the `pendusmc` anaconda environment (unless you have globally installed the project's dependencies). This can be done by typing the following command in your terminal	
```
source activate pendusmc
```

### 2.1 The PenduSMC notebook
The `PenduSMC` notebook was created to accompany the paper. It is recommended to keep it open while reading the paper and running the relevant cell as you progress in the paper to better understand the concept presented.

It will take you through the different steps of the pendulum inverse problem: modelling the dynamics of the pendulum, modelling the error and prior knowledge and finally approximating the solution of the inverse problem.

To run it, you can simply start Jupyter by typing in your terminal from the root of the directory
```
jupyter notebook
```

and in your browser open the `PenduSMC.ipynb` file.

### 2.2 Plots
The `plots/` directory contains the code used by the authors to create benchmarks of the numerical methods and the plots present in the paper.

The first step is to generate the benchmark data. To do this, you can run the `plots/benchmark.sh` script by typing
```
bash plots/benchmark.sh
```

You can edit the script to configure the following benchmark options:
* `CORES` (default=`4`) : number of parallel processors to use for the benchmarks.
* `TRIALS` (default=`50`) : number of trials to run per method per number of particles.
* `PARTICLES` (default=`(16 32 64 128 256 512 1024 2048 4096)`) : list of number of particles to benchmarks.
  
Warning: the script may take a while to run.

Once you are done with running the benchmarks, you can reproduce the plots presented in the paper by exploring the `plots/Plots.ipynb` notebook. To do this, start Jupyter from the root of the project by typing in your terminal
```
jupyter notebook
```
and then simply open the `Plots.ipynb` notebook in the `plots/` folder.

## Citation
```
@ARTICLE{BLU2018,
   author = {{Bult\'e}, M and {Latz}, J and {Ullmann}, E},
    title = "{A practical example of non-linear Bayesian filtering of model parameters}",
  journal = {Submitted},
     year = 2018,
}
```
