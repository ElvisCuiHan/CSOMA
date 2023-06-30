[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg )](https://github.com/ElvisCuiHan/CSOMA/blob/main/LICENSE.md)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/ambv/black)

![CSOMA Logo](https://github.com/ElvisCuiHan/CSOMA/blob/main/csoma-main.png?width="300")
---

# Competitive Swarm Optimizer with Mutated Agents (CSOMA)
The `CSOMA` repository offers the competitive swarm optimizer with mutated agents in both **Python** and **Matlab**.

### General Information

The `Python` codes are forked from the [`pyswarms` package](https://github.com/ljvmiranda921/pyswarms) with some minor modification.
The `Matlab` codes are original.

We also provide two additional applications (*in folders section 3.5 and section 3.6*) that use CSOMA algorithm.

### Introduction to CSO-MA

[Competitive Swarm Optimizer (CSO)](https://ieeexplore.ieee.org/document/6819057) is a relatively novel swarm-based algorithm that has been proven to be very effective in solving different types of optimization problems. CSO has had successful applications to solve large and hard optimization problems. For example, [Gu et al.](https://link.springer.com/article/10.1007/s00500-016-2385-6) applied CSO to select variables for high-dimensional classification models, and [Xiong et al.](https://www.sciencedirect.com/science/article/abs/pii/S1568494618300784) used CSO to study a power system economic dispatch, which is typically a complex nonlinear multivariable strongly coupled optimization problem with equality and inequality constraints.

#### CSO Algorithm Description

The Competitive Swarm Optimizer algorithm, or CSO for short, minimizes a given function $\mathbf{x}$ over a user-specified compact space $\boldsymbol{\Omega}$ by first generating a set of candidate solutions. In our case, they take the form of a swarm of `n` particles at positions $\mathbf{x}_1, \cdots, \mathbf{x}_n$, along with their corresponding random velocities $\mathbf{v}_1, \cdots, \mathbf{v}_n$.

After the initial swarm is generated, at each iteration we randomly divide the swarm into $\left \lfloor \frac{n}{2} \right \rfloor$ pairs and compare their objective function values. We identify $\mathbf{x}^t_i$ as the winner and \(\mathbf{x}^t_j\) as the loser if these two are competed at the iteration $t$ and $\mathbf{x}^t_i < \mathbf{x}^t_j$. The winner retains the status quo, and the loser learns from the winner. The two defining equations for CSO are:

```math
v^{t+1}_{j} = R_1 ⊙ v^t_{j} + R_2 ⊙ (x^t_{i} - x^t_{j}) + φR_3 ⊙ (x̄^t - x^t_{j})
```
and
```math
x^{t+1}_{j} = x^t_{j} + v^{t+1}_{j}
```

where $\mathbf{R}_1, \mathbf{R}_2, \mathbf{R}_3$ are all random vectors whose elements are drawn from $U(0, 1)$; operation $\odot$ represents element-wise multiplication; vector $\bar{\mathbf{x}}^t$ is simply the swarm center at iteration $t$; social factor $\phi$ controls the influence of the neighboring particles to the loser, and a large value is helpful for enhancing swarm diversity (but possibly impacts convergence rate). This process iterates until some stopping criteria are met.


where $\mathbf{R}_1, \;\mathbf{R}_2, \;\mathbf{R}_3$ are all random vectors whose elements are drawn from $U(0, 1)$; operation $\otimes$ also represents element-wise multiplication; vector $\bar{\mathbf{x}}^t$ is simply the swarm center at iteration $t$; social factor $\phi$ controls the influence of the neighboring particles to the loser and a large value is helpful for enhancing swarm diversity (but possibly impacts convergence rate).  This process iterates until some stopping criteria are met.

#### CSOMA Description

An improvement on CSO and we call it the enhanced version, [**Competitive Swarm Optimizer with Mutated Agents**](https://link.springer.com/article/10.1007/s12293-020-00305-6) or, in short, **CSO-MA**. After pairing up the swarm in groups of two at each iteration, we randomly choose a loser particle $p$ as an agent, randomly pick a variable indexed as $q$ and then randomly change the value of $\textbf{x}_{pq}$ to either $\textbf{xmax}_{q}$ or $\textbf{xmin}_q$, where $\textbf{xmax}_q$ and $\textbf{xmin}_q$ represent, respectively, the upper bound and lower bound of the $q$-th variable. If the current optimal value is already close to the global optimum, this change will not hurt since we implement this experiment on a loser particle, which is not leading the movement for the whole swarm; otherwise, this chosen agent restarts a journey from the boundary and has a chance to escape from a local optimum.  

The computational complexity of CSO is $\mathcal{O}(nD)$, where $n$ is the swarm size and $D$ is the dimension of the problem. Since our modification only adds one coordinate mutation operation to each particle, its computational complexity is the same as that of CSO. The improved performance of CSO-MA over CSO-MA to find the optimum for many complex multi-dimensional benchmark functions has been validated.

### Usage

#### Matlab Codes

In Matlab, we run the following codes to obtain a high dimensional [D-optimal design](https://en.wikipedia.org/wiki/Optimal_design). The function `glm_fisher` computes the log-determinant of the Fisher information in a [generalized linear model](https://en.wikipedia.org/wiki/Generalized_linear_model).

```matlab
n = 210;
lb = -ones(1, n);
ub = ones(1, n);
phi = 0.25;
maxiter = 5000;
swarmsize = 64;
b = 2 * rand(1, n) - 1;

theta = [0.72, -0.25, 0.11, 0.91, 0.47, 0.63, -0.80, 0.86, ...
0.22, 0.19, -0.82, -0.31, 0.33, -0.12, 0.10, 0.41]'; % 8-1
theta = [-0.50, -0.10, -0.18, -0.48, 0.74, -0.63, -0.96, 0.90, ...
    0.36, -0.03, -0.93, -0.21, -0.84, -0.30, -0.67, 0.97]'; % 8-2
%theta = [0.54, -2.70, 0.37, 1.60, 2.47, -2.44, 2.42, -0.23, ...
%    -0.29, 3.00, -2.03, 1.26, -2.04, -1.86, -2.79, 0.21]'; % 9-1
%theta = [0.17, -1.01, -0.88, -2.53, 0.34, -2.01, -1.23, 2.04, ...
%    -0.82, -0.96, 1.26, -2.81, -0.17, 1.39, 1.64, -1.55]'; % 9-2
%theta = 0.1 * ones(16, 1);

obj_fun = @(b)glm_fisher(b, theta);
[value, design] = csoma(obj_fun, lb, ub, swarmsize, phi, maxiter);
```

#### Python Codes

For Python, we run the following chunk to obtain an optimal design via **CSOMA** in a trinomial dose-response model.

```python
import CSOMA as ps
import pandas as pd
import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt
from scipy.linalg import pinvh

def info(dose, theta, delete_last=True):
    """
    This function computes the information matrix in a proportional odds model.
    """
    C_trans = np.array([[1, 0, -1, 0, 0], [0, 1, 0, -1, 0], [0, 0, 0, 0, 1]])
    L = np.array([[1, 0, 0], [1, 1, 0], [0, 1, 1], [0, 0, 1], [1, 1, 1]])
    #theta = np.array([2.5062, 7.8042, -0.9795, 0])
    X = np.array([[1, 0, dose, 0], [0, 1, dose, 0], [0, 0, 0, 1]])
    A = np.array([[1, 0, 0], [-1, 1, 0], [0, -1, 2]])
    pi = A.dot(1 / (1 + np.exp(-X.dot(theta))))
    D_inv = np.diag(1 / L.dot(pi))
    M_inv = np.diag(1 / pi)
    G = np.linalg.inv(C_trans.dot(D_inv).dot(L)).dot(X)
    
    Fi = G.T.dot(M_inv).dot(G)
    
    if delete_last:
        return Fi[:-1, :-1]
    else:
        return Fi

def compute_pi(dose):
    theta = np.array([2.5062, 7.8004, -0.9791, 0])
    X = np.array([[1, 0, dose, 0], [0, 1, dose, 0], [0, 0, 0, 1]])
    A = np.array([[1, 0, 0], [-1, 1, 0], [0, -1, 2]])
    pi = A.dot(1 / (1 + np.exp(-X.dot(theta))))

    return pi

def D_optim(b, **kwargs):
    """D-optim design

    Parameters
    ----------
    b : numpy.ndarray
        sets of inputs shape :code:'(n_particles, dimensions)'
        usually for a simple logistic model, dimension is 8.

    Returns
    ----------
    numpy.ndarray
        computed cost of size :code:'(n_particles, )'
    """
    theta, = kwargs.values()
    #print(theta)

    n, d = b.shape
    loss = np.zeros(n)
    
    for i in range(n):
        m = 1e-7 * np.eye(3) #np.zeros((3, 3))
        x = b[i, :(d//2)]
        p = b[i, (d//2):]
        p = p / np.sum(p)
        
        for j in range((d//2)):
            m += p[j] * info(x[j], theta) #p[j] * (ca.dot(ca.T))#
            
        #m = np.linalg.inv(m)
        
        loss[i] = np.linalg.det(m)
        if p[-1] < 0:
            loss[i] -= 1e200
        
    return -loss

n = 10 # number of particles
d = 6  # dimension of the problem
b = np.random.random((n, d))
low = 0
upp = 12
theta = np.array([2.5062, 7.8004, -0.9791, 0])

bounds = [tuple(np.concatenate([[low] * (d//2), [0]* (d//2)])),
          tuple(np.concatenate([[upp] * (d//2), [1]* (d//2)]))]

options = {'c1': 0.5, 'c2': 0.3, 'w': .9, 'phi': .2}
optimizer = ps.single.CSOMA(n_particles=n, dimensions=d, options=options, bounds=bounds)
best_cost, best_pos = optimizer.optimize(D_optim, iters=500, theta=theta)
best_pos[(d//2):] /= sum(best_pos[(d//2):])

print("Design points:")
print(np.round(best_pos[:(d//2)], 5))
print("Design weights:")
print(np.round(best_pos[(d//2):], 5))
```
