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

### Usage

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
