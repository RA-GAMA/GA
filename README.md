# Genetic Algorithm

This code shows a basic-level genetic algorithm applied to the salesman trip problem, using a few coordinates and considering direct routes only.

This code is intended to learn how GA is applied to permutation cases like this one.

The code considers applications on Selection | Crossover | Mutation | Survival methods:

# Parent Selection:
- Roulette Wheel.
- Tournament.
- Rank.
- Stochastic Universal Sampling (SUS).

# Crossover:
- Swap Tails
- Swap Section
- Uniform Distribution

# Mutation:
- Insertion
- Swap
- Scramble
- Inversion

# Survival:
- Generational
- μ, λ
- μ + λ
- Steady State Replacement (SSR)
- Elitism

# RANDOMNESS:
To ensure randomness, the code considers a selection of each possible method for each process (Parent Selection, Crossover, Mutation, and Survival) in every iteration.

Additionally, mutation and crossover rates are applied.
