# libposiform

A header-only C/C++ library for representing and manipulating solutions to pseudo-Boolean problems.

This library can be used for implementing methods for solving pseudo-Boolean optimization problem instances, where the objective function is represented as a posiform.

# Installation

Simply download the `libposiform.h` file and include it in your source.

# API

## Instance manipulation

The following structures and functions can be used to initialize and manipulate a pseudo-Boolean optimization problem instance.

### `struct posiform_instance`
A pseudo-Boolean optimization problem instance, represented as a multi-linear posiform.

### `void posiform_instance_init(struct posiform_instance *ins)`
Initializes a pseudo-Boolean optimization problem instance.

Parameters:
* `ins`: the pseudo-Boolean optimization problem instance. 

### `void posiform_instance_free(struct posiform_instance *ins)`
Frees the memory associated with a pseudo-Boolean optimization problem instance.

Parameters:
* `ins`: the pseudo-Boolean optimization problem instance. 

### `void posiform_instance_set_clause_coefficient(struct posiform_instance *ins, size_t i, double c)`
Associates a coefficient to a clause of the posiform.

If the clause does not exist, it is created.

Any previously set coefficient for the clause is overwritten.

Parameters:
* `ins`: the pseudo-Boolean optimization problem instance. 
* `i`: the index of the clause.
* `c`: the coefficient.

### `void posiform_instance_add_literal(struct posiform_instance *ins, size_t i, size_t j, bool sv)`
Adds a literal to a clause of the posiform.

If the clause does not exist, it is created.

Parameters:
* `ins`: the pseudo-Boolean optimization problem instance. 
* `i`: the index of the clause.
* `j`: the index of the variable.
* `sv`: the value of the variable which satisfies the literal.

## Solution manipulation

The following structures and functions can be used to initialize and manipulate solutions to a pseudo-Boolean optimization problem instance.

### `struct posiform_solution`
A solution to a pseudo-Boolean optimization problem instance.

This structure allows fast re-evaluation of the objective function.

### `void posiform_solution_swap(struct posiform_solution *sol, size_t k, size_t k2)`
Swaps the position of two variables in a solution.

Parameters:
* `sol`: the solution.
* `k1`: a position of x_s.
* `k2`: another position of x_s.

### `void posiform_solution_init(struct posiform_solution *sol, const struct posiform_instance *ins)`
Initializes a solution to a pseudo-Boolean optimization problem instance.

The time complexity of this function is linear on the size of the posiform.

Parameters:
* `sol`: the solution.
* `ins`: the pseudo-Boolean optimization problem instance.

### `void posiform_solution_free(struct posiform_solution *sol)`
Frees memory associated with a solution.

Parameters:
* `sol`: the solution.

### `void posiform_solution_flip(struct posiform_solution *sol, size_t j)`
Changes the value of a variable in a solution.

The time complexity of this function is linear on the number of variables whose dx value will change.

Parameters:
* `sol`: the solution.
* `j`: the index of the variable.

### `void posiform_solution_copy(const struct posiform_solution *sol, struct posiform_solution *other)`
Copies a solution to another.

Parameters:
* `sol`: the solution to copy from.
* `other`: the solution to copy to.
