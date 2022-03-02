#pragma once

#include <stdlib.h>
#include <stdbool.h>

/**
 * @brief A pseudo-Boolean optimization problem instance, represented as a multi-linear posiform.
 */
struct posiform_instance
{
    /** The number of variables. */
    size_t n;
    /** The number of clauses in the posiform. */
    size_t m;
    /** The array of clauses in the posiform. */
    struct
    {
        /** The coefficient of the clause in the posiform. */
        double c;
        /** The array of literals in the clause. */
        struct
        {
            /** The index of the variable of the literal. */
            size_t j;
            /** The value the variable needs to assume to satisfy the literal. */
            bool sv;
        } * literals;
        /** The number of literals in the clause. */
        size_t size;
    } * clauses;
    /** The total number of literals throughout all clauses of the posiform. */
    size_t size;
};

/**
 * @brief Initializes a pseudo-Boolean optimization problem instance.
 *
 * @param ins the pseudo-Boolean optimization problem instance.
 */
void posiform_instance_init(struct posiform_instance *ins)
{
    ins->n = ins->m = ins->size = 0;
    ins->clauses = NULL;
}

/**
 * @brief Frees the memory associated with a pseudo-Boolean optimization problem instance.
 *
 * @param ins the pseudo-Boolean optimization problem instance.
 */
void posiform_instance_free(struct posiform_instance *ins)
{
    for (size_t i = 0; i < ins->m; i++)
        free(ins->clauses[i].literals);
    free(ins->clauses);
}

/**
 * @brief Ensures that the pseudo-Boolean optimization problem instance has the specified number of clauses.
 *
 * @param ins the pseudo-Boolean optimization problem instance.
 * @param m the number of clauses.
 */
void posiform_instance_ensure_m(struct posiform_instance *ins, size_t m)
{
    if (m <= ins->m)
        return;
    __typeof__(ins->clauses) prev = ins->clauses;
    ins->clauses = (__typeof__(ins->clauses))calloc(m, sizeof(ins->clauses[0]));
    for (size_t k = 0; k < ins->m; k++)
        ins->clauses[k] = prev[k];
    free(prev);
    for (size_t k = ins->m; k < m; k++)
    {
        ins->clauses[k].c = 0;
        ins->clauses[k].size = 0;
        ins->clauses[k].literals = NULL;
    }
    ins->m = m;
}

/**
 * @brief Associates a coefficient to a clause of the posiform.
 *
 * If the clause does not exist, it is created.
 *
 * Any previously set coefficient for the clause is overwritten.
 *
 * @param ins the pseudo-Boolean optimization problem instance.
 * @param i the index of the clause.
 * @param c the coefficient.
 */
void posiform_instance_set_clause_coefficient(struct posiform_instance *ins, size_t i, double c)
{
    posiform_instance_ensure_m(ins, i + 1);
    ins->clauses[i].c = c;
}

/**
 * @brief Adds a literal to a clause of the posiform.
 *
 * If the clause does not exist, it is created.
 *
 * @param ins the pseudo-Boolean optimization problem instance.
 * @param i the index of the clause.
 * @param j the index of the variable.
 * @param sv the value of the variable which satisfies the literal.
 */
void posiform_instance_add_literal(struct posiform_instance *ins, size_t i, size_t j, bool sv)
{
    posiform_instance_ensure_m(ins, i + 1);
    __typeof__(ins->clauses[i].literals) prev = ins->clauses[i].literals;
    ins->clauses[i].literals = (__typeof__(ins->clauses[i].literals))calloc(ins->clauses[i].size + 1, sizeof(ins->clauses[i].literals[0]));
    for (size_t k = 0; k < ins->clauses[i].size; k++)
        ins->clauses[i].literals[k] = prev[k];
    free(prev);
    ins->clauses[i].literals[ins->clauses[i].size].j = j;
    ins->clauses[i].literals[ins->clauses[i].size].sv = sv;
    ins->clauses[i].size++;
    if (j >= ins->n)
        ins->n = j + 1;
    ins->size++;
}

/**
 * @brief A base solution to a pseudo-Boolean optimization problem instance.
 *
 * This structure does not allow efficient re-evaluation of the objective function.
 */
struct posiform_base_solution
{
    /** The number of variables. */
    size_t n;
    /** The solution vector. */
    bool *x;
    /** The objective function value. */
    double fx;
};

/**
 * @brief Initializes and evaluates a base solution to the pseudo-Boolean optimization problem instance.
 *
 * The solution vector is zero-initialized.
 *
 * The time complexity of this function is linear on the size of the posiform.
 *
 * @param sol the solution.
 * @param ins the pseudo-Boolean optimization problem instance.
 */
void posiform_base_solution_init(struct posiform_base_solution *sol, const struct posiform_instance *ins)
{
    sol->n = ins->n;
    sol->fx = 0;
    sol->x = (bool *)calloc(ins->n, sizeof(bool));
    for (size_t i = 0; i < ins->m; i++)
    {
        sol->fx += ins->clauses[i].c;
        for (size_t k = 0; k < ins->clauses[i].size; k++)
            if (ins->clauses[i].literals[k].sv)
            {
                sol->fx -= ins->clauses[i].c;
                break;
            }
    }
}

/**
 * @brief Frees memory associated with a base solution.
 *
 * @param sol the solution.
 */
void posiform_base_solution_free(struct posiform_base_solution *s)
{
    free(s->x);
}

/**
 * @brief A solution to a pseudo-Boolean optimization problem instance.
 *
 * This structure allows fast re-evaluation of the objective function.
 */
struct posiform_solution
{
    /** The base solution structure. */
    struct posiform_base_solution base;
    /** The changes in objective function value upon changing each variable. */
    double *dx;

    /** The coefficients of each clause of the posiform. */
    double *s_c;
    /** The numbers of unsatisfied literals in each clause of the posiform. */
    size_t *s_num_unsats;

    /** The indices where each clause's variables starts in s_x. */
    size_t *s_x_inits;
    /** The number of variables which appear in each clause. */
    size_t *s_x_sizes;
    /** The variables appearing in each clause. */
    size_t *s_x;
    /** The positions in x_s. */
    size_t *s_x_pos;

    /** The indices where each variable's clauses starts in x_s. */
    size_t *x_s_inits;
    /** The number of clauses where each variable appears in. */
    size_t *x_s_sizes;
    /** The clauses where each variable appears in. */
    size_t *x_s;
    /** The positions in s_x. */
    size_t *x_s_pos;
};

/**
 * @brief Swaps the position of two variables in a solution.
 *
 * @param sol the solution.
 * @param k1 a position of x_s.
 * @param k2 another position of x_s.
 */
void posiform_solution_swap(struct posiform_solution *sol, size_t k, size_t k2)
{
    size_t temp;
    temp = sol->s_x[sol->x_s_pos[k]];
    sol->s_x[sol->x_s_pos[k]] = sol->s_x[sol->x_s_pos[k2]];
    sol->s_x[sol->x_s_pos[k2]] = temp;
    temp = sol->s_x_pos[sol->x_s_pos[k]];
    sol->s_x_pos[sol->x_s_pos[k]] = sol->s_x_pos[sol->x_s_pos[k2]];
    sol->s_x_pos[sol->x_s_pos[k2]] = temp;
    temp = sol->x_s_pos[k];
    sol->x_s_pos[k] = sol->x_s_pos[k2];
    sol->x_s_pos[k2] = temp;
}

/**
 * @brief Initializes a solution to a pseudo-Boolean optimization problem instance.
 *
 * The time complexity of this function is linear on the size of the posiform.
 *
 * @param sol the solution.
 * @param ins the pseudo-Boolean optimization problem instance.
 */
void posiform_solution_init(struct posiform_solution *sol, const struct posiform_instance *ins)
{
    posiform_base_solution_init(&sol->base, ins);
    sol->dx = (double *)calloc(ins->n, sizeof(double));
    sol->s_c = (double *)calloc(ins->m, sizeof(double));
    sol->s_num_unsats = (size_t *)calloc(ins->m, sizeof(size_t));
    sol->s_x_inits = (size_t *)calloc(ins->m, sizeof(size_t));
    sol->s_x_sizes = (size_t *)calloc(ins->m, sizeof(size_t));
    sol->s_x = (size_t *)calloc(ins->size, sizeof(size_t));
    sol->s_x_pos = (size_t *)calloc(ins->size, sizeof(size_t));
    sol->x_s_inits = (size_t *)calloc(ins->n, sizeof(size_t));
    sol->x_s_sizes = (size_t *)calloc(ins->n, sizeof(size_t));
    sol->x_s = (size_t *)calloc(ins->size, sizeof(size_t));
    sol->x_s_pos = (size_t *)calloc(ins->size, sizeof(size_t));
    for (size_t i = 0; i < ins->m; i++)
        sol->s_c[i] = ins->clauses[i].c;

    for (size_t k = 0, i = 0; i < ins->m; k += ins->clauses[i++].size)
        sol->s_x_inits[i] = k;
    for (size_t i = 0; i < ins->m; i++)
        for (size_t k = 0; k < ins->clauses[i].size; k++)
            sol->x_s_sizes[ins->clauses[i].literals[k].j]++;
    for (size_t k = 0, j = 0; j < ins->n; k += sol->x_s_sizes[j++])
        sol->x_s_inits[j] = k;
    for (size_t j = 0; j < ins->n; j++)
        sol->x_s_sizes[j] = 0;

    for (size_t i = 0; i < ins->m; i++)
        for (size_t k = 0; k < ins->clauses[i].size; k++)
        {
            size_t j = ins->clauses[i].literals[k].j;
            sol->s_x[sol->s_x_inits[i] + sol->s_x_sizes[i]] = j;
            sol->x_s[sol->x_s_inits[j] + sol->x_s_sizes[j]] = i;
            sol->s_x_pos[sol->s_x_inits[i] + sol->s_x_sizes[i]] = sol->x_s_inits[j] + sol->x_s_sizes[j];
            sol->x_s_pos[sol->x_s_inits[j] + sol->x_s_sizes[j]] = sol->s_x_inits[i] + sol->s_x_sizes[i];
            sol->s_x_sizes[i]++;
            sol->x_s_sizes[j]++;
        }

    for (size_t i = 0; i < ins->m; i++)
    {
        size_t init = sol->s_x_inits[i], size = sol->s_x_sizes[i];
        for (size_t k = 0; k < size; k++)
            if (ins->clauses[i].literals[k].sv)
                posiform_solution_swap(sol, sol->s_x_pos[init + k], sol->s_x_pos[init + (sol->s_num_unsats[i]++)]);
        double c = sol->s_c[i];
        if (sol->s_num_unsats[i] == 0)
            for (size_t k = 0; k < size; k++)
                sol->dx[sol->s_x[init + k]] -= c;
        else if (sol->s_num_unsats[i] == 1)
            sol->dx[sol->s_x[init]] += c;
    }
}

/**
 * @brief Frees memory associated with a solution.
 *
 * @param sol the solution.
 */
void posiform_solution_free(struct posiform_solution *sol)
{
    posiform_base_solution_free(&sol->base);
    free(sol->dx);
    free(sol->s_c);
    free(sol->s_num_unsats);
    free(sol->s_x_inits);
    free(sol->s_x_sizes);
    free(sol->s_x);
    free(sol->s_x_pos);
    free(sol->x_s_inits);
    free(sol->x_s_sizes);
    free(sol->x_s);
    free(sol->x_s_pos);
}

/**
 * @brief Changes the value of a variable in a solution.
 *
 * The time complexity of this function is linear on the number of variables whose dx value will change.
 *
 * @param sol the solution.
 * @param j the index of the variable.
 */
void posiform_solution_flip(struct posiform_solution *sol, size_t j)
{
    sol->base.x[j] = 1 - sol->base.x[j];
    sol->base.fx += sol->dx[j];
    for (size_t k = sol->x_s_inits[j]; k < sol->x_s_inits[j] + sol->x_s_sizes[j]; k++)
    {
        size_t i = sol->x_s[k];
        double c = sol->s_c[i];
        size_t init = sol->s_x_inits[i], size = sol->s_x_sizes[i];
        size_t pos_in_s_x = sol->x_s_pos[k];
        size_t num_unsats_before = sol->s_num_unsats[i];
        bool is_unsat = pos_in_s_x - init < num_unsats_before;
        size_t num_unsats_after = num_unsats_before + 1 - 2 * is_unsat;
        posiform_solution_swap(sol, k, sol->s_x_pos[init + num_unsats_before - is_unsat]);
        sol->s_num_unsats[i] = num_unsats_after;
        if (num_unsats_after == 0)
        {
            for (size_t l = init; l < init + size; l++)
                sol->dx[sol->s_x[l]] -= c;
            sol->dx[j] -= c;
        }
        else if (num_unsats_after == 1 && num_unsats_before == 0)
        {
            for (size_t l = init; l < init + size; l++)
                sol->dx[sol->s_x[l]] += c;
            sol->dx[j] += c;
        }
        else if (num_unsats_after == 1 && num_unsats_before == 2)
        {
            sol->dx[sol->s_x[init]] += c;
        }
        else if (num_unsats_after == 2 && num_unsats_before == 1)
        {
            sol->dx[sol->s_x[init]] -= c;
        }
    }
}

/**
 * @brief Copies a solution to another.
 *
 * Both must be initialized with the same value of `n`.
 *
 * @param sol the solution to copy from.
 * @param other the solution to copy to.
 */
void posiform_base_solution_copy(const struct posiform_base_solution *sol, struct posiform_base_solution *other)
{
    other->fx = sol->fx;
    for (size_t i = 0; i < sol->n; i++)
        other->x[i] = sol->x[i];
}

/**
 * @brief Copies a solution to another.
 *
 * @param sol the solution to copy from.
 * @param other the solution to copy to.
 */
void posiform_solution_copy_to_base_solution(const struct posiform_solution *sol, struct posiform_base_solution *other)
{
    posiform_base_solution_copy(&sol->base, other);
}

/**
 * @brief Copies a solution to another.
 *
 * @param sol the solution to copy from.
 * @param other the solution to copy to.
 */
void posiform_solution_copy(const struct posiform_solution *sol, struct posiform_solution *other)
{
    for (size_t i = 0; i < sol->base.n; i++)
        if (sol->base.x[i] != other->base.x[i])
            posiform_solution_flip(other, i);
}