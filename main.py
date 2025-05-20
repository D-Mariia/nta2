# main.py
import argparse
from sympy import factorint, mod_inverse
from functools import reduce
from operator import mul


def prod(iterable):
    return reduce(mul, iterable, 1)


def factorization(n):
    return factorint(n)


def build_tables(alpha, n, factors, m):
    tables = {}
    for p, l in factors.items():
        table = []
        for j in range(p):
            exponent = (n * j) // p
            value = pow(alpha, exponent, m)
            table.append((j, value))
        tables[p] = table
    return tables


def find_x_i(alpha, beta, factors, tables, m, n):
    x_i = {}
    for p, l in factors.items():
        x_vals = []

        exp0 = n // p
        beta_exp = pow(beta, exp0, m)
        for j, val in tables[p]:
            if val == beta_exp:
                x_vals.append(j)
                break
        else:
            x_vals.append(None)

        for i in range(1, l):
            if None in x_vals:
                x_vals.append(None)
                continue

            exp_sum = sum(x_vals[j] * (p ** j) for j in range(i))
            alpha_inv = mod_inverse(alpha, m)
            inv_pow = pow(alpha_inv, exp_sum, m)
            temp = (beta * inv_pow) % m
            exp = n // (p ** (i + 1))
            left = pow(temp, exp, m)

            for j, val in tables[p]:
                if val == left:
                    x_vals.append(j)
                    break
            else:
                x_vals.append(None)

        x_i[p] = x_vals
    return x_i


def congruences(x_i, factors):
    congruences = []
    for p, x_vals in x_i.items():
        l = factors[p]
        modulus = p ** l
        x_p = sum(x_vals[i] * (p ** i) for i in range(l)) % modulus
        congruences.append((x_p, modulus))
    return congruences


def ktl(congruences):
    x = 0
    M = prod(mod for _, mod in congruences)

    for ai, mi in congruences:
        Mi = M // mi
        Ni = mod_inverse(Mi, mi)
        term = ai * Mi * Ni
        x += term

    return x % M


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--alpha", type=int, required=True)
    parser.add_argument("--beta", type=int, required=True)
    parser.add_argument("--modulo", type=int, required=True)

    args = parser.parse_args()

    alpha = args.alpha
    beta = args.beta
    m = args.modulo
    n = m - 1

    factors = factorization(n)
    tables = build_tables(alpha, n, factors, m)
    x_i = find_x_i(alpha, beta, factors, tables, m, n)
    congr = congruences(x_i, factors)
    result = ktl(congr)

    print(f"Результат: x ≡ {result} (mod {n})")


if __name__ == "__main__":
    main()
