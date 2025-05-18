import math
from math import prod
from sympy import factorint, mod_inverse

alpha = 101036261437008927733827582235
beta = 37215769967503060481752795474
n = 156012931668612832596880843366
m = 156012931668612832596880843367


def factorization(n):
    factors = factorint(n)
    print("Канонічний розклад n = {}:".format(n))
    for p, l in factors.items():
        print(f"  {p}^{l}")
    print()
    return factors
factors = factorization(n)


def build_tables(alpha, n, factors):
    tables = {}
    for p, l in factors.items():
        table = []
        for j in range(p):
            exponent = (n * j) // p
            value = pow(alpha, exponent, m)
            table.append((j, value))
        tables[p] = table
    return tables
tables = build_tables(alpha, n, factors)


def find_x_i(alpha, beta, factors, tables):

    print("\nОбчислення всіх x_i для кожного p^l:")
    x_i = {}

    for p, l in factors.items():
        print(f"\nПростий дільник p = {p}, ступінь l = {l}")
        x_vals = []

        exp0 = n // p
        beta_exp = pow(beta, exp0, m)
        print(f"  Обчислюємо x0: {beta}^{n}//{p} ≡ {beta_exp} (mod {m})")

        found = False
        for j, val in tables[p]:
            if val == beta_exp:
                x_vals.append(j)
                print(f" З таблиці: x0 = {j}")
                found = True
                break
        if not found:
            print(" Не знайдено x0 у таблиці.")
            x_vals.append(None)

        for i in range(1, l):
            if None in x_vals:
                x_vals.append(None)
                print(f"Помилка!")
                continue

            exp_sum = sum(x_vals[j] * (p ** j) for j in range(i))
            alpha_inv = mod_inverse(alpha, m)
            inv_pow = pow(alpha_inv, exp_sum, m)
            temp = (beta * inv_pow) % m
            exp = n // (p ** (i + 1))
            left = pow(temp, exp, m)

            print(f"\n  Обчислюємо x_{i}:")
            print(f"    {beta} * {alpha}^(-{exp_sum}) ≡ {temp} (mod {m})")
            print(f"    ({temp})^({n}//{p**(i+1)}) ≡ {left} (mod {m})")

            found = False
            for j, val in tables[p]:
                if val == left:
                    x_vals.append(j)
                    print(f"  З таблиці: x_{i} = {j}")
                    found = True
                    break
            if not found:
                x_vals.append(None)
                print(f"  Не знайдено x_{i} у таблиці.")

        x_i[p] = x_vals
        print(f"\n Всі значення x_i для p = {p}: {x_vals}")

    return x_i
x_i = find_x_i(alpha, beta, factors, tables)


def congruences(x_i, factors):
    congruences = []

    print("\n Система:")
    for p, x_vals in x_i.items():
        l = factors[p]
        modulus = p ** l
        x_p = sum(x_vals[i] * (p ** i) for i in range(l)) % modulus
        print(f"  x =  {x_p} (mod {modulus})")
        congruences.append((x_p, modulus))

    return congruences
congruences = congruences(x_i, factors)


def ktl(congruences):
    x = 0
    M = prod(mod for _, mod in congruences)

    for ai, mi in congruences:
        Mi = M // mi
        Ni = mod_inverse(Mi, mi)
        term = ai * Mi * Ni
        x += term

    result = x % M
    print(f"\n Результат КТЛ: x = {result} (mod {M})")
    return result

x_final = ktl(congruences)

print(f"\n Відповідь {alpha} ^ {x_final} = {beta} (mod {m})")