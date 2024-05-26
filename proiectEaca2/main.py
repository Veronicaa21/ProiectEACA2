import numpy as np
import cmath
import math
import matplotlib.pyplot as plt

a = complex(input("Dati valoarea coeficientului a: "))
b = complex(input("Dati valoarea coeficientului b: "))
c = complex(input("Dati valoarea coeficientului c: "))
d = complex(input("Dati valoarea coeficientului d: "))


def ComplexToString(z):
    x = round(z.real, 3)
    y = round(z.imag, 3)
    if y == 0:
        return str(x)
    if x == 0:
        return str(y) + "j"
    return str(x) + ("+" if y > 0 else "") + str(y) + "j"


def RezolvareEcuatie(a, b, c, d):
    assert a != 0, "Coeficientul a nu poate fi zero."

    p = c / a - b * b / (3.0 * a * a)
    q = 2.0 * b * 3 / (27.0 * a * 3) - b * c / (3.0 * a * a) + d / a
    Delta = q * q + 4.0 * p ** 3 / 27.0
    w = (-q + cmath.sqrt(Delta)) / 2

    rho, theta = cmath.polar(w)
    radacini = [cmath.rect(math.pow(rho, 1 / 3), (theta + 2 * k * math.pi) / 3) for k in range(3)]
    r = [u - p / (3 * u) - b / (3.0 * a) for u in radacini]

    return r


def FunctiaDeGrad3(a, b, c, d, x):
    return a * x ** 3 + b * x ** 2 + c * x + d


def F(a, b, c, d, r):
    x = np.linspace(-10, 10, 400)
    y = FunctiaDeGrad3(a, b, c, d, x)

    plt.figure(figsize=(10, 6))
    plt.plot(x, y, label=f'{a}x^3 + {b}x^2 + {c}x + {d}')

    real_r = [r1.real for r1 in r if abs(r1.imag) < pow(10, -5)]
    for r1 in real_r:
        plt.plot(r1, 0, 'ro')
        plt.text(r1, 0, ComplexToString(complex(r1, 0)), fontsize=12, ha='right')

    for r1 in r:
        plt.plot(r1.real, r1.imag, 'bo')
        plt.text(r1.real, r1.imag, ComplexToString(r1), fontsize=12, ha='right')

    plt.axhline(0, color='black', linewidth=0.5)
    plt.axvline(0, color='black', linewidth=0.5)
    plt.grid(color='gray', linestyle='--', linewidth=0.5)
    plt.title('Graficul functiei de grad 3')
    plt.xlabel('Partea reala / x')
    plt.ylabel(' Partea imaginara / f(x)')
    plt.legend()
    plt.show()


r = RezolvareEcuatie(a, b, c, d)
F(a, b, c, d, r)