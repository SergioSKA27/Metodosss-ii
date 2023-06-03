import sympy as sy
import numpy as np
import pandas as pd
from sympy.plotting.plot import MatplotlibBackend, Plot
import platform
import os
from sympy.plotting import plot3d,plot3d_parametric_line





if platform.system() == 'Linux':
    CLEARW = 'clear'
elif platform.system() ==  'Windows':
    CLEARW = 'cls'


def integrate_simpson_13(f, a, b, n):
    h = (b - a) / n
    x = sy.symbols('x')
    # Calcular los valores de la función en los puntos a, b y los puntos intermedios
    x_values = [a + i * h for i in range(n+1)]
    y_values = [f.subs({x: t}).evalf() for t in x_values]

    # Verificar si el número de puntos es impar y ajustar n si es necesario
    if n % 2 != 0:
        n += 1
        h = (b - a) / n

    # Aplicar la fórmula de Simpson 1/3
    result = (h / 3) * sum(
        y_values[i] + 4 * y_values[i + 1] + y_values[i + 2]
        for i in range(0, n, 2)
    )

    print('El resultado es :', result)
    input()

def main_menu():
    os.system(CLEARW)
    print('__________________________________________________')
    print('|                                                 |')
    print('|           INTEGRACION NUMERICA                  |')
    print('|_________________________________________________|')
    print('|                                                 |')
    print('|             [1] Método                          |')
    print('|             [2] Salir                           |')
    print('|                                                 |')
    print('|_________________________________________________|')
    while True:
        try:
            op = int(input('Seleccione una opcion del menu: '))
            while op != 1 and op != 2:
                op = int(input('Seleccione una opcion del menu: '))
            break
        except:
            print('Un error ha ocurrido (Ingrese un valor numerico) :(')
            continue

    return op

def metodo():
    os.system(CLEARW)
    print(" _________________________________________________")
    print("|                                                 |")
    print("|              INTEGRACIÓN NUMÉRICA               |")
    print("|_________________________________________________|")
    print("|                                                 |")
    print("|        Selecciona una función:                  |")
    print("|                                                 |")
    print("|   1. x^4*((3+2x^2)^(1/2)/3)                     |")
    print("|   2. x^2/(x^2+4)^(1/3)                          |")
    print("|                                                 |")
    print("|_________________________________________________|")

    funcs = [sy.parse_expr('x**4*((3+2*x**2)**(1/2)/3)'),
             sy.parse_expr('x**2/(x**2+4)**(1/3)')]
    while True:
        try:
            func = int(input(': '))
            break
        except:
            print('Un error ha ocurrido (Ingrese un valor numerico) :(')
            continue


    print('La funcion seleccionada es: ')
    sy.pprint(funcs[func])
    print('')
    sy.textplot(funcs[func],-10,10)

    input('Presione enter para continuar....')


    os.system(CLEARW)
    print(" _________________________________________________")
    print("|                                                 |")
    print("|              INTEGRACIÓN NUMÉRICA               |")
    print("|_________________________________________________|")
    print("|                                                 |")
    print("|      Por favor, ingrese el intervalo [a, b]:    |")
    print("|      Separado por coma:                         |")
    print("|                                                 |")
    print("|_________________________________________________|")

    while True:
        try:
            inter = str(input(': '))
            interval = list(map(float,inter.split(',')))
            break
        except:
            print('Un error ha ocurrido (Ingrese valores numericos separados por coma) :(')
            continue

    print(" _________________________________________________")
    print("|                                                 |")
    print("|              INTEGRACIÓN NUMÉRICA               |")
    print("|_________________________________________________|")
    print("|                                                 |")
    print("| Por favor, ingrese el número de subintervalos:  |")
    print("|                                                 |")
    print("|_________________________________________________|")

    while True:
        try:
            numinter = int(input(': '))
            break
        except:
            print('Un error ha ocurrido (Ingrese un valor numerico) :(')
            continue


    os.system(CLEARW)

    print('---------------------')

    integrate_simpson_13(funcs[func], interval[0], interval[1],numinter)


    os.system(CLEARW)
    print(" _________________________________________________")
    print("|                                                 |")
    print("|              INTEGRACIÓN NUMÉRICA               |")
    print("|_________________________________________________|")
    print("|                                                 |")
    print("| ¿Desea realizar otro ajuste con otra tabla      |")
    print("| de valores?                                     |")
    print("| Por favor, responda con 's' (sí) o 'n' (no):    |")
    print("|                                                 |")
    print("|_________________________________________________|")

    try:
        newtable= str(input(': '))

        while newtable != 's' and newtable != 'n':
            newtable= str(input(': '))
        if newtable.lower() == 's':
            metodo()
        else:
            return
    except:
        return


if __name__ == "__main__":
    while True:
        o = main_menu()
        if o == 2:
            break
        else:
            metodo()
