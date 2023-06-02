import sympy as sy
import numpy as np
import pandas as pd
from sympy.plotting.plot import MatplotlibBackend, Plot
from sympy.plotting import plot3d,plot3d_parametric_line
import platform
import os
from tabulate import tabulate
import streamlit as st

if platform.system() == 'Linux':
    CLEARW = 'clear'
elif platform.system() ==  'Windows':
    CLEARW = 'cls'

def jacobian(ff,symb):
    """
    It takes a vector of functions and a vector of symbols and returns the Jacobian matrix of the functions with respect to
    the symbols
    :param ff: the function
    :param symb: the symbols that are used in the function
    :return: A matrix of the partial derivatives of the function with respect to the variables.
    """
    m = []

    for i in range(0,len(ff)):
        aux  = []
        for j in range(0,len(symb)):
            aux.append(sy.diff(ff[i],symb[j]))
        m.append(aux)

    return np.array(m)

def eval_matrix(matrix , v,symb):
    """
    It takes a matrix, a list of symbols and a list of values, and returns the matrix with the symbols substituted by the
    values

    :param matrix: the matrix of the system of equations
    :param v: the vector of values for the variables
    :param symb: the symbols that will be used in the matrix
    :return: the matrix with the values of the variables substituted by the values of the vector v.
    """
    e = 0
    mm = []
    for i in range(0,len(matrix)):
        aux = []
        ev = []
        for k in range(0,len(symb)):
            ev.append((symb[k],v[k]))
        for j in range(len(matrix[i])):
            aux.append(matrix[i][j].subs(ev).evalf())
        mm.append(aux)
    return np.array(mm)

def evalVector(ff, x0,symb):
    """
    > Given a list of symbolic expressions, a list of values for the symbols, and a list of the symbols, evaluate the
    symbolic expressions at the given values

    :param ff: the vector of functions
    :param x0: initial guess
    :param symb: the symbols that are used in the symbolic expression
    :return: the value of the function at the point x0.
    """
    v = []
    for i in range(0,len(ff)):
        ev = []

        for k in range(0,len(x0)):
            ev.append((symb[k],x0[k]))

        v.append(ff[i].subs(ev).evalf())
    return np.array(v)

def NewtonMethod( ff, x0,symb ):
    """
    The function takes in a vector of functions, a vector of initial guesses, and a vector of symbols. It then calculates
    the Jacobian matrix, the Jacobian matrix evaluated at the initial guess, the inverse of the Jacobian matrix evaluated at
    the initial guess, the vector of functions evaluated at the initial guess, and then the Newton step.

    The function returns the Newton step.

    :param ff: the function we want to find the root of
    :param x0: initial guess
    :param symb: the symbols used in the function
    :return: The return value is the x_np1 value.
    """
    j = jacobian(ff,symb)
    #print("Jacobian Matrix")
    #pprint(Matrix(j))
    jev = sy.Matrix( eval_matrix(j,x0,symb))
    #print("J(",x0,")")
    #pprint(jev)

    jinv = jev.inv()
    #print("F(",x0,")")
    ffev = sy.Matrix(evalVector(np.transpose(ff),x0,symb))
    #print("J^-1(",x0,")*","F(",x0,")")
    mm = sy.Matrix(jinv)*ffev
    #pprint(mm)
    x_np1 = sy.Matrix(np.transpose(np.array(x0)))
    #pprint(x_np1-mm)
    return list(x_np1-mm)

def norm_inf(x_0,x_1):
    """
    > The function `norm_inf` takes two vectors `x_0` and `x_1` and returns the maximum absolute difference between the two
    vectors

    :param x_0: the initial guess
    :param x_1: the vector of the current iteration
    :return: The maximum difference between the two vectors.
    """
    a = [abs(x_1[i]-x_0[i]) for i in range(len(x_0))]
    return max(a)

def newton_method(ff,x_0,symbs,error,maxiter):
    """
    Given a function (x,y)$, a starting point $, and a list of symbols,
    the function will return the next point $ in the Newton's method sequence

    :param ff: the function to be minimized
    :param x_0: initial guess
    :param symbs: the symbols that we're using in the function
    :return: the final value of x_0, the list of x values, and the list of y values.
    """
    #pprint(Matrix(x_0))
    xs = []
    ys = []
    xns = [x_0]
    erros = []
    iterr = 0
    while True and iterr < maxiter:

        x_1 = NewtonMethod(ff,x_0,symbs)
        #print(x_1)
        ninf = norm_inf(x_0,x_1)
        erros.append(ninf)
        #print(ninf)

        x_0 = list(x_1)
        xns.append(tuple(x_0))
        xs.append(x_0[0])
        ys.append(x_0[1])
        if ninf < error:
            #print("Iteraciones: ",iterr)
            break
        iterr = iterr+1

    #print(x_0)
    return xns,xs,ys,erros

def get_sympy_subplots(plot:Plot):
    """
    It takes a plot object and returns a matplotlib figure object

    :param plot: The plot object to be rendered
    :type plot: Plot
    :return: A matplotlib figure object.
    """
    backend = MatplotlibBackend(plot)

    backend.process_series()
    backend.fig.tight_layout()
    return backend.plt


def main_menu():
    os.system(CLEARW)
    print('__________________________________________________')
    print('|                                                 |')
    print('|              METODO DE NEWTON                   |')
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
    print("|                                                  |")
    print("|                MÉTODO DE NEWTON                  |")
    print("|__________________________________________________|")
    print("|                                                  |")
    print("|    Por favor, seleccione un sistema a resolver:  |")
    print("|                                                  |")
    print("| 1. f(x,y) = x^2 + xy - 10 = 0                    |")
    print("|    f(x,y) = y + 3xy^2 - 50 = 0                   |")
    print("|                                                  |")
    print("| 2. f(x,y) = x^2 + y^2 - 9 = 0                    |")
    print("|    f(x,y) = -exp(x) - 2y - 3 = 0                 |")
    print("|                                                  |")
    print("| 3. f(x,y,z) = 2x^2 - 4x + y^2 + 3z^2 + 6z + 2 = 0|")
    print("|    f(x,y,z) = x^2 + y^2 - 2y + 2z^2 - 5 = 0      |")
    print("|    f(x,y,z) = 3x^2 - 12x + y^2 - 3z^2 + 8 = 0    |")
    print("|                                                  |")
    print("| 4. f(x,y,z) = x^2 - 4x + y^2 = 0                 |")
    print("|    f(x,y,z) = x^2 - x - 12y + 1 = 0              |")
    print("|    f(x,y,z) = 3x^2 - 12x + y^2 - 3z^2 + 8 = 0    |")
    print("|                                                  |")
    print("| 5. Otro sistema                                  |")
    print("|                                                  |")
    print("|__________________________________________________|")


    while True:
        try:
            sistema = int(input(': '))
            break
        except:
            print('Un error ha ocurrido (Ingrese un valor numerico) :(')
            continue


    if sistema == 1:
        x,y = sy.symbols('x,y')
        symbs = [x,y]
        system = [sy.parse_expr('x**2+x*y-10',transformations='all'),
                  sy.parse_expr('y+3*x*y**2-50',transformations='all')]

    if sistema == 2:
        x,y = sy.symbols('x,y')
        symbs = [x,y]
        system = [sy.parse_expr('x**2+y**2-9',transformations='all'),
                  sy.parse_expr('-exp(x)-2*y-3',transformations='all')]

    if sistema == 3:
        x,y,z = sy.symbols('x,y,z')
        symbs = [x,y,z]
        system = [sy.parse_expr('2*x**2-4*x+y**2+3*z**2+6*z+2',transformations='all'),
                  sy.parse_expr('x**2+y**2-2*y+2*z**2-5',transformations='all'),
                  sy.parse_expr('3x**2-12*x+y**2-3z**2+8',transformations='all')]

    if sistema == 4:
        x,y,z = sy.symbols('x,y,z')
        symbs = [x,y,z]
        system = [sy.parse_expr('x**2-4*x+y**2',transformations='all'),
                  sy.parse_expr('x**2-x-12*y+1',transformations='all'),
                  sy.parse_expr('3*x**2-12*x+y**2-3*z**2+8',transformations='all')]

    print('El sistema seleccionado es: ')

    sy.pprint(sy.Matrix(system),use_unicode=False)

    input('Presione enter para continuar....')


    os.system(CLEARW)
    print(" _________________________________________________")
    print("|                                                 |")
    print("|                MÉTODO DE NEWTON                 |")
    print("|_________________________________________________|")
    print("|                                                 |")
    print("|    Por favor, ingrese la aproximación inicial:  |")
    print("|    x_1, ... , x_n                               |")
    print("|_________________________________________________|")

    while True:
        try:
            initaprox = str(input(': '))
            aprox0 = list(map(float,initaprox.split(',')))
            print(aprox0)
            break
        except:
            print('Un error ha ocurrido (Ingrese '+str(len(system))+' valores numericos separados por coma) :(')
            continue


    os.system(CLEARW)
    print(" _________________________________________________")
    print("|                                                 |")
    print("|                MÉTODO DE NEWTON                 |")
    print("|_________________________________________________|")
    print("|                                                 |")
    print("|    Por favor, ingrese el número máximo de       |")
    print("|    iteraciones:                                 |")
    print("|                                                 |")
    print("|_________________________________________________|")

    while True:
        try:
            maxiteraciones = int(input(': '))
            break
        except:
            print('Un error ha ocurrido (Ingrese un valor numerico) :(')
            continue

    os.system(CLEARW)
    print(" _________________________________________________")
    print("|                                                 |")
    print("|                MÉTODO DE NEWTON                 |")
    print("|_________________________________________________|")
    print("|                                                 |")
    print("|    Por favor, ingrese el error de tolerancia:   |")
    print("|                                                 |")
    print("|_________________________________________________|")

    while True:
        try:
            errort = float(input(': '))
            break
        except:
            print('Un error ha ocurrido (Ingrese un valor numerico) :(')
            continue


    try:
        method = newton_method(system,aprox0,symbs,errort,maxiteraciones)
        tabb = []

        for i in range(0,len(method[1])):
            aux = list(method[0][i])
            aux.append(method[3][i])
            tabb.append(aux)

        cols = list(map(str, list(symbs)))
        cols.append('Error')
        #st.write(tabb)
        table = pd.DataFrame(tabb,columns=cols)

    except:
        print('Error al tratar de resolver el sistema :(')

    os.system(CLEARW)
    print(" _________________________________________________")
    print("|                                                 |")
    print("|                MÉTODO DE NEWTON                 |")
    print("|_________________________________________________|")
    print("|                                                 |")
    print("|   La solución está dada por:                    |")
    print("|                                                 |")
    print("|                                                 |")
    print("|_________________________________________________|")


    print(table)
    print('')
    print('La solucion esta dada por: ', (method[0][-1]))
    input('Presione enter para continuar....')

    os.system(CLEARW)
    print(" _________________________________________________")
    print("|                                                 |")
    print('|           APROXIMACION POLINOMIAL               |')
    print("|_________________________________________________|")
    print("|                                                 |")
    print("| ¿Desea graficar el sistema con solucion?        |")
    print("| Por favor, responda con 's' (sí) o 'n' (no):    |")
    print("|                                                 |")
    print("|_________________________________________________|")

    try:
        plots= str(input(': '))

        while plots != 's' and plots != 'n':
            plots= str(input(': '))
        if plots.lower() == 's':
            print('-Para continuar cierre la ventana con el grafico-')
            auxpd = plot3d(system[0],(x,-1+method[1][-1],method[1][-1]+1),(y,-1+method[2][-1],method[2][-1]+1),show=False,alpha=0.3,title='Grafica de la Solución')
            auxpd.append(plot3d(system[1],(x,-1+method[1][-1],method[1][-1]+1),(y,-1+method[2][-1],method[2][-1]+1),show=False,alpha=0.5)[0])
            pda = get_sympy_subplots(auxpd)
            zs = []
            for i in range(0,len(method[1])):
                zs.append(system[0].subs({x : method[1][i],y: method[2][i]}).evalf())
            pda.plot(method[1],method[2],zs,'o',markersize=30,color='red')
            pda.show()
    except:
        print('Ha  ocurrido un error al graficar')

    os.system(CLEARW)
    print(" _________________________________________________")
    print("|                                                 |")
    print('|           APROXIMACION POLINOMIAL               |')
    print("|_________________________________________________|")
    print("|                                                 |")
    print("| ¿Desea realizar otra aproximacion con otra tabla|")
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
