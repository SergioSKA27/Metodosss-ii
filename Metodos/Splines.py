import sympy as sy
import numpy as np
import pandas as pd
from sympy.plotting.plot import MatplotlibBackend, Plot
import platform
import os
from tabulate import tabulate
import streamlit as st

if platform.system() == 'Linux':
    CLEARW = 'clear'
elif platform.system() ==  'Windows':
    CLEARW = 'cls'

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

def spline_sujeto(fx,v,fpx0,fpx1 ):
    """
    It takes a list of x values, a list of y values, and the first and second derivatives of the first and last points, and
    returns a plot of the cubic spline interpolation

    :param fx: the function values
    :param v: the x values of the points
    :param fpx0: the first derivative of the function at the first point
    :param fpx1: the derivative of the function at the last point
    """

    inter = []
    fxinter = []
    prettysplines =[]
    hi =[]
    for i in range(0,len(v)-1):
        inter.append((v[i],v[i+1]))
        fxinter.append((fx[i],fx[i+1]))

    #print(inter)
    for i in range(0,len(inter)):
        hi.append(inter[i][1]-inter[i][0])

    m = np.zeros(len(v)**2).reshape(len(fx),len(fx))
    #print(hi)
    #print(m)
    for i in range(0,len(v)):
        for j in range(0,len(v)):
            if (i == j and i == 0 and j == 0) :
                m[i][j] = 2*hi[i]
                m[i][j+1] = hi[i]
                continue
            elif (j == i and i == len(v)-1 and j == len(v)-1):
                m[i][j] = 2*hi[-1]
                m[i][j-1] = hi[-1]
                continue
            else:
                if (i == j):
                    m[i][j] = 2*(hi[i-1]+hi[i])
                    m[i][j-1] = hi[i-1]
                    m[i][j+1] = hi[i]

    b = np.zeros(len(v))
    b[0] = ((3/hi[0])*(fx[1]-fx[0]))- (3*fpx0)
    b[-1] = (3*fpx1)-((3/hi[-1])*(fx[-1]-fx[len(fx)-2]))

    for i in range(1,len(v)-1):
        b[i] = ((3/hi[i])*(fx[i+1]-fx[i]))-((3/hi[i-1])*(fx[i]-fx[i-1]))

    #print(m)
    #pprint(Matrix(b.transpose()))

    c = (sy.Matrix(m).inv())*sy.Matrix(b.transpose())
    #pprint(c)
    b = []

    for i in range(0,len(hi)):
        b.append(((fx[i+1]-fx[i])/hi[i])-((((2*c[i])+c[i+1])*hi[i])/3))

    #pprint(Matrix(b))

    d = []

    for i in range(0,len(hi)):
        d.append((c[i+1]-c[i])/(3*hi[i]))

    #pprint(Matrix(d))


    x = sy.symbols('x')
    spl = []
    for i in range(0,len(inter)):
        spl.append(sy.expand(fx[i]+ (b[i]*(x-v[i]))+(c[i]*((x-v[i])**2)) + (d[i]*((x-v[i])**3))))
        prettysplines.append(sy.latex(sy.re(sy.expand(fx[i]+ (b[i]*(x-v[i]))+(c[i]*((x-v[i])**2)) + (d[i]*((x-v[i])**3))))))

    #pprint(Matrix(spl))



    p = sy.plot(spl[0], (x,inter[0][0],inter[0][1]),show=False)

    for i in range(1, len(spl)):
        paux = sy.plot(spl[i],(x,inter[i][0],inter[i][1]),show=False)
        p.append(paux[0])


    p2 = get_sympy_subplots(p)
    p2.plot(v,fx,"o")
    return spl,p2,prettysplines,hi

def main_menu():
    os.system(CLEARW)
    print('__________________________________________________')
    print('|                                                 |')
    print('|         INTERPOLACIÓN POR SPLINES CÚBICOS       |')
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
    print("|         INTERPOLACIÓN POR SPLINES CÚBICOS       |")
    print("|_________________________________________________|")
    print("|                                                 |")
    print("|      Por favor, ingrese el número de puntos:    |")
    print("|                                                 |")
    print("|_________________________________________________|")

    while True:
        try:
            num_puntos = int(input(': '))
            break
        except:
            print('Un error ha ocurrido (Ingrese un valor numerico) :(')
            continue

    os.system(CLEARW)
    print(" _________________________________________________")
    print("|                                                 |")
    print("|        INTERPOLACIÓN POR SPLINES CÚBICOS        |")
    print("|_________________________________________________|")
    print("|                                                 |")
    print("|  Por favor, ingrese los puntos (x_i, y_i)       |")
    print("|  separados por coma 'x,y':                      |")
    print("|                                                 |")
    print("|_________________________________________________|")
    xs = []
    fxs = []
    for i in range(num_puntos):
        while True:
            try:
                ite = 'Ingrese el punto (x_'+str(i)+',y_'+str(i)+'): '
                s =str(input(ite))
                x,y = s.split(',')
                xs.append(float(x))
                fxs.append(float(y))
                break
            except Exception as e:
                print('Un error ha ocurrido (Ingrese un valores numericos separados por coma) :(')
                print(e)

    os.system(CLEARW)
    points = pd.DataFrame({'x_i': xs , 'f(x_i)': fxs})
    print(" _________________________________________________")
    print("|                                                 |")
    print("|        INTERPOLACIÓN POR SPLINES CÚBICOS        |")
    print("|_________________________________________________|")
    print("|                                                 |")
    print("|  ¿Los datos ingresados son correctos?           |")
    print("|  Por favor, responda con 's' (sí) o 'n' (no):   |")
    print("|                                                 |")
    print("|_________________________________________________|")

    with pd.option_context('display.max_rows', None,
                       'display.max_columns', None,
                       'display.precision', 3,
                       ):
        print(points)

    fix = input(': ')

    while fix != 's' and fix != 'n':
        fix = input(': ')

    if fix == 'n':
        while True:
            try:
                index = int(input('Ingrese el indice que desea modificar: '))
                if index >= 0 and index < num_puntos:
                    try:
                        nidex = str(input('Ingrese el nuevo valor del punto (x_'+str(index)+',y_'+str(index)+'): '))
                        x,y = nidex.split(',')
                        xs[index] = float(x)
                        fxs[index] = float(y)
                    except:
                        print('Un error ha ocurrido (Ingrese un valores numericos separados por coma) :(')
                        continue
                else:
                    print('Ingrese un indice valido')
                    continue

                nit = str(input('''Desea modificar otro elemento responda con 's' (sí) o 'n' (no): '''))
                if nit.lower() == 's':
                    continue




                break
            except:
                continue


    method = spline_sujeto(fxs,xs,-1,-1)
    intervals = [[xs[i-1],xs[i]] for i in range(1,num_puntos)]
    os.system(CLEARW)
    print(" _________________________________________________")
    print("|                                                 |")
    print("|        INTERPOLACIÓN POR SPLINES CÚBICOS        |")
    print("|_________________________________________________|")
    print("|                                                 |")
    print("|   Los splines están dados por:                  |")
    print("|                                                 |")
    print("|                                                 |")
    print("|_________________________________________________|")


    splines = pd.DataFrame({'Spline': method[0], 'Intervalos': intervals,'h_i': method[3]})

    with pd.option_context('display.max_rows', None,
                       'display.max_columns', None,
                       'display.precision', 3,
                       ):
        print(points)

    with pd.option_context('display.max_rows', None,
                       'display.max_columns', None,
                       'display.precision', 3,
                       ):
        print(splines)

    print('Los splines son: ')
    for j in range(len(method[0])):
        print(str(j)+': ')
        sy.pprint(method[0][j],use_unicode=False)

    input('Presione enter para continuar')
    os.system(CLEARW)
    print(" _________________________________________________")
    print("|                                                 |")
    print("|        INTERPOLACIÓN POR SPLINES CÚBICOS        |")
    print("|_________________________________________________|")
    print("|                                                 |")
    print("| ¿Desea graficar los splines y los puntos dados? |")
    print("| Por favor, responda con 's' (sí) o 'n' (no):    |")
    print("|                                                 |")
    print("|_________________________________________________|")

    try:
        plots= str(input(': '))

        while plots != 's' and plots != 'n':
            plots= str(input(': '))
        if plots.lower() == 's':
            print('-Para continuar cierre la ventana con el grafico-')
            method[1].show()

    except:
        print('Algo salio mal al graficar :(')


    os.system(CLEARW)
    print(" _________________________________________________")
    print("|                                                 |")
    print("|        INTERPOLACIÓN POR SPLINES CÚBICOS        |")
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

