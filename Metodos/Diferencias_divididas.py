import sympy as sy
import numpy as np
import pandas as pd
from sympy.plotting.plot import MatplotlibBackend, Plot
import platform
import os


if platform.system() == 'Linux':
    CLEARW = 'clear'
elif platform.system() ==  'Windows':
    CLEARW = 'cls'


def calcular_diferencias_divididas(x, fx):
    n = len(x)
    tabla_diferencias = [[0] * n for _ in range(n)]

    # Inicializar la primera columna de la tabla con los valores de fx
    for i in range(n):
        tabla_diferencias[i][0] = fx[i]

    # Calcular las diferencias divididas
    for j in range(1, n):
        for i in range(n - j):
            tabla_diferencias[i][j] = (tabla_diferencias[i+1][j-1] - tabla_diferencias[i][j-1]) / (x[i+j] - x[i])

    return tabla_diferencias[0]



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

def diff_div(v,fx,order):
    """
    > The function takes in a list of values, a list of function values, and an order, and returns a list of divided
    differences

    :param v: the list of x values
    :param fx: the function you want to differentiate
    :param order: the order of the derivative you want to take
    :return: the difference quotient of the function f(x)
    """

    m = []

    for i in range(0,len(fx)):
        #print(fx[i])
        if i + 1 < len(fx) and i +order < len(v):
            #print(v[i+1],v[i],"/",fx[i+order]," ",fx[i])
            m.append((fx[i+1]-fx[i])/(v[i+order]-v[i]))
    return m

def divided_diff(fx,v):
    """
    The function takes in a list of x values and a list of f(x) values, and returns a list of lists of divided differences

    :param fx: the function to be interpolated
    :param v: list of x values
    :return: The divided difference table is being returned.
    """
    x = v
    nfx = fx
    m = []
    for i in range(0,len(v)-1):
        nx = diff_div(v,nfx,i+1)
        #print(nx)
        m.append(nx)
        nfx = nx

    #print(m)
    return m

def table_divided_dif(x,fx,diff):
    table = []
    for i in range(0,len(fx)):
        table.append([x[i],fx[i]])
        table.append([' ',' '])

    for i in range(1,len(table)):
        for j in range(0,len(diff)):
            for k in range(0,len(diff[j])):
                table[i].append(diff[j][k])


    return table


def print_divdiff(v):
    """
    The function `print_divdiff` prints a divided difference table given a list of values.

    :param v: The input parameter "v" is expected to be a list of lists, where the first list contains the x values and the
    second list contains the corresponding f(x) values. The remaining lists contain the divided differences of the function
    """
    table = []
    s = ''
    for i in range(0,len(v[0])):
        table.append([str(v[0][i] ),str(v[1][i]) ])
        table.append([])
    #print(table)
    for k in range(2,len(v)):
        aux = k-1
        for j in range(0,len(v[k])):
            for t in range(0,k):
                table[aux+j].append(' ')
            table[aux+j].append( str(v[k][j]))
            aux = aux + 1
    m = 0
    for i in table:
        if len(i) > m:
            m = len(i)

    for t in table:
        while len(t) != m:
            t.append('')

    return table
def Newton_interpolation(fx,v):
    """
    It takes in a list of x values and a list of f(x) values, and returns a polynomial that interpolates the points

    :param fx: a list of the function values
    :param v: list of x values
    :return: The function is being returned.
    """
    diff = divided_diff(fx,v)
    x = sy.symbols('x')

    expr = v[0]
    expr_reg = v[0]

    for i in range(0,len(diff)):
        s = diff[i][0]
        p = 1
        for k in range(0,len(v)):

            p = p*(x-v[k])
            #print(p, "p",k)
            if k == i:
                break
        s = s * p
        expr = expr + s

    for i in range(0,len(diff)):
        s = diff[i][-1]
        p = 1
        for k in range(len(v)-1,0,-1):

            p = p*(x-v[k])
            #print(p, "p",k)
            if k == i:
                break
        s = s * p
        expr_reg = expr_reg + s

    #pprint(expr)

    p = sy.plot(expr,(x,-10,10),show=False)
    #p.append(sy.plot(expr_reg,(x,-10,10),show=False)[0])
    p2 = get_sympy_subplots(p)
    p2.plot(v,fx,"o")



    #p2.show()
    difdiv = [v,fx]
    for i in diff:
        difdiv.append(i)

    return sy.expand(expr),p2,(table_divided_dif(v,fx,difdiv)),sy.expand(expr_reg)


def main_menu():
    os.system(CLEARW)
    print('__________________________________________________')
    print('|                                                 |')
    print('|           APROXIMACION POLINOMIAL               |')
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
    print('|           APROXIMACION POLINOMIAL               |')
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
    print('|           APROXIMACION POLINOMIAL               |')
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
    print('|           APROXIMACION POLINOMIAL               |')
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


    method = Newton_interpolation(fxs,xs)

    os.system(CLEARW)
    print(" _________________________________________________")
    print("|                                                 |")
    print('|           APROXIMACION POLINOMIAL               |')
    print("|_________________________________________________|")
    print("|                                                 |")
    print("|   El polinomio esta dado por:                   |")
    print("|                                                 |")
    print("|                                                 |")
    print("|_________________________________________________|")


    #divdiff = pd.DataFrame(method[2])

    with pd.option_context('display.max_rows', None,
                       'display.max_columns', None,
                       'display.precision', 3,
                       ):
        print(points)

    #with pd.option_context('display.max_rows', None,
    #                   'display.max_columns', None,
    #                   'display.precision', 3,
    #                   ):
    #    print(divdiff)
    print(calcular_diferencias_divididas(xs,fxs))
    print('El polinomio P(x) esta dado por: ')
    sy.pprint(method[0],use_unicode=False)


    input('Presione enter para continuar')
    os.system(CLEARW)
    print(" _________________________________________________")
    print("|                                                 |")
    print('|           APROXIMACION POLINOMIAL               |')
    print("|_________________________________________________|")
    print("|                                                 |")
    print("| ¿Desea graficar el polinomio y los puntos dados?|")
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




    while True:
        os.system(CLEARW)
        print(" _________________________________________________")
        print("|                                                 |")
        print("|              APROXIMACIÓN POLINOMIAL            |")
        print("|_________________________________________________|")
        print("|                                                 |")
        print("|   Por favor, ingrese un punto para interpolar:  |")
        print("|                                                 |")
        print("|_________________________________________________|")


        try:
            inter_value = float(input(': '))
            x = sy.symbols('x')
            ff = sy.lambdify(x,method[0])

            if not (inter_value >= min(xs) and inter_value <= max(xs)):
                print('Ingrese un punto que se encuentre en '+ str(min(xs)) + '<= x <=' +str(max(xs)))

            print("El valor de P("+str((inter_value))+') = '+str(ff(inter_value)))

            niter = str(input('''Desea interpolar otro valor responda con 's' (sí) o 'n' (no): '''))
            if niter.lower() == 's':
                continue
            break
        except:
            print('Ha ocurrido un error(Ingrese un valor numerico')
            continue




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

