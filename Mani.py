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


def main_menu_newton():
    """
    This function displays a menu for the Newton method and prompts the user to select an option, returning the selected
    option.
    :return: the user's choice from the main menu for the Newton method, either 1 or 2.
    """
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

def metodo_newton():
    """
    This function implements the Newton's method for solving systems of equations and allows the user to select a system,
    input initial approximations, set maximum iterations and error tolerance, and optionally plot the solution.
    :return: The function does not have a return statement, so it does not return anything.
    """
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
    print("| ¿Desea realizar otra aproximacion con otro      |")
    print("| sistema?                                        |")
    print("| Por favor, responda con 's' (sí) o 'n' (no):    |")
    print("|                                                 |")
    print("|_________________________________________________|")

    try:
        newtable= str(input(': '))

        while newtable != 's' and newtable != 'n':
            newtable= str(input(': '))
        if newtable.lower() == 's':
            metodo_newton()
        else:
            return
    except:
        return



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


def main_menu_interpolation():
    """
    This function displays a menu for polynomial approximation and prompts the user to select an option, returning the
    selected option.
    :return: the user's choice from the main menu, either 1 or 2.
    """
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

def metodo_interpolation():
    """
    This function performs polynomial interpolation using the Newton interpolation method, allowing the user to input the
    number of points and their values, modify them if necessary, and interpolate new values.
    :return: The function does not have a return statement, so it returns None by default.
    """
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


    input('Presione enter para continuar....')
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
            #print('-Para continuar cierre la ventana con el grafico-')
            sy.textplot(method[0],min(xs),max(xs))
            input('Presione enter para continuar....')


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
            metodo_interpolation()
        else:
            return
    except:
        return



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

def main_menu_splines():
    """
    This function displays a menu for selecting a method or exiting, and returns the user's choice.
    :return: the user's choice from the main menu for cubic spline fitting, either 1 or 2.
    """
    os.system(CLEARW)
    print('__________________________________________________')
    print('|                                                 |')
    print('|           AJUSTE POR SPLINES CÚBICOS            |')
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

def metodo_splines():
    """
    This function implements the cubic spline interpolation method to fit a curve to a set of given points and allows the
    user to modify the input points, display the resulting splines, and optionally plot the points and splines.
    :return: The function does not have a return statement, so it returns None by default.
    """
    os.system(CLEARW)
    print(" _________________________________________________")
    print("|                                                 |")
    print("|           AJUSTE POR SPLINES CÚBICOS            |")
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
    print("|          AJUSTE POR SPLINES CÚBICOS             |")
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
    print("|          AJUSTE POR SPLINES CÚBICOS             |")
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
    print("|          AJUSTE POR SPLINES CÚBICOS             |")
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

    input('Presione enter para continuar....')
    os.system(CLEARW)
    print(" _________________________________________________")
    print("|                                                 |")
    print("|          AJUSTE POR SPLINES CÚBICOS             |")
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
            #print('-Para continuar cierre la ventana con el grafico-')
            for s in range(len(method[0])):
                print('Spline ',s,' :')
                sy.textplot(method[0][s],intervals[s][0],intervals[s][1])

        input('Presione enter para continuar....')

    except:
        print('Algo salio mal al graficar :(')


    os.system(CLEARW)
    print(" _________________________________________________")
    print("|                                                 |")
    print("|          AJUSTE POR SPLINES CÚBICOS             |")
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
            metodo_splines()
        else:
            return
    except:
        return



def integracion_simpson_3_8(a, b, f, n):
    """
    This function calculates the definite integral of a function using Simpson's 3/8 rule with a given number of
    subintervals.

    :param a: The lower limit of integration
    :param b: The upper limit of integration
    :param f: a function to be integrated
    :param n: the number of subintervals used in the Simpson's 3/8 rule for numerical integration
    """

    if n % 3 != 0:
        raise ValueError(
            "El número de subintervalos debe ser divisible por 3.")

    h = (b - a) / n
    xs = a

    x = sy.symbols('x')
    suma = f.subs({x: a}).evalf()

    for i in range(1, n):
        xs += h
        coeficiente = 3 if i % 3 == 0 else 2
        suma += coeficiente * f.subs({x: xs})

    suma += f.subs({x: b}).evalf()

    integral = (3 * h / 8) * suma
    print('El resultado es :', integral)
    input()

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

def main_menu_integral():
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

def metodo_integral():
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
            if func != 1 and func != 2:
                continue
            break
        except:
            print('Un error ha ocurrido (Ingrese un valor numerico) :(')
            continue


    print('La funcion seleccionada es: ')
    sy.pprint(funcs[func-1])
    print('')
    sy.textplot(funcs[func-1],-10,10)

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
    try:
        integrate_simpson_13(funcs[func-1], interval[0], interval[1],numinter)
    except:
        try:
            integracion_simpson_3_8(interval[0],interval[1],funcs[func-1],numinter)
        except:
            print('Algo Salio Mal :(')
        pass

    input('Presione enter para continuar....')



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
            metodo_integral()
        else:
            return
    except:
        return












def mostrar_menu():
    os.system(CLEARW)
    print(" _________________________________________________________")
    print("|                                                         |")
    print("|                  MÉTODOS NUMÉRICOS II                   |")
    print("|_________________________________________________________|")
    print("|                                                         |")
    print("|              Presentación: Desarrolladores              |")
    print("|                                                         |")
    print("|  [1] Sistemas de ecuaciones no lineales                 |")
    print("|  [2] Interpolación y ajuste de curvas                   |")
    print("|      - [2.1] Interpolación (Diferencias divididas)      |")
    print("|      - [2.2] Ajuste por spline cúbico                   |")
    print("|  [3] Diferenciación e integración (datos igualmente     |")
    print("|      espaciados)                                        |")
    print("|      - [3.1] Integración                                |")
    print("|                                                         |")
    print("|  [4] Salir                                              |")
    print("|                                                         |")
    print("|_________________________________________________________|")

    while True:
        try:
            op = float(input('Seleccione una opcion del menu: '))
            while op != 1 and op != 2 and op == 2.1 and op == 2.2 and op == 2.3 and op ==3  and op == 3.1 and op == 3.2 and op ==3.3 and op ==4:
                op = float(input('Seleccione una opcion del menu: '))
            break
        except:
            print('Un error ha ocurrido (Ingrese un valor numerico) :(')
            continue
    return op





if __name__ == '__main__':
    while True:
        o = mostrar_menu()
        if o == 1:
            op = main_menu_newton()
            if op == 1:
                metodo_newton()
            else:
                continue
        elif o == 2.1:
            op = main_menu_interpolation()
            if op == 1:
                metodo_interpolation()
            else:
                continue
        elif o == 2.2:
            op = main_menu_splines()
            if op == 1:
                metodo_splines()
            else:
                continue
        elif o == 3.1:
            op = main_menu_integral()
            if op == 1:
                metodo_integral()
            else:
                continue
        elif o == 4:
            break

        else:
            print('Ingrese el numero que aparece al lado de la opcion')

