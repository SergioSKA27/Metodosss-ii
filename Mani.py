import sympy as sy
import numpy as np
import pandas as pd
from sympy.plotting.plot import MatplotlibBackend, Plot
import platform
import os
from tabulate import tabulate
import streamlit as st

def mostrar_menu():
    print(" _________________________________________________________")
    print("|                                                         |")
    print("|                  MÉTODOS NUMÉRICOS II                    |")
    print("|_________________________________________________________|")
    print("|                                                         |")
    print("|              Presentación: Desarrolladores              |")
    print("|                                                         |")
    print("|  [1] Sistemas de ecuaciones no lineales                 |")
    print("|  [2] Interpolación y ajuste de curvas                   |")
    print("|      - [2.1] Leer la tabla                              |")
    print("|      - [2.2] Interpolación (Diferencias divididas)      |")
    print("|      - [2.3] Ajuste por spline cúbico                   |")
    print("|  [3] Diferenciación e integración (datos igualmente     |")
    print("|      espaciados)                                        |")
    print("|      - [3.1] Leer tabla                                 |")
    print("|      - [3.2] Diferenciación centrada                    |")
    print("|      - [3.3] Integración                                |")
    print("|                                                         |")
    print("|  [4] Salir                                              |")
    print("|                                                         |")
    print("|_________________________________________________________|")





if __name__ == '__main__':
    mostrar_menu()

