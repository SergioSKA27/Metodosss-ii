import sympy as sy
import numpy as np
import pandas as pd
from sympy.plotting.plot import MatplotlibBackend, Plot
import platform
import os
from tabulate import tabulate
import streamlit as st



def main():
    st.title("MÉTODOS NUMÉRICOS II")
    st.write("Presentación: Desarrolladores")
    st.write("")

    option = st.sidebar.selectbox("Seleccione una opción",
                                  ["Sistemas de ecuaciones no lineales",
                                   "Interpolación y ajuste de curvas",
                                   "Diferenciación e integración (datos igualmente espaciados)",
                                   "Salir"])

    if option == "Sistemas de ecuaciones no lineales":
        st.write("Opción seleccionada: Sistemas de ecuaciones no lineales")
        # Aquí puedes agregar el código correspondiente a esa opción

    elif option == "Interpolación y ajuste de curvas":
        submenu_option = st.sidebar.selectbox("Seleccione una subopción",
                                              ["Leer la tabla",
                                               "Interpolación (Diferencias divididas)",
                                               "Ajuste por spline cúbico"])

        if submenu_option == "Leer la tabla":
            st.write("Opción seleccionada: Leer la tabla")
            # Código correspondiente a Leer la tabla

        elif submenu_option == "Interpolación (Diferencias divididas)":
            st.write("Opción seleccionada: Interpolación (Diferencias divididas)")
            # Código correspondiente a Interpolación (Diferencias divididas)

        elif submenu_option == "Ajuste por spline cúbico":
            st.write("Opción seleccionada: Ajuste por spline cúbico")
            # Código correspondiente a Ajuste por spline cúbico

    elif option == "Diferenciación e integración (datos igualmente espaciados)":
        submenu_option = st.sidebar.selectbox("Seleccione una subopción",
                                              ["Leer tabla",
                                               "Diferenciación centrada",
                                               "Integración"])

        if submenu_option == "Leer tabla":
            st.write("Opción seleccionada: Leer tabla")
            # Código correspondiente a Leer tabla

        elif submenu_option == "Diferenciación centrada":
            st.write("Opción seleccionada: Diferenciación centrada")
            # Código correspondiente a Diferenciación centrada

        elif submenu_option == "Integración":
            st.write("Opción seleccionada: Integración")
            # Código correspondiente a Integración

    elif option == "Salir":
        st.warning("¡Gracias por utilizar la aplicación!")
        st.stop()

if __name__ == "__main__":
    main()
