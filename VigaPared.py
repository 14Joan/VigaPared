import sys
import csv
import datetime
import numpy as np
import matplotlib
matplotlib.use('QtAgg')  
import matplotlib.pyplot as plt
import scienceplots
from numpy.linalg import LinAlgError
from PyQt6.QtWidgets import (QApplication, QWidget, QHBoxLayout, QVBoxLayout,
                             QLabel,QLineEdit, QPushButton, QTableWidget,
                             QTableWidgetItem, QComboBox, QMessageBox)
from PyQt6.QtGui import QFont
from PyQt6.QtCore import Qt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas 
from matplotlib.figure import Figure
from Graficar import Graficar
from Calcular import Calcular
from Auxiliar import Auxiliar

class Ventana(QWidget):
    def __init__(self):
        super().__init__()
        self.initializeUI()
        # Almacenamiento de Valores
        self.AsH = None  # Acedro de refuerzo horizontal [cm2]
        self.AsV = None  # Acedro de refuerzo vertical [cm2]
        self.U = None   # Desplazamientos [mm]
        self.F = None  # Niveles de carga [kN]
        self.R = None  # Reacciones y fuerzas en nudos [kN]
        self.eN = None  # Deformaciones en el hormigón [adimensional]
        # Deformaciones principales en el hormigón [adimensional]
        self.eP = None
        # Deformaciones en el refuerzo horizontal [adimensional]
        self.esh = None
        self.esv = None
        # Deformaciones en el refuerzo vertical [adimensional]
        self.NF = None  # Número de fisuras
        self.SN = None  # Tensiones en el hormigón [MPa]
        self.SP = None  # Tensiones principales en el hormigón [MPa]
        self.Sn = None  # Tensiones en nudos [MPa]
        self.Ssh = None  # Tensiones en el refuerzo horizontal [MPa]
        self.Ssv = None  # Tensiones en el refuerzo vertical [MPa]
        self.Carga_Def = None 
        # Carga, deflexión, número de iteraciones y precisión

    def initializeUI(self):
        self.setFixedSize(1100, 720)
        self.setWindowTitle("ANÁLISIS Y DISEÑO DE VIGAS PARED")
        self.setUpMainWindow()
        self.show()

    def setUpMainWindow(self):
        # Título
        label_titulo = QLabel("<b>ANÁLISIS Y DISEÑO DE VIGAS PARED<b>",
                              self)
        label_titulo.setFont(QFont("Arial", 18))
        label_titulo.setAlignment(Qt.AlignmentFlag.AlignHCenter)
        # DIMENSIONES
        label_dimensiones = QLabel("<b>DIMENSIONES<b>", self)
        label_altura = QLabel("Altura (m)", self)
        label_altura.setFixedWidth(70)
        self.edit_altura = QLineEdit(self)
        label_Ydivisiones = QLabel("Divisiones Vert.", self)
        label_Ydivisiones.setFixedWidth(90)
        self.edit_Ydivisiones = QLineEdit(self)
        label_longitud = QLabel("Longitud (m)", self)
        label_longitud.setFixedWidth(70)
        self.edit_longitud = QLineEdit(self)
        label_Xdivisiones = QLabel("Divisiones Horiz.", self)
        label_Xdivisiones.setFixedWidth(90)
        self.edit_Xdivisiones = QLineEdit(self)
        label_espesor = QLabel("Espesor (cm)", self)
        label_espesor.setFixedWidth(70)
        self.edit_espesor = QLineEdit(self)
        # MATERIALES
        label_materiales = QLabel("<b>PROPIEDADES DE LOS MATERIALES<b>",
                                  self)
        label_hormigon = QLabel("Tipo de Hormigón (MPa)", self)
        label_hormigon.setFixedWidth(135)
        self.edit_hormigon = QLineEdit(self)
        label_acero = QLabel("Tipo de Acero (MPa)", self)
        label_acero.setFixedWidth(135)
        self.edit_acero = QLineEdit(self)
        # CARGAS
        label_cargas = QLabel("<b>CARGAS (kN/m)<b>", self)
        self.tabla_cargas = QTableWidget(2, 5)
        self.tabla_cargas.setHorizontalHeaderLabels(
            ["q1", "q2", "X1 (m)", "X2 (m)", "Y (m)"])
        self.tabla_cargas.setFixedHeight(110)
        for i in range(0, 6):
            self.tabla_cargas.setColumnWidth(i, 60)
        # APOYOS
        label_apoyos = QLabel("<b>APOYOS<b>", self)
        self.tabla_apoyos = QTableWidget(2, 4)
        self.tabla_apoyos.setHorizontalHeaderLabels(
            ["X1 (m)", "X2 (m)", "Y (m)", "Restricción"])
        self.tabla_apoyos.setFixedHeight(110)
        for i in range(0, 3):
            self.tabla_apoyos.setColumnWidth(i, 65)
        combo_apoyos1 = QComboBox()
        combo_apoyos1.addItems(["Ambas", "Vertical", "Horizontal"])
        self.tabla_apoyos.setCellWidget(0, 3, combo_apoyos1)
        combo_apoyos2 = QComboBox()
        combo_apoyos2.addItems(["Ambas", "Vertical", "Horizontal"])
        self.tabla_apoyos.setCellWidget(1, 3, combo_apoyos2)
        # CUANTÍAS
        label_cuantias = QLabel("<b>CUANTÍAS (0/00)<b>", self)
        self.tabla_cuantias = QTableWidget(5, 1)
        self.tabla_cuantias.setVerticalHeaderLabels(["Horizontal (0-15)%H",
                                                "Horizontal (15-55)%H",
                                                "Horizontal (55-100)%H",
                                                "Vertical (0-15,85-100)%L",
                                                "Vertical (15-85)%L"])
        self.tabla_cuantias.horizontalHeader().setVisible(False)
        self.tabla_cuantias.setColumnWidth(0, 165)
        # BOTONES
        boton_graficar = QPushButton("GRAFICAR", self)
        boton_graficar.clicked.connect(self.graficar_viga)
        boton_calcular = QPushButton("CALCULAR", self)
        boton_calcular.clicked.connect(self.calcular_viga)
        boton_salir = QPushButton("CERRAR", self)
        boton_salir.clicked.connect(self.close)
        # RESULTADOS GRAFICOS
        label_graficos = QLabel("<b>GRÁFICOS DE: <b>", self)
        combo_graficos = QComboBox()
        combo_graficos.addItems(["Esquema Estructural",
                                 "Estructura Deformada",
                                 "Tensiones X Hormigón (MPa)",
                                 "Tensiones Y Hormigón (MPa)",
                                 "Tensiones XY Hormigón (MPa)",
                                 "Tensiones Principales Hormigón (MPa)",
                                 "Patrón de Fisuración",
                                 "Carga vs Deflexión"])
        combo_graficos.currentIndexChanged.connect(self.mostrar_graficos)
        label_etapas = QLabel("Etapa de cargado: ", self)
        self.combo_etapas = QComboBox()
        self.combo_etapas.addItems([])
        self.figura = Figure()
        self.lienzo = FigureCanvas(self.figura)
        # RESULTADOS NUMERICOS
        label_resultados_num = QLabel("Mostrar resultados de: ", self)
        combo_resultados_num = QComboBox()
        combo_resultados_num.addItems(["Desplazamientos X (mm)",
                                    "Desplazamientos Y (mm)",
                                    "Fuerzas X (kN)", "Fuerzas Y (kN)",
                                    "Tensiones X Hormigón (MPa)",
                                    "Tensiones Y Hormigón (MPa)",
                                    "Tensiones XY Hormigón (MPa)",
                                    "Tensiones Máximas Hormigón (MPa)",
                                    "Tensiones Mínimas Hormigón (MPa)",
                                    "Ángulo Principal Hormmigón (º)",
                                    "Deformaciones Máximas Hormigón (0/00)",
                                    "Deformaciones Mínimas Hormigón (0/00)",
                                    "Acero de Refuerzo Horizontal (cm2)",
                                    "Acero de Refuerzo Vertical (cm2)",
                                    "Deformaciones X Acero (0/00)",
                                    "Deformaciones Y Acero (0/00)",
                                    "Tensiones X Acero (MPa)",
                                    "Tensiones Y Acero (MPa)",
                                    "Carga vs Deflexión"])
        combo_resultados_num.currentIndexChanged.connect(
            self.mostrar_datos_tabla)
        self.tabla_resultados_num = QTableWidget(1, 1)
        self.tabla_resultados_num.setFixedWidth(750)
        self.tabla_resultados_num.setFixedHeight(200)
        # ASIGNANDO POSICIONES EN LA INTERFAZ
        # LAYOUT DE DATOS
        self.QHBox_Altura = QHBoxLayout()
        self.QHBox_Altura.addWidget(label_altura)
        self.QHBox_Altura.addWidget(self.edit_altura)
        self.QHBox_Altura.addWidget(label_Ydivisiones)
        self.QHBox_Altura.addWidget(self.edit_Ydivisiones)
        self.QHBox_Longitud = QHBoxLayout()
        self.QHBox_Longitud.addWidget(label_longitud)
        self.QHBox_Longitud.addWidget(self.edit_longitud)
        self.QHBox_Longitud.addWidget(label_Xdivisiones)
        self.QHBox_Longitud.addWidget(self.edit_Xdivisiones)
        self.QHBox_Espesor = QHBoxLayout()
        self.QHBox_Espesor.addWidget(label_espesor)
        self.QHBox_Espesor.addWidget(self.edit_espesor)
        self.QHBox_Espesor.addSpacing(175)
        self.QHBox_Hormigon = QHBoxLayout()
        self.QHBox_Hormigon.addWidget(label_hormigon)
        self.QHBox_Hormigon.addWidget(self.edit_hormigon)
        self.QHBox_Acero = QHBoxLayout()
        self.QHBox_Acero.addWidget(label_acero)
        self.QHBox_Acero.addWidget(self.edit_acero)
        self.QHBox_Botones = QHBoxLayout()
        self.QHBox_Botones.addWidget(boton_graficar)
        self.QHBox_Botones.addWidget(boton_calcular)
        self.QHBox_Botones.addWidget(boton_salir)
        # ASIGNANDO POSICIONES EN EL LAYOUT DE DATOS
        self.QVBox_Datos = QVBoxLayout()
        self.QVBox_Datos.addWidget(label_dimensiones)
        self.QVBox_Datos.addLayout(self.QHBox_Altura)
        self.QVBox_Datos.addLayout(self.QHBox_Longitud)
        self.QVBox_Datos.addLayout(self.QHBox_Espesor)
        self.QVBox_Datos.addWidget(label_materiales)
        self.QVBox_Datos.addLayout(self.QHBox_Hormigon)
        self.QVBox_Datos.addLayout(self.QHBox_Acero)
        self.QVBox_Datos.addWidget(label_cargas)
        self.QVBox_Datos.addWidget(self.tabla_cargas)
        self.QVBox_Datos.addWidget(label_apoyos)
        self.QVBox_Datos.addWidget(self.tabla_apoyos)
        self.QVBox_Datos.addWidget(label_cuantias)
        self.QVBox_Datos.addWidget(self.tabla_cuantias)
        self.QVBox_Datos.addLayout(self.QHBox_Botones)
        # LAYOUT DE RESULTADOS
        self.QHBox_Graficos = QHBoxLayout()
        self.QHBox_Graficos.addSpacing(80)
        self.QHBox_Graficos.addWidget(label_graficos)
        self.QHBox_Graficos.addWidget(combo_graficos)
        self.QHBox_Graficos.addSpacing(40)
        self.QHBox_Graficos.addWidget(label_etapas)
        self.QHBox_Graficos.addWidget(self.combo_etapas)
        self.QHBox_Graficos.addSpacing(80)
        self.QHBox_Resultados_num = QHBoxLayout()
        self.QHBox_Resultados_num.addSpacing(150)
        self.QHBox_Resultados_num.addWidget(label_resultados_num)
        self.QHBox_Resultados_num.addWidget(combo_resultados_num)
        self.QHBox_Resultados_num.addSpacing(150)
        # ASIGNANDO POSICIONES EN EL LAYOUT DE RESULTADOS
        self.QVBox_Resultados = QVBoxLayout()
        self.QVBox_Resultados.addLayout(self.QHBox_Graficos)
        self.QVBox_Resultados.addWidget(self.lienzo)
        self.QVBox_Resultados.addLayout(self.QHBox_Resultados_num)
        self.QVBox_Resultados.addWidget(self.tabla_resultados_num)
        # ASIGNANDO POSICIONES EN EL LAYOUT PRINCIPAL
        self.QHBox_11 = QHBoxLayout()
        self.QHBox_11.addLayout(self.QVBox_Datos)
        self.QHBox_11.addLayout(self.QVBox_Resultados)
        self.QVBox_1 = QVBoxLayout()
        self.QVBox_1.addWidget(label_titulo)
        self.QVBox_1.addLayout(self.QHBox_11)
        self.setLayout(self.QVBox_1)

    def validar_datos(self):
        self.figura.clear()
        ax = self.figura.add_subplot(111)
        # Datos de Cargas
        matriz_cargas = []
        for fila in range(self.tabla_cargas.rowCount()):
            fila_cargas = []
            if self.tabla_cargas.item(fila, 1) == None:
                break
            for columna in range(self.tabla_cargas.columnCount()):
                item = self.tabla_cargas.item(fila, columna)
                fila_cargas.append(item.text())
            matriz_cargas.append(fila_cargas)
        M_cargas = matriz_cargas
        # Datos apoyos
        matriz_apoyos = []
        vector_apoyos = []
        for fila in range(self.tabla_apoyos.rowCount()):
            fila_apoyos = []
            if self.tabla_apoyos.item(fila, 0) == None:
                break
            else:
                celda = self.tabla_apoyos.cellWidget(fila, 3)
                vector_apoyos.append(celda.currentText())
                for columna in range(self.tabla_apoyos.columnCount()-1):
                    item = self.tabla_apoyos.item(fila, columna)
                    fila_apoyos.append(item.text())
                matriz_apoyos.append(fila_apoyos)
        M_apoyos = matriz_apoyos
        V_apoyos = vector_apoyos
        # Datos cuantias
        vector_cuantias = []
        for fila in range(self.tabla_cuantias.rowCount()):
            if self.tabla_cuantias.item(fila, 0) == None:
                break
            else:
                item = self.tabla_cuantias.item(fila, 0)
                vector_cuantias.append(item.text())
        V_cuantias = vector_cuantias
        # Verificando validez de datos
        try:
            # Graficando la viga
            H = float(self.edit_altura.text())
            L = float(self.edit_longitud.text())
            NX = float(self.edit_Xdivisiones.text())
            NY = float(self.edit_Ydivisiones.text())
            t = float(self.edit_espesor.text())
            fc = 0.53*float(self.edit_hormigon.text())
            fy = float(self.edit_acero.text())
            if M_cargas == []:
                QMessageBox.information(
                    self, "Datos incompletos", "Introduzca datos de cargas",
                    QMessageBox.StandardButton.Ok)
            elif (M_apoyos or V_apoyos) == []:
                QMessageBox.information(
                    self, "Datos incompletos", "Introduzca datos de apoyos",
                    QMessageBox.StandardButton.Ok)
            elif V_cuantias == []:
                QMessageBox.information(
                    self, "Datos incompletos", "Introduzca datos del refuerzo",
                    QMessageBox.StandardButton.Ok)
            else:
                return H, L, NX, NY, t, fc, fy, M_cargas, M_apoyos, V_apoyos, \
                        V_cuantias, ax
        except ValueError:
            QMessageBox.information(
                self, "Datos incompletos", "Introduzca todas las dimensiones",
                QMessageBox.StandardButton.Ok)

    def graficar_viga(self):
        graficar = Graficar()
        plt.style.use(['science', 'notebook', 'grid'])
        try:
            H, L, NX, NY, t, fc, fy, M_cargas, M_apoyos, V_apoyos, \
                V_cuantias, ax = self.validar_datos()
            graficar.GraficarViga(H, L, NX, NY, ax, self.lienzo)
            graficar.GraficarCargas(H, L, NX, NY, M_cargas, ax, self.lienzo)
            graficar.GraficarApoyos(
                H, L, NX, NY, M_apoyos, V_apoyos, ax, self.lienzo)
        except TypeError or ValueError:
            QMessageBox.information(
                self, "Datos incompletos", "Introduzca todos los datos",
                QMessageBox.StandardButton.Ok)

    def calcular_viga(self):
        graficar = Graficar()
        calcular = Calcular()
        auxiliar = Auxiliar()
        # Configuracion de grafico
        figura, ejes = plt.subplots()
        plt.xlabel('Deflexión máxima (mm)')  # Etiqueta del eje x
        plt.ylabel('Carga aplicada')  # Etiqueta del eje y
        # Validando datos
        H, L, NX, NY, t, fc, fy, M_cargas, M_apoyos, V_apoyos,\
            V_cuantias,ax = self.validar_datos()
        # Cantidad de nudos y elementos
        n = int(NX*NY + NX + NY + 1)     # Numero de nudos
        N = int(NX*NY)  # Número de elementos
        # Número de elementos de acero de refuerzo horizontal
        Nsh = int((NY+1)*NX)
        # Número de elementos de acero de refuerzo vertical
        Nsv = int(NY*(NX+1))
        # Áreas de refuerzo
        self.AsH, self.AsV = calcular.CalcularAreas(
            H, L, NY, NX, t, V_cuantias)
        # Grados de libertad
        GLF = auxiliar.GradosLibertadFijos(H, L, NY, NX, M_apoyos, V_apoyos)
        # Cargas
        PPropio = auxiliar.CargasVolumetricas(H, L, NY, NX, t)
        CExterna = auxiliar.CargasSuperficales(H, L, NY, NX, M_cargas)
        NivCS = np.round(np.arange(0.01, 1.21, 0.01), 2)
        # Almacenamiento de Valores
        self.F = np.array([np.zeros(2*n)])
        self.U = np.array([np.zeros(2*n)])
        self.R = np.array([np.zeros(2*n)])
        self.eN = np.array([np.zeros((N, 3))])
        self.eP = np.array([np.zeros((N, 3))])
        self.esh = np.array([np.zeros(Nsh)])
        self.esv = np.array([np.zeros(Nsv)])
        self.NF = np.array([np.zeros((N, 2))])
        self.SN = np.array([np.zeros((N, 3))])
        self.SP = np.array([np.zeros((N, 3))])
        self.Sn = np.array([np.zeros((n, 3))])
        self.Ssh = np.array([np.zeros(Nsh)])
        self.Ssv = np.array([np.zeros(Nsv)])
        # Carga, deflexión, número de iteraciones y precisión
        self.Carga_Def = np.zeros((1, 4))
        # Graficando
        plt.ion()  # Modo interactivo
        line, = ejes.plot([], [], 'o--', color='purple', lw=0.8, ms=8)
        x, y = self.Carga_Def[-1, 0], self.Carga_Def[-1, 1]
        ejes.set_xlim(0, 1)
        ejes.set_ylim(0, 1)
        # Añadir x al conjunto de datos
        line.set_xdata(np.append(line.get_xdata(), x))
        # Añadir y al conjunto de datos
        line.set_ydata(np.append(line.get_ydata(), y))
        plt.pause(0.01)
        # Actualizando lista desplegable
        self.combo_etapas.clear()
        self.ListaCombo = np.array([])
        for j in range(NivCS.shape[0]-1):
            try:
                self.F = np.vstack(
                    (self.F, np.floor(1/(j+1))*PPropio + NivCS[j]*CExterna))
                # Cargando
                # Condiciones iniciales
                DFj = self.F[j+1] - self.F[j]
                Uj = self.U[j]
                Rj = self.R[j]
                eNj = self.eN[j]
                eshj = self.esh[j]
                esvj = self.esv[j]
                NFj = self.NF[j]
                SNj = self.SN[j]
                SPj = self.SP[j]
                Snj = self.Sn[j]
                Sshj = self.Ssh[j]
                Ssvj = self.Ssv[j]
                # Condiciones iniciales para iteración
                i = 0
                F00 = np.zeros(2*n)
                # Iterando
                while True:
                    DF0_i = DFj - F00
                    K0_i = calcular.CalcularRigidez(
                        H, L, NY, NX, t, fc, SPj, NFj, 'tan', fy, eshj, esvj,
                            self.AsH, self.AsV)
                    DU0_i = calcular.CalcularDesplazamientosFinal(
                        H, L, NY, NX, M_apoyos, V_apoyos, DF0_i, K0_i)
                    U0_i = Uj + DU0_i
                    DeN0_i, Desh0_i, Desv0_i = calcular.CalcularDeformaciones(
                        H, L, NY, NX, DU0_i, self.AsH, self.AsV)
                    eN0_i = eNj + DeN0_i
                    eP0_i = calcular.CalcularDeformacionesPrincipales(
                        NY, NX, eN0_i)
                    esh0_i = eshj + Desh0_i
                    esv0_i = esvj + Desv0_i
                    DSN0_i, NF0_i, DSn0_i = calcular.CalcularTensionesHormigon(
                        H, L, NY, NX, fc, DU0_i, SPj, NFj, 'sec')
                    Sshj, Ssvj = calcular.CalcularTensionesAcero(
                        NY, NX, fy, eshj, esvj)
                    Ssh0_i, Ssv0_i = calcular.CalcularTensionesAcero(
                        NY, NX, fy, esh0_i, esv0_i)
                    SN0_i = SNj + DSN0_i
                    Sn0_i = Snj + DSn0_i
                    SP0_i = calcular.CalcularTensionesPrincipales(
                        NY, NX, SN0_i)
                    DSsh0_i = Ssh0_i - Sshj
                    DSsv0_i = Ssv0_i - Ssvj
                    F0_i = calcular.CalcularFuerzas(
                        H, L, NY, NX, t, DSN0_i, DSsh0_i, DSsv0_i, self.AsH,
                        self.AsV)
                    # Criterio de convergencia
                    DU_Ref = np.max(
                        np.abs(np.delete(DU0_i, GLF)/np.delete(U0_i, GLF)))
                    if i > 1000 or DU_Ref < 0.005:
                        break
                    # Valores para la siguiente iteración
                    DFj = DF0_i
                    F00 = F0_i
                    Uj = U0_i
                    Rj = Rj + F0_i
                    eNj = eN0_i
                    eshj = esh0_i
                    esvj = esv0_i
                    NFj = NF0_i
                    SNj = SN0_i
                    SPj = SP0_i
                    Snj = Sn0_i
                    Sshj = Ssh0_i
                    Ssvj = Ssv0_i
                    i += 1
                if DU_Ref > 0.005:
                    QMessageBox.information(
                        self, "Cálculo terminado",
                        f'Límite de resistencia alcanzado a {NivCS[j-1]} '
                            f'veces la sobrecarga, limite de iteraciones',
                              QMessageBox.StandardButton.Ok)
                    print(
                        f'Límite de resistencia alcanzado a {NivCS[j-1]} '
                            f'veces la sobrecarga, limite de iteraciones')
                    break
                # Almacenando resultados
                self.U = np.vstack((self.U, [U0_i]))
                self.R = np.vstack((self.R, [Rj]))
                self.eN = np.vstack((self.eN, [eN0_i]))
                self.eP = np.vstack((self.eP, [eP0_i]))
                self.esh = np.vstack((self.esh, [esh0_i]))
                self.esv = np.vstack((self.esv, [esv0_i]))
                self.NF = np.vstack((self.NF, [NF0_i]))
                self.SN = np.vstack((self.SN, [SN0_i]))
                self.SP = np.vstack((self.SP, [SP0_i]))
                self.Sn = np.vstack((self.Sn, [Sn0_i]))
                self.Ssh = np.vstack((self.Ssh, [Ssh0_i]))
                self.Ssv = np.vstack((self.Ssv, [Ssv0_i]))
                self.Carga_Def = np.vstack((self.Carga_Def, np.array(
                    [-min(U0_i[1+np.arange(0, len(U0_i), 2)]), NivCS[j],
                     i-1, DU_Ref])))
                self.ListaCombo = np.append(
                    self.ListaCombo, f'PP + {NivCS[j]}*CV')
            except LinAlgError:
                QMessageBox.information(
                    self, "Cálculo terminado",
                    f'Límite de resistencia alcanzado a {NivCS[j-1]} '
                        f'veces la sobrecarga, matriz singular',
                          QMessageBox.StandardButton.Ok)
                print(
                    f'Límite de resistencia alcanzado a {NivCS[j-1]} '
                        f'veces la sobrecarga, matriz singular')
                # Límite de resistencia alcanzado
                break
            # Graficando
            x, y = self.Carga_Def[-1, 0], self.Carga_Def[-1, 1]
            ejes.set_xlim(0, 1.1*x)
            ejes.set_ylim(0, 1.1*y)
            # Añadir x al conjunto de datos
            line.set_xdata(np.append(line.get_xdata(), x))
            # Añadir y al conjunto de datos
            line.set_ydata(np.append(line.get_ydata(), y))
            plt.pause(0.01)
        plt.ioff()  # Desactivar el modo interactivo
        plt.show()
        plt.close()
        graficar.GraficarCurva(self.Carga_Def, ax, self.lienzo)
        QMessageBox.information(
            self, "Fin del análisis", f'Análisis terminado a {NivCS[j]} '
                f'veces la sobrecarga', QMessageBox.StandardButton.Ok)
        self.combo_etapas.addItems(self.ListaCombo)
        # Exportando resultados
        encabezados = ['Deflexión (mm)', 'Nivel de Carga',
                       'No. de iteraciones', 'Error']
        datos = self.Carga_Def
        marca_temporal = datetime.datetime.now().strftime('%d%m%Y%H%M%S')
        # Creando el archivo CSV en modo escritura
        with open(f'Resultados_{marca_temporal}.csv', 'w', newline='') as file:
            writer = csv.writer(file)
            writer.writerow(encabezados)
            for row in self.Carga_Def:
                writer.writerow(row)
        self.Enc_Col_Nudo, self.Enc_Fila_Nudo, self.Enc_Col_Elem, \
            self.Enc_Fila_Elem = calcular.EncabezadosTablas(H, L, NY, NX)
        return self.F, self.U, self.R, self.eN, self.esh, self.esv,\
                self.NF, self.SN, self.SP, self.Sn, self.Ssh, self.Ssv, self.AsH,\
                self.AsV, self.Carga_Def, self.Enc_Col_Nudo, \
                self.Enc_Fila_Nudo, self.Enc_Col_Elem, self.Enc_Fila_Elem

    def mostrar_graficos(self, index):
        graficar = Graficar()
        auxiliar = Auxiliar()
        H, L, NX, NY, t, fc, fy, M_cargas, M_apoyos, V_apoyos, \
            V_cuantias, ax = self.validar_datos()
        etapa = self.combo_etapas.currentIndex()+1
        if index == 0:
            graficar.GraficarViga(H, L, NX, NY, ax, self.lienzo)
            graficar.GraficarCargas(H, L, NX, NY, M_cargas, ax, self.lienzo)
            graficar.GraficarApoyos(
                H, L, NX, NY, M_apoyos, V_apoyos, ax, self.lienzo)
        elif index == 1:
            graficar.GraficarDeformada(
                H, L, NY, NX, self.U[etapa], self.Carga_Def[-1, 0], ax,
                self.lienzo)
        elif index == 2:
            Tabla_SX = auxiliar.TablaNudos(H, L, NY, NX,
                                           self.Sn[etapa, :, 0])
            graficar.GraficarTensionesX(
                H, L, NY, NX, fc, Tabla_SX, ax, self.lienzo)
        elif index == 3:
            Tabla_SY = auxiliar.TablaNudos(H, L, NY, NX,self.Sn[etapa, :, 1])
            graficar.GraficarTensionesY(
                H, L, NY, NX, fc, Tabla_SY, ax, self.lienzo)
        elif index == 4:
            Tabla_SXY = auxiliar.TablaNudos(H, L, NY, NX,
                                            self.Sn[etapa, :, 2])
            graficar.GraficarTensionesXY(
                H, L, NY, NX, fc, Tabla_SXY, ax, self.lienzo)
        elif index == 5:
            graficar.GraficarTensionesP(
                H, L, NY, NX, self.SP[etapa], ax, self.lienzo)
        elif index == 6:
            graficar.GraficarFisuras(
                H, L, NY, NX, fc, self.SP[etapa], self.NF[etapa],
                ax, self.lienzo)
        elif index == 7:
            graficar.GraficarCurva(self.Carga_Def, ax, self.lienzo)

    def mostrar_datos_tabla(self, index):
        self.tabla_resultados_num.verticalHeader().setVisible(True)
        auxiliar = Auxiliar()
        H, L, NX, NY, t, fc, fy, M_cargas, M_apoyos, V_apoyos, \
            V_cuantias, ax = self.validar_datos()
        etapa = self.combo_etapas.currentIndex()+1
        if index == 0:
            self.tabla_resultados_num.setRowCount(len(self.Enc_Fila_Nudo))
            self.tabla_resultados_num.setColumnCount(len(self.Enc_Col_Nudo))
            self.tabla_resultados_num.setVerticalHeaderLabels(
                self.Enc_Fila_Nudo)
            self.tabla_resultados_num.setHorizontalHeaderLabels(
                self.Enc_Col_Nudo)
            Tabla_DespV, Tabla_DespH = auxiliar.TablaGradosLibertad(
                H, L, NY, NX, self.U[etapa])
            for fila in range(len(self.Enc_Fila_Nudo)):
                for columna in range(len(self.Enc_Col_Nudo)):
                    item = QTableWidgetItem(Tabla_DespH[fila][columna])
                    self.tabla_resultados_num.setItem(fila, columna, item)
        elif index == 1:
            self.tabla_resultados_num.setRowCount(len(self.Enc_Fila_Nudo))
            self.tabla_resultados_num.setColumnCount(len(self.Enc_Col_Nudo))
            self.tabla_resultados_num.setVerticalHeaderLabels(
                self.Enc_Fila_Nudo)
            self.tabla_resultados_num.setHorizontalHeaderLabels(
                self.Enc_Col_Nudo)
            Tabla_DespV, Tabla_DespH = auxiliar.TablaGradosLibertad(
                H, L, NY, NX, self.U[etapa])
            for fila in range(len(self.Enc_Fila_Nudo)):
                for columna in range(len(self.Enc_Col_Nudo)):
                    item = QTableWidgetItem(Tabla_DespV[fila][columna])
                    self.tabla_resultados_num.setItem(fila, columna, item)
        elif index == 2:
            self.tabla_resultados_num.setRowCount(len(self.Enc_Fila_Nudo))
            self.tabla_resultados_num.setColumnCount(len(self.Enc_Col_Nudo))
            self.tabla_resultados_num.setVerticalHeaderLabels(
                self.Enc_Fila_Nudo)
            self.tabla_resultados_num.setHorizontalHeaderLabels(
                self.Enc_Col_Nudo)
            Tabla_FzasV, Tabla_FzasH = auxiliar.TablaGradosLibertad(
                H, L, NY, NX, self.R[etapa])
            for fila in range(len(self.Enc_Fila_Nudo)):
                for columna in range(len(self.Enc_Col_Nudo)):
                    item = QTableWidgetItem(Tabla_FzasH[fila][columna])
                    self.tabla_resultados_num.setItem(fila, columna, item)
        elif index == 3:
            self.tabla_resultados_num.setRowCount(len(self.Enc_Fila_Nudo))
            self.tabla_resultados_num.setColumnCount(len(self.Enc_Col_Nudo))
            self.tabla_resultados_num.setVerticalHeaderLabels(
                self.Enc_Fila_Nudo)
            self.tabla_resultados_num.setHorizontalHeaderLabels(
                self.Enc_Col_Nudo)
            Tabla_FzasV, Tabla_FzasH = auxiliar.TablaGradosLibertad(
                H, L, NY, NX, self.R[etapa])
            for fila in range(len(self.Enc_Fila_Nudo)):
                for columna in range(len(self.Enc_Col_Nudo)):
                    item = QTableWidgetItem(Tabla_FzasV[fila][columna])
                    self.tabla_resultados_num.setItem(fila, columna, item)
        elif index == 4:
            self.tabla_resultados_num.setRowCount(len(self.Enc_Fila_Nudo))
            self.tabla_resultados_num.setColumnCount(len(self.Enc_Col_Nudo))
            self.tabla_resultados_num.setVerticalHeaderLabels(
                self.Enc_Fila_Nudo)
            self.tabla_resultados_num.setHorizontalHeaderLabels(
                self.Enc_Col_Nudo)
            Tabla_SX = auxiliar.TablaNudos(H, L, NY, NX, self.Sn[etapa,:,0])
            for fila in range(len(self.Enc_Fila_Nudo)):
                for columna in range(len(self.Enc_Col_Nudo)):
                    item = QTableWidgetItem(Tabla_SX[fila][columna])
                    self.tabla_resultados_num.setItem(fila, columna, item)
        elif index == 5:
            self.tabla_resultados_num.setRowCount(len(self.Enc_Fila_Nudo))
            self.tabla_resultados_num.setColumnCount(len(self.Enc_Col_Nudo))
            self.tabla_resultados_num.setVerticalHeaderLabels(
                self.Enc_Fila_Nudo)
            self.tabla_resultados_num.setHorizontalHeaderLabels(
                self.Enc_Col_Nudo)
            Tabla_SY = auxiliar.TablaNudos(H, L, NY, NX, self.Sn[etapa,:,1])
            for fila in range(len(self.Enc_Fila_Nudo)):
                for columna in range(len(self.Enc_Col_Nudo)):
                    item = QTableWidgetItem(Tabla_SY[fila][columna])
                    self.tabla_resultados_num.setItem(fila, columna, item)
        elif index == 6:
            self.tabla_resultados_num.setRowCount(len(self.Enc_Fila_Nudo))
            self.tabla_resultados_num.setColumnCount(len(self.Enc_Col_Nudo))
            self.tabla_resultados_num.setVerticalHeaderLabels(
                self.Enc_Fila_Nudo)
            self.tabla_resultados_num.setHorizontalHeaderLabels(
                self.Enc_Col_Nudo)
            Tabla_SXY = auxiliar.TablaNudos(H, L, NY, NX, self.Sn[etapa,:,2])
            for fila in range(len(self.Enc_Fila_Nudo)):
                for columna in range(len(self.Enc_Col_Nudo)):
                    item = QTableWidgetItem(Tabla_SXY[fila][columna])
                    self.tabla_resultados_num.setItem(fila, columna, item)
        elif index == 7:
            self.tabla_resultados_num.setRowCount(len(self.Enc_Fila_Elem))
            self.tabla_resultados_num.setColumnCount(len(self.Enc_Col_Elem))
            self.tabla_resultados_num.setVerticalHeaderLabels(
                self.Enc_Fila_Elem)
            self.tabla_resultados_num.setHorizontalHeaderLabels(
                self.Enc_Col_Elem)
            Tabla_SMax = auxiliar.TablaElementos(
                H, L, NY, NX, self.SP[etapa, :, 0])
            for fila in range(len(self.Enc_Fila_Elem)):
                for columna in range(len(self.Enc_Col_Elem)):
                    item = QTableWidgetItem(Tabla_SMax[fila][columna])
                    self.tabla_resultados_num.setItem(fila, columna, item)
        elif index == 8:
            self.tabla_resultados_num.setRowCount(len(self.Enc_Fila_Elem))
            self.tabla_resultados_num.setColumnCount(len(self.Enc_Col_Elem))
            self.tabla_resultados_num.setVerticalHeaderLabels(
                self.Enc_Fila_Elem)
            self.tabla_resultados_num.setHorizontalHeaderLabels(
                self.Enc_Col_Elem)
            Tabla_SMin = auxiliar.TablaElementos(
                H, L, NY, NX, self.SP[etapa, :, 1])
            for fila in range(len(self.Enc_Fila_Elem)):
                for columna in range(len(self.Enc_Col_Elem)):
                    item = QTableWidgetItem(Tabla_SMin[fila][columna])
                    self.tabla_resultados_num.setItem(fila, columna, item)
        elif index == 9:
            self.tabla_resultados_num.setRowCount(len(self.Enc_Fila_Elem))
            self.tabla_resultados_num.setColumnCount(len(self.Enc_Col_Elem))
            self.tabla_resultados_num.setVerticalHeaderLabels(
                self.Enc_Fila_Elem)
            self.tabla_resultados_num.setHorizontalHeaderLabels(
                self.Enc_Col_Elem)
            Tabla_SAng = auxiliar.TablaElementos(
                H, L, NY, NX, self.SP[etapa, :, 2])
            for fila in range(len(self.Enc_Fila_Elem)):
                for columna in range(len(self.Enc_Col_Elem)):
                    item = QTableWidgetItem(Tabla_SAng[fila][columna])
                    self.tabla_resultados_num.setItem(fila, columna, item)
        elif index == 10:
            self.tabla_resultados_num.setRowCount(len(self.Enc_Fila_Elem))
            self.tabla_resultados_num.setColumnCount(len(self.Enc_Col_Elem))
            self.tabla_resultados_num.setVerticalHeaderLabels(
                self.Enc_Fila_Elem)
            self.tabla_resultados_num.setHorizontalHeaderLabels(
                self.Enc_Col_Elem)
            Tabla_DefMax = auxiliar.TablaElementos(
                H, L, NY, NX, 1000*self.eP[etapa, :, 0])
            for fila in range(len(self.Enc_Fila_Elem)):
                for columna in range(len(self.Enc_Col_Elem)):
                    item = QTableWidgetItem(Tabla_DefMax[fila][columna])
                    self.tabla_resultados_num.setItem(fila, columna, item)
        elif index == 11:
            self.tabla_resultados_num.setRowCount(len(self.Enc_Fila_Elem))
            self.tabla_resultados_num.setColumnCount(len(self.Enc_Col_Elem))
            self.tabla_resultados_num.setVerticalHeaderLabels(
                self.Enc_Fila_Elem)
            self.tabla_resultados_num.setHorizontalHeaderLabels(
                self.Enc_Col_Elem)
            Tabla_DefMin = auxiliar.TablaElementos(
                H, L, NY, NX, 1000*self.eP[etapa, :, 1])
            for fila in range(len(self.Enc_Fila_Elem)):
                for columna in range(len(self.Enc_Col_Elem)):
                    item = QTableWidgetItem(Tabla_DefMin[fila][columna])
                    self.tabla_resultados_num.setItem(fila, columna, item)
        elif index == 12:
            self.tabla_resultados_num.setRowCount(len(self.Enc_Fila_Nudo))
            self.tabla_resultados_num.setColumnCount(len(self.Enc_Col_Elem))
            self.tabla_resultados_num.setVerticalHeaderLabels(
                self.Enc_Fila_Nudo)
            self.tabla_resultados_num.setHorizontalHeaderLabels(
                self.Enc_Col_Elem)
            Tabla_AsH = auxiliar.TablaAsHorizontal(H, L, NY, NX, self.AsH)
            for fila in range(len(self.Enc_Fila_Nudo)):
                for columna in range(len(self.Enc_Col_Elem)):
                    item = QTableWidgetItem(Tabla_AsH[fila][columna])
                    self.tabla_resultados_num.setItem(fila, columna, item)
        elif index == 13:
            self.tabla_resultados_num.setRowCount(len(self.Enc_Fila_Elem))
            self.tabla_resultados_num.setColumnCount(len(self.Enc_Col_Nudo))
            self.tabla_resultados_num.setVerticalHeaderLabels(
                self.Enc_Fila_Elem)
            self.tabla_resultados_num.setHorizontalHeaderLabels(
                self.Enc_Col_Nudo)
            Tabla_AsV = auxiliar.TablaAsVertical(H, L, NY, NX, self.AsV)
            for fila in range(len(self.Enc_Fila_Elem)):
                for columna in range(len(self.Enc_Col_Nudo)):
                    item = QTableWidgetItem(Tabla_AsV[fila][columna])
                    self.tabla_resultados_num.setItem(fila, columna, item)
        elif index == 14:
            self.tabla_resultados_num.setRowCount(len(self.Enc_Fila_Nudo))
            self.tabla_resultados_num.setColumnCount(len(self.Enc_Col_Elem))
            self.tabla_resultados_num.setVerticalHeaderLabels(
                self.Enc_Fila_Nudo)
            self.tabla_resultados_num.setHorizontalHeaderLabels(
                self.Enc_Col_Elem)
            Tabla_DefSh = auxiliar.TablaAsHorizontal(
                H, L, NY, NX, 1000*self.esh[etapa])
            for fila in range(len(self.Enc_Fila_Nudo)):
                for columna in range(len(self.Enc_Col_Elem)):
                    item = QTableWidgetItem(Tabla_DefSh[fila][columna])
                    self.tabla_resultados_num.setItem(fila, columna, item)
        elif index == 15:
            self.tabla_resultados_num.setRowCount(len(self.Enc_Fila_Elem))
            self.tabla_resultados_num.setColumnCount(len(self.Enc_Col_Nudo))
            self.tabla_resultados_num.setVerticalHeaderLabels(
                self.Enc_Fila_Elem)
            self.tabla_resultados_num.setHorizontalHeaderLabels(
                self.Enc_Col_Nudo)
            Tabla_DefSv = auxiliar.TablaAsVertical(
                H, L, NY, NX, 1000*self.esv[etapa])
            for fila in range(len(self.Enc_Fila_Elem)):
                for columna in range(len(self.Enc_Col_Nudo)):
                    item = QTableWidgetItem(Tabla_DefSv[fila][columna])
                    self.tabla_resultados_num.setItem(fila, columna, item)
        elif index == 16:
            self.tabla_resultados_num.setRowCount(len(self.Enc_Fila_Nudo))
            self.tabla_resultados_num.setColumnCount(len(self.Enc_Col_Elem))
            self.tabla_resultados_num.setVerticalHeaderLabels(
                self.Enc_Fila_Nudo)
            self.tabla_resultados_num.setHorizontalHeaderLabels(
                self.Enc_Col_Elem)
            Tabla_Ssh = auxiliar.TablaAsHorizontal(
                H, L, NY, NX, self.Ssh[etapa])
            for fila in range(len(self.Enc_Fila_Nudo)):
                for columna in range(len(self.Enc_Col_Elem)):
                    item = QTableWidgetItem(Tabla_Ssh[fila][columna])
                    self.tabla_resultados_num.setItem(fila, columna, item)
        elif index == 17:
            self.tabla_resultados_num.setRowCount(len(self.Enc_Fila_Elem))
            self.tabla_resultados_num.setColumnCount(len(self.Enc_Col_Nudo))
            self.tabla_resultados_num.setVerticalHeaderLabels(
                self.Enc_Fila_Elem)
            self.tabla_resultados_num.setHorizontalHeaderLabels(
                self.Enc_Col_Nudo)
            Tabla_Ssv = auxiliar.TablaAsVertical(H, L, NY,NX,self.Ssv[etapa])
            for fila in range(len(self.Enc_Fila_Elem)):
                for columna in range(len(self.Enc_Col_Nudo)):
                    item = QTableWidgetItem(Tabla_Ssv[fila][columna])
                    self.tabla_resultados_num.setItem(fila, columna, item)
        elif index == 18:
            self.tabla_resultados_num.verticalHeader().setVisible(False)
            Carga_Def_str = np.round(self.Carga_Def, 4).astype(str)
            self.tabla_resultados_num.setRowCount(len(self.Carga_Def[:, 0]))
            self.tabla_resultados_num.setColumnCount(4)
            self.tabla_resultados_num.setHorizontalHeaderLabels(
                ["Def. máx. (mm)", "Niv. Carga", "No. Iteraciones", "Error"])
            for fila in range(len(self.Carga_Def[:, 0])):
                for columna in ([1, 0, 2, 3]):
                    item = QTableWidgetItem(Carga_Def_str[fila][columna])
                    self.tabla_resultados_num.setItem(fila, columna, item)

if __name__ == '__main__':
    app = QApplication(sys.argv)
    ventana = Ventana()
    sys.exit(app.exec())
