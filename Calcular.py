import numpy as np
from Auxiliar import Auxiliar


class Calcular():
    def CalcularRigidez(self, H, L, NY, NX, t, fc, SP, NF, T_mod, fy, esh,
                        esv, AsH, AsV):
        # Determina la rigidez global KN/m
        DX = L/NX
        DY = H/NY
        auxiliar = Auxiliar()
        n = int(NX*NY + NX + NY + 1)     # Numero de nudos
        N = int(NX*NY)  # Número de elementos
        # Número de elementos de acero de refuerzo horizontal
        Nsh = int((NY+1)*NX)
        # Número de elementos de acero de refuerzo vertical
        Nsv = int(NY*(NX+1))
        # Matriz de rigidez global
        Rigidez = np.zeros((2*n, 2*n))
        # Aporte del hormigón

        def int_k(C, DX, DY, r, s):
            B = 1/2*np.array([[(1+s)/DX, 0, -(1+s)/DX, 0, -(1-s)/DX, 0,
                               (1-s)/DX, 0],
                              [0, (1+r)/DY, 0, (1-r)/DY, 0, -
                               (1-r)/DY, 0, -(1+r)/DY],
                              [(1+r)/DY, (1+s)/DX, (1-r)/DY, -(1+s)/DX,
                                -(1-r)/DY, -(1-s)/DX, -(1+r)/DY, (1-s)/DX]])
            k = B.T @ C @ B
            return k
        for i in range(N):
            Ct, Cs, N_NF = auxiliar.MatrizElasticaHormigon(
                fc, SP[i, :], NF[i, :])
            if T_mod == 'tan':
                C = 1000*Ct    # Por 1000 para convertir a KN/m2
            elif T_mod == 'sec':
                C = 1000*Cs    # Por 1000 para convertir a KN/m2
            # Rigidez local  hormigón
            K_loc = DX*DY*t/4/100*(int_k(C, DX, DY, 3**-0.5, 3**-0.5)+
                                   int_k(C, DX, DY, -3**-0.5,3**-0.5)+
                                   int_k(C, DX, DY, -3**-0.5, -3**-0.5)+
                                   int_k(C, DX, DY, 3**-0.5, -3**-0.5))
            NoCol = np.floor(i/NY+1).astype(int)
            p = np.array([2*(i+NoCol+NY+1), 2*(i+NoCol+NY+1)+1, 2*(i+NoCol),
                          2*(i+NoCol)+1,2*(i+NoCol-1), 2*(i+NoCol-1)+1,
                          2*(i+NoCol+NY), 2*(i+NoCol+NY)+1], dtype=int)
            Rigidez[np.ix_(p, p)] = Rigidez[np.ix_(p, p)] + K_loc
        # Aporte del acero
        AsH = np.array(AsH).astype(float)
        AsV = np.array(AsV).astype(float)
        for i in range(Nsh):
            # Rigidez local del Acero horizontal
            Eh = auxiliar.ModuloElasticoAcero(fy, esh[i])
            ksh = AsH[i]*Eh/DX*np.array([[1, -1], [-1, 1]])/10    # En KN/m
            p = np.array([2*(i+NY+1), 2*i], dtype=int)
            Rigidez[np.ix_(p, p)] = Rigidez[np.ix_(p, p)] + ksh
        for i in range(Nsv):
            # Rigidez local del Acero horizontal
            Ev = auxiliar.ModuloElasticoAcero(fy, esv[i])
            ksv = AsV[i]*Ev/DY*np.array([[1, -1], [-1, 1]])/10    # En KN/m
            NoCol = np.floor(i/NY+1).astype(int)
            p = np.array([2*(i+NoCol+1)-1, 2*(i+NoCol)-1], dtype=int)
            Rigidez[np.ix_(p, p)] = Rigidez[np.ix_(p, p)] + ksv
        return Rigidez

    def CalcularVectorCargas(self, H, L, NY, NX, t, M_cargas, NCV, NCS):
        # Calcula el vector de cargas en kN
        # Nivel de carga, a multiplicar por el vetor de cargas superficiales
        auxiliar = Auxiliar()
        # Cargas superficiales
        Rs = auxiliar.CargasSuperficales(H, L, NY, NX, M_cargas)
        # Cargas volumétricas
        Rv = auxiliar.CargasVolumetricas(H, L, NY, NX, t)
        # Cargas totales
        R = NCV*Rv + NCS*Rs
        return R

    def CalcularDesplazamientosFinal(self, H, L, NY, NX, M_apoyos, V_apoyos,
                                     R, M_Rigidez):
        # Calcula los desplazamientos en mm
        auxiliar = Auxiliar()
        # Apoyos
        GL_Fijos = auxiliar.GradosLibertadFijos(
            H, L, NY, NX, M_apoyos, V_apoyos)
        # Reduciendo matriz de rigidez
        Rigidez_red = np.delete(M_Rigidez, GL_Fijos, axis=0)
        Rigidez_red = np.delete(Rigidez_red, GL_Fijos, axis=1)
        R_red = np.delete(R, GL_Fijos)
        U = np.linalg.inv(Rigidez_red)@R_red
        # Completando el vector de desplazamientos
        for i in range(GL_Fijos.shape[0]):
            U = np.insert(U, GL_Fijos[i], 0)
        U = U*1000  # en mm
        return U

    def CalcularReacciones(self, H, L, NY, NX, M_apoyos, V_apoyos, U,
                           M_Rigidez):
        # Calcula las reacciones en kN
        n = int(NX*NY + NX + NY + 1)     # Numero de nudos
        auxiliar = Auxiliar()
        GL_Fijos = auxiliar.GradosLibertadFijos(
            H, L, NY, NX, M_apoyos, V_apoyos)
        Aux = np.zeros(2*n)
        Reac = M_Rigidez @ U/1000
        # Completando el vector de desplazamientos
        Aux[GL_Fijos] = Reac[GL_Fijos]
        Reac = Aux
        return Reac

    def CalcularDesplazamientos(self, H, L, NY, NX, t, fc, M_cargas, M_apoyos,
                                V_apoyos, AsH, AsV):
        # Código para el cálculo de los desplazamientos mm
        DX = L/NX
        DY = H/NY
        v = 0.2
        Ec = 1.45*fc*1000/3.5
        Es = 200000
        n = int(NX*NY + NX + NY + 1)     # Numero de nudos
        N = int(NX*NY)  # Número de elementos
        M_cargas = np.array(M_cargas).astype(float)
        M_apoyos = np.array(M_apoyos).astype(float)
        V_apoyos = np.array(V_apoyos)
        AsH = np.array(AsH).astype(float)
        AsV = np.array(AsV).astype(float)
        # Número de elementos de acero de refuerzo horizontal
        Nsh = int((NY+1)*NX)
        # Número de elementos de acero de refuerzo vertical
        Nsv = int(NY*(NX+1))
        # Cálculo del vector de cargas
        auxiliar = Auxiliar()
        # Cargas superficiales
        Rs = auxiliar.CargasSuperficales(H, L, NY, NX, M_cargas)
        # Cargas volumétricas
        Rv = auxiliar.CargasVolumetricas(H, L, NY, NX, t)
        # Cargas totales
        R = Rs + Rv
        # Apoyos
        GL_Fijos = auxiliar.GradosLibertadFijos(
            H, L, NY, NX, M_apoyos, V_apoyos)
        # Cálculo de las rigideces
        # Rigidez local  hormigón
        C = Ec/(1-v**2)*np.array([[1, v, 0], [v, 1, 0], [0, 0, (1-v)/2]])*1000

        def int_rigidez(DX, DY, r, s):
            B = 1/2*np.array([[(1+s)/DX, 0, -(1+s)/DX, 0, -(1-s)/DX, 0,
                               (1-s)/DX, 0],
                              [0, (1+r)/DY, 0, (1-r)/DY, 0, -
                               (1-r)/DY, 0, -(1+r)/DY],
                              [(1+r)/DY, (1+s)/DX, (1-r)/DY, -(1+s)/DX,
                               -(1-r)/DY, -(1-s)/DX, -(1+r)/DY, (1-s)/DX]])
            k = B.T @ C @ B
            return k
        rigidez_loc = DX*DY*t/4/100*(int_rigidez(DX, DY, 3**-0.5, 3**-0.5)+
                                     int_rigidez(DX, DY, -3**-0.5, 3**-0.5)+
                                     int_rigidez(DX, DY, -3**-0.5, -3**-0.5)+
                                     int_rigidez(DX, DY, 3**-0.5, -3**-0.5))
        Rigidez = np.zeros((2*n, 2*n))
        # Rigidez local del Acero
        ksh = Es/DX*np.array([[1, -1], [-1, 1]])*1000
        ksv = Es/DY*np.array([[1, -1], [-1, 1]])*1000
        # Ensamblando la matriz global
        # Aporte del acero
        for i in range(Nsh):
            p = np.array([2*(i+NY+1), 2*i], dtype=int)
            Rigidez[np.ix_(p, p)] = Rigidez[np.ix_(p, p)] + ksh/10000*AsH[i]
        for i in range(Nsv):
            NoCol = np.floor(i/NY+1).astype(int)
            p = np.array([2*(i+NoCol+1)-1, 2*(i+NoCol)-1], dtype=int)
            Rigidez[np.ix_(p, p)] = Rigidez[np.ix_(p, p)] + ksv/10000*AsV[i]
        # Aporte del hormigón
        for i in range(N):
            NoCol = np.floor(i/NY+1).astype(int)
            p = np.array([2*(i+NoCol+NY+1), 2*(i+NoCol+NY+1)+1, 2*(i+NoCol),
                          2*(i+NoCol)+1,2*(i+NoCol-1), 2*(i+NoCol-1)+1,
                          2*(i+NoCol+NY), 2*(i+NoCol+NY)+1], dtype=int)
            Rigidez[np.ix_(p, p)] = Rigidez[np.ix_(p, p)] + rigidez_loc
        # Calculando los desplazamientos
        Rigidez_red = np.delete(Rigidez, GL_Fijos, axis=0)
        Rigidez_red = np.delete(Rigidez_red, GL_Fijos, axis=1)
        R_red = np.delete(R, GL_Fijos)
        U = np.linalg.inv(Rigidez_red)@R_red
        # Completando el vector de desplazamientos
        for i in range(GL_Fijos.shape[0]):
            U = np.insert(U, GL_Fijos[i], 0)
        # Creando matriz para la tabla de resultados
        U = np.round(U*1000, 4)
        Fzas = np.round(Rigidez@U, 4)
        return U, Fzas

    def CalcularDeformaciones(self, H, L, NY, NX, U, AsH, AsV):
        # Calcula las deformaciones ex, ey, exy en el centro de cada elemento
        # de hormigón y ex o ey en el centro de cada refuerzo, adimensional
        DY = H/NY
        DX = L/NX
        N = int(NX*NY)
        U = U/1000  # U en m
        # Calculo de las deformaciones en el centro de cada elemento de hormigón
        B = 1/2*np.array([[1/DX, 0, -1/DX, 0, -1/DX, 0, 1/DX, 0],
                          [0, 1/DY, 0, 1/DY, 0, -1/DY, 0, -1/DY],
                          [1/DY, 1/DX, 1/DY, -1/DX, -1/DY, -1/DX, -1/DY,
                           1/DX]])
        eN = np.zeros((N, 3))  # 3 columnas para ex, ey, exy
        for i in range(N):
            NoCol = np.floor(i/NY+1).astype(int)
            p = np.array([2*(i+NoCol+NY+1), 2*(i+NoCol+NY+1)+1, 2*(i+NoCol),
                          2*(i+NoCol)+1, 2*(i+NoCol-1), 2*(i+NoCol-1)+1,
                          2*(i+NoCol+NY), 2*(i+NoCol+NY)+1], dtype=int)
            Ui = U[p]
            eN[i, :] = B@Ui     # Deformaciones ex, ey, exy, adimensionales
        # Cálculo de deformaciones en el acero de refuerzo
        # Número de elementos de acero de refuerzo horizontal
        Nsh = int((NY+1)*NX)
        # Número de elementos de acero de refuerzo vertical
        Nsv = int(NY*(NX+1))
        # Deformaciones refuerzo horizontal
        Bh = np.array([1, -1])/DX
        esh = np.zeros(Nsh)
        for i in range(Nsh):
            p = np.array([2*(i+NY+1), 2*i], dtype=int)
            if AsH[i] == 0:
                esh[i] = 0
            else:
                esh[i] = Bh@U[p]    # Adimensional
        # Deformaciones refuerzo vertical
        Bv = np.array([1, -1])/DY
        esv = np.zeros(Nsv)
        for i in range(Nsv):
            NoCol = np.floor(i/NY+1).astype(int)
            p = np.array([2*(i+NoCol+1)-1, 2*(i+NoCol)-1], dtype=int)
            if AsV[i] == 0:
                esv[i] = 0
            else:
                esv[i] = Bv@U[p]    # Adimensional
        return eN, esh, esv

    def CalcularDeformacionesPrincipales(self, NY, NX, eN):
        N = int(NX*NY)
        eP = np.zeros((N, 3))  # 3 columnas para emax, emin, angulo
        for i in range(N):
            eP[i, 0] = (eN[i, 0] + eN[i, 1])/2 + \
                (((eN[i, 0] - eN[i, 1])/2)**2 + (eN[i, 2]/2)**2)**0.5
            eP[i, 1] = (eN[i, 0] + eN[i, 1])/2 - \
                (((eN[i, 0] - eN[i, 1])/2)**2 + (eN[i, 2]/2)**2)**0.5
            eP[i, 2] = np.arctan2(eN[i, 2], eN[i, 0] - eN[i, 1])*180/np.pi/2
        return eP

    def CalcularFuerzas(self, H, L, NY, NX, t, SN, Ssh, Ssv, AsH, AsV):
        # Calcula las fuerzas internas  apartir del estado tensional en el
        # centro de cada elemento de hormigón y acero kN
        DY = H/NY
        DX = L/NX
        N = int(NX*NY)
        n = int(NX*NY + NX + NY + 1)     # Numero de nudos
        F = np.zeros(2*n)
        # Cálculo de las fuerzas en cada elemento de hormigón

        def int_fuerzas(DX, DY, r, s, Si):
            B = 1/2*np.array([[(1+s)/DX, 0, -(1+s)/DX, 0, -(1-s)/DX, 0,
                               (1-s)/DX, 0],
                              [0, (1+r)/DY, 0, (1-r)/DY, 0, -
                               (1-r)/DY, 0, -(1+r)/DY],
                              [(1+r)/DY, (1+s)/DX, (1-r)/DY, -(1+s)/DX,
                               -(1-r)/DY, -(1-s)/DX, -(1+r)/DY, (1-s)/DX]])
            Fi = B.T @ Si
            return Fi
        for i in range(N):
            NoCol = np.floor(i/NY+1).astype(int)
            F_i = DX*DY/4*t*10*(int_fuerzas(DX, DY, 3**-0.5, 3**-0.5, SN[i,:])+
                            int_fuerzas(DX, DY, 3**-0.5, -3**-0.5, SN[i,:])+
                            int_fuerzas(DX, DY,-3**-0.5,3**-0.5, SN[i,:])+
                            int_fuerzas(DX, DY,-3**-0.5,-3**-0.5, SN[i,:]))
            p = np.array([2*(i+NoCol+NY+1), 2*(i+NoCol+NY+1)+1, 2*(i+NoCol),
                          2*(i+NoCol)+1,2*(i+NoCol-1), 2*(i+NoCol-1)+1,
                          2*(i+NoCol+NY), 2*(i+NoCol+NY)+1], dtype=int)
            F[p] = F[p] + F_i
        # Cálculo de las fuerzas en cada elemento de acero
        # Número de elementos de acero de refuerzo horizontal
        Nsh = int((NY+1)*NX)
        # Número de elementos de acero de refuerzo vertical
        Nsv = int(NY*(NX+1))
        # Fuerzas refuerzo horizontal
        for i in range(Nsh):
            Fsh = np.array([1, -1])*Ssh[i]*AsH[i]/10  # En kN
            p = np.array([2*(i+NY+1), 2*i], dtype=int)
            F[p] = F[p] + Fsh
        # Fuerzas refuerzo vertical
        for i in range(Nsv):
            Fsv = np.array([1, -1])*Ssv[i]*AsV[i]/10  # En kN
            NoCol = np.floor(i/NY+1).astype(int)
            p = np.array([2*(i+NoCol+1)-1, 2*(i+NoCol)-1], dtype=int)
            F[p] = F[p] + Fsv
        return F

    def CalcularAreas(self, H, L, NY, NX, t, V_cuantias):
        # Calcula las áreas de cada elemento del refuerzo en cm2
        NY = int(NY)
        NX = int(NX)
        V_cuantias = np.array(V_cuantias).astype(float)
        # Vector refuerzo horizontal
        AsH = np.zeros((NY+1)*NX)
        # Vector refuerzo  verticales
        AsV = np.zeros(NY*(NX+1))
        # Determinando que filas y columnas de barras corresponden a
        # cada cuantía
        Fila_Ch1 = np.floor((NY+1)*0.15+1)
        Fila_Ch2 = np.floor((NY+1)*0.55+1) - Fila_Ch1
        Fila_Ch3 = NY+1 - Fila_Ch1 - Fila_Ch2
        Col_Cv1 = np.floor((NX+1)*0.15+1)
        Col_Cv2 = NX+1 - 2*Col_Cv1
        NoBarras = np.array([Fila_Ch1, Fila_Ch2, Fila_Ch3, Col_Cv1, Col_Cv2])
        Franjas = np.array([0.15*H, 0.4*H, 0.45*H, 0.15*L, 0.7*L])
        Areas_Franja = t*Franjas*V_cuantias/10
        Area_barra = Areas_Franja/NoBarras
        # Asignando áreas a cada elemento horizontal
        for fila in range(NY+1):
            if fila < NoBarras[0]:
                AsH[np.arange(0, NX)*(NY+1)+fila] = Area_barra[0]
            elif fila < NoBarras[0]+NoBarras[1]:
                AsH[np.arange(0, NX)*(NY+1)+fila] = Area_barra[1]
            else:
                AsH[np.arange(0, NX)*(NY+1)+fila] = Area_barra[2]
        # Asignando áreas a cada elemento vertical
        for columna in range(NX+1):
            if columna < NoBarras[3]:
                AsV[np.arange(0, NY) + NY*columna] = Area_barra[3]
            elif columna < NoBarras[3] + NoBarras[4]:
                AsV[np.arange(0, NY) + NY*columna] = Area_barra[4]
            else:
                AsV[np.arange(0, NY) + NY*columna] = Area_barra[3]
        return AsH, AsV

    def CalcularTensionesHormigon(self, H, L, NY, NX, fc, U, SP0, NF, T_mod):
        # Calcular las tensiones X, Y, XY y en el centro de cada elemento,
        # las tensiones promedio X, Y, XY en cada nudo del hormigón en MPa
        DY = H/NY
        DX = L/NX
        n = int(NX*NY + NX + NY + 1)     # Numero de nudos
        N = int(NX*NY)  # Número de elementos
        U = U/1000  # desplazamientos en m
        auxiliar = Auxiliar()
        # Calculo de esfuerzos en el hormigón
        # Coordendas de las 4 esquians y el centro
        Coord = np.array([[1, -1, -1, 1, 0],
                          [1, 1, -1, -1, 0]])
        # Matrices de almacenamiento de tensiones en los 5 puntos
        SX = np.zeros((N, 5))
        SY = np.zeros((N, 5))
        SXY = np.zeros((N, 5))
        # Matriz de almacenamiento de tensiones SX, SY, SXY en el centro de
        # cada elemento
        SN = np.zeros((N, 3))
        # Matriz de almacenamiento de tensiones X, Y, XY promedio en cada nudo
        Sn = np.zeros((n, 3))
        # Número de fisuras y su orientacion en cada elemento
        N_NF = np.zeros((N, 2))
        for i in range(N):
            # Matriz de elasticidad hormigón
            Ct, Cs, N_NF[i, :] = auxiliar.MatrizElasticaHormigon(
                fc, SP0[i, :], NF[i, :])
            if T_mod == 'tan':
                C = Ct
            elif T_mod == 'sec':
                C = Cs
            NoCol = np.floor(i/NY+1).astype(int)
            p = np.array([2*(i+NoCol+NY+1), 2*(i+NoCol+NY+1)+1, 2*(i+NoCol),
                          2*(i+NoCol)+1, 2*(i+NoCol-1), 2*(i+NoCol-1)+1,
                          2*(i+NoCol+NY), 2*(i+NoCol+NY)+1], dtype=int)
            Ui = U[p]
            for j in range(5):
                B = 1/2*np.array([[(1+Coord[1, j])/DX, 0, -(1+Coord[1, j])/DX,
                            0, -(1-Coord[1, j])/DX, 0, (1-Coord[1, j])/DX,0],
                            [0, (1+Coord[0, j])/DY, 0, (1-Coord[0, j])/DY,
                            0, -(1-Coord[0, j])/DY, 0, -(1+Coord[0, j])/DY],
                            [(1+Coord[0, j])/DY, (1+Coord[1, j])/DX,
                            (1-Coord[0, j])/DY, -(1+Coord[1, j])/DX,
                            -(1-Coord[0, j])/DY, -(1-Coord[1, j])/DX,
                            -(1+Coord[0, j])/DY, (1-Coord[1, j])/DX]])
                SNE = C @ B @ Ui
                SX[i, j] = SNE[0]
                SY[i, j] = SNE[1]
                SXY[i, j] = SNE[2]
        # Tensiones en el centro de cada elemento
        SN = np.vstack((SX[:, 4], SY[:, 4], SXY[:, 4]))
        SN = SN.T
        # Actualizando numero de fisuras
        NF = N_NF
        # Cálculo de tensiones promedio en cada nudo
        # Nudos de esquina
        Aux1 = np.array([[0, NY, (NY+1)*NX, (NY+1)*(NX+1)-1],
                        [0, NY-1, (NX-1)*NY, NX*NY-1],
                        [2, 1, 3, 0]]).astype(int)
        for m in range(4):
            Sn[Aux1[0, m], 0] = SX[Aux1[1, m], Aux1[2, m]]
            Sn[Aux1[0, m], 1] = SY[Aux1[1, m], Aux1[2, m]]
            Sn[Aux1[0, m], 2] = SXY[Aux1[1, m], Aux1[2, m]]
        # Nudos borde izquierdo
        n2 = np.arange(1, NY)
        NCol2 = 1
        Aux2 = np.array([n2, n2-NCol2, n2-NCol2+1]).astype(int)
        Sn[Aux2[0, :], 0] = (SX[Aux2[1, :], 1] + SX[Aux2[2, :], 2])/2
        Sn[Aux2[0, :], 1] = (SY[Aux2[1, :], 1] + SY[Aux2[2, :], 2])/2
        Sn[Aux2[0, :], 2] = (SXY[Aux2[1, :], 1] + SXY[Aux2[2, :], 2])/2
        # Nudos borde derecho
        n3 = NX*(NY+1) + np.arange(1, NY)
        NCol3 = NX + 1
        Aux3 = np.array([n3, n3-NCol3-NY, n3-NCol3-NY+1]).astype(int)
        Sn[Aux3[0, :], 0] = (SX[Aux3[1, :], 0] + SX[Aux3[2, :], 3])/2
        Sn[Aux3[0, :], 1] = (SY[Aux3[1, :], 0] + SY[Aux3[2, :], 3])/2
        Sn[Aux3[0, :], 2] = (SXY[Aux3[1, :], 0] + SXY[Aux3[2, :], 3])/2
        # Nudos borde superior
        n4 = np.arange(2, NX+1)*(NY+1)-1
        NCol4 = np.arange(2, NX+1)
        Aux4 = np.array([n4, n4-NCol4-NY, n4-NCol4]).astype(int)
        Sn[Aux4[0, :], 0] = (SX[Aux4[1, :], 0] + SX[Aux4[2, :], 1])/2
        Sn[Aux4[0, :], 1] = (SY[Aux4[1, :], 0] + SY[Aux4[2, :], 1])/2
        Sn[Aux4[0, :], 2] = (SXY[Aux4[1, :], 0] + SXY[Aux4[2, :], 1])/2
        # Nudos borde inferior
        n5 = np.arange(1, NX)*(NY+1)
        NCol5 = np.arange(2, NX+1)
        Aux5 = np.array([n5, n5-NCol5-NY+1, n5-NCol5+1]).astype(int)
        Sn[Aux5[0, :], 0] = (SX[Aux5[1, :], 3] + SX[Aux5[2, :], 2])/2
        Sn[Aux5[0, :], 1] = (SY[Aux5[1, :], 3] + SY[Aux5[2, :], 2])/2
        Sn[Aux5[0, :], 2] = (SXY[Aux5[1, :], 3] + SXY[Aux5[2, :], 2])/2
        # Nudos interiores
        for k in range(1, int(NY)):
            n6 = k + np.arange(1, NX)*(NY+1)
            NCol6 = np.arange(2, NX+1)
            Aux6 = np.array([n6, n6-NCol6-NY, n6-NCol6-NY+1,
                            n6-NCol6, n6-NCol6+1]).astype(int)
            Sn[Aux6[0, :], 0] = (
                SX[Aux6[1, :], 0] + SX[Aux6[2, :], 3] + SX[Aux6[3, :], 1]
                + SX[Aux6[4, :], 2])/4
            Sn[Aux6[0, :], 1] = (
                SY[Aux6[1, :], 0] + SY[Aux6[2, :], 3] + SY[Aux6[3, :], 1]
                + SY[Aux6[4, :], 2])/4
            Sn[Aux6[0, :], 2] = (
                SXY[Aux6[1, :], 0] + SXY[Aux6[2, :], 3] + SXY[Aux6[3, :], 1]
                + SXY[Aux6[4, :], 2])/4
        return SN, NF, Sn

    def CalcularTensionesPrincipales(self, NY, NX, SN):
        N = int(NX*NY)  # Número de elementos
        # Matriz de almacenamiento de tensiones principales Smax, Smin,
        # Angulo en el centro de cada elemento
        SP = np.zeros((N, 3))
        for i in range(N):
            SP[i, 0] = (SN[i, 0] + SN[i, 1])/2 + \
                (((SN[i, 0] - SN[i, 1])/2)**2 + (SN[i, 2])**2)**0.5
            SP[i, 1] = (SN[i, 0] + SN[i, 1])/2 - \
                (((SN[i, 0] - SN[i, 1])/2)**2 + (SN[i, 2])**2)**0.5
            SP[i, 2] = np.arctan2(
                SN[i, 2], (SN[i, 0] - SN[i, 1])/2)*180/np.pi/2
        return SP

    def CalcularTensionesAcero(self, NY, NX, fy, esh, esv):
        # Calcular las tensiones en cada elemento de acero en MPa
        auxiliar = Auxiliar()
        # Número de elementos de acero de refuerzo horizontal
        Nsh = int((NY+1)*NX)
        # Número de elementos de acero de refuerzo vertical
        Nsv = int(NY*(NX+1))
        # Tensiones refuerzo horizontal
        Ssh = np.zeros(Nsh)
        for i in range(Nsh):
            Ssh[i] = auxiliar.TensionAcero(fy, esh[i])
        # Esfuerzos refuerzo vertical
        Ssv = np.zeros(Nsv)
        for i in range(Nsv):
            Ssv[i] = auxiliar.TensionAcero(fy, esv[i])
        return Ssh, Ssv

    def EncabezadosTablas(self, H, L, NY, NX):
        DY = H/NY
        DX = L/NX
        # Encabezados para tablas de datos de nudos
        Enc_col_nudo = np.linspace(0, L, int(NX+1))
        Enc_fila_nudo = np.linspace(H, 0, int(NY+1))
        Enc_col_nudo = np.round(Enc_col_nudo, 3).astype(str)
        Enc_col_nudo = np.core.defchararray.add(Enc_col_nudo, ' m')
        Enc_fila_nudo = np.round(Enc_fila_nudo, 3).astype(str)
        Enc_fila_nudo = np.core.defchararray.add(Enc_fila_nudo, ' m')
        # Encabezados para tablas de datos de elementos
        Enc_col_elem = np.linspace(DX/2, L-DX/2, int(NX))
        Enc_fila_elem = np.linspace(H-DY/2, DY/2, int(NY))
        Enc_col_elem = np.round(Enc_col_elem, 3).astype(str)
        Enc_col_elem = np.core.defchararray.add(Enc_col_elem, ' m')
        Enc_fila_elem = np.round(Enc_fila_elem, 3).astype(str)
        Enc_fila_elem = np.core.defchararray.add(Enc_fila_elem, ' m')
        return Enc_col_nudo, Enc_fila_nudo, Enc_col_elem, Enc_fila_elem
