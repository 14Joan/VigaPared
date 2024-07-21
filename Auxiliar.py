import numpy as np
class Auxiliar():
    def CargasSuperficales(self, H, L, NY, NX, M_cargas):
        # Código para el álculo de del vector de cargas superficiales [kN]
        DX = L/NX
        DY = H/NY
        n = int(NX*NY + NX + NY + 1)     # Numero de nudos
        M_cargas = np.array(M_cargas).astype(float)
        # Cargas superficiales
        Rs = np.zeros(2*n)  # Vector de cargas superficiales
        Dq = np.zeros(M_cargas.shape[0])
        No_nudos = np.zeros(M_cargas.shape[0], dtype=int)
        Coord_cargas = np.array([])
        for fila in range(M_cargas.shape[0]):
            Dq[fila] = (M_cargas[fila, 1]-M_cargas[fila, 0]) / \
                (M_cargas[fila, 3]-M_cargas[fila, 2])*DX
            No_nudos[fila] = np.round(
                (M_cargas[fila, 3]-M_cargas[fila, 2])/DX+1).astype(int)
            for i in range(No_nudos[fila]):
                if i == 0:
                    Fila_coord_cargas = np.array(
                        [fila, M_cargas[fila, 2]+i*DX, M_cargas[fila, 4],
                         (1/2*M_cargas[fila, 0]+1/6*Dq[fila])*DX])
                elif i == No_nudos[fila]-1:
                    Fila_coord_cargas = np.array(
                        [fila, M_cargas[fila, 2]+i*DX, M_cargas[fila, 4],
                         (1/2*M_cargas[fila, 0]+
                          1/6*(3*No_nudos[fila]-4)*Dq[fila])*DX])
                else:
                    Fila_coord_cargas = np.array(
                        [fila, M_cargas[fila, 2]+i*DX, M_cargas[fila, 4],
                         (M_cargas[fila, 0]+i*Dq[fila])*DX])
                if Coord_cargas.size == 0:
                    Coord_cargas = Fila_coord_cargas
                else:
                    Coord_cargas = np.vstack((Coord_cargas, Fila_coord_cargas))
        # Ordenando el vector de cargas
        for j in range(Coord_cargas.shape[0]):
            NCol = np.round(Coord_cargas[j, 1]/DX + 1).astype(int)
            NFila = np.round(Coord_cargas[j, 2]/DY + 1).astype(int)
            NNudo = int((NCol-1)*(NY+1)+NFila-1)
            Rs[2*NNudo+1] = Rs[2*NNudo+1] + Coord_cargas[j, 3]
        return Rs

    def CargasVolumetricas(self, H, L, NY, NX, t):
        # Codigo para el cálculo del vector de cargas volumétricas [kN]
        dens = 25   # KN/m3
        DX = L/NX
        DY = H/NY
        n = int(NX*NY + NX + NY + 1)     # Numero de nudos
        N = int(NX*NY)  # Número de elementos
        # Cargas volumétricas
        Rv = np.zeros(2*n)
        Rvi = -t*DX*DY*dens/4/100
        # Formando el vector de cargas volumétricas
        for i in range(N):
            NoCol = np.floor(i/NY+1).astype(int)
            Rv[int(2*(i+NoCol-1)+1)] = Rv[int(2*(i+NoCol-1)+1)] + Rvi
            Rv[int(2*(i+NoCol)+1)] = Rv[int(2*(i+NoCol)+1)] + Rvi
            Rv[int(2*(i+NoCol+NY)+1)] = Rv[int(2*(i+NoCol+NY)+1)] + Rvi
            Rv[int(2*(i+NoCol+NY+1)+1)] = Rv[int(2*(i+NoCol+NY+1)+1)] + Rvi
        return Rv

    def GradosLibertadFijos(Self, H, L, NY, NX, M_apoyos, V_apoyos):
        # Código que determina los grados de libertad fijos acorde a la
        # distribucion y tipo de apoyos
        DX = L/NX
        DY = H/NY
        M_apoyos = np.array(np.float64(M_apoyos))
        No_nudos_A = np.zeros(M_apoyos.shape[0], dtype=int)
        Coord_apoyos = np.array([])
        for fila in range(M_apoyos.shape[0]):
            No_nudos_A[fila] = np.round(
                (M_apoyos[fila, 1]-M_apoyos[fila, 0])/DX+1).astype(int)
            for i in range(No_nudos_A[fila]):
                if V_apoyos[fila] == 'Ambas':
                    T_apoyo = 0
                elif V_apoyos[fila] == 'Vertical':
                    T_apoyo = 1
                elif V_apoyos[fila] == 'Horizontal':
                    T_apoyo = 2
                Fila_coord_apoyos = np.array(
                    [T_apoyo, M_apoyos[fila, 0]+i*DX, M_apoyos[fila, 2]])
                if Coord_apoyos.size == 0:
                    Coord_apoyos = Fila_coord_apoyos
                else:
                    Coord_apoyos = np.vstack((Coord_apoyos,
                                              Fila_coord_apoyos))
        # Lista para grados de libertad fijos
        GL_Fijos = np.array([], dtype=int)
        for j in range(Coord_apoyos.shape[0]):
            NCol = np.round(Coord_apoyos[j, 1]/DX + 1).astype(int)
            NFila = np.round(Coord_apoyos[j, 2]/DY + 1).astype(int)
            NNudo = int((NCol-1)*(NY+1)+NFila-1)
            if Coord_apoyos[j, 0] == 0:
                GL_Fijos = np.append(GL_Fijos, [2*NNudo])
                GL_Fijos = np.append(GL_Fijos, [2*NNudo+1])
            elif Coord_apoyos[j, 0] == 1:
                GL_Fijos = np.append(GL_Fijos, [2*NNudo+1])
            elif Coord_apoyos[j, 0] == 2:
                GL_Fijos = np.append(GL_Fijos, [2*NNudo])
        return GL_Fijos

    def ModuloElasticoAcero(self, fy, es):
        # Código que determina el modulo de elasticidad del acero acorde
        # a la deformacion [MPa]
        Es1 = 200000
        Es2 = 200*fy/(2+0.001*fy)
        Es3 = 2000
        e1 = 0.8*fy/Es1
        e2 = e1 + 0.2*fy/Es2
        e3 = 0.002
        # eu = 0.12
        eu = 0.025
        if np.abs(es) <= e1:
            Es = Es1
        elif es > e1 and es <= e2:
            Es = Es2
        elif es > e2 and es <= eu:
            Es = Es3
        elif es >= -e3 and es < -e1:
            Es = Es2
        elif es >= -eu and es < -e3:
            Es = Es3
        else:
            Es = 0
        return Es

    def TensionAcero(self, fy, es):
        # Código que determina el nivel de tensión del acero acorde a la
        # deformacion [MPa]
        Es1 = 200000
        Es2 = 200*fy/(2+0.001*fy)
        Es3 = 2000
        e1 = 0.8*fy/Es1
        e2 = e1 + 0.2*fy/Es2
        e3 = 0.002
        # eu = 0.12
        eu = 0.025
        if np.abs(es) <= e1:
            fs = Es1*es
        elif es > e1 and es <= e2:
            fs = Es2*(es-e1) + 0.8*fy
        elif es > e2 and es <= eu:
            fs = Es3*(es-e2) + fy
        elif es >= -e3 and es < -e1:
            fs = Es2*(es+e1) - 0.8*fy
        elif es >= -eu and es < -e3:
            fs = Es3*(es+e3) - 0.8*fy - Es2*(e1-e3)
        else:
            fs = 0
        return fs

    def MatrizElasticaHormigon(self, fc, SP, NF):
        # Código que determina la matriz elástica del hormigón dado un
        # estado tensional y el nivel de fisuracion previo [MPa]
        # Ajuste de signos y unidades
        ang = SP[2]*np.pi/180  # Ángulo de tensión principal máxima Rad
        # Transformación a compresión positiva, tracción negativa
        SPT = np.array([-SP[1], -SP[0], SP[2]])
        # Tensiones octahedricas [MPa] y su direccion [º]
        So = (SPT[0] + SPT[1])/3
        To = (SPT[0]**2 + SPT[1]**2 + (SPT[0]-SPT[1])**2)**0.5/3
        # Orientación de la tensión desviadora
        if So == 0 and To == 0:
            Ango = 0
        else:
            Ango = np.arccos(So/(2**0.5*To))*180/np.pi
        SO = np.array([So, To, Ango])
        # Determinacion de suficiencia tensional
        if So/fc > -0.05:
            Toc = fc*0.944*(So/fc + 0.05)**0.724
            Toe = fc*0.633*(So/fc + 0.05)**0.857
            if Ango >= 0 and Ango < 60:
                AngR = Ango*np.pi/180
            elif Ango >= 60 and Ango < 120:
                AngR = (120-Ango)*np.pi/180
            elif Ango >= 120 and Ango < 180:
                AngR = (Ango-120)*np.pi/180
            elif Ango >= 180 and Ango < 240:
                AngR = (240-Ango)*np.pi/180
            elif Ango >= 240 and Ango < 300:
                AngR = (Ango-240)*np.pi/180
            else:
                AngR = (360-Ango)*np.pi/180
            # Interpolando la capacidad a corte del elemento
            Tou = (2*Toc*(Toc**2-Toe**2)*np.cos(AngR) +
                   Toc*(2*Toe-Toc)*(4*(Toc**2-Toe**2)*(np.cos(AngR))**2
                                    + 5*Toe**2 - 4*Toc*Toe)**0.5) \
                / (4*(Toc**2-Toe**2)*(np.cos(AngR))**2 + (2*Toe-Toc)**2)
        else:
            Toc = 0
            Toe = 0
            Tou = 0
        # Determinacion del estado de fisuracion
        if So/fc > -0.05:
            if To <= Tou:
                NF[0] = NF[0]
            # if To <= Tou:
            # else:
                # NF[0] = np.min([NF[0] + 1, 2])
            elif To > Tou and np.abs(NF[1]-90-SP[2]) < 45:
                NF[0] = 1
            elif To > Tou and np.abs(NF[1]-90-SP[2]) >= 45:
                NF[0] = np.min([NF[0] + 1, 2])
        else:
            NF[0] = 2
        NF[1] = SP[2]+90
        # Determinacion de la matriz de elasticidad
        # Constantes elásticass iniciales
        if fc <= 15:
            fc = 15
        elif fc >= 65:
            fc = 65
        Ke = 11000 + 3.2*fc**2
        Ge = 9224 + 136*fc + 3296*10**-15*fc**8.273
        Ee = 9*Ke*Ge/(3*Ke + Ge)
        ve = (3*Ke - 2*Ge)/(6*Ke + 2*Ge)
        # Constantes elásticass secantes y tangentes
        if fc <= 31.7:
            A = 0.516
            C = 3.573
            d = 2.12 + 0.0183*fc
        else:
            A = 0.516/(1 + 0.0027*(fc-31.7)**2.397)
            C = 3.573/(1 + 0.0134*(fc-31.7)**1.414)
            d = 2.7
        b = 2 + 1.81*10**-8*fc**4.461
        if So/fc <= -0.05:
            Ks = 0
            Kt = 0
        elif So/fc <= 0:
            Ks = Ke/(1 + A*(np.abs(So/fc))**(b-1))
            Kt = Ke/(1 + b*A*(np.abs(So/fc))**(b-1))
        elif So/fc <= 2:
            Ks = Ke/(1 + A*(So/fc)**(b-1))
            Kt = Ke/(1 + b*A*(So/fc)**(b-1))
        else:
            Ks = Ke/(1 + 2**(b-1)*A*b - 2**b*(b-1)*A*(So/fc)**-1)
            Kt = Ke/(1 + 2**(b-1)*A*b)
        Gs = Ge/(1 + C*(To/fc)**(d-1))
        Gt = Ge/(1 + d*C*(To/fc)**(d-1))
        # Coeficientes secantes
        Es = 9*Ks*Gs/(3*Ks + Gs)
        vs = (3*Ks - 2*Gs)/(6*Ks + 2*Gs)
        # coeficientes tangentes
        Et = 9*Kt*Gt/(3*Kt + Gt)
        vt = (3*Kt - 2*Gt)/(6*Kt + 2*Gt)
        # Factor de retencion de corte
        SRF = 0.3
        if NF[0] == 0:
            if So == 0:
                Cs = Ee/(1-ve**2)*np.array([[1, ve, 0],
                                            [ve, 1, 0],
                                            [0, 0, (1-ve)/2]])
                Ct = Cs
            else:
                Ct = Et/(1-vt**2)*np.array([[1, vt, 0],
                                            [vt, 1, 0],
                                            [0, 0, (1-vt)/2]])
                Cs = Es/(1-vs**2)*np.array([[1, vs, 0],
                                            [vs, 1, 0],
                                            [0, 0, (1-vs)/2]])
        elif NF[0] == 1:
            MRe = np.array([[(np.cos(ang))**2, (np.sin(ang))**2,
                             0.5*np.sin(2*ang)],
                            [(np.sin(ang))**2, (np.cos(ang))
                             ** 2, -0.5*np.sin(2*ang)],
                            [-np.sin(2*ang), np.sin(2*ang), np.cos(2*ang)]])
            Ct = np.array([[0, 0, 0],
                           [0, Et, 0],
                           [0, 0, SRF*Gt]])
            Ct = MRe.T @ Ct @ MRe
            Cs = np.array([[0, 0, 0],
                           [0, Es, 0],
                           [0, 0, SRF*Gs]])
            Cs = MRe.T @ Cs @ MRe
        else:
            Cs = np.zeros((3, 3))
            Ct = Cs
        return Ct, Cs, NF

    def TablaGradosLibertad(self, H, L, NY, NX, Vect_GL):
        # Organiza lo vectores de desplazamientos, fuerzas, etc para
        # presentarse en forma de tablas
        DX = L/NX
        DY = H/NY
        Enc_col_nudo = np.linspace(0, L, int(NX+1))
        Enc_fila_nudo = np.linspace(H, 0, int(NY+1))
        Tabla_GL_V = np.array([])
        Tabla_GL_H = np.array([])
        for fila in Enc_fila_nudo:
            Fila_GL_V = np.array([])
            Fila_GL_H = np.array([])
            NFila = np.round(fila/DY+1).astype(int)
            for col in Enc_col_nudo:
                NCol = np.round(col/DX+1).astype(int)
                NNudo = int((NCol-1)*(NY+1)+NFila-1)
                Fila_GL_V = np.append(Fila_GL_V, Vect_GL[2*NNudo+1])
                Fila_GL_H = np.append(Fila_GL_H, Vect_GL[2*NNudo])
            if Tabla_GL_V.size == 0 and Tabla_GL_H.size == 0:
                Tabla_GL_V = Fila_GL_V
                Tabla_GL_H = Fila_GL_H
            else:
                Tabla_GL_V = np.vstack((Tabla_GL_V, Fila_GL_V))
                Tabla_GL_H = np.vstack((Tabla_GL_H, Fila_GL_H))
        Tabla_GL_V = np.round(Tabla_GL_V, 4).astype(str)
        Tabla_GL_H = np.round(Tabla_GL_H, 4).astype(str)
        return Tabla_GL_V, Tabla_GL_H

    def TablaNudos(self, H, L, NY, NX, Vect_Nudos):
        # Organiza lo vectores de desplazamientos, fuerzas, etc para
        # presentarse en forma de tablas
        DX = L/NX
        DY = H/NY
        Enc_col_nudo = np.linspace(0, L, int(NX+1))
        Enc_fila_nudo = np.linspace(H, 0, int(NY+1))
        Tabla_Nudos = np.array([])
        for fila in Enc_fila_nudo:
            Fila_Nudos = np.array([])
            NFila = np.round(fila/DY+1).astype(int)
            for col in Enc_col_nudo:
                NCol = np.round(col/DX+1).astype(int)
                NNudo = int((NCol-1)*(NY+1)+NFila-1)
                Fila_Nudos = np.append(Fila_Nudos, Vect_Nudos[NNudo])
            if Tabla_Nudos.size == 0:
                Tabla_Nudos = Fila_Nudos
            else:
                Tabla_Nudos = np.vstack((Tabla_Nudos, Fila_Nudos))
        Tabla_Nudos = np.round(Tabla_Nudos, 4).astype(str)
        return Tabla_Nudos

    def TablaElementos(self, H, L, NY, NX, Vect_Elem):
        # Organiza lo vectores de desplazamientos, fuerzas, etc para 
        # presentarse en forma de tablas
        DY = H/NY
        DX = L/NX
        Enc_col_elem = np.linspace(DX/2, L-DX/2, int(NX))
        Enc_fila_elem = np.linspace(H-DY/2, DY/2, int(NY))
        Tabla_Elem = np.array([])
        for fila in Enc_fila_elem:
            Fila_Elem = np.array([])
            NFila = np.round((fila+DY/2)/DY).astype(int)
            for col in Enc_col_elem:
                NCol = np.round((col+DX/2)/DX).astype(int)
                NElem = int((NCol-1)*NY+NFila-1)
                Fila_Elem = np.append(Fila_Elem, Vect_Elem[NElem])
            if Tabla_Elem.size == 0:
                Tabla_Elem = Fila_Elem
            else:
                Tabla_Elem = np.vstack((Tabla_Elem, Fila_Elem))
        Tabla_Elem = np.round(Tabla_Elem, 4).astype(str)
        return Tabla_Elem

    def TablaAsHorizontal(self, H, L, NY, NX, Vect_RefH):
        # Organiza lo vectores del refuerzo horizontal para 
        # presentarse en forma de tablas
        DY = H/NY
        DX = L/NX
        Enc_fila_nudo = np.linspace(H, 0, int(NY+1))
        Enc_col_elem = np.linspace(DX/2, L-DX/2, int(NX))
        Tabla_RefH = np.array([])
        for fila in Enc_fila_nudo:
            Fila_RefH = np.array([])
            NFila = np.round(fila/DY+1).astype(int)
            for col in Enc_col_elem:
                NCol = np.round((col+DX/2)/DX).astype(int)
                NElem = int((NCol-1)*(NY+1)+NFila-1)
                Fila_RefH = np.append(Fila_RefH, Vect_RefH[NElem])
            if Tabla_RefH.size == 0:
                Tabla_RefH = Fila_RefH
            else:
                Tabla_RefH = np.vstack((Tabla_RefH, Fila_RefH))
        Tabla_RefH = np.round(Tabla_RefH, 4).astype(str)
        return Tabla_RefH

    def TablaAsVertical(self, H, L, NY, NX, Vect_RefV):
        # Organiza lo vectores del refuerzo vertical para presentarse en
        # forma de tablas
        DY = H/NY
        DX = L/NX
        Enc_col_nudo = np.linspace(0, L, int(NX+1))
        Enc_fila_elem = np.linspace(H-DY/2, DY/2, int(NY))
        Tabla_RefV = np.array([])
        for fila in Enc_fila_elem:
            Fila_RefV = np.array([])
            NFila = np.round((fila+DY/2)/DY).astype(int)
            for col in Enc_col_nudo:
                NCol = np.round(col/DX+1).astype(int)
                NElem = int((NCol-1)*NY+NFila-1)
                Fila_RefV = np.append(Fila_RefV, Vect_RefV[NElem])
            if Tabla_RefV.size == 0:
                Tabla_RefV = Fila_RefV
            else:
                Tabla_RefV = np.vstack((Tabla_RefV, Fila_RefV))
        Tabla_RefV = np.round(Tabla_RefV, 4).astype(str)
        return Tabla_RefV
