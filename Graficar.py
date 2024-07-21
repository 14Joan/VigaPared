import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import numpy as np
class Graficar():
    def GraficarViga(self, H, L, NX, NY, ax, Lienzo):
        DX = L/NX
        DY = H/NY
        N = int(NX*NY)
        # Dibujando
        ax.set_aspect('equal')
        for i in range(N):
            NC = np.floor(i/NY+1).astype(int)
            NF = np.floor(i%NY+1).astype(int)
            poligono_x = [DX*(NC-1), DX*NC, DX*NC, DX*(NC-1), DX*(NC-1)] 
            # Coordenadas x de los vértices del polígono
            poligono_y = [DY*(NF-1), DY*(NF-1), DY*NF, DY*NF, DY*(NF-1)] 
            # Coordenadas y de los vértices del polígono
            ax.plot(poligono_x, poligono_y, color='cyan', lw=2)
        Lienzo.draw()
    def GraficarCargas(self, H, L, NX, NY, M_cargas, ax, Lienzo):
        DX = L/NX
        DY = H/NY
        D = (DX*DY)**0.5
        cargas = np.array(np.float64(M_cargas))
        #Referencias para graficar
        flecha_x = D*np.array([0, -1/4, 0, 0, 0, 1/4, 0])
        flecha_y = -D*np.array([0, 1/4, 0, 1, 0, 1/4, 0])
        qmax = np.max(np.abs(cargas[:,:2]))
        #Graficando
        ax.set_aspect('equal')
        for fila in range(cargas.shape[0]):
            NoFlechas = np.round((cargas[fila,3]
                                  -cargas[fila,2])/(DX/2)+1).astype(int)
            for i in range(NoFlechas):
                X = cargas[fila,2] + (DX/2)*i
                q = ((cargas[fila,1]-cargas[fila,0])/
                     (cargas[fila,3]-cargas[fila,2])*
                     (X-cargas[fila,2]) + cargas[fila,0])/qmax
                ax.plot(flecha_x + X , q*flecha_y + cargas[fila,4],
                        color='red')
        Lienzo.draw()
    def GraficarApoyos(self, H, L, NX, NY, M_apoyos, V_apoyos, ax, Lienzo):
        DX = L/NX
        DY = H/NY
        D = (DX*DY)**0.5
        apoyos = np.array(np.float64(M_apoyos))
        tipo_apoyo = np.array(V_apoyos)
        #Referencias para graficar
        triangulo_x = D/2*np.array([0, -0.25*3**0.5, 0.25*3**0.5, 0])
        triangulo_y = D/2*np.array([0, -0.5*3**0.5, -0.5*3**0.5, 0])
        #Graficando
        ax.set_aspect('equal')
        for fila in range(apoyos.shape[0]):
            NoApoyos = np.round((apoyos[fila,1]-apoyos[fila,0])
                                /DX+1).astype(int)
            for i in range(NoApoyos):
                X = apoyos[fila,0] + DX*i
                if tipo_apoyo[fila] == 'Vertical':
                    ax.add_patch(Circle((X, -D/2*0.25*3**0.5 + apoyos[fila,2]),
                                        -D/2*0.25*3**0.5, fill=False,
                                        color='green', lw=2))
                elif tipo_apoyo[fila] == 'Horizontal':
                    ax.add_patch(Circle((-D/2*0.25*3**0.5 + X, apoyos[fila,2]),
                                    -D/2*0.25*3**0.5, fill=False,
                                    color='green', lw=2))
                elif tipo_apoyo[fila] == 'Ambas':
                    ax.plot(triangulo_x + X , triangulo_y + apoyos[fila,2],
                            color='green')
        Lienzo.draw()
    def GraficarDeformada(self, H, L, NY, NX, U, DefMax, ax, Lienzo):
        DX = L/NX
        DY = H/NY
        N = int(NX*NY)
        #Graficando
        ax.set_aspect('equal')
        Dref = 0.95*(DX*DY)**0.5/np.abs(DefMax)
        # Estructura sin deformar
        for i in range(N):
            NC = np.floor(i/NY+1).astype(int)
            NF = np.floor(i%NY+1).astype(int)
            poligono_x = [DX*(NC-1), DX*NC, DX*NC, DX*(NC-1), DX*(NC-1)]
            poligono_y = [DY*(NF-1), DY*(NF-1), DY*NF, DY*NF, DY*(NF-1)]
            ax.plot(poligono_x, poligono_y, color='cyan', linewidth=2)
        # Estructura deformada
        for i in range(N):
            NC = np.floor(i/NY+1).astype(int)
            NF = np.floor(i%NY+1).astype(int)
            cuad_x = np.array([DX*(NC-1), DX*(NC), DX*(NC), DX*(NC-1),
                               DX*(NC-1)])
            cuad_y = np.array([DY*(NF-1), DY*(NF-1), DY*(NF), DY*(NF),
                               DY*(NF-1)] )
            coord_x = Dref*np.array([U[int(2*(i+NC-1))], U[int(2*(i+NC+NY))],
                                     U[int(2*(i+NC+NY+1))], U[int(2*(i+NC))],
                                     U[int(2*(i+NC-1))]])
            coord_y = Dref*np.array([U[int(2*(i+NC-1)+1)],
                                     U[int(2*(i+NC+NY)+1)],
                                     U[int(2*(i+NC+NY+1)+1)],
                                     U[int(2*(i+NC)+1)],
                                     U[int(2*(i+NC-1)+1)]])
            ax.plot(cuad_x + coord_x, cuad_y + coord_y, color='blue')
        Lienzo.draw()
    def GraficarTensionesX(self, H, L, NY, NX, fc, Tabla_SX, ax, Lienzo):
        ax.set_aspect('equal')
        Coord_X = np.linspace(0, L, int(NX+1))
        Coord_Y = np.linspace(0, H, int(NY+1))
        x, y = np.meshgrid(Coord_X, Coord_Y)
        z = np.array(Tabla_SX).astype(float)
        z = z[::-1]
        smin = np.min(z)
        smax = np.max(z)
        grafico = ax.contourf(x,y,z, levels=100, vmin = np.max([-4*fc, smin]),
                              vmax = np.min([4*fc, smax]), cmap='rainbow')
        cbar = plt.colorbar(grafico)
        cbar.set_label('Tensiones Normales X (MPa)', rotation=270,
                       labelpad=30)
        # Muestra el gráfico
        Lienzo.draw()
    def GraficarTensionesY(self, H, L, NY, NX, fc, Tabla_SY, ax, Lienzo):
        ax.set_aspect('equal')
        Coord_X = np.linspace(0, L, int(NX+1))
        Coord_Y = np.linspace(0, H, int(NY+1))
        x, y = np.meshgrid(Coord_X, Coord_Y)
        z = np.array(Tabla_SY).astype(float)
        z = z[::-1]
        smin = np.min(z)
        smax = np.max(z)
        grafico = ax.contourf(x,y,z, levels=100, vmin = np.max([-4*fc, smin]),
                              vmax = np.min([4*fc, smax]), cmap='rainbow')
        cbar = plt.colorbar(grafico)
        cbar.set_label('Tensiones Normales Y (MPa)', rotation=270,
                       labelpad=30)
        Lienzo.draw()
    def GraficarTensionesXY(self, H, L, NY, NX, fc, Tabla_SXY, ax, Lienzo):
        ax.set_aspect('equal')
        Coord_X = np.linspace(0, L, int(NX+1))
        Coord_Y = np.linspace(0, H, int(NY+1))
        x, y = np.meshgrid(Coord_X, Coord_Y)
        z = np.array(Tabla_SXY).astype(float)
        z = z[::-1]
        smin = np.min(z)
        smax = np.max(z)
        grafico = ax.contourf(x,y,z, levels=100, vmin = np.max([-2*fc, smin]),
                              vmax = np.min([2*fc, smax]), cmap='rainbow')
        cbar = plt.colorbar(grafico)
        cbar.set_label('Tensiones Cortantes XY (MPa)', rotation=270,
                       labelpad=30)
        Lienzo.draw()
    def GraficarTensionesP(self, H, L, NY, NX, SP, ax, Lienzo):
        DY = H/NY
        DX = L/NX
        D = 0.8*(DX*DY)**0.5
        N = int(NX*NY)
        SP = np.array(SP).astype(float)
        SMax = np.max(np.abs(SP[:,:2]))
        Flecha_T = np.array([[ 0, -1, -0.7, -1, -0.7, -1, 1, 0.7, 1, 0.7,
                              1, 0],
                            [0, 0, -0.2, 0, 0.2, 0, -0, 0.2, -0, -0.2,
                             -0, 0,]])
        Flecha_C = np.array([[0, -0.3, 0, -1, 0, -0.3, 0, 0.3, 0, 1, 0,
                              0.3, 0],
                            [0, -0.2, 0, 0, 0, 0.2, 0, 0.2, 0, -0, 0,
                             -0.2, 0]])
        #Graficando
        ax.set_aspect('equal')
        def MRot(Ang):
            MRot = np.array([[np.cos(Ang*np.pi/180), -np.sin(Ang*np.pi/180)],
                             [np.sin(Ang*np.pi/180), np.cos(Ang*np.pi/180)]])
            return MRot
        for i in range(N):
            NC = np.floor(i/NY+1).astype(int)
            NF = np.floor(i%NY+1).astype(int)
            # Graficando estructura sin deformar
            poligono_x = [DX*(NC-1), DX*NC, DX*NC, DX*(NC-1), DX*(NC-1)]
            poligono_y = [DY*(NF-1), DY*(NF-1), DY*NF, DY*NF, DY*(NF-1)]
            ax.plot(poligono_x, poligono_y, color='cyan', linewidth=2)
            # Graficando flechas
            Centro_X = DX/2*(2*NC-1)
            Centro_Y = DY/2*(2*NF-1)
            CoefSMax = D*np.abs(SP[i,0])/SMax
            CoefSMin = D*np.abs(SP[i,1])/SMax
            FlechaTR_max = MRot(SP[i,2])@Flecha_T
            FlechaCR_max = MRot(SP[i,2])@Flecha_C
            FlechaTR_min = MRot(SP[i,2]+90)@Flecha_T
            FlechaCR_min = MRot(SP[i,2]+90)@Flecha_C
            if SP[i,0] > 0:
                ax.plot(Centro_X + CoefSMax*FlechaTR_max[0],
                        Centro_Y + CoefSMax*FlechaTR_max[1], color='red')
            else:
                ax.plot(Centro_X + CoefSMax*FlechaCR_max[0],
                        Centro_Y + CoefSMax*FlechaCR_max[1], color='blue')
            if SP[i,1] > 0:
                ax.plot(Centro_X + CoefSMin*FlechaTR_min[0],
                        Centro_Y + CoefSMin*FlechaTR_min[1], color='red')
            else:
                ax.plot(Centro_X + CoefSMin*FlechaCR_min[0],
                        Centro_Y + CoefSMin*FlechaCR_min[1], color='blue')
        Lienzo.draw()
    def GraficarFisuras(self, H, L, NY, NX, fc, SP, NF, ax, Lienzo):
        DY = H/NY
        DX = L/NX
        D = 0.5*(DX*DY)**0.5
        N = int(NX*NY)   #Número de elementos
        SP = np.array(SP).astype(float)
        NF = np.array(NF).astype(float)
        Fisura = np.array([[ -1, 1], [0, 0]])
        #Graficando
        ax.set_aspect('equal')
        def MRot(Ang):
            MRot = np.array([[np.cos(Ang*np.pi/180), -np.sin(Ang*np.pi/180)],
                             [np.sin(Ang*np.pi/180), np.cos(Ang*np.pi/180)]])
            return MRot
        for i in range(N):
            NCol = np.floor(i/NY+1).astype(int)
            NFila = np.floor(i%NY+1).astype(int)
            # Graficando estructura sin deformar
            poligono_x = [DX*(NCol-1), DX*NCol, DX*NCol, DX*(NCol-1),
                          DX*(NCol-1)]
            poligono_y = [DY*(NFila-1), DY*(NFila-1), DY*NFila, DY*NFila,
                          DY*(NFila-1)]
            ax.plot(poligono_x, poligono_y, color='cyan', linewidth=2)
            # Graficando fisuras
            Centro_X = DX/2*(2*NCol-1)
            Centro_Y = DY/2*(2*NFila-1)
            Fisura_R1 = MRot(NF[i,1])@Fisura
            Fisura_R2 = MRot(NF[i,1]+90)@Fisura
            if NF[i,0] == 1:
                ax.plot(Centro_X + D*Fisura_R1[0], Centro_Y + D*Fisura_R1[1],
                        color='black')
            elif NF[i,0] >= 2:
                ax.plot(Centro_X + D*Fisura_R1[0], Centro_Y + D*Fisura_R1[1],
                        color='black')
                ax.plot(Centro_X + D*Fisura_R2[0], Centro_Y + D*Fisura_R2[1],
                        color='black')
        Lienzo.draw()
    def GraficarCurva(self, Carga_Def, ax, Lienzo):
        Carga_Def = np.array(np.float64(Carga_Def))
        ax.set_xlabel('Deformación máxima [mm]')
        ax.set_ylabel('Nivel de carga aplicada')
        # Dibujando
        for i in range(Carga_Def.shape[0]):
            ax.plot(Carga_Def[:,0], Carga_Def[:,1], 'o--', color='purple',
                    lw=1, ms=5) 
        Lienzo.draw()