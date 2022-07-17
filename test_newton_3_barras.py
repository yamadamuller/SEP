#Algoritmo Newton-Raphson para cálculo do fluxo de potência
#Aluno: Mateus Yamada Muller

import numpy as np
import sympy as sym
#-------------------------------------------------------------
#definindo funções
def y_barra(n_barras,impedancias):
    Y_barra = np.zeros((n_barras, n_barras), dtype=complex)
    for linha in range(n_barras):
        for coluna in range(n_barras):
            var = 'Y{}{}'.format(linha + 1, coluna + 1)
            Y_barra[linha, coluna] += impedancias[var]  #itera a partir dos valores do dicionário
    return Y_barra

def ativa(k,barra):
    P = list()
    for N in range(k):
        if barra==(N+1): #sin(thetaii) = 0 / cos(thetaii) = 1
            P.append(sym.sympify('V{}*({}*V{})'.format(
                barra, G_barra[barra - 1, N], barra)))
        else:
            P.append(sym.sympify('V{}*(({}*cos(theta{}{})+{}*sin(theta{}{}))*V{})'.format(
                barra, G_barra[barra - 1, N], barra, N + 1, B_barra[barra - 1, N], barra, N + 1, N + 1)))

    P = sym.sympify(sum(P))
    return P

def reativa(k,barra):
    Q = list()
    for N2 in range(k):
        if barra==(N2+1): #sin(thetaii) = 0 / cos(thetaii) = 1
            Q.append(sym.sympify('V{}*(-{}*V{})'.format(
                barra, B_barra[barra - 1, N2], barra)))
        else:
            Q.append(sym.sympify('V{}*(({}*sin(theta{}{})-{}*cos(theta{}{}))*V{})'.format(
                barra, G_barra[barra - 1, N2], barra, N2 + 1, B_barra[barra - 1, N2], barra, N2 + 1, N2 + 1)))
    Q = sym.sympify(sum(Q))
    return Q

def jacobiano(vars,func):
    J = sym.zeros(len(func),len(vars)) #Matriz do jacobiano zerada para iterações
    for linha in range(len(func)):
        for coluna in range(len(vars)):
            if linha == 0:
                if coluna == 1:
                    J[linha, coluna] = sym.sympify('-16.25906104*V2*V3*cos(theta2) - 2.57435133*V2*V3*sin(theta2)')
                else:
                    J[linha,coluna] = sym.diff(func[linha],vars[coluna]) #realiza a derivada de uma função em relação a uma variável
            if linha == 1:
                if coluna == 0:
                    J[linha, coluna] = sym.sympify('-16.25906104*V2*V3*cos(theta3) - 2.57435133*V2*V3*sin(theta3)')
                else:
                    J[linha, coluna] = sym.diff(func[linha], vars[coluna])  # realiza a derivada de uma função em relação a uma variável
            if linha == 2:
                if coluna == 0:
                    J[linha, coluna] = sym.sympify('-16.25906104*V2*V3*sin(theta3) + 2.57435133*V2*V3*cos(theta3)')
                else:
                    J[linha, coluna] = sym.diff(func[linha], vars[coluna])  # realiza a derivada de uma função em relação a uma variável
            if linha == 3:
                if coluna == 1:
                    J[linha, coluna] = sym.sympify('-16.25906104*V2*V3*sin(theta2) + 2.57435133*V2*V3*cos(theta2)')
                else:
                    J[linha, coluna] = sym.diff(func[linha], vars[coluna])  # realiza a derivada de uma função em relação a uma variável
    return J
#-------------------------------------------------------------
#Valores de entrada

k = 3 #número de barras do sistema

admitancias = {
"Y11" : 1/complex(0,0.03),
"Y21": -1/(complex(0,0.03)),
"Y31": 0,
"Y12" : -1/(complex(0,0.03)),
"Y22" : (1/(complex(0.0095,0.06)) + complex(0,0.0201)) + (1/complex(0,0.03)),
"Y32": -1/(complex(0.0095,0.06)),
"Y13": 0,
"Y23": -1/(complex(0.0095,0.06)),
"Y33": (1/(complex(0.0095,0.06)) + complex(0,0.0201))
}
#-------------------------------------------------------------
#Matriz Admitância Y-Barra
Y_barra = y_barra(k,admitancias) #matriz admitância Y_barra
#-------------------------------------------------------------
#Fluxo de Potência
G_barra = Y_barra.real #parte real da matriz admitancia
B_barra = Y_barra.imag #parte complexa da matriz

P2 = ativa(k,2) #P2 com V e theta desconhecidos
Q2 = reativa(k,2) #Q2 com V e theta desconhecidos
P3 = ativa(k,3) #P3 com V e theta desconhecidos
Q3 = reativa(k,3) #Q2 com V e theta desconhecidos
#-------------------------------------------------------------
#Newton-Raphson
#Dados de exercício e chute inicial
P2_dados = 0
Q2_dados = 0
P3_dados = -0.92
Q3_dados = -0.25

V1 = 1.05 #módulo v1 dado de exercício
theta1 = np.radians(0) #ângulo theta1 dado de exercício em radianos
V2 = 1 #chute do módulo inicial em p.u
theta2 = np.radians(0) #chute inicial do ângulo em radianos
V3 = 1 #chute do módulo inicial em p.u
theta3 = np.radians(0) #chute inicial do ângulo em radianos


#Dados para controle das iterações
residuos = 0 #resíduos iniciais
iter_max = 15 #número de iterações permitidas
tol = 1e-4 #tolerância

for i in range(iter_max):
    print('iteração {}'.format(i+1))

    list_var = 'theta2 theta3 V3 V2'
    list_func = [P2.subs('theta21', 'theta2').subs('theta23', 'theta2').subs('theta32','theta3'),
                 P3.subs('theta21', 'theta2').subs('theta23', 'theta2').subs('theta32','theta3'),
                 Q3.subs('theta21', 'theta2').subs('theta23', 'theta2').subs('theta32','theta3'),
                 Q2.subs('theta21', 'theta2').subs('theta23', 'theta2').subs('theta32','theta3')]
    vars = sym.symbols(list_var)  # variáveis das funções P e Q
    func = sym.sympify(list_func)  # transforma string das funções em operadores

    residuos = np.zeros((k+1, 1))
    residuos[0, 0] = P2_dados - P2.subs('V1', V1).subs('V2', V2).subs('V3', V3).subs('theta21',theta2-theta1).subs('theta23',theta2-theta3)
    residuos[1, 0] = P3_dados - P3.subs('V1', V1).subs('V2', V2).subs('V3', V3).subs('theta32',theta3-theta2)
    residuos[2, 0] = Q3_dados - Q3.subs('V1', V1).subs('V2', V2).subs('V3', V3).subs('theta32',theta3-theta2)
    residuos[3, 0] = Q2_dados - Q2.subs('V1', V1).subs('V2', V2).subs('V3', V3).subs('theta21',theta2-theta1).subs('theta23',theta2-theta3)
    if np.any(np.abs(residuos)>tol):
        print('Resíduos maiores que a tolerância')

        J = jacobiano(vars, func)  # calcula o jacobiano chamando a função

        J[0,0] = (
            '33.333333333333336*V1*V2*cos(theta21) + V2*V3*(2.5743513312106225*sin(theta2) + 16.259061039224985*cos(theta2))')
        J[0,3] = (
            '33.333333333333336*V1*sin(theta21) + 5.148702662421245*V2 + V3*(16.259061039224985*sin(theta2) - 2.5743513312106225*cos(theta2))')
        J[3,0] = (
            '33.333333333333336*V1*V2*sin(theta21) + V2*V3*(16.259061039224985*sin(theta2) - 2.5743513312106225*cos(theta2))')
        J[3,3] = (
            '-33.333333333333336*V1*cos(theta21) + 99.14458874511664*V2 + V3*(-2.5743513312106225*sin(theta2) - 16.259061039224985*cos(theta2))')

        J = np.array(J.subs('V1', V1).subs('V2', V2).subs('V3', V3).subs('theta21', theta2 - theta1).subs('theta2',theta2 - theta3).subs('theta3', theta3 - theta2)) \
            .astype(np.float64)  # substitui o valor de chute e dados na matriz
        J_inv = np.linalg.inv(J)  # calcula a matriz inversa do jacobiano

        delta_x = np.matmul(J_inv, residuos)  # calculo dos próximos valores da iteração
        theta2 = theta2 + delta_x[0, 0]  # novo valor candidato para theta2
        theta3 = theta3 + delta_x[1, 0]  # novo valor candidato para theta3
        V3 = V3 + delta_x[2, 0]  # novo valor candidatos para V3
        V2 = V2 + delta_x[3, 0]  # novo valor candidatos para V2

        i = i + 1  # adiciona um ao contador
    else:
        print('Valores V2 = {}[p.u], Theta2 = {}[rad], V3 = {}[p.u] e Theta3 = {}[rad] são a solução'.format(V2,theta2,V3,theta3))
        break
#-------------------------------------------------------------
#Injeção de Potência na Barra 1
#Análise númerica com os resultados do Newton

P1 = ativa(k,1).subs('V1',V1).subs('V2',V2).subs('theta12',theta1-theta2)
Q1 = reativa(k,1).subs('V1',V1).subs('V2',V2).subs('theta12',theta1-theta2)

print('Sg = {} + j{} [p.u]'.format(P1,Q1))

