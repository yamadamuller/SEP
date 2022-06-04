#Algoritmo Newton-Raphson para cálculo do fluxo de potência
#Aluno: Mateus Yamada Muller

import numpy as np
import sympy as sym
from sympy import sin,cos
from sympy.abc import theta

#-------------------------------------------------------------
#definindo funções
def y_barra(n_barras,impedancias):
    Y_barra = np.zeros((n_barras, n_barras), dtype=complex)
    for linha in range(n_barras):
        for coluna in range(n_barras):
            var = 'Y{}{}'.format(linha + 1, coluna + 1)
            Y_barra[linha, coluna] += impedancias[var]  # itera para montar a matriz de impedancias Z_barra
    return Y_barra

def ativa(k):
    P = list()
    for N in range(k):
        P.append(sym.sympify('V2*(({}*cos(theta{}{})+{}*sin(theta{}{}))*V{})'.format(
            G_barra[k - 1, N], k, N + 1, B_barra[k - 1, N], k, N + 1, N + 1)))

    P_soma = sym.sympify(sum(P))
    P = P_soma.subs('theta22', np.radians(theta1)).subs('theta21', 'theta2')
    return P

def reativa(k):
    Q = list()
    for N2 in range(k):
        Q.append(sym.sympify('V2*(({}*sin(theta{}{})-{}*cos(theta{}{}))*V{})'.format(
            G_barra[k - 1, N2], k, N2 + 1, B_barra[k - 1, N2], k, N2 + 1, N2 + 1)))

    Q_soma = sym.sympify(sum(Q))
    Q = Q_soma.subs('theta22', np.radians(theta1)).subs('theta21', 'theta2')
    return Q

def jacobiano(vars,func): #definindo função que gera o Jacobiano
    J = sym.zeros(len(func),len(vars)) #Matriz do jacobiano zerada para iterações
    for linha in range(len(func)):
        for coluna in range(len(vars)):
            J[linha,coluna] = sym.diff(func[linha],vars[coluna]) #realiza a derivada de uma função em relação a uma variável
    return J
#-------------------------------------------------------------
#Valores de entrada

n_barras = 2 #número de barras do sistema
impedancias = {
"Y11" : 1/complex(0.01,0.05),
"Y21": -1/(complex(0.01,0.05)),
"Y12" : -1/(complex(0.01,0.05)),
"Y22" : 1/complex(0.01,0.05)
}
#dicionário contendo as impedâncias
# das barras e entre barras do sistema
#-------------------------------------------------------------
#Matriz Admitância Y-Barra
Y_barra = y_barra(n_barras,impedancias) #matriz admitância Y_barra
#-------------------------------------------------------------
#Fluxo de Potência
G_barra = Y_barra.real #parte real da matriz admitancia
B_barra = Y_barra.imag #parte complexa da matriz admitancia

k = 2 #barras PQ que deseja-se calcular a tensão por Newton-Raphson
theta1 = 0 #ângulo da tensão na barra de referência

#P = list()
#for N in range(k):
#    P.append(sym.sympify('V2*(({}*cos(theta{}{})+{}*sin(theta{}{}))*V{})'.format(
#        G_barra[k-1,N],k,N+1,B_barra[k-1,N],k,N+1,N+1)))

#P_soma = sym.sympify(sum(P))
#P = P_soma.subs('theta22',np.radians(theta1)).subs('theta21','theta2')
P = ativa(k)

#Q = list()
#for N2 in range(k):
#    Q.append(sym.sympify('V2*(({}*sin(theta{}{})-{}*cos(theta{}{}))*V{})'.format(
#        G_barra[k-1,N2],k,N2+1,B_barra[k-1,N2],k,N+1,N2+1)))

#Q_soma = sym.sympify(sum(Q))
#Q = Q_soma.subs('theta22',np.radians(theta1)).subs('theta21','theta2')
Q = reativa(k)
#-------------------------------------------------------------
#Newton-Raphson

#Dados de exercício e chute inicial
P2_dados = -1
Q2_dados = 0
V1 = 1.0112 #módulo v1 dado de exercício
theta1 = np.radians(0) #ângulo theta1 dado de exercício em radianos
V2 = 1 #chute do módulo inicial em p.u
theta2 = np.radians(0) #chute inicial do ângulo em radianos

#Dados para controle das iterações
residuos = 0 #residuos inicial
iter_max = 10 #iteração
tol = 1e-4 #tolerância

for i in range(iter_max):
    print('iteração {}'.format(i))

    list_var = 'theta2 V2'
    list_func = [P,Q]
    vars = sym.symbols(list_var)  # variáveis das funções P e Q
    func = sym.sympify(list_func)  # transforma string das funções em operadores

    residuos = np.zeros((2, 1))
    residuos[0, 0] = P2_dados - func[0].subs('V1', V1).subs('V2', V2).subs('theta2', theta2)
    residuos[1, 0] = Q2_dados - func[1].subs('V1', V1).subs('V2', V2).subs('theta2', theta2)
    if ((np.abs(residuos[0,0])>tol) or (np.abs(residuos[1,0])>tol)):
        print('Valores V2 = {} e Theta2 = {} maiores que a tolerância'.format(V2,np.degrees(theta2)))

        J = jacobiano(vars, func)  # calcula o jacobiano chamando a função definida anteriormente
        J = np.array(J.subs('V1', V1).subs('V2', V2).subs('theta2', theta2)).astype(
            np.float64)  # substitui o valor de chute e dados na matriz
        J_inv = np.linalg.inv(J)  # calcula a matriz inversa do jacobiano

        delta_x = np.matmul(J_inv, residuos)  # calculo dos próximos valores da iteração
        V2 = V2 + delta_x[1, 0]  # novo valor candidatos para V2
        theta2 = theta2 + delta_x[0, 0]  # novo valor candidato para theta2

        i = i + 1  # adiciona um ao contador
    else:
        print('Valores V2 = {} e Theta2 = {} são a solução'.format(V2,np.degrees(theta2)))
        break
