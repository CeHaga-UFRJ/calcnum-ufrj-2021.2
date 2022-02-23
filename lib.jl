#=
Aplica o metodo de Newton-Raphson
f: Funcao a ser aproximada
df: Derivada de f
x: Chute inicial
n: Numero de iteracoes

Retorno: A aproximacao apos n iteracoes do metodo
=#
function newton_raphson(f, df, x, n)
    # x0 = chute inicial
    xn = x
    for i = 1:n
        # Aplica o metodo n vezes
        xn = xn - f(xn)/df(xn)
    end
    return xn
end

#=
Aplica o metodo de Taylor para calcular ln centrado em a = 1 com erro definido
x: Parametro de ln
erro: Erro aceitavel

Retorno: A aproximacao apos n iteracoes do metodo
=#
function taylor_ln(x, erro)
    # Calcula o primeiro termo da serie
    res = x-1.0
    n = 2
    while(true)
        # Calcula o erro
        xn = (-1)^(n+1)*(x-1)^n/n
        
        # Se for aceitavel, para o programa
        if(abs(xn) < erro) break end
        
        # Senao, o erro ja eh o proximo termo
        res += xn
        n += 1
    end
    return res
end

#=
Aplica o metodo da bissecao para um funcao ate atingir um certo erro
f: Funcao a ser encontrada a raiz
a: Limite inferior do intervalo
b: Limite superior do intervalo
error: Erro maximo

Retorno: A aproximacao do metodo abaixo de um erro
=#
function bisection(f, a, b, error)
    # Se a >= b o intervalo esta ao contrario, erro de entrada
    if (a >= b) return "Intervalo invalido por limites", nothing end

    # Se nas "pontas" do intervalo tem uma raiz, retorna o proprio valor com 0 passos
    if (f(a) == 0) return a, 0 end
    if (f(b) == 0) return b, 0 end
    
    # Se possuem o mesmo sinal, erro de entrada
    if (f(a)/f(b) >= 0) return "Intervalo invalido por sinais", nothing end

    # Calcula o numero de passos, sabendo que o erro eh metade do intervalo e usando a formula descrita anteriormente
    steps = floor(log2((b-a)/(2*error)))+1
    for step = 1:steps
        # Encontra o meio
        m = (a+b)/2.0

        # Se o meio for raiz, retorna junto com o passo atual
        if (f(m) == 0) return m, step end

        # Troca o limite com mesmo sinal pro meio
        if (f(m)/f(b) >= 0)
            b = m
        else
            a = m
        end
    end

    return (a+b)/2.0, steps
end

#=
Realiza o metodo da bissecao ate um certo erro, depois realiza o metodo de Newton por determinados passos
f: Funcao a ser encontrada a raiz
df: Sua derivada
a: Limite inferior do intervalo
b: Limite superior do intervalo
range_size: Tamanho do intervalo
steps: Passos no metodo de Newton

Retorno: A aproximacao do metodo da bissecao abaixo de um erro seguido de n iteracoes do metodo de Newton
=#
function bisection_newton(f, df, a, b, range_size, steps)
    # Realiza o metodo da bisecao ate o intervalo determinado
    result, steps_bisection = bisection(f, a, b, range_size/2.0)
    
    # Se deu erro na bisecao retorna o erro
    if(isnothing(steps_bisection)) return result end
    
    # Continua a aproximacao com metodo de Newton
    result = newton_raphson(f, df, result, steps)
    
    return result
end

#=
Monta uma matriz de Vandermonde de dado grau
x: Vetor usado como base
rows: Quantidade de pontos
degree: Grau do polinomio

Retorno: Matriz de Vandermonde
=#
function vandermonde(x, rows, degree) 
    # Cria uma matriz vazia 
    V = zeros(rows, degree+1)
    
    # Para cada coluna adiciona uma potencia de u
    for i = 1:degree+1
         V[:, i] = x.^(i-1)
    end
    
    return V
end

#=
Realiza a interpolacao e retorna os coeficientes
x: Pontos do dominio
y: Pontos avaliados
degree: Grau do polinomio

Retorno: Coeficientes da interpolacao
=#
function interpolation(x, y, degree)
    rows = length(y)
    
    # Cria a matriz de Vandermonde
    V = vandermonde(x, rows, degree)
    
    # E resolve o sistema linear
    c = V \ y
    
    return c
end

#=
Avalia um polinomio em dados pontos a partir de sua lista de coeficientes
c: Coeficientes no formato (ntermos, 1)
x: Pontos a serem avaliados no formato (npontos, 1)

Retorno: Valor do polinomio em n pontos
=#
function evaluate(c, x, degree)
    columns = length(c)
    
    # Cria a matriz de Vandermonde
    V = vandermonde(x, columns)
    
    # E monta o sistema linear
    y = V * c
    
    return y
end

#=
Realiza a interpolacao por partes em 5 pontos usando 2 cubicas 
x: Pontos do dominio ordenados no formato (5, 1)
y: Pontos da imagem correspondentes no formato (5, 1)

Retorno: Coeficiente das cubicas
=#
function spline_interpolation(x, y)
    # Cria uma matriz vazia
    V = zeros(8,8)
    
    # O primeiro "retangulo" 3x4 eh vandermonde com x1, x2 e x3
    # O segundo "retangulo" 3x4 eh 0
    V[1:3, 1:4] = vandermonde(x[1:3], 4, 3)
    
    # O terceiro "retangulo" 3x4 eh 0
    # o quarto "retangulo" 3x4 eh vandermonde com x3, x4 e x5
    V[4:6, 5:8] = vandermonde(x[3:5], 4, 3)
    
    # A setima linha eh igualar as derivadas
    V[7,:] = [0 1 2x[3] 3x[3]^2 0 -1 -2x[3] -3x[3]^2]
    
    # A oitava linha eh igualar as segundas derivadas
    V[8,:] = [0 0 2 6x[3] 0 0 -2 -6x[3]]
    
    # Monta as saidas de cada equacao
    out_y = [y[1], y[2], y[3], y[3], y[4], y[5], 0, 0]
    
    # E resolve o sistema linear
    c = V \ out_y
    
    # Retorna os coeficientes da primeira e da segunda cubica
    return c[1:4], c[5:8]
end

#=
Realiza a interpolacao bilinear dado 4 pontos e seus valores respectivos
points: Uma lista com os valores [x1, x2, y1, y2]
z: Uma lista com as alturas correspondentes [z11, z12, z21, z22]

# Retorno: Coeficientes da funcao
=#
function bilinear_interpolation(points, z)
    x1, x2, y1, y2 = points
    # Para cada altura, multiplica pelo polinomio que da 1 no ponto e 0 nos outros
    # A funcao eh normalizada depois para que o resultado seja realmente 1
    f(x, y) = (z[1]*(x2 - x)*(y2 - y) + z[2]*(x2 - x)*(y - y1) + z[3]*(x - x1)*(y2 - y) + z[4]*(x - x1)*(y - y1))/((x2 - x1)*(y2 - y1))
    # Como o polinomio esta no formato a + bx + cy + dxy, podemos fazer algumas substituições para encontrar esses coeficientes
    a = f(0, 0) # = a + 0 + 0 + 0
    b = f(1, 0) - a # = a + b + 0 + 0
    c = f(0, 1) - a # = a + 0 + c + 0
    d = f(1, 1) - a - b - c # = a + b + c + d
    return [a, b, c, d]
end

#=
Calcula o Mean Squared Error para um determinado modelo
predict: Pontos previstos
real: Pontos reais

Retorno: Mean Squared Error
=#
function mse(predict, real)
    # Faz a soma do quadrado de cada diferenca e retorna a raiz
    error = sum((predict .- real).^2)
    return sqrt(error)
end

#=
Calcula a integral usando o metodo dos trapezios
f: Funcao a ser calculada a integral
a: Limite inferior de integracao
b: Limite superior de integracao
n: Quantidade de trapezios

Retorno: Aproximacao da integral com o metodo dos trapezios apos n iteracoes
=#
function integral_trapezoid(f, a, b, n = 1000)
    # Calcula distancia entre cada ponto
    h = (b-a)/n
    
    # Calcula a soma usando a formula extendida das areas
    S = 0
    xi = a
    # O meio sera somado duas vezes, as pontas uma
    for i = 0:n-1
        S += f(xi)
        xi += h
        S += f(xi)
    end
    S *= h/2
    
    return S
end

#=
Calcula a integral dupla usando o metodo dos trapezios
f: Funcao a ser calculada a integral
a: Limite inferior da integracao dy
b: Limite superior da integracao dy
h: Limite inferior da integracao dx
g: Limite superior da integracao dx
n: Quantidade de trapezios

Retorno: Aproximacao da integral dupla com o metodo dos trapezios apos n iteracoes
=#
function double_integral_trapezoid(f, a, b, h, g, n = 1000)
    # Calcula a integral para um y especifico
    function outer_integral(y)
        return integral_trapezoid(x -> f(x,y), h(y), g(y), n)
    end
    # Calcula a integral para cada y
    return integral_trapezoid(outer_integral, a, b)
end

#=
Resolve um sistema Ax = y onde A eh uma matriz diagonal
A: Matrix diagonal no formato (n,n)
y: Lado direito da equacao no formato (n,1)

Retorno: Solucao x
=#
function resolve_diagonal(A, y)
    # Podemos pegar o tamanho de y, pois eh o mesmo que A
    n = length(y)
    
    # Divide o lado direito pelo coeficiente de cada variavel
    x = [y[i]/A[i,i] for i = 1:n]
    
    return x
end

#=
Resolve um sistema Ax = y onde A eh uma matriz triangular superior
A: Matrix triangular superior no formato (n,n)
y: Lado direito da equacao no formato (n,1)

Retorno: Solucao x
=#
function resolve_triangular_superior(A, y)
    # Podemos pegar o tamanho de y, pois eh o mesmo que A
    n = length(y)   
    x = zeros(n)
    
    # Na triangular superior, comecamos de baixo para cima, substituindo as variaveis anteriores nas equacoes
    for i = n:(-1):1
        x[i] = (y[i] - sum([A[i, k] * x[k] for k = i+1:n])) / A[i,i] 
    end
    
    return x
end

#=
Resolve um sistema Ax = y onde A eh uma matriz triangular inferior
A: Matrix triangular inferior no formato (n,n)
y: Lado direito da equacao no formato (n,1)

Retorno: Solucao x
=#
function resolve_triangular_inferior(A, y)
    # Podemos pegar o tamanho de y, pois eh o mesmo que A
    n = length(y)
    x = zeros(n)
    
    # Na triangular inferior, comecamos de cima para baixo, substituindo as variaveis anteriores nas equacoes
    for i = 1:n
        x[i] = (y[i] - sum([A[i, k] * x[k] for k = 1:i-1])) / A[i,i] 
    end
    
    return x
end

#=
Realiza a decomposicao LU de uma matriz quadrada A
A: Matrix no formato (n,n)

Retorno: Matrizes L & U da decomposicao LU
=#
function decomposicao_lu(A)
    # Podemos salvar apenas um tamanho pois sao o mesmo
    n, = size(A)
    
    # U comeca como uma copia de A, enquanto L comeca como uma matriz identidade
    U = copy(A)
    L = Matrix(1.0I, n, n)
    
    for i = 1:n
        for j = i+1:n
            # Eh calculado o coeficiente dividindo o numero da matriz pelo pivot
            l = U[j,i] / U[i,i]
            # O coeficiente eh o elemento de L
            L[j,i] = l
            # E o coeficiente eh usado para alterar a linha de U
            U[j,:] += -l * U[i,:]
        end
    end
    return L, U
end

#=
Calcula a matriz inversa usando decomposicao LU
A: Matriz no formato (n,n) para encontrar a inversa

Retorno: Inversa de A
=#
function inverse_LU(A)
    # Decompoe em LU
    L, U = decomposicao_lu(A)
    n, = size(A)
    
    # Inicializa a inversa
    inv_A = zeros(n,n)
    
    # Para cada coluna
    for i = 1:n
        # Cria um vetor one-hot (Identidade final)
        y = zeros(n)
        y[i] = 1
        
        # Resolve o sistema para a coluna i
        Y = resolve_triangular_inferior(L, y)
        x = resolve_triangular_superior(U, Y)
        
        # Substitui a coluna i da inversa
        inv_A[:,i] = x
    end
    
    return inv_A
end

#=
Calcula a matriz usada para o metodo das diferencas finitas
n: Tamanho da matriz

Retorno: Matriz tridiagonal [1 -2 1]
=#
function finite_difference_matrix(n)
    A = zeros(n,n)
    
    # Comeca com -2 1
    A[1,1] = -2
    A[1,2] = 1
    
    # Termina com 1 -2
    A[n,n-1] = 1
    A[n,n] = -2
    
    # Forma uma tridiagonal com 1 -2 1
    for i = 2:n-1
        A[i,i-1] = 1
        A[i,i] = -2
        A[i,i+1] = 1
    end
    
    return A
end

#=
Realiza o metodo das diferencas finitas e retorna os pontos no meio
n: Numero de pontos
sd: Funcao da segunda derivada
yi & yf: Intervalo de y
xi & xf: Intervalo de x

Retorno: Pontos do meio apos o metodo de diferencas finitas
=#
function finite_difference(n, sd, yi, yf, xi, xf)
    # Calcula o tamanho do intervalo
    h = (xf - xi)/(n-1)
    
    A = finite_difference_matrix(n-2)
    b = zeros(n-2)
    
    # Cada equação vai possuir esse termo
    for i = 1:n-2
        b[i] = sd(xi + i*h) * 2 * h^2
    end
    
    # Com excecao do inicio e final, que possuem os limites ja conhecidos
    b[1] -= yi
    b[n-2] -= yf
    
    # Resolve e retorna
    y = A \ b
    return y
end

#=
Calcula a temperatura de cada ponto no problema do lago e resolve usando LU
n: Tamanho do lado do lago
tempN: Temperatura norte
tempW: Temperatura oeste
tempS: Temperatura sul
tempE: Temperatura leste

Retorno: Temperatura em cada ponto do lago
=#
function lake_degree(n, tempN, tempW, tempS, tempE)
    A = zeros(n^2, n^2)
    for i = 1:n^2
        A[i,i] = 4
        # Olha norte
        if(i > n) A[i, i-n] = -1 end
        
        # Olha sul
        if(i <= n*(n-1)) A[i, i+n] = -1 end
        
        # Olha oeste
        if(i % n != 1) A[i, i-1] = -1 end
        
        # Olha lest
        if(i % n != 0) A[i, i+1] = -1 end
    end
    
    y = zeros(n^2)
    
    # Parte de cima soma norte
    y[1 : n] += tempN * ones(n)
    
    # Parte de baixo soma sul
    y[n * (n-1) + 1 : n^2] += tempS * ones(n)
    
    # Parte da esquerda soma oeste
    for i = 1:n:n^2 y[i] += tempW end
    
    # Parte da direita soma leste
    for i = n:n:n^2 y[i] += tempE end
    
    # Resolve usando decomposicao LU
    L, U = decomposicao_lu(A)
    Y = resolve_triangular_inferior(L, y)
    x = resolve_triangular_superior(U, Y)
    
    return x
end