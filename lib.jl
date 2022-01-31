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
Calcula o erro para um determinado modelo
x: Pontos no dominio
y: Pontos da imagem
model: Modelo a ser avaliado
=#
function calculate_error(x, y, model)
    # Faz a soma do quadrado de cada diferenca e retorna a raiz
    error = sum((y .- model.(x)).^2)
    return sqrt(error)
end