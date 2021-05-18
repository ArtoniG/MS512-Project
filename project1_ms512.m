
%% Projeto 1 - MS512
%
% Guilherme Artoni                  RA 160318
% Julia Martins Guardia             RA 176979
% Thiago Henrique Da Cruz           RA 206333
% Willian Rodrigues de Abreu Silva  RA 178739
% 
%%

# CRIA UMA FUNCAO QUE EXECUTA A DECOMPOSICAO QR
# EM UMA MATRIZ DE POSTO COMPLETO 
# RETORNANDO A MATRIZ R E VETOR RESPOSTA c
function [R, c] = QR_decomposition (A,b)
  
  # ARMAZENA O NUMERO DE LINHAS DA MATRIZ A
  rowsA = rows(A);
  
  # ARMAZENA O NUMERO DE COLUNAS DA MATRIZ A
  colsA = columns(A);
  
  # CALCULA AS MATRIZES Q
  for i = 1:colsA
 
     # CENARIO EXCLUSIVO DA PRIMEIRA ITERACAO
     if (i == 1)
    
      # ARMAZENA A PRIMEIRA COLUNA DA MATRIZ A EM UM OBJETO
      u = A(i:rowsA,i);
      
      # CALCULA SIGMA
      sigma = sign(u(1)) * sqrt(sum(u.**2))
      
      # SOMA SIGMA A PRIMEIRA COORDENADA DE u
      u(1) = u(1) + sigma;
  
      # CALCULA GAMMA
      gamma = 1/(sigma*u(1));
      
      # CALCULA MATRIZ Q ORTOGONAL
      Q = eye(rows(u)) - gamma.*u*u';
      
      # CALCULA MATRIZ R TRIANGULAR SUPERIOR
      R = Q*A;
      
      # CALCULA O NOVO VETOR RESPOSTA 
      c = Q*b;
    
    # CENARIO APOS PRIMEIRA ITERACAO
    else

      # ARMAZENA A PRIMEIRA COLUNA DA MATRIZ R EM UM OBJETO 
      u = R(i:rowsA,i);

      # CALCULA SIGMA
      sigma = sign(u(1)) * sqrt(sum(u.**2))
  
      # SOMA SIGMA A PRIMEIRA COORDENADA DE u
      u(1) = u(1) + sigma;
 
      # CALCULA GAMMA
      gamma = 1/(sigma*u(1));    
    
      # CRIA UMA MATRIZ IDENTIDADE COM MESMO NUMERO DE LINHAS DE A
      I = eye(rowsA);
      
      # CALCULA MATRIZ Q
      Q = eye(rows(u)) - gamma.*u*u';
      
      # INSERE A MATRIZ Q NA MATRIZ IDENTIDADE
      I(i:rowsA,i:rowsA) =  Q;
      
      # ATUALIZA A MATRIZ R
      R = I*R;
      
      # ATUALIZA O VETOR RESPOSTA
      c = I*c;
      
    endif 
    
  endfor
  
endfunction

%% Questão 1

# CRIA OBJETOS COM OS DADOS DA QUESTAO
t = [-1 -0.75 -0.5 0 0.25 0.5 0.75];
y = [1.00 0.8125 0.75 1.00 1.3125 1.75 2.3125];

# MOSTRA CURVA GERADA PELOS DADOS
plot(t,y);
xlabel("t");
ylabel("y");
title("Curva gerada pelos dados fornecidos no problema");


# MATRIZ DE DESENHO 
r = length(t);
A = [ones(r,1) t' t'.**2];

# EXECUTA A DECOMPOSICAO QR
[R, c] = QR_decomposition(A,y');

R_hat = R(1:columns(A),:);
c_hat = c(1:columns(A));

solution = inverse(R_hat)*c_hat;

%% Questão 2

# CRIA OBJETOS COM OS DADOS DA QUESTAO
t = [1000 1050 1060 1080 1110 1130];
y = [6010 6153 6421 6399 6726 6701];
z = [9422 9300 9220 9150 9042 8800];

# MATRIZ DE DESENHO 
r = length(t);
A = [1000*ones(r,1) t'-1065];

# EXECUTA A DECOMPOSICAO QR EM y
[Ry, cy] = QR_decomposition(A,y');

# SOLUCAO EM y
solution_y = inverse(Ry(1:columns(Ry),:))*cy(1:columns(Ry));

# CALCULA OS VALORES PREDITOS EM y
y_hat = 1000*solution_y(1) + solution_y(2).*A(:,2);

# MOSTRA CURVA GERADA PELOS DADOS
plot(t-1065,y,t-1065,y_hat);
xlabel('x (unidades arbitrárias)');
ylabel('y (unidades arbitrárias)');
title('Valores Observados x Valores Preditos através do método de Quadrados Mínimos');
legend('Observado', 'Predito');

# EXECUTA A DECOMPOSICAO QR EM z
[Rz, cz] = QR_decomposition(A,z');

# SOLUCAO EM z
solution_z = inverse(Rz(1:columns(Rz),:))*cz(1:columns(Rz));

# CALCULA OS VALORES PREDITOS EM z
z_hat = 1000*solution_z(1) + solution_z(2).*A(:,2);

# MOSTRA CURVA GERADA PELOS DADOS
plot(t-1065,z,t-1065,z_hat);
xlabel('x (unidades arbitrárias)');
ylabel('z (unidades arbitrárias)');
title('Valores Observados x Valores Preditos através do método de Quadrados Mínimos');
legend('Observado', 'Predito');

%% Questão 3

# CRIA OBJETOS COM OS DADOS DA QUESTAO
A = [1 2 3 4; 
     5 6 7 8; 
     9 10 11 12; 
     1 1 1 1; 
     3 2 1 0]
b = [10 26 42 4 6];

# EXECUTA A DECOMPOSICAO QR
[R, c] = QR_decomposition(A,b');

% A matriz A como foi definida possui posto 2 pois há combinações lineares 

A(:,2) + A(:,2) - A(:,1) 
A(:,3)

A(:,4) - A(:,2) + A(:,1)
A(:,3)
