# PROYECTO 2

# PROBLEMA 1
# Implement the Jacobi symmetric eigenvalue algorithm. In particular, write a computer function that takes as input a symmetric n x n matrix A and outputs a diagonal matrix D
# with the eigenvalues the diagonal. Keep track of the J matrices so that you also output a matrix of eigenvectors.

# PROGRAMAS QUE OCUPAMOS PARA PODER SACAR LOS AUTOVALORES

function fuedia(A::Array) # Tiene de estrada la matriz y regresa el "tamaño" de lo que está afuera de la diagonal
 n=size(A,1)
 h=0
 for i=1:n-1
  for j=i+1:n
   h=h+2*A[i,j]^2
  end
 end
 return sqrt(h)
end

function indigüento(A::Array) # Tiene como entrar la matriz y busca el elemento más grande fuera de la diagonal, y regresa los índices
 n=size(A,1)
 k,l,amax=0,0,0
 for i=1:(n-1)                # Checa los valores que está arriba de la diagonal porque es simétrica
  for j=(i+1):n
   if abs(A[i,j])>=amax
    amax=abs(A[i,j])
    k=i
    l=j
   end
  end
 end
 return [k,l]
end

function valorotagüento(A::Array,coo)   #Tiene de entrada la matriz y regresa los valores de nuestra matriz rotacional que son c y s.
 if A[coo[1],coo[2]]!=0
  tau=(A[coo[2],coo[2]]-A[coo[1],coo[1]])/(2*A[coo[1],coo[2]])
  if tau>=0
   t=1/(tau+sqrt(1+tau^2))
  else
   t=-1/(-tau+sqrt(1+tau^2))
  end
  c=1/sqrt(1+t^2)
  s=t*c
 else
  c=1
  s=0
 end
 return [c,s]
end

function matrizrotgüento(n,c::Array,d) # (n,el vector de las matrices rotacionales, los valores de c y s)
 M=eye(n)                               # Regresa la matriz de n x n rotacional definida en el reporte
 M[c[1],c[1]]=d[1]
 M[c[1],c[2]]=d[2]
 M[c[2],c[1]]=-d[2]
 M[c[2],c[2]]=d[1]
 return M
end


# MÉTODO 1-PARA SACAR LOS AUTOVALORES Y AUTOVECTORES

function eignocla(A::Array)           # Tiene de entrada la matriz y regresa los autovalores y autovectores
 n=size(A,1)
 B=eye(n)
 tol=10.0^(-10)
 while fuedia(A)>tol                  # Checa si el "tamaño" del tamaño fuera de la diagonal es menor a nuestra tolerancia
  for i=1:((n^2-n)/2)
  coo=indigüento(A)                   # Obtiene los índices del valor mas grande
  num=valorotagüento(A,coo)
  mat=matrizrotgüento(n,coo,num)
  B=B*mat                             # Para obtener la matriz de autovectores
  A=transpose(mat)*A*mat
  end
 end
 return A,B                           # A siendo la matriz de autovectores de la diagonal y B la matriz de autovectores
end

# MÉTODO 2-PARA SACAR LOS AUTOVALORES Y AUTOVECTORES

function eigclasi(A::Array)        # Tiene de entrada la matriz y regresa los autovalores y autovectores
 n=size(A,2)
 B=eye(n)
 tol=10.0^(-10)
 while fuedia(A)>tol               # Checa si el "tamaño" del tamaño fuera de la diagonal es menor a nuestra tolerancia
  for i=1:n-1                      # A este tiempo de ciclo que checa lo que está afuera de la diagonal se llama "sweep"
   for j=i+1:n
    num=valorotagüento(A,[i,j])
    mat=matrizrotgüento(n,[i,j],num)
    B=B*mat                        # Para obtener la matriz de autovectores
    A=transpose(mat)*A*mat
   end
  end
 end
 return A,B                        # A siendo la matriz de autovectores de la diagonal y B la matriz de autovectores
end

# PROBLEMA 2
# Analyze the number of sweeps needed until convergences as a function of n. That is, for n=10,20,50,100,150,..., report the average number of sweeps needed.
# Produce a plot that shows the number of sweeps versus n.

function eigclascon2(A::Array)    # Lo mismo el programa de arriba, solo que ahora cuenta los "sweeps"
 n=size(A,2)
 B=eye(n)
 contador=0
 tol=10.0^(-10)
 while fuedia(A)>tol
  for i=1:(n^2/2-n/2)
   coo=indigüento(A)                   # Obtiene los índices del valor mas grande
   num=valorotagüento(A,coo)
   mat=matrizrotgüento(n,coo,num)
   B=B*mat                             # Para obtener la matriz de autovectores
   A=transpose(mat)*A*mat
  end
  contador+=1
  println("$contador-$(fuedia(A))") # Imprime para ver como va convergiendo
 end
 return A,B,contador
end

function aleagüento(n)    # Crea matrices de n x n aleatorias simétricas
 M=eye(n)
 for i=1:n
  for j=1:n
   m=rand()*10
   M[j,i]=m
   M[i,j]=m
  end
 end
 return M
end

function promegüento(n,p)   # ( Tamaño de la matriz, número de repeticiones)
 h=0
 for i=1:p                  # Para obtener el promedio de "sweeps" para llegar a la convergencia
  Alea=aleagüento(n)
  sweep=eigclascon2(Alea)
  h=h+sweep[3]
 end
 return h/p
end

promegüento(20,10)

# También puedo hacerlo con el método de -Cyclic-by-row-

function eigclasi2(A::Array)        # Tiene de entrada la matriz y regresa los autovalores y autovectores
 n=size(A,2)
 B=eye(n)
 contador=0
 tol=10.0^(-10)
 while fuedia(A)>tol               # Checa si el "tamaño" del tamaño fuera de la diagonal es menor a nuestra tolerancia
  for i=1:n-1                      # A este tiempo de ciclo que checa lo que está afuera de la diagonal se llama "sweep"
   for j=i+1:n
    num=valorotagüento(A,[i,j])
    mat=matrizrotgüento(n,[i,j],num)
    B=B*mat                        # Para obtener la matriz de autovectores
    A=transpose(mat)*A*mat
   end
  end
  contador+=1
  println("$contador-$(fuedia(A))")
 end
 return A,B,contador                        # A siendo la matriz de autovectores de la diagonal y B la matriz de autovectores
end

function promegüento2(n,p)   # ( Tamaño de la matriz, número de repeticiones)
 h=0
 for i=1:p                  # Para obtener el promedio de "sweeps" para llegar a la convergencia
  Alea=aleagüento(n)
  sweep=eigclasi(Alea)
  h=h+sweep[3]
 end
 return h/p
end

promegüento2(100,2)


# Desde 10 hasta 150 me dio que el promedio de "sweeps" son 4

using PyPlot

plot([10,20,50, 100, 150],[4, 4, 4, 4, 5],"r.",markersize=20,label="Sweeps-Classic");plot([10:2:150],map(x->log(x),[10:2:150]);label="\$log(n)\$");legend(loc="upper left",fancybox="true")

plot([10,20,50,100,150],[6,7,8,8,9],"g.",markersize=20,label="Sweeps-Cyclic-by-row-");legend(loc="upper left",fancybox="true")


plot([10,20,50, 100, 150,150],[4, 4, 4, 4, 4,5],"r.",markersize=20,label="Sweeps-Classic");plot([10:2:150],map(x->log(x),[10:2:150]);label="\$log(n)\$");legend(loc="upper left",fancybox="true")

plot([10,10,20,20,20,50,50,100,150],[6,5,5,6,7,8,7,8,9],"g.",markersize=20,label="Sweeps-Cyclic-by-row-");legend(loc="upper left",fancybox="true")

# PROBLEMA 4
# The Jacobi method for the symmetric eigenvalue problem can be adapted to a Jacobi algorithm to comput the SVD. Instead of solving a sequence 2x2 symmetric eigenvalue problems,
# solve a sequence of 2x2 SVD problems. Thus, for a given index pair (p,q), compute a pair of rotations such that...

# Para matrices n x n

function fuediasvd(A::Array) # Tiene de estrada la matriz y regresa el "tamaño" de lo que está afuera de la diagonal
 n=size(A,2)
 h=0
 for i=1:n
  for j=1:n
   if j!=i
    h=h+A[i,j]^2
   end
  end
 end
 return sqrt(h)
end

function indigüentosvd(A::Array) # Tiene como entrar la matriz y busca el elemento más grande fuera de la diagonal, y regresa los índices
 n=size(A,2)
 k,l,amax=0,0,0
 for i=1:n
  for j=1:n
   if j!=i
    if abs(A[i,j])>=amax
     amax=abs(A[i,j])
     k=i
     l=j
    end
   end
  end
 end
 return [k,l]
end

function indimatrizsvd(A::Array,a::Array) # Regresa un matriz de 2 x 2 para poder aplicar el logaritmo de SVD
 M=eye(2)
 if a[1]<a[2]
  M[1,1]=A[a[1],a[1]]
  M[1,2]=A[a[1],a[2]]
  M[2,1]=A[a[2],a[1]]
  M[2,2]=A[a[2],a[2]]
 elseif a[1]>a[2]
  M[1,1]=A[a[2],a[2]]
  M[1,2]=A[a[2],a[1]]
  M[2,1]=A[a[1],a[2]]
  M[2,2]=A[a[1],a[1]]
 end
 return M
end

function hacersimsvd(A::Array)            # Obtiene la matriz para que B=MN par N dado M hace simétrica
 M=eye(2)
 c1=(A[1,1]+A[2,2])/(A[1,1]*A[1,2]+A[2,1]*A[2,2])
 s1=(-A[1,2]+A[2,1])/(A[1,1]*A[1,2]+A[2,1]*A[2,2])
 m=1/sqrt(s1^2+c1^2)
 c=c1*m
 s=s1*m
 M[1,1]=c
 M[1,2]=s
 M[2,1]=-s
 M[2,2]=c
 return M*A,M
end

function svd1(A::Array)                     # Toma una matriz y aplica el SVD
 n=size(A,1)
 U=eye(n)
 V=eye(n)
 tol=10.0^(-6)
 while fuedia(A)>tol                         # Checa si el "tamaño" del tamaño fuera de la diagonal es menor a nuestra tolerancia
  for i=1:((n^2-n)/2)
   coo=indigüentosvd(A)                      # Obtiene los índices del valor mas grande
   num=indimatrizsvd(A,coo)                  # Hace la matriz de 2x2
   CS=hacersimsvd(num)                       # Saca s c y d1 d2
   C2S2=eigclasi(CS[1])                      # Obtiene c2 y s2
   C1S1=CS[2]*transpose(C2S2[2])             # Obtiene c1 y s1
   mat1=matrizrotgüento(n,coo,C1S1[:,1])     # Obtiene la matriz de c1,s2
   mat2=matrizrotgüento(n,coo,C2S2[2][:,1])  # Obtiene la matriz de c2,s2
   U=U*transpose(mat1)                       # Para U y V
   V=V*mat2
   A=mat1*A*mat2                             # Para sacar SIGMA
  end
 end
 return (A,V,U)
end

function svd2(A::Array)                     # Toma una matriz y aplica el SVD
 n=size(A,1)
 U=eye(n)
 V=eye(n)
 for j=1:10
  for i=1:((n^2-n)/2)
   coo=indigüentosvd(A)                      # Obtiene los índices del valor mas grande
   num=indimatrizsvd(A,coo)                  # Hace la matriz de 2x2
   CS=hacersimsvd(num)                       # Saca s c y d1 d2
   C2S2=eigclasi(transpose(CS[1]))                      # Obtiene c2 y s2
   C1S1=transpose(C2S2[2]*transpose(CS[2]))           # Obtiene c1 y s1
   mat1=matrizrotgüento(n,coo,C1S1[:,1])     # Obtiene la matriz de c1,s2
   mat2=matrizrotgüento(n,coo,C2S2[2][:,1])  # Obtiene la matriz de c2,s2
   U=U*mat1                       # Para U y V
   V=V*mat2
   A=mat1*A*mat2                             # Para sacar SIGMA
  end
 end
 return (A,V,U)
end

AA=[4 11 5; 14 8 6;7 -2 5]

Mat=svd2(AA)

Mat[3]*Mat[1]*transpose(Mat[2])             # Nos regresará la matriz original

# PROBLEMA 8
# Resolver soluciones de sistemas lineales
(
3x3 Array{Float64,2}:
 21.1747        0.0           1.67433e-8
 -8.88178e-16   8.95013      -1.8875e-13
  2.38034e-7   -3.79677e-12  -2.74383   ,

3x3 Array{Float64,2}:
 0.725021  -0.547344  -0.418042
 0.546998   0.826439  -0.133388
 0.418495  -0.131959   0.898581,

3x3 Array{Float64,2}:
 0.519939   0.697381  -0.493278
 0.804605  -0.205926   0.55696
 0.286834  -0.68648   -0.668185)
function jacobilin(A,b,x0,n)                 # Tiene de entrada la matriz, el vector de soluciones, un vector cualquiera y el número de iteraciones
 D=eye(length(b),length(b))
 R=zeros(length(b))                          # Regresa el vector de soluciones
 for i=1:length(b)
  D[i,i]=A[i,i]
 end
 R=A-D
 for i=1:n
  x0=inv(D)*(transpose(transpose(b))-R*transpose(transpose(x0)))
 end
 return x0
end
