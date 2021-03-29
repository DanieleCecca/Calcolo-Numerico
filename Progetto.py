# -*- coding: utf-8 -*-
"""
Created on Mon Feb  1 09:59:35 2021

@author: Daniele
"""

import numpy as np
import sympy as sp




def RiduzioneScalini(A):
    
    """
    Algoritmo di riduzione a gradini.

    INPUT
    A: matrice appartenente alli'insieme delle matrici R^m*n
    
    OUTPUT
    B:matrice ridotta
    """
    B=np.copy(A)
    righe,colonne=A.shape
    tol=1e-15
    
    j=0
    i=0
    while i < righe-1:
                
        #TECNICA DEL PIVOT
        zero=False
        if abs(B[i][j]) < tol:
            zero=True
            for g in range(i+1,righe):
                if abs(B[g][j]) > tol:
                     B[[i, g]]=B[[g, i]]
                     zero=False
                     break
                 
          #se sono tutti zero controlla a destra      
            if zero:
                j=j+1
                while j < colonne:
                    if abs(B[i][j]) > tol:
                        zero=False
                        break
                    j=j+1
        if zero:
            break
        
     #calcolo del moltiplicatore di Gauss   
        for k in range(i+1,righe):
            m = -B[k][j]/B[i][j]

            
           #calcolo nuova riga
            for h in range (j,colonne):
                somma = B[i][h]*m + B[k][h]
                B[k][h]=somma

     

        j+=1
        i+=1

    return B


def rank(A):
    """
    Algoritmo per poter calcolare il rango di una matrice

    INPUT
    A: matrice appartenente alli'insieme delle matrici R^m*n
    
    OUTPUT
    rank : rango di A
    """
    B=RiduzioneScalini(A)
    rank=0
    tol=1e-15
    righe,colonne = B.shape
    j=0
    for i in range(righe):
        while j < colonne:
            if abs(B[i][j]) > tol:
                rank+=1
                break
            j+=1
    return rank


def base(A):
    """
    Algoritmo che permette di calcolare la base dell'immagine di A.

    INPUT
    A: matrice appartenente alli'insieme delle matrici R^m*n
    
    OUTPUT
    base : base dell'immagine di A
    """
    
    B=RiduzioneScalini(A)
    dim=rank(A)
    pos=np.zeros(dim)
    righe,colonne = B.shape
    Base=np.zeros((righe,dim))
    
    tol=1e-15
    
    j=0
    for i in range(dim):
        while j < colonne:
            if abs(B[i][j]) > tol:
                pos[i]=j
                break
            j+=1
            
    for i in range (dim):
        Base[: , i] = A[:, int(pos[i])]
    return Base


def libere(A):
    
    B=RiduzioneScalini(A)
    dim=rank(A)
    righe,colonne = B.shape
    pos=np.zeros(dim)
    posLibere=np.zeros(colonne-dim)
    
    Libere=np.zeros((righe,colonne-dim))
    
    tol=1e-15
    
    j=0
    for i in range(dim):
        while j < colonne:
            if abs(B[i][j]) > tol:
                pos[i]=j
                break
            j+=1
    k=0
    valore=0  
    while k < colonne - dim:
        verifica=False
        for i in range(dim):
           if valore!=pos[i]:  
               verifica=True
           else:
               verifica=False
               break
        if verifica==True:
           posLibere[k]=valore
           
           k+=1
           valore+=1
        else:
            valore+=1
        
            
    for i in range (colonne - dim):
        Libere[: , i] = A[:, int(posLibere[i])]
    return Libere


def rimuoviRigaZeri(A):
    righe,colonne=A.shape
    
    verifica=False
    for i in range (righe-1,-1,-1):
        for j in range (colonne):
            if A[i][j]==0:
                verifica=True
            else:
                verifica=False
                break
        if verifica==True:
            A=np.delete(A,i,0)
    return A
    

#funziona solo se b iniziale =0
def solSottoDim(A):
    tol=1e-15
    
    
    righe,colonne=A.shape
    dim=rank(A)
    numParametri=(colonne-dim)
    
    x=np.zeros((colonne, numParametri))
    pos=np.zeros(dim)
    posLibere=np.zeros(colonne-dim)
    
    C=RiduzioneScalini(A)

    B=base(rimuoviRigaZeri(RiduzioneScalini(A)))
    c=(libere(rimuoviRigaZeri(RiduzioneScalini(A))))*-1
    
    sol=np.linalg.solve(B,c)

    

    
      
    j=0
    for i in range(dim):
        while j < colonne:
            if abs(C[i][j]) > tol:
                pos[i]=j
                break
            j+=1
    
    k=0
    valore=0  
    while k < numParametri:
        verifica=False
        for i in range(dim):
           if valore!=pos[i]:  
               verifica=True
           else:
               verifica=False
               break
        if verifica==True:
           posLibere[k]=valore
           
           k+=1
           valore+=1
        else:
            valore+=1
    
    
    for j in range(numParametri):
        for i in range (colonne):
            for k in range(dim):
                  if i==pos[k]:

                      x[i][j]=sol[k][j]
                      break
    j=0            
    while j <numParametri:
        for i in range (colonne):             
            for k in range(numParametri):
                  if i==posLibere[k]:
                      x[i][j]=1
                      j+=1
                      break
                      
   
    return x  


def RisolSistema(A,b):
    righe,colonne=A.shape
   
    B=RiduzioneScalini(A)
    c=np.copy(b)
    
    Ab=np.hstack([A,b])
    Ab=RiduzioneScalini(Ab)
    
    c[:,0]=Ab[:,colonne]
    
    print("\nAb:\n",Ab,"\n")


    
    if rank(B) != rank(Ab):
        print("sistema non compatibile")
        return 0
    
   
    if colonne > righe or colonne < righe:
        return solSottoDim(A)
    
    else:
        print(B)
        print(c)
        return np.linalg.solve(B,c)


def kern(A):
    """
    Algoritmo che permette di calcolare la base del kernel.

    INPUT
    A: matrice appartenente alli'insieme delle matrici R^m*n
    
    OUTPUT
    kern : base del kernel di A
    """
    righe,colonne=A.shape
    b=np.zeros((righe,1))
    return RisolSistema(A,b)
 
   
def scomposizione(A,x):
    """
    Algoritmo che permette di scomporre un vettore come proiezioni 
    di due vettori in Im(A) e nel suo ortogonale.

    INPUT
    A: matrice appartenente alli'insieme delle matrici R^m*n
    x: vettore da scomporre
    
    OUTPUT
    kern : base del kernel di A
    """
    
    b1=np.zeros((len(x),1))
    b2=np.zeros((len(x),1))
    
    w=kern(A.T)
    v=base(A)
    C=np.hstack([v,w])
    
    b=RisolSistema(C,x)
    print("\nb:\n",b,"\n")
    
    righeBase,colonneBase=v.shape
    righeOrt,colonneOrt=w.shape
   
    j=0
    i=0
    while j <colonneBase:
        b1[:,0]=(b[j][0]*v[:,j])+b1[:,0]
        j+=1
        
    while i <colonneOrt:
        b2[:,0]=(b[j][0]*w[:,i])+b2[:,0]
        j+=1
        i+=1
        
        
    print("b1:\n",b1,"\n")
    print("b2:\n",b2,"\n")
                
                


A=np.array([[-1,-3,2,-3],
            [1,3,1,0],
            [-1,-3,5,-6]])

AT=A.T

x=np.array([[0],
            [6],
            [0]])

B=RiduzioneScalini(A)
rango=rank(A)
Base=base(A)
Libere=libere(A)






print("A:\n",A,"\n")
print("B:\n",B,"\n")
print("rank:",rango,"\n")
print("base:\n",Base,"\n")
print("libere:\n",Libere,"\n")
print("AT:\n",AT,"\n")

#C=rimuoviRigaZeri(B)
#print("C:\n",C,"\n")


Kern=kern(A)
print("base kern:\n",Kern,"\n")

Ortogonale=kern(AT)
print("base ortogonale:\n",Ortogonale,"\n")


scomposizione(A,x)




