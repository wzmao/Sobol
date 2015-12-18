# -*- coding: utf-8 -*-
"""This module contains some Sobol functions.
"""

__author__ = 'Fei Xie'
__all__ = ["Sobolpts"]

def Sobolmat(c,m,r): #c(q+1维列表)是q次多项式系数，m(q维列表)是m初始值，r是所需方向数的个数，输出V是r阶上三角阵
	from numpy import arange,zeros,array
	q=c.size-1
	c*=(2**arange(q+1))
	mm=m
	m=zeros(r,dtype=int)
	if len(mm)>r:
		m=array(mm)
	else:
		m[:len(mm)]=array(mm)
	if q==0:
		V=eye(r)
	elif q>0:
		V=zeros((r,r))
		for i in range (q,r):
			n=0
			for j in c[1:]*m[i-1::-1][:q]:
				n=n^j
			n=n^m[i-q]
			m[i]=n
		for i in range (0,r):
			V[:i+1,i]=Binary(m[i],i+1)
	return V

def Binary(k,r=0): #自然数k的二进制表示，输出列表
	from numpy import array,binary_repr,int8
	return (array(list(binary_repr(k,width=r)))=='1').astype(int8)

def Sobolpts(n0,N,dim,M=1,**kwarg): 
	'''n0(正整数)点列的起始点，d维数，p(d*1)多项式(每个数的二进制表示是多项式系数)，m(q*d)初始值矩阵，每一维对应的m的维数q不同，输出P(M*N*dim)Sobol点列'''

	from .Constant import getconstant
	# from numpy import array,zeros,ones,empty,eye,exp,sqrt,log,log2,diag,transpose,arange,int8,binary_repr,zeros_like,hstack,inf
	from numpy import log2,zeros,int8
	m = getconstant("m0")[:dim]
	p = getconstant("p0")[:dim]
	nmax = n0+N-1
	rmax = 1+int(log2(nmax)) #nmax的二进制位数，方向阵的阶数
	r = 1 if n0==1 else 1+int(log2(n0-1))
	#生成V，每一维有对应的不同的方向阵，都是rmax阶的上三角阵
	V=zeros((dim,rmax,rmax),dtype=int8)
	for d in range(dim):
		q=int(log2(p[d])) #对应维多项式的阶数，m的维数
		c=Binary(p[d]) #q维向量
		V[d]=Sobolmat(c,m[d],rmax)
	return 0
	#计算第0个点的y值
	a=Binary(n0-1)
	g=Graycode(n0-1)
	Y=zeros((dim,rmax),dtype=int8)
	for d in range (0,dim):
		for i in range (0,r):
			temp=0
			for j in range (i,r):
				temp=temp^(V[d,i,j]*g[r-i-1])
			Y[d,i]=temp
	#Scrambling
	L=zeros((M,dim,rmax,rmax),dtype=int8)
	e=randint(0,2,size=(M,dim,rmax))
	if M==1:
		for i in range(dim):
			L[0,i]=eye(rmax)
		e[0,:]=0
	elif M>1:
		for j in range(rmax):
			for m in range(M):
				for i in range(dim):
					L[m,i,j,j]=1
		for m in range(M):
			for i in range(dim):
				for j in range(1,rmax):
					L[m,i,j,:j]=randint(0,2,size=j)	
	# 由前一个点的y值计算当前点的y值，并得到最终的Sobol点
	qnext=2**r #记录下一次r需要加1的位置
	P=zeros((M,N,dim))
	pwrs=(0.5)**(arange(rmax)+1)
	for k in range(n0,nmax+1):
		if k==qnext:
			r+=1
			l=r-1
			qnext*=2
		else:
			for i in range(len(a)-1,-1,-1):
				if a[i]==0:
					l=r-1-i
					break
		a=Nextbinary(a)
		for d in range (0,dim):
			# Y[d]^=V[d,:,l]
			for i in range(0,l+1):
				Y[d,i]^=V[d,i,l]
			for m in range(M):
				P[m,k-n0,d]=((L[m,d].dot(Y[d])+e[m,d])%2).dot(pwrs)
	return P
