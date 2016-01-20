      SUBROUTINE tridag(a,b,c,r,u,n,np)
      INTEGER n,NMAX,np
      REAL a(np),b(np),c(np),r(np),u(np)
      PARAMETER (NMAX=500)
      INTEGER j
      REAL bet,gam(NMAX)
      if(b(1).eq.0)pause'tridag: rewrite equations'
      bet=b(1)
      u(1)=r(1)/bet
      do 11 j=2,n
        gam(j)=c(j-1)/bet
        bet=b(j)-a(j)*gam(j)
        if(bet.eq.0)pause 'tridag failed'
        u(j)=(r(j)-a(j)*u(j-1))/bet
11    continue
      do 12 j=n-1,1,-1
        u(j)=u(j)-gam(j+1)*u(j+1)
12    continue
      return
      END
c
      SUBROUTINE tridagc(a,b,c,r,u,n,np)
      INTEGER n,NMAX,np
      complex a(np),b(np),c(np),r(np),u(np)
      PARAMETER (NMAX=500)
      INTEGER j
      complex bet,gam(NMAX)
      if(b(1).eq.(0.,0.))pause 'tridag: rewrite equations'
      bet=b(1)
      u(1)=r(1)/bet
      do 11 j=2,n
        gam(j)=c(j-1)/bet
        bet=b(j)-a(j)*gam(j)
        if(bet.eq.(0.,0.))pause 'tridag failed'
        u(j)=(r(j)-a(j)*u(j-1))/bet
11    continue
      do 12 j=n-1,1,-1
        u(j)=u(j)-gam(j+1)*u(j+1)
12    continue
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software %5'KH!9.
