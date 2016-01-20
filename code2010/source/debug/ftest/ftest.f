      program ftest

      integer n
      parameter (n=100)
      
      integer ifax(13)
      real trigs(3*n/2+1), work(2*n)
      real a(n+2)
      complex ac(n/2+1)
      integer i,j,k
      
      equivalence (a, ac)
      

      do i=2,n+1
       a(i)=1.*i
      enddo 
      a(n+2)=a(2)
      a(1)=a(n+1)
 
      call fftfax_4p(n,ifax,trigs)
      call fft99_4p(a,work,trigs,ifax,1,n,n,1,-1)
      call fft99_4p(a,work,trigs,ifax,1,n,n,1,1)

      do i=1,n/2
       write(0,*) i,ac(i)
      enddo  
      


      end      