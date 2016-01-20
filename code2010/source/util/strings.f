      module strings
!
!     routines to manipulate output and character type conversion
!     

      public atoi
      public atof
      public isempty

      contains


      character*3 index(0:999)
      character*5 index2(0:99999)


!...  fill index array for various output purposes
      m=0
      m2=0
      do i=0,9
         do j=0,9
            do k=0,9
               index3(m)=CHAR(i+48)//CHAR(j+48)//CHAR(k+48)
               m=m+1
      	       do jj=0,9
      	         do kk=0,9
      		   index5(m2)=char(i+48)//char(j+48)//char(k+48)//char(jj+48)//char(kk+48)
      		   m2=m2+1
                 enddo
      	       enddo
            enddo
         enddo
      enddo

!
       function atoi(a)
!    convert string to integer
       implicit none
!
       int :: atoi
       character*(*) :: a
!
       end function atoi
! ****************************
       function atof(a)
!    convert string to float
       implicit none
!
       float :: atof
       character*(*) :: a
!
       end function atof
! ****************************
      end module strings
!