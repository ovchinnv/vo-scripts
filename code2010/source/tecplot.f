C
C-----------------------------------------------------------------------
C                  ***************************                         
C                  *        tecplot.f        *                        
C                  ***************************                       
C----------------------------------------------------------------------- 
C
C This file contains all necessary routines to create a data-file
C for tecplot (ascii-format). 
C
C	- TECiniA: open data file
C	- TECzneA: write zone header
C	- TECdatA: write zone data
C	- TECnodA: write node data
C	- TECendA: close data file
C
C routine TECnodA needs some massage before using it!
C
C----------------------------------------------------------------------- 
C
C    INTEGER FUNCTION TECiniA(Title, Variables, FileName, FileUnit)
C    INTEGER FUNCTION TECzneA(ZoneTitle, imax, jmax, kmax, ZFormat, FileUnit)
C    INTEGER FUNCTION TECdatA(n, Data, FileUnit)
C    INTEGER FUNCTION TECnodA(EType, NData, FileUnit)
C    INTEGER FUNCTION TECendA(FileUnit)
C
C
C----------------------------------------------------------------------- 
C 
C-----FUNCTION-TECiniA--------------------------P. FLOHR--20/02/1994----
C
      INTEGER FUNCTION TECiniA(Title, Variables, FileName, FileUnit) 
C
C-----------------------------------------------------------------------
C
      CHARACTER*(*)  Title, Variables, FileName
      INTEGER        FileUnit
      CHARACTER*(72) Text
      INTEGER        ioc
C
      ioc = 0
C
      OPEN (UNIT=FileUnit,FILE=FileName,IOSTAT=ioc,STATUS='UNKNOWN')
      TECiniA = ioc
      REWIND (FileUnit)
      Text = 'TITLE = "'//Title//'"'
      WRITE(FileUnit,' (A) ') Text
      Text = 'VARIABLES = '//Variables
      WRITE(FileUnit,' (A) ') Text
      RETURN
      END
C
C 
C-----FUNCTION-TECzneA--------------------------P. FLOHR--20/02/1994----
C
      INTEGER FUNCTION TECzneA(ZoneTitle, imax, jmax, kmax, 
     &                         ZoneFormat, FileUnit) 
C
C-----------------------------------------------------------------------
C
      CHARACTER*(*)  ZoneTitle, ZoneFormat
      INTEGER        imax, jmax, kmax, FileUnit
      INTEGER        ioc
C
      ioc = 0
C
      WRITE(UNIT=FileUnit,FMT=1000,IOSTAT=ioc) ZoneTitle
      WRITE(UNIT=FileUnit,FMT=2000,IOSTAT=ioc) imax,jmax,kmax,ZoneFormat
 1000 FORMAT ('ZONE T = "',(A),'"')
 2000 FORMAT ('I = ',I4,', J = ',I4,', K = ',I4,', F = ',(A))
C
      TECzneA = ioc
      RETURN
      END
C
C 
C-----FUNCTION-TECdatA--------------------------P. FLOHR--20/02/1994----
C
      INTEGER FUNCTION TECdatA(n, Data, FileUnit)
C
C-----------------------------------------------------------------------
C
      INTEGER   n
      DIMENSION Data(n)
      INTEGER   FileUnit
      INTEGER   ioc, i
C
      ioc = 0     
C
      WRITE(UNIT=FileUnit,FMT=*,IOSTAT=ioc) (Data(i), i = 1,n)
      TECdatA = ioc
      RETURN
      END
C
C 
C-----FUNCTION-TECnodA--------------------------P. FLOHR--20/02/1994----
C
      INTEGER FUNCTION TECnodA(ElementType, n, NData, FileUnit) 
C
C-----------------------------------------------------------------------
C
      CHARACTER*(*) ElementType
      INTEGER       n, FileUnit
      DIMENSION     NData(n)
      INTEGER       ioc
C
      ioc = 0
C
      TECnodA = ioc
      RETURN
      END
C
C 
C-----FUNCTION-TECendA--------------------------P. FLOHR--20/02/1994----
C
      INTEGER FUNCTION TECendA(FileUnit) 
C
C-----------------------------------------------------------------------
C
      INTEGER FileUnit
      INTEGER ioc
C
      ioc = 0
C
      CLOSE (UNIT=FileUnit,IOSTAT=ioc)
C
      TECendA = ioc
      RETURN
      END
C
C-----------------------------------------------------------------------
C
