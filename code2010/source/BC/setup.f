





! boundary numbering:        1 - 'west'            2 - 'east'
!                            3 - 'front'           4 - 'back'
!                            5 - 'bottom'          6 - 'top'
!
       ibegin(west)=ix1
       iend(west)  =ix1
       jbegin(west)=jy1+1
       jend(west)  =jy2
       kbegin(west)=kz1+1
       kend(west)  =kz2
!
       ibegin(east)=ix2+1
       iend(east)  =ix2+1
       jbegin(east)=jy1+1
       jend(east)  =jy2
       kbegin(east)=kz1+1
       kend(east)  =kz2
!
       ibegin(front)=ix1+1
       iend(front)  =ix2
       jbegin(front)=jy1
       jend(front)  =jy1
       kbegin(front)=kz1+1
       kend(front)  =kz2
!
       ibegin(back)=ix1+1
       iend(back)  =ix2
       jbegin(back)=jy2+1
       jend(back)  =jy2+1
       kbegin(back)=kz1+1
       kend(back)  =kz2
!
       ibegin(bottom)=ix1+1
       iend(bottom)  =ix2
       jbegin(bottom)=jy1+1
       jend(bottom)  =jy2
       kbegin(bottom)=kz1
       kend(bottom)  =kz1
!
       ibegin(top)=ix1+1
       iend(top)  =ix2
       jbegin(top)=jy1+1
       jend(top)  =jy2
       kbegin(top)=kz2+1
       kend(top)  =kz2+1
!
!
!*** Starting and ending indexes for u,v,w momentum
!*** colocation points
!
!    u - momentum
       IF(ITYPE(west).EQ.500) THEN
        IBU = IX1+1
        IEU = IX2
        IPER= 0
       ELSE
        IBU = IX1+1
        IEU = IX2-1
        IPER= 1
       ENDIF
*    v - momentum
       IF(ITYPE(front).EQ.500) THEN
        JBV = JY1+1
        JEV = JY2
        JPER= 0
       ELSE
        JBV = JY1+1
        JEV = JY2-1
        JPER= 1
       ENDIF
*    w - momentum
       IF((ITYPE(bottom).EQ.500).or.(itype(top).eq.500).or.(ITYPE(top).EQ.0)) THEN
        KBW = KZ1+1
        KEW = KZ2
        KPER= 0
       ELSE
        KBW = KZ1+1
        KEW = KZ2-1
        KPER= 1
       ENDIF
!
       END subroutine boundary_initialize
!
       