  ! IMP UP <--> BATH UP
  do iup=1,dimdw
     mup  = H%up(iup)
     nup  = bdecomp(mup,Ns)
     impi_up = iup - ishift_dw

     do iorb=1,Norb
        do kp=1,Nbath
           alfa=getBathStride(iorb,kp)
           !
           if( (diag_hybr(1,iorb,kp)/=0d0) .AND. (nup(iorb)==1) .AND. (nup(alfa)==0) )then              
              call c(iorb,m,k1,sg1)
              call cdg(alfa,k1,k2,sg2)
              jup = binary_search(H%up,k2)
              htmp = diag_hybr(1,iorb,kp)*sg1*sg2
              !
              call sp_insert_element(spH0up,htmp,impi_up,jup)
              !
           endif
           if( (diag_hybr(1,iorb,kp)/=0d0) .AND. (ib(iorb)==0) .AND. (ib(alfa)==1) )then
              call c(alfa,m,k1,sg1)
              call cdg(iorb,k1,k2,sg2)
              jup=binary_search(H%up,k2)
              htmp = diag_hybr(1,iorb,kp)*sg1*sg2
              !
              call sp_insert_element(spH0up,htmp,impi_up,jup)
              !
           endif
        enddo
     enddo
     !
  enddo

  !IMP DW <--> BATH DW
  do idw=1,dimdw
     mdw  = H%dw(idw)
     ndw  = bdecomp(mdw,Ns)
     impi_dw = idw - ishift_dw

     do iorb=1,Norb
        do kp=1,Nbath
           alfa=getBathStride(iorb,kp)
           !
           if( (diag_hybr(Nspin,iorb,kp)/=0d0) .AND. (ib(iorb)==1) .AND. (ib(alfa)==0) )then
              call c(iorb,m,k1,sg1)
              call cdg(alfa,k1,k2,sg2)
              jdw=binary_search(H%dw,k2)
              htmp=diag_hybr(Nspin,iorb,kp)*sg1*sg2
              !
              call sp_insert_element(spH0dw,htmp,impi_dw,jdw)
              !
           endif
           if( (diag_hybr(Nspin,iorb,kp)/=0d0) .AND. (ib(iorb)==0) .AND. (ib(alfa)==1) )then
              call c(alfa,m,k1,sg1)
              call cdg(iorb,k1,k2,sg2)
              jdw=binary_search(H%dw,k2)
              htmp=diag_hybr(Nspin,iorb,kp)*sg1*sg2
              !
              call sp_insert_element(spH0dw,htmp,impi_dw,jdw)
              !
           endif
        enddo
     enddo
     !
  enddo
