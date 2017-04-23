  if(bath_type/="replica") then
     !
     !diagonal bath hamiltonian: +energy of the bath=\sum_a=1,Norb\sum_{l=1,Nbath}\e^a_l n^a_l
     htmp=zero
     do iorb=1,size(dmft_bath%e,2)
        do kp=1,Nbath
           alfa=getBathStride(iorb,kp)
           htmp =htmp + &
                dmft_bath%e(1,iorb,kp)*ib(alfa) +& !UP
                dmft_bath%e(Nspin,iorb,kp)*ib(alfa+Ns) !DW
        enddo
     enddo
     !
     call sp_insert_element(spH0,htmp,impi,i)
     !
  else
     !`
     !diagonal bath hamiltonian: +energy of the bath=\sum_a=1,Norb\sum_{l=1,Nbath}\e^a_l n^a_l
     htmp=zero
     do kp=1,Nbath
        do iorb=1,Norb
           alfa = getBathStride(iorb,kp) !iorb + kp*Norb
           htmp = htmp + &
                dmft_bath%h(1,1,iorb,iorb,kp)*ib(alfa) + & !UP
                dmft_bath%h(Nspin,Nspin,iorb,iorb,kp)*ib(alfa+Ns) !DW
        enddo
     enddo
     !
     call sp_insert_element(spH0,htmp,impi,i)
     !
     !off-diagonal elements
     do kp=1,Nbath
        do iorb=1,Norb
           do jorb=1,Norb
              if(iorb==jorb)cycle
              !UP
              alfa = getBathStride(iorb,kp)
              beta = getBathStride(jorb,kp)
              Jcondition = &
                   (dmft_bath%h(1,1,iorb,jorb,kp)/=zero) .AND. &
                   (ib(beta)==1)                         .AND. &
                   (ib(alfa)==0)
              if (Jcondition)then
                 call c(beta,m,k1,sg1)
                 call cdg(alfa,k1,k2,sg2)
                 j = binary_search(H%map,k2)
                 htmp = dmft_bath%h(1,1,iorb,jorb,kp)*sg1*sg2
                 !
                 call sp_insert_element(spH0,htmp,impi,j)
                 !
              endif
              !DW
              alfa = getBathStride(iorb,kp)
              beta = getBathStride(jorb,kp)
              Jcondition = &
                   (dmft_bath%h(Nspin,Nspin,iorb,jorb,kp)/=zero) .AND. &
                   (ib(beta+Ns)==1)                         .AND. &
                   (ib(alfa+Ns)==0)
              if (Jcondition)then
                 call c(beta+Ns,m,k1,sg1)
                 call cdg(alfa+Ns,k1,k2,sg2)
                 j = binary_search(H%map,k2)
                 htmp = dmft_bath%h(Nspin,Nspin,iorb,jorb,kp)*sg1*sg2
                 !
                 call sp_insert_element(spH0,htmp,impi,j)
                 !
              endif
              !
           enddo
        enddo
     enddo
     !
  endif



  !anomalous pair-creation/destruction
  if(ed_mode=="superc")then
     do iorb=1,size(dmft_bath%e,2)
        do kp=1,Nbath
           ms=getBathStride(iorb,kp)
           !\Delta_l c_{\up,ms} c_{\dw,ms}
           if( (dmft_bath%d(1,iorb,kp)/=0d0) .AND. (ib(ms)==1) .AND. (ib(ms+Ns)==1) )then
              call c(ms,m,k1,sg1)
              call c(ms+Ns,k1,k2,sg2)
              j=binary_search(H%map,k2)
              htmp=one*dmft_bath%d(1,iorb,kp)*sg1*sg2
              !
              call sp_insert_element(spH0,htmp,impi,j)
              !
           endif
           !\Delta_l cdg_{\up,ms} cdg_{\dw,ms}
           if( (dmft_bath%d(1,iorb,kp)/=0d0) .AND. (ib(ms)==0) .AND. (ib(ms+Ns)==0) )then
              call cdg(ms+Ns,m,k1,sg1)
              call cdg(ms,k1,k2,sg2)
              j=binary_search(H%map,k2)
              htmp=one*dmft_bath%d(1,iorb,kp)*sg1*sg2 !
              !
              call sp_insert_element(spH0,htmp,impi,j)
              !
           endif
        enddo
     enddo
  endif
