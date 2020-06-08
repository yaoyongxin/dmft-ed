  do i=MpiIstart,MpiIend
     m = H%map(i)
     ib = bdecomp(m,Ns)
     !
     htmp = zero
     !
     ! LONG RANGE DENSITY-DENSITY INTERACTION TERMS
     do isite=1,Norb*Nbath
        !
        n_i = dble(ib(isite))
        !
        do distance=1,size(Neigh)
           !
           do jsite=1,Neigh(distance)
              n_j = dble(ib(Vstride(isite,distance,jsite)))
              htmp = htmp + Vnn(distance) * n_i * n_j
           enddo
           !
        enddo
        !
     enddo
     !
     select case(MpiStatus)
     case (.true.)
        call sp_insert_element(MpiComm,spH0,htmp,i,i)
     case (.false.)
        call sp_insert_element(spH0,htmp,i,i)
     end select
     !
     !
     !
     ! HOPPING TERMS
     do isite=1,Norb*Nbath
        do jsite=1,4
           !
           hopndx = Vstride(isite,1,jsite)
           !
           Jcondition = (Thopping/=zero) .AND. (ib(isite)==1) .AND. (ib(hopndx)==0 )
           if (Jcondition) then
              call c(isite,m,k1,sg1)
              call cdg(hopndx,k1,k2,sg2)
              j = binary_search(H%map,k2)
              htmp = Thopping * sg1 * sg2
              !
              select case(MpiStatus)
              case (.true.)
                 call sp_insert_element(MpiComm,spH0,htmp,i,j)
              case (.false.)
                 call sp_insert_element(spH0,htmp,i,j)
              end select
              !
           endif
           !
           Jcondition = (Thopping/=zero) .AND. (ib(hopndx)==1) .AND. (ib(isite)==0 )
           if (Jcondition) then
              call c(hopndx,m,k1,sg1)
              call cdg(isite,k1,k2,sg2)
              j = binary_search(H%map,k2)
              htmp = Thopping * sg1 * sg2
              !
              select case(MpiStatus)
              case (.true.)
                 call sp_insert_element(MpiComm,spH0,htmp,i,j)
              case (.false.)
                 call sp_insert_element(spH0,htmp,i,j)
              end select
              !
           endif
           !
        enddo
     enddo
     !
     !
     !
  enddo
