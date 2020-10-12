  do i=MpiIstart,MpiIend
     m_8 = H8%map(i)
     ib = bdecomp8(m_8,Ns)
     !
     htmp = zero
     !
     ! LONG RANGE DENSITY-DENSITY INTERACTION TERMS
     do isite=1,Norb*Nbath
        !
        n_i = dble(ib(isite))
        !
        do distance=1,size(Neigh_int)
           do ineig=firstNeig,Neigh_int(distance)
              !
              jsite = Vstride(isite,distance,ineig)
              !
              n_j = dble(ib(jsite))
              htmp = htmp + Vmat(distance) * n_i * n_j
              !
           enddo
           !
        enddo
        !
        !symmetry breaking term
        htmp = htmp + Ust * Radius(vec2lat(isite,1),vec2lat(isite,2)) * n_i
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
              call c8(isite,m_8,k1_8,sg1)
              call cdg8(hopndx,k1_8,k2_8,sg2)
              j = binary_search8(H8%map,k2_8)
              htmp = Thopping * (sg1*sg2)**BosonExp
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
