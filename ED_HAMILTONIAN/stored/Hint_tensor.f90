  if(Utensor)then
     !
     !
     do i=MpiIstart,MpiIend
        m = H%map(i)
        ib = bdecomp(m,2*Ns)
        !
        do iorb=1,Norb
           nup(iorb)=dble(ib(iorb))
           ndw(iorb)=dble(ib(iorb+Ns))
        enddo


        ! FULL INTERACTION IN THE FORM:
        !    ( c^+_iorb_s  c_jorb_s )( c^+_korb_s'  c_lorb_s' )
        if(Norb>1)then
           do iorb=1,Norb
              do jorb=1,Norb
                 do korb=1,Norb
                    do lorb=1,Norb
                       !
                       !Hartree block {s,s'}={up,dw}
                       Jcondition=(&
                       (ib(lorb+Ns)==1).AND.&
                       (ib(korb+Ns)==0).AND.&
                       (ib(jorb)==1).AND.&
                       (ib(iorb)==0))
                       !
                       if(Jcondition)then
                          call   c(lorb+Ns, m,k1,sg1)
                          call cdg(korb+Ns,k1,k2,sg2)
                          call   c(jorb,k2,k3,sg3)
                          call cdg(iorb,k3,k4,sg4)
                          j=binary_search(H%map,k4)
                          htmp = one*Umat(iorb,jorb,korb,lorb,1)*sg1*sg2*sg3*sg4
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
                       !Hartree block {s,s'}={up,up}/{dw,dw}
                       do ispin=1,2
                          Jcondition=(&
                         (ib(lorb+Ns*(ispin-1))==1).AND.&
                         (ib(korb+Ns*(ispin-1))==0).AND.&
                         (ib(jorb+Ns*(ispin-1))==1).AND.&
                         (ib(iorb+Ns*(ispin-1))==0))
                         !
                         if(Jcondition)then
                            call   c(lorb+Ns*(ispin-1), m,k1,sg1)
                            call cdg(korb+Ns*(ispin-1),k1,k2,sg2)
                            call   c(jorb+Ns*(ispin-1),k2,k3,sg3)
                            call cdg(iorb+Ns*(ispin-1),k3,k4,sg4)
                            j=binary_search(H%map,k4)
                            htmp = one*0.5d0*Umat(iorb,jorb,korb,lorb,2)*sg1*sg2*sg3*sg4
                            !
                            select case(MpiStatus)
                            case (.true.)
                               call sp_insert_element(MpiComm,spH0,htmp,i,j)
                            case (.false.)
                               call sp_insert_element(spH0,htmp,i,j)
                            end select
                            !
                         endif
                       enddo
                       !
                    enddo
                 enddo
              enddo
           enddo
        endif


     enddo
     !
     !
  endif
