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
                  if((iorb==jorb).and.(korb==lorb))then
                    Jcondition=((ib(jorb)==1).AND.(ib(lorb+Ns)==1))
                  elseif((iorb==jorb).and.(korb/=lorb))then
                    Jcondition=((ib(jorb)==1).AND.(ib(korb+Ns)==0).AND.(ib(lorb+Ns)==1))
                  elseif((iorb/=jorb).and.(korb==lorb))then
                    Jcondition=((ib(iorb)==0).AND.(ib(jorb)==1).AND.(ib(lorb+Ns)==1))
                  elseif((iorb/=jorb).and.(korb/=lorb))then
                    Jcondition=((ib(iorb)==0).AND.(ib(jorb)==1).AND.(ib(korb+Ns)==0).AND.(ib(lorb+Ns)==1))
                  else
                    Jcondition=.false.
                  endif
                  !
                  if(Jcondition)then
                    call   c(lorb+Ns, m,k1,sg1)
                    call cdg(korb+Ns,k1,k2,sg2)
                    call   c(jorb,k2,k3,sg3)
                    call cdg(iorb,k3,k4,sg4)
                    !
                    htmp = one*Umat(iorb,jorb,korb,lorb,1)*sg1*sg2*sg3*sg4
                    j=binary_search(H%map,k4)
                    !
                    if((abs(htmp).ne.0d0).and.(j.gt.0))then
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
                  endif
                  !
                  !Hartree block {s,s'}={up,up}/{dw,dw}
                  do ispin=1,2
                    if((iorb==jorb).and.(korb==lorb))then
                       Jcondition=((ib(jorb+Ns*(ispin-1))==1).AND.(ib(lorb+Ns*(ispin-1))==1))
                    elseif((iorb==jorb).and.(korb/=lorb).and.(iorb/=korb).and.(jorb/=lorb))then
                       Jcondition=((ib(jorb+Ns*(ispin-1))==1).AND.(ib(korb+Ns*(ispin-1))==0).AND.(ib(lorb+Ns*(ispin-1))==1))
                    elseif((iorb/=jorb).and.(korb==lorb).and.(iorb/=korb).and.(jorb/=lorb))then
                       Jcondition=((ib(iorb+Ns*(ispin-1))==0).AND.(ib(jorb+Ns*(ispin-1))==1).AND.(ib(lorb+Ns*(ispin-1))==1))
                    elseif((iorb/=jorb).and.(korb/=lorb).and.(iorb/=korb).and.(jorb/=lorb))then
                       Jcondition=((ib(iorb+Ns*(ispin-1))==0).AND.(ib(jorb+Ns*(ispin-1))==1).AND.(ib(korb+Ns*(ispin-1))==0).AND.(ib(lorb+Ns*(ispin-1))==1))
                    else
                       Jcondition=.false.
                    endif
                    !
                    if(Jcondition)then
                       call   c(lorb+Ns*(ispin-1), m,k1,sg1)
                       call cdg(korb+Ns*(ispin-1),k1,k2,sg2)
                       call   c(jorb+Ns*(ispin-1),k2,k3,sg3)
                       call cdg(iorb+Ns*(ispin-1),k3,k4,sg4)
                       !
                       htmp = one*Umat(iorb,jorb,korb,lorb,2)*sg1*sg2*sg3*sg4
                       j=binary_search(H%map,k4)
                       !
                       if((abs(htmp).ne.0d0).and.(j.gt.0))then
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
                    endif
                  enddo
                  !
               enddo
            enddo
         enddo
     enddo
     !
     !
  endif


enddo
