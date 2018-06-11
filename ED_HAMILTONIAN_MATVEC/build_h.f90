  states: do i=first_state,last_state
     m = H%map(i)
     impi = i-ishift
     ib = bdecomp(m,2*Ns)
     !
     do iorb=1,Norb
        nup(iorb)=dble(ib(iorb))
        ndw(iorb)=dble(ib(iorb+Ns))
     enddo
     !
     !
     !IMPURITY  HAMILTONIAN
     include "ED_HAMILTONIAN_MATVEC/Himp.f90"
     !
     !LOCAL INTERACTION
     select case(kanamori)
     case default
        include "ED_HAMILTONIAN_MATVEC/Hint.f90"
     case (.false.)
        !>DEBUG
        U0=Uloc(1)/2.d0
        U1=Ust/2.d0
        Nflavors=2*Norb
        allocate(nvec(Nflavors))      !1:Norb = n_up, Norb+1:2Norb=n_dw
        allocate(Umatrix(Nflavors,Nflavors))
!       Umatrix(1,:)=[0d0,U1,U0,U1,U0,U1,U0,U1]
!       Umatrix(2,:)=[U1,0d0,U1,U0,U1,U0,U1,U0]
!       Umatrix(3,:)=[U0,U1,0d0,U1,U0,U1,U0,U1]
!       Umatrix(4,:)=[U1,U0,U1,0d0,U1,U0,U1,U0]
!       Umatrix(5,:)=[U0,U1,U0,U1,0d0,U1,U0,U1]
!       Umatrix(6,:)=[U1,U0,U1,U0,U1,0d0,U1,U0]
!       Umatrix(7,:)=[U0,U1,U0,U1,U0,U1,0d0,U1]
!       Umatrix(8,:)=[U1,U0,U1,U0,U1,U0,U1,0d0]

!       >TEST: 2 orbitals square lattice
!       Umatrix(1,:)=[0.d0,0.d0,U0,0.d0]
!       Umatrix(2,:)=[0.d0,0.d0,0.d0,U0]
!       Umatrix(3,:)=[U0,0.d0,0.d0,0.d0]
!       Umatrix(4,:)=[0.d0,U0,0.d0,0.d0]
!       <TEST
!       >TEST: 4 orbitals square lattice
        Umatrix(1,:)=[0.d0,0.d0,0.d0,0.d0,U0,0.d0,0.d0,0.d0]
        Umatrix(2,:)=[0.d0,0.d0,0.d0,0.d0,0.d0,U0,0.d0,0.d0]
        Umatrix(3,:)=[0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,U0,0.d0]
        Umatrix(4,:)=[0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,U0]
        Umatrix(5,:)=[U0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0]
        Umatrix(6,:)=[0.d0,U0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0]
        Umatrix(7,:)=[0.d0,0.d0,U0,0.d0,0.d0,0.d0,0.d0,0.d0]
        Umatrix(8,:)=[0.d0,0.d0,0.d0,U0,0.d0,0.d0,0.d0,0.d0]
!       <TEST
        nvec = [nup,ndw]
 
        include "ED_HAMILTONIAN_MATVEC/Hint_Umatrix.f90"

        deallocate(nvec,Umatrix)

     end select
     !
     !BATH HAMILTONIAN
     include "ED_HAMILTONIAN_MATVEC/Hbath.f90"
     !
     !IMPURITY- BATH HYBRIDIZATION
     include "ED_HAMILTONIAN_MATVEC/Himp_bath.f90"
     !
     !
  enddo states



