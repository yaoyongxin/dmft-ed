MODULE ED_INPUT_VARS
  USE SF_VERSION
  USE SF_PARSE_INPUT
  USE SF_IOTOOLS, only:str
  implicit none

  !GIT VERSION
  include "revision.inc"  !this file is generated at compilation time in the Makefile


  !input variables
  !=========================================================
  integer              :: Nbath               !Nbath=# of bath sites (per orbital or not depending on bath_type)
  integer              :: Norb                !Norb =# of impurity orbitals
  integer              :: Nspin               !Nspin=# spin degeneracy (max 2)
  integer              :: nloop               !max dmft loop variables
  real(8),dimension(3) :: Uloc                !local interactions
  real(8)              :: Ust                 !intra-orbitals interactions
  real(8)              :: Jh                  !J_Hund: Hunds' coupling constant
  real(8)              :: Jx                  !J_X: coupling constant for the spin-eXchange interaction term
  real(8)              :: Jp                  !J_P: coupling constant for the Pair-hopping interaction term
  real(8)              :: xmu                 !chemical potential
  real(8)              :: deltasc             !breaking symmetry field
  real(8)              :: beta                !inverse temperature
  real(8)              :: eps                 !broadening
  real(8)              :: wini,wfin           !
  integer              :: Nsuccess            !
  logical              :: chispin_flag             !
  logical              :: chidens_flag             !
  logical              :: chipair_flag             !
  logical              :: HFmode              !flag for HF interaction form U(n-1/2)(n-1/2) VS Unn
  real(8)              :: cutoff              !cutoff for spectral summation
  real(8)              :: gs_threshold        !Energy threshold for ground state degeneracy loop up
  real(8)              :: dmft_error          !dmft convergence threshold
  real(8)              :: sb_field            !symmetry breaking field
  real(8)              :: hwband              !half-bandwidth for the bath initialization: flat in -hwband:hwband
  !
  integer              :: ed_verbose          !
  logical              :: ed_sparse_H         !flag to select  storage of sparse matrix H (mem--, cpu++) if TRUE, or direct on-the-fly H*v product (mem++, cpu--) if FALSE
  logical              :: ed_solve_offdiag_gf !flag to select the calculation of the off-diagonal impurity GF. this is T by default if bath_type/=normal
  logical              :: ed_print_Sigma      !flag to print impurity Self-energies
  logical              :: ed_print_G          !flag to print impurity Green`s functions
  logical              :: ed_print_G0         !flag to print impurity non-interacting Green`s functions
  logical              :: ed_twin             !flag to reduce (T) or not (F,default) the number of visited sector using twin symmetry.
  logical              :: ed_sectors          !flag to reduce sector scan for the spectrum to specific sectors +/- ed_sectors_shift
  integer              :: ed_sectors_shift    !shift to the ed_sectors scan
  character(len=7)     :: ed_mode             !flag to set ed symmetry type: normal=normal (default), superc=superconductive, nonsu2=broken SU(2)
  logical              :: ed_para             !flag to enforce non-magnetic solution
  real(8)              :: ed_vsf_ratio        !
  real(8)              :: ed_bath_noise_thr   !
  !
  character(len=12)    :: lanc_method         !select the lanczos method to be used in the determination of the spectrum. ARPACK (default), LANCZOS (T=0 only)
  real(8)              :: lanc_tolerance      !Tolerance for the Lanczos iterations as used in Arpack and plain lanczos.
  integer              :: lanc_niter          !Max number of Lanczos iterations
  integer              :: lanc_ngfiter        !Max number of iteration in resolvant tri-diagonalization
  integer              :: lanc_ncv_factor     !Set the size of the block used in Lanczos-Arpack by multiplying the required Neigen (Ncv=lanc_ncv_factor*Neigen+lanc_ncv_add)
  integer              :: lanc_ncv_add        !Adds up to the size of the block to prevent it to become too small (Ncv=lanc_ncv_factor*Neigen+lanc_ncv_add)
  integer              :: lanc_nstates_sector !Max number of required eigenvalues per sector
  integer              :: lanc_nstates_total  !Max number of states hold in the finite T calculation
  integer              :: lanc_nstates_step   !Number of states added at each step to determine the optimal spectrum size at finite T
  integer              :: lanc_dim_threshold  !Min dimension threshold to use Lanczos determination of the spectrum rather than Lapack based exact diagonalization.
  !
  integer              :: cg_Niter            !Max number of iteration in the fit
  real(8)              :: cg_Ftol             !Tolerance in the cg fit
  integer              :: cg_Weight           !CGfit mode 0=normal,1=1/n weight, 2=1/w weight
  character(len=5)     :: cg_Scheme           !fit scheme: delta (default), weiss for G0^
  integer              :: cg_method           !fit routine type:0=CGnr (default), 1=minimize (old f77), 2=CG+
  integer              :: cg_stop             !fit stop condition:0-3, 0=default
  integer              :: cg_grad             !gradient evaluation: 0=analytic, 1=numeric
  real(8)              :: cg_eps              !fit eps tolerance
  !
  logical              :: finiteT             !flag for finite temperature calculation
  character(len=7)     :: bath_type           !flag to set bath type: normal (1bath/imp), hybrid(1bath)
  real(8)              :: nread               !fixed density. if 0.d0 fixed chemical potential calculation.
  real(8)              :: nerr                !fix density threshold. a loop over from 1.d-1 to required nerr is performed
  real(8)              :: ndelta              !initial chemical potential step
  real(8)              :: ncoeff              !multiplier for the initial ndelta read from a file (ndelta-->ndelta*ncoeff)
  integer              :: niter               !
  logical              :: Jz_basis            !switch to spherical harmonics basis (Interaction need to be rotated)
  logical              :: Jz_max              !limit the number of sectorn spanned by the solver (like ed_sectors_shift)
  real(8)              :: Jz_max_value        !maximum Jz reached by sector scan
  !
  logical              :: replica_operators   !flag to set the bath representation via user provided fixed operators
  integer              :: Nop                 !number of user provided fixed operators
  !
  logical              :: Utensor             !flag to use tensorial representation of the interaction
  !
  logical              :: plaquette           !flag to perform ED calculation on the lattice
  logical              :: HardCoreBoson       !flag to switch to a bosonic model
  integer              :: filling             !filling of the lattice
  real(8),dimension(5) :: Vnn                 !non-local interactions
  real(8)              :: Thopping            !hopping integral
  logical              :: ladder              !flag to remove PBC in the lattice y direction


  !Some parameters for function dimension:
  !=========================================================
  integer              :: Lmats
  integer              :: Lreal
  integer              :: Lfit
  integer              :: Ltau

  !LOG AND Hamiltonian UNITS
  !=========================================================
  character(len=100)   :: Hfile,HLOCfile,UTENSfile
  integer,save         :: LOGfile



contains


  !+-------------------------------------------------------------------+
  !PURPOSE  : READ THE INPUT FILE AND SETUP GLOBAL VARIABLES
  !+-------------------------------------------------------------------+
  subroutine ed_read_input(INPUTunit,comm)
#ifdef _MPI
    USE MPI
    USE SF_MPI
#endif
    character(len=*) :: INPUTunit
    integer,optional :: comm
    logical          :: master=.true.
    integer          :: i,rank=0
#ifdef _MPI
    if(present(comm))then
       master=get_Master_MPI(comm)
       rank  =get_Rank_MPI(comm)
    endif
#endif


    !DEFAULT VALUES OF THE PARAMETERS:
    call parse_input_variable(Norb,"NORB",INPUTunit,default=1,comment="Number of impurity orbitals.")
    call parse_input_variable(Nbath,"NBATH",INPUTunit,default=6,comment="Number of bath sites:(normal=>Nbath per orb)(hybrid=>Nbath total)(replica=>Nbath=Nreplica)")
    call parse_input_variable(Nspin,"NSPIN",INPUTunit,default=1,comment="Number of spin degeneracy (max 2)")
    call parse_input_variable(uloc,"ULOC",INPUTunit,default=[2.d0,0.d0,0.d0],comment="Values of the local interaction per orbital (max 3)")
    call parse_input_variable(ust,"UST",INPUTunit,default=0.d0,comment="Value of the inter-orbital interaction term")
    call parse_input_variable(Jh,"JH",INPUTunit,default=0.d0,comment="Hunds coupling")
    call parse_input_variable(Jx,"JX",INPUTunit,default=0.d0,comment="S-E coupling")
    call parse_input_variable(Jp,"JP",INPUTunit,default=0.d0,comment="P-H coupling")
    call parse_input_variable(beta,"BETA",INPUTunit,default=1000.d0,comment="Inverse temperature, at T=0 is used as a IR cut-off.")
    call parse_input_variable(xmu,"XMU",INPUTunit,default=0.d0,comment="Chemical potential. If HFMODE=T, xmu=0 indicates half-filling condition.")
    call parse_input_variable(deltasc,"DELTASC",INPUTunit,default=0.02d0,comment="Value of the SC symmetry breaking term.")
    call parse_input_variable(nloop,"NLOOP",INPUTunit,default=100,comment="Max number of DMFT iterations.")
    call parse_input_variable(dmft_error,"DMFT_ERROR",INPUTunit,default=0.00001d0,comment="Error threshold for DMFT convergence")
    call parse_input_variable(sb_field,"SB_FIELD",INPUTunit,default=0.1d0,comment="Value of a symmetry breaking field for magnetic solutions.")
    !
    call parse_input_variable(ed_twin,"ED_TWIN",INPUTunit,default=.false.,comment="flag to reduce (T) or not (F,default) the number of visited sector using twin symmetry.")
    call parse_input_variable(ed_sectors,"ED_SECTORS",INPUTunit,default=.false.,comment="flag to reduce sector scan for the spectrum to specific sectors +/- ed_sectors_shift.")
    call parse_input_variable(ed_sectors_shift,"ED_SECTORS_SHIFT",INPUTunit,1,comment="shift to ed_sectors")
    call parse_input_variable(ed_sparse_H,"ED_SPARSE_H",INPUTunit,default=.true.,comment="flag to select  storage of sparse matrix H (mem--, cpu++) if TRUE, or direct on-the-fly H*v product (mem++, cpu--) if FALSE ")
    call parse_input_variable(ed_solve_offdiag_gf,"ED_SOLVE_OFFDIAG_GF",INPUTunit,default=.false.,comment="flag to select the calculation of the off-diagonal impurity GF. this is T by default if bath_type/=normal")
    call parse_input_variable(ed_print_Sigma,"ED_PRINT_SIGMA",INPUTunit,default=.true.,comment="flag to print impurity Self-energies")
    call parse_input_variable(ed_print_G,"ED_PRINT_G",INPUTunit,default=.true.,comment="flag to print impurity Greens function")
    call parse_input_variable(ed_print_G0,"ED_PRINT_G0",INPUTunit,default=.true.,comment="flag to print non-interacting impurity Greens function")
    !
    call parse_input_variable(nsuccess,"NSUCCESS",INPUTunit,default=1,comment="Number of successive iterations below threshold for convergence")
    call parse_input_variable(Lmats,"LMATS",INPUTunit,default=5000,comment="Number of Matsubara frequencies.")
    call parse_input_variable(Lreal,"LREAL",INPUTunit,default=5000,comment="Number of real-axis frequencies.")
    call parse_input_variable(Ltau,"LTAU",INPUTunit,default=1000,comment="Number of imaginary time points.")
    call parse_input_variable(Lfit,"LFIT",INPUTunit,default=1000,comment="Number of Matsubara frequencies used in the \Chi2 fit.")
    call parse_input_variable(nread,"NREAD",INPUTunit,default=0.d0,comment="Objective density for fixed density calculations.")
    call parse_input_variable(nerr,"NERR",INPUTunit,default=1.d-4,comment="Error threshold for fixed density calculations.")
    call parse_input_variable(ndelta,"NDELTA",INPUTunit,default=0.1d0,comment="Initial step for fixed density calculations.")
    call parse_input_variable(ncoeff,"NCOEFF",INPUTunit,default=1d0,comment="multiplier for the initial ndelta read from a file (ndelta-->ndelta*ncoeff).")
    call parse_input_variable(wini,"WINI",INPUTunit,default=-5.d0,comment="Smallest real-axis frequency")
    call parse_input_variable(wfin,"WFIN",INPUTunit,default=5.d0,comment="Largest real-axis frequency")
    call parse_input_variable(chispin_flag,"CHISPIN_FLAG",INPUTunit,default=.false.,comment="Flag to activate spin susceptibility calculation.")
    call parse_input_variable(chidens_flag,"CHIDENS_FLAG",INPUTunit,default=.false.,comment="Flag to activate density susceptibility calculation.")
    call parse_input_variable(chipair_flag,"CHIPAIR_FLAG",INPUTunit,default=.false.,comment="Flag to activate pair susceptibility calculation.")

    call parse_input_variable(hfmode,"HFMODE",INPUTunit,default=.true.,comment="Flag to set the Hartree form of the interaction (n-1/2). see xmu.")
    call parse_input_variable(eps,"EPS",INPUTunit,default=0.01d0,comment="Broadening on the real-axis.")
    call parse_input_variable(cutoff,"CUTOFF",INPUTunit,default=1.d-9,comment="Spectrum cut-off, used to determine the number states to be retained.")
    call parse_input_variable(gs_threshold,"GS_THRESHOLD",INPUTunit,default=1.d-9,comment="Energy threshold for ground state degeneracy loop up")
    call parse_input_variable(hwband,"HWBAND",INPUTunit,default=2d0,comment="half-bandwidth for the bath initialization: flat in -hwband:hwband")
    !
    call parse_input_variable(lanc_method,"LANC_METHOD",INPUTunit,default="arpack",comment="select the lanczos method to be used in the determination of the spectrum. ARPACK (default), LANCZOS (T=0 only)")
    call parse_input_variable(lanc_nstates_sector,"LANC_NSTATES_SECTOR",INPUTunit,default=6,comment="Initial number of states per sector to be determined.")
    call parse_input_variable(lanc_nstates_total,"LANC_NSTATES_TOTAL",INPUTunit,default=1,comment="Initial number of total states to be determined.")
    call parse_input_variable(lanc_nstates_step,"LANC_NSTATES_STEP",INPUTunit,default=2,comment="Number of states added to the spectrum at each step.")
    call parse_input_variable(lanc_ncv_factor,"LANC_NCV_FACTOR",INPUTunit,default=3,comment="Set the size of the block used in Lanczos-Arpack by multiplying the required Neigen (Ncv=lanc_ncv_factor*Neigen+lanc_ncv_add)")
    call parse_input_variable(lanc_ncv_add,"LANC_NCV_ADD",INPUTunit,default=5,comment="Adds up to the size of the block to prevent it to become too small (Ncv=lanc_ncv_factor*Neigen+lanc_ncv_add)")
    call parse_input_variable(lanc_niter,"LANC_NITER",INPUTunit,default=512,comment="Number of Lanczos iteration in spectrum determination.")
    call parse_input_variable(lanc_ngfiter,"LANC_NGFITER",INPUTunit,default=200,comment="Number of Lanczos iteration in GF determination. Number of momenta.")
    call parse_input_variable(lanc_tolerance,"LANC_TOLERANCE",INPUTunit,default=0.000000000001d0,comment="Tolerance for the Lanczos iterations as used in Arpack and plain lanczos.")
    call parse_input_variable(lanc_dim_threshold,"LANC_DIM_THRESHOLD",INPUTunit,default=256,comment="Min dimension threshold to use Lanczos determination of the spectrum rather than Lapack based exact diagonalization.")
    !
    call parse_input_variable(cg_niter,"CG_NITER",INPUTunit,default=500,comment="Max. number of Conjugate-Gradient iterations.")
    call parse_input_variable(cg_scheme,"CG_SCHEME",INPUTunit,default='weiss',comment="Conjugate-Gradient fit scheme: delta or weiss.")
    call parse_input_variable(cg_ftol,"CG_FTOL",INPUTunit,default=0.00001d0,comment="Conjugate-Gradient tolerance.")
    call parse_input_variable(cg_method,"CG_METHOD",INPUTunit,default=0,comment="Conjugate-Gradient method: 0=NR, 1=minimize, 2=CG+.")
    call parse_input_variable(cg_stop,"CG_STOP",INPUTunit,default=0,comment="Conjugate-Gradient stopping condition.")
    call parse_input_variable(cg_eps,"CG_EPS",INPUTunit,default=0.000001d0,comment="Conjugate-Gradient eps tolerance.")
    call parse_input_variable(cg_grad,"CG_GRAD",INPUTunit,default=0,comment="Gradient evaluation method: 0=analytic (default), 1=numeric.")
    call parse_input_variable(cg_weight,"CG_WEIGHT",INPUTunit,default=0,comment="Conjugate-Gradient weight form: 0=1.0 ,1=1/n , 2=1/w.")
    call parse_input_variable(ed_mode,"ED_MODE",INPUTunit,default='normal',comment="Flag to set ED type: normal=normal, superc=superconductive, nonsu2=broken SU(2)")
    call parse_input_variable(ed_para,"ED_PARA",INPUTunit,default=.false.,comment="Flag to force paramagnetic solution (only used in ed_mode=nonsu2 now).")
    call parse_input_variable(ed_vsf_ratio,"ED_VSF_RATIO",INPUTunit,default=0.1d0,comment="Ration of the spin-flip hopping to spin-hold ones in the nonSU2 channel.")
    call parse_input_variable(ed_bath_noise_thr,"ED_BATH_NOISE_THR",INPUTunit,default=0.d0,comment="Noise added to the impurity hybridization")
    !
    call parse_input_variable(bath_type,"BATH_TYPE",INPUTunit,default='normal',comment="flag to set bath type: normal (1bath/imp), hybrid(1bath), replica(1replica/imp)")
    call parse_input_variable(Hfile,"HFILE",INPUTunit,default="hamiltonian",comment="File where to retrieve/store the bath parameters.")
    call parse_input_variable(HLOCfile,"impHfile",INPUTunit,default="inputHLOC.in",comment="File read the input local H.")
    call parse_input_variable(LOGfile,"LOGFILE",INPUTunit,default=6,comment="LOG unit.")
    call parse_input_variable(Jz_basis,"JZ_BASIS",INPUTunit,default=.false.,comment="switch to spherical harmonics basis (Interaction need to be rotated)")
    call parse_input_variable(Jz_max,"JZ_MAX",INPUTunit,default=.false.,comment="limit the number of sectorn spanned by the solver (like ed_sectors_shift)")
    call parse_input_variable(Jz_max_value,"JZ_MAX_VALUE",INPUTunit,default=-1000.d0,comment="if (Jz_basis.and.Jz_max):maximum Jz reached by sector scan elseif(Jz_basis.and..not.Jz_max): minimum Ntot from where sector scan starts")
    call parse_input_variable(ed_verbose,"ED_VERBOSE",INPUTunit,default=3,comment="Verbosity level: 0=almost nothing --> 3:all.")
    !
    call parse_input_variable(replica_operators,"REPLICA_OPERATORS",INPUTunit,default=.false.,comment="Flag to set the bath representation via user provided fixed operators")
    call parse_input_variable(Nop,"NOP",INPUTunit,default=0,comment="Number of user provided fixed operators")
    !
    call parse_input_variable(Utensor,"UTENSOR",INPUTunit,default=.false.,comment="Flag to use tensorial representation of the interaction")
    call parse_input_variable(UTENSfile,"UTENSfile",INPUTunit,default="Utensor.in",comment="File read the input tensorial representation of the interaction")
    !
    call parse_input_variable(plaquette,"PLAQUETTE",INPUTunit,default=.false.,comment="Flag to perform ED calculation on the lattice")
    call parse_input_variable(HardCoreBoson,"HARDCOREBOSON",INPUTunit,default=.false.,comment="Flag to switch to a bosonic model")
    call parse_input_variable(filling,"FILLING",INPUTunit,default=4,comment="Filling of the lattice")
    call parse_input_variable(Vnn,"VNN",INPUTunit,default=[0.d0,0.d0,0.d0,0.d0,0.d0,0.d0,0.d0],comment="Values of the non-local interaction per site (max 7)")
    call parse_input_variable(Thopping,"THOPPING",INPUTunit,default=0.d0,comment="Hopping integral")
    !
#ifdef _MPI
    if(present(comm))then
       if(.not.master)then
          LOGfile=9999-rank
          open(LOGfile,file="stdOUT.rank"//str(rank)//".ed")
          do i=1,get_Size_MPI(comm)
             if(i==rank)write(*,"(A,I0,A,I0)")"Rank ",rank," writing to unit: ",LOGfile
          enddo
       endif
    endif
#endif
    !
    !
    Ltau=max(int(beta),Ltau)
    if(master)then
       call save_input_file(INPUTunit)
       call scifor_version()
       call code_version(revision)
    endif
    !Act on the input variable only after printing.
    !In the new parser variables are hard-linked into the list:
    !any change to the variable is immediately copied into the list... (if you delete .ed it won't be printed out)
    call substring_delete(Hfile,".restart")
    call substring_delete(Hfile,".ed")
  end subroutine ed_read_input




  subroutine substring_delete (s,sub)
    !! S_S_DELETE2 recursively removes a substring from a string.
    !    The remainder is left justified and padded with blanks.
    !    The substitution is recursive, so
    !    that, for example, removing all occurrences of "ab" from
    !    "aaaaabbbbbQ" results in "Q".
    !  Parameters:
    !    Input/output, character ( len = * ) S, the string to be transformed.
    !    Input, character ( len = * ) SUB, the substring to be removed.
    !    Output, integer ( kind = 4 ) IREP, the number of occurrences of
    !    the substring.
    integer          :: ihi
    integer          :: irep
    integer          :: loc
    integer          :: nsub
    character(len=*) ::  s
    integer          :: s_length
    character(len=*) :: sub
    s_length = len ( s )
    nsub = len ( sub )
    irep = 0
    ihi = s_length
    do while ( 0 < ihi )
       loc = index ( s(1:ihi), sub )
       if ( loc == 0 ) then
          return
       end if
       irep = irep + 1
       call s_chop ( s, loc, loc+nsub-1 )
       ihi = ihi - nsub
    end do
    return
  end subroutine substring_delete

  subroutine s_chop ( s, ilo, ihi )
    !! S_CHOP "chops out" a portion of a string, and closes up the hole.
    !  Example:
    !    S = 'Fred is not a jerk!'
    !    call s_chop ( S, 9, 12 )
    !    S = 'Fred is a jerk!    '
    !  Parameters:
    !    Input/output, character ( len = * ) S, the string to be transformed.
    !    Input, integer ( kind = 4 ) ILO, IHI, the locations of the first and last
    !    characters to be removed.
    integer               ::ihi
    integer               ::ihi2
    integer               ::ilo
    integer               ::ilo2
    character ( len = * ) :: s
    integer               ::s_length
    s_length = len ( s )
    ilo2 = max ( ilo, 1 )
    ihi2 = min ( ihi, s_length )
    if ( ihi2 < ilo2 ) then
       return
    end if
    s(ilo2:s_length+ilo2-ihi2-1) = s(ihi2+1:s_length)
    s(s_length+ilo2-ihi2:s_length) = ' '
    return
  end subroutine s_chop


END MODULE ED_INPUT_VARS
