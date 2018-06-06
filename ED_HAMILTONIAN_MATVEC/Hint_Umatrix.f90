  htmp = zero
  do iflav=1,Nflavors
     do jflav=1,Nflavors
        htmp = htmp + nvec(iflav)*Umatrix(iflav,jflav)*nvec(jflav)
     enddo
  enddo

  !if using the Hartree-shifted chemical potential: mu=0 for half-filling
  !sum up the contributions of hartree terms:
  if(hfmode)then
     do iflav=1,Nflavors
        do jflav=1,Nflavors
           htmp = htmp - 0.5d0*Umatrix(iflav,jflav)*(nvec(iflav)+nvec(jflav))
        enddo
     enddo
  endif
  !
  call sp_insert_element(spH0,htmp,impi,i)
  !
