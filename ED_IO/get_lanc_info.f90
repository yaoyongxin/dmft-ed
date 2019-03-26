subroutine ed_get_neigen_total(nlii,Nlat) 
  integer                      :: Nlat
  integer,dimension(Nlat) :: nlii
  nlii=0d0
  if(allocated(neigen_totalii))then
     if(Nlat>size(neigen_totalii)) stop "ed_get_docc error: required N_sites > evaluated N_sites"
     nlii=neigen_totalii
  endif
end subroutine ed_get_neigen_total
