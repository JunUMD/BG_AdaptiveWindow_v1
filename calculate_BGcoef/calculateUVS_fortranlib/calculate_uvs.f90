module  calculateuvs_single_mod

!****************************************************************************************
!...BG enhancement algorithm
!...reference:
!...Author   date       ver     modify
!Tiger@UMD, June,2012   2.0    Modified from  original 
!                              version of NSMC created by Tiger in 2007
!JunZhou@UMD Oct,2020   3.0    (1) Modified the integration of U
!                              (2) Adapt the code for adaptive window
!--------------------------------------------------------------------------
!----COPYRIGHT
!
!-------------------------------------------------------------------------
!--- Reference:
!--- A.Stogryn,1978,"Estimates of brightness temperatures from scanning radiometer data"
!--- , IEEE Transactions on antennas and propagation, Vol.AP-26, No.5, PP.720-726 

!****************************************************************************************
contains

subroutine calculateuvs_single(save_file,tgt_grd_dir,halfwin,atm_grd_dir,window_fname,ifr, &
                              numxgrd,numygrd)
implicit none
integer,intent(in)::ifr,numxgrd,numygrd,halfwin
character(len=*),intent(in)::save_file,tgt_grd_dir,atm_grd_dir,window_fname

character(len=2)::str_ifr,str_isc,str_ifv
character(len=200)::file_source_antgrid,file_target_antgrid
integer(kind=2),parameter :: fid_win=40,fid_antsrc=10,fid_anttgt=20,fid_uvs=30
integer(kind=2),parameter :: nparam=6
integer(kind=2) i,idim,jdim,ipix,jpix,isc,ifv,idx_src
integer(kind=2) num,idim_source,jdim_source,idim_target,jdim_target
integer(kind=2) idim_source_ix,jdim_source_iy,idim_source_jx,jdim_source_jy
integer(kind=2) min_nx_target,min_ny_target,min_nx_source,min_ny_source
integer(kind=2) min_nx_ix,max_nx_ix,min_ny_iy,max_ny_iy,min_nx_jx
integer(kind=2) max_nx_jx,min_ny_jy,max_ny_jy,max_nx_target
integer(kind=2) max_ny_target,max_nx_source,max_ny_source
integer(kind=4),allocatable :: antpos_target(:)
integer(kind=2),allocatable :: reshp_antpos_source(:,:) 
integer(kind=4), allocatable  :: temppos_source(:)  
real(kind=8) zsum
real(kind=8),allocatable,dimension(:,:) :: u,v   
real(kind=8),allocatable,dimension(:,:) :: s             
integer(kind=2),allocatable,dimension(:,:)     :: antpos_target_nx,antpos_target_ny   
integer(kind=2),allocatable,dimension(:,:,:)   :: reshp_antpos_source_nx,reshp_antpos_source_ny 
real(kind=8)   ,allocatable,dimension(:,:)     :: antwgt_target,antwgt_target_ave 
real(kind=8)   ,allocatable,dimension(:,:,:)   :: reshp_antwgt_source,reshp_antwgt_source_ave 
real(kind=8)   ,allocatable,dimension(:,:) :: tempwgt_source,tempwgt_source_ave

integer(kind=4)::windowsize
integer(kind=4),allocatable,dimension(:,:)::windowindx

integer(kind=4) :: foot_nx_target,foot_ny_target,foot_nx_source,foot_ny_source


open(fid_win,file = trim(window_fname) ,form='unformatted',access='stream',convert='big_endian',status="old")
read(fid_win) windowsize
allocate(windowindx(windowsize,2))
read(fid_win) windowindx
read(fid_win) idx_src
close(fid_win)

allocate(antpos_target(6))
allocate(antwgt_target(numxgrd,numygrd))
allocate(antwgt_target_ave(numxgrd,numygrd))
allocate(antpos_target_nx(numxgrd,numygrd))
allocate(antpos_target_ny(numxgrd,numygrd))


allocate(temppos_source(6))
allocate(tempwgt_source(numxgrd,numygrd))
allocate(tempwgt_source_ave(numxgrd,numygrd))
allocate(reshp_antpos_source(windowsize,6))
allocate(reshp_antwgt_source(windowsize,numxgrd,numygrd))
allocate(reshp_antwgt_source_ave(windowsize,numxgrd,numygrd))
allocate(reshp_antpos_source_nx(windowsize,numxgrd,numygrd))
allocate(reshp_antpos_source_ny(windowsize,numxgrd,numygrd))
allocate(u(windowsize,1))
allocate(v(windowsize,1))
allocate(s(windowsize,windowsize))

   antpos_target_nx       = 0
   antpos_target_ny       = 0
   antwgt_target          = 0
   antwgt_target_ave      = 0
   antpos_target          = 0
   reshp_antpos_source    = 0
   reshp_antpos_source_nx = 0
   reshp_antpos_source_ny = 0
   reshp_antwgt_source    = 0
   reshp_antwgt_source_ave= 0
   u = 0
   v = 0
   s = 0

   write(str_ifr,"(I2)") ifr
   write(str_isc,"(I2)") halfwin
   file_target_antgrid=trim(tgt_grd_dir)//'isc'//str_isc//'_ifv'//trim(adjustl(str_ifr))//'.bin'
   open(fid_anttgt,file = trim(file_target_antgrid) ,form='unformatted',access='stream',convert='big_endian',status="old")
   read(fid_anttgt) antpos_target  
   read(fid_anttgt) antwgt_target  
   read(fid_anttgt) antwgt_target_ave  
   close(fid_anttgt)

   min_nx_target = antpos_target(1)  
   max_nx_target = antpos_target(2) 
   min_ny_target = antpos_target(3) 
   max_ny_target = antpos_target(4) 
   foot_nx_target = antpos_target(5) 
   foot_ny_target = antpos_target(6)
   do idim = 1,numxgrd
     do jdim = 1,numygrd
        antpos_target_nx(idim,jdim) = idim + min_nx_target-1   
        antpos_target_ny(idim,jdim) = jdim + min_ny_target-1
     end do
   end do
 

   do i = 1,windowsize
         isc=windowindx(i,1) 
         ifv=windowindx(i,2)


         write(str_isc,"(I2)") isc 
         write(str_ifv,"(I2)") ifv 
         file_source_antgrid=trim(atm_grd_dir)//&
                             'isc'//trim(adjustl(str_isc))//'_ifv'//trim(adjustl(str_ifv))//'.bin'      
         open(fid_antsrc,file = trim(file_source_antgrid),form='unformatted',access='stream',&
                                convert='big_endian',status="old")
         read(fid_antsrc) temppos_source
         read(fid_antsrc) tempwgt_source
         read(fid_antsrc) tempwgt_source_ave
         close(fid_antsrc)
          
         reshp_antpos_source(i,:)   = temppos_source
         reshp_antwgt_source(i,:,:) = tempwgt_source
         reshp_antwgt_source_ave(i,:,:) = tempwgt_source_ave
         foot_nx_source=temppos_source(5)
         foot_ny_source=temppos_source(6)
         
   end do
   
   do i = 1,windowsize
         min_nx_source  = reshp_antpos_source(i,1) 
         min_ny_source  = reshp_antpos_source(i,3)
         max_nx_source  = reshp_antpos_source(i,2)
         max_ny_source  = reshp_antpos_source(i,4)
         do idim = 1,numxgrd
           do jdim = 1,numygrd
              reshp_antpos_source_nx(i,idim,jdim) = idim + min_nx_source-1  
              reshp_antpos_source_ny(i,idim,jdim) = jdim + min_ny_source-1
           end do
         end do
   end do 
  
 
   num = 0 
   do ipix = 1,windowsize

      min_nx_ix = reshp_antpos_source(ipix,1) 
      max_nx_ix = reshp_antpos_source(ipix,2)
      min_ny_iy = reshp_antpos_source(ipix,3)
      max_ny_iy = reshp_antpos_source(ipix,4)
      do idim_source = 1,numxgrd
        do jdim_source = 1,numygrd

          if( reshp_antwgt_source(ipix,idim_source,jdim_source) .ne. 0) then
             U(ipix,1) = U(ipix,1) + reshp_antwgt_source(ipix,idim_source,jdim_source)
          endif

          if(   reshp_antpos_source_nx(ipix,idim_source,jdim_source) .le. max_nx_target &
          .and. reshp_antpos_source_nx(ipix,idim_source,jdim_source) .ge. min_nx_target &
          .and. reshp_antpos_source_ny(ipix,idim_source,jdim_source) .le. max_ny_target &
          .and. reshp_antpos_source_ny(ipix,idim_source,jdim_source) .ge. min_ny_target ) then

             idim_target = reshp_antpos_source_nx(ipix,idim_source,jdim_source) - min_nx_target + 1  
             jdim_target = reshp_antpos_source_ny(ipix,idim_source,jdim_source) - min_ny_target + 1
             if(   idim_target .gt. 0 .and. idim_target .le. numxgrd &
             .and. jdim_target .gt. 0 .and. jdim_target .le. numygrd ) then
             
                if( reshp_antwgt_source(ipix,idim_source,jdim_source) .ne. 0       &
                     .and. antwgt_target(idim_target,jdim_target) .ne. 0)then

                   V(ipix,1) = V(ipix,1) + reshp_antwgt_source(ipix,idim_source,jdim_source)* &
                                           antwgt_target(idim_target,jdim_target) 
                endif
             endif
          endif

       end do 
      end do

    
      do jpix = 1,windowsize        
         min_nx_jx = reshp_antpos_source(jpix,1)
         max_nx_jx = reshp_antpos_source(jpix,2)
         min_ny_jy = reshp_antpos_source(jpix,3)
         max_ny_jy = reshp_antpos_source(jpix,4)
        
         Zsum = 0
         do idim_source_ix = 1,numxgrd
            do jdim_source_iy = 1,numygrd
                  idim_source_jx = reshp_antpos_source_nx(ipix,idim_source_ix,jdim_source_iy) - min_nx_jx + 1
                  jdim_source_jy = reshp_antpos_source_ny(ipix,idim_source_ix,jdim_source_iy) - min_ny_jy + 1
                  if(   idim_source_jx .ge. 1 .and. idim_source_jx .le. numxgrd &
                  .and. jdim_source_jy .ge. 1 .and. jdim_source_jy .le. numygrd) then 
                     if(      reshp_antwgt_source(jpix,idim_source_jx,jdim_source_jy) .ne. 0 &
                        .and. reshp_antwgt_source(ipix,idim_source_ix,jdim_source_iy) .ne. 0) then
                       
                           zsum = zsum + reshp_antwgt_source(ipix,idim_source_ix,jdim_source_iy) &
                                        *reshp_antwgt_source(jpix,idim_source_jx,jdim_source_jy) 
                     end if
                  end if
            end do 
         end do 
         s(ipix,jpix)=zsum 
         s(jpix,ipix) = s(ipix,jpix)
      end do 
   end do 

   open(fid_uvs,file = trim(save_file),form='unformatted',access='stream',convert='big_endian',status="replace") 
   write(fid_uvs) u
   write(fid_uvs) v
   write(fid_uvs) s
   close(fid_uvs)
deallocate(antpos_target_nx)
deallocate(antpos_target_ny)
deallocate(antwgt_target)
deallocate(antpos_target)
deallocate(reshp_antpos_source)
deallocate(reshp_antwgt_source)
deallocate(reshp_antwgt_source_ave)
deallocate(reshp_antpos_source_nx)
deallocate(reshp_antpos_source_ny)
deallocate(windowindx)
deallocate(tempwgt_source)
deallocate(tempwgt_source_ave)
deallocate(u,v,s)

return 
end subroutine calculateuvs_single
end module  calculateuvs_single_mod
