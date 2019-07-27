      program driver_hilb
       !-------------------------------------------
       !  Driver to test Hilbert transform modules
       !
       !  Author: Jonathas Maciel
       !  e-mail: jonathassilvamaciel@gmail.com
       !-------------------------------------------
       implicit none
       integer(4),parameter :: nt=1000
       integer(4),parameter :: FID=10
       real(4)              :: ttotal
       real(4)              :: f0=15
       real(4)              :: t0=0.5
       real(4)              :: dt=0.001
       real(4)              :: phase=45
       real(4)              :: t(nt)
       real(4)              :: h(nt)
       character(len=100)   :: filename="hilb.dat"
       integer(4)           :: ierr,it,reclen


       call Ricker(nt,dt,t0,f0,t)
       call Hilbert_Transform(nt,t,h)
       !call FHT(nt,phase,t,h)

       reclen = 2*nt*4
       open(unit=FID,             &
            file=trim(filename),  &
            form='unformatted',   &
            access='direct',      &
            recl=reclen,          &
            status='replace',     &
            action='write',       &
            iostat=ierr)
       if(ierr/=0) STOP 'ERROR: File hilb.dat can not open'
       write(unit=FID,rec=1) t(1:nt),h(1:nt)
       close(unit=FID)



       call system("xgraph nplot=2 linecolor=0,4 linewidth=2 &
                   windowtitle='Hilbert T. project' &
                   label1='Time [s]' label2='Amplitude' &
                   n=1000 f1=0.0 d1=0.001 style=normal < hilb.dat")


 
      end program driver_hilb




      subroutine Ricker(nt,dt,t0,f0,t)
       !----------------------------------------
       ! Ricker pulse
       !----------------------------------------
       implicit none
       integer(4),intent(in) :: nt
       real(4),intent(in)    :: dt
       real(4),intent(in)    :: t0
       real(4),intent(in)    :: f0
       real(4),intent(out)   :: t(nt)
       real(4),parameter     :: pi=3.14159265359
       real(4)               :: arg
       integer(4)            :: it

       do it = 1,nt        
          arg = (pi*f0*(real(it-1)*dt - t0))**2
          t(it) = (1.0-2.0*arg)*exp(-arg)
       end do
 
      end subroutine Ricker



      subroutine Hilbert_Transform(nt,t,h)
       !----------------------------------------
       ! Hilbert Transform 
       ! 
       ! Author: Jonathas Maciel
       ! E-mail: jonathassilvamaciel@gmail.com
       !----------------------------------------
       implicit none
       integer(4),intent(in) :: nt
       real(4),intent(in)    :: t(nt)
       real(4),intent(out)   :: h(nt)
       real(4)               :: buffer(-nt:nt)
       real(4),parameter     :: pi=3.14159265359
       real(4)               :: soma
       integer(4)            :: it,i

       buffer(-nt:0)=0.0
       buffer(1:nt)=t(1:nt)
       do it = 1,nt       
          soma = 0.0
          do i=-nt,nt
             if (it == i) cycle
             soma = soma + buffer(i)/real(i-it)
          end do
          h(it) = soma/pi
       end do
 
      end subroutine Hilbert_Transform



      subroutine FHT(nt,phase,t,fh)
       !----------------------------------------
       ! Fractional Hilbert Transform (FHT)
       ! 
       ! Author: Jonathas Maciel
       ! E-mail: jonathassilvamaciel@gmail.com
       !----------------------------------------
       implicit none
       integer(4),intent(in) :: nt
       real(4),intent(in)    :: phase
       real(4),intent(in)    :: t(nt)
       real(4),intent(out)   :: fh(nt)
       real(4)               :: buffer(-nt:nt)
       real(4)               :: h(1:nt)
       real(4),parameter     :: pi=3.14159265359
       real(4)               :: soma
       real(4)               :: rad
       integer(4)            :: it,i

       rad = phase*pi/180.0

       buffer(-nt:0)=0.0
       buffer(1:nt)=t(1:nt)
       do it = 1,nt       
          soma = 0.0
          do i=-nt,nt
             if (it == i) cycle
             soma = soma + buffer(i)/real(i-it)
          end do
          h(it) = soma/pi
       end do

       do it=1,nt
          fh(it) = cos(rad)*t(it) + sin(rad)*h(it)
       end do

      end subroutine FHT


