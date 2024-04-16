module solver
!$ use omp_lib
implicit none
contains

subroutine  solve_matrix_x (p, ap, ae, aw, an, as, bb, m, n, relux_factor, iter_max)
        use global_2d
        implicit none
        real,intent(inout),dimension(0:md,0:nd)::p
        real,intent(in),dimension(0:md,0:nd)::ap, ae, aw, an, as, bb
        integer,intent(in)::m, n
        real,intent(in)::relux_factor
        integer,intent(in)::iter_max 
        ! local variables
        real::error, correct
        real,dimension(0:md/2+1,0:nd/2+1)::f00,f01,f10,f11
        real,dimension(0:md/2+1,0:nd/2+1)::ap00,ap01,ap10,ap11, &
                                                ae00,ae01,ae10,ae11, &
                                                aw00,aw01,aw10,aw11, &
                                                an00,an01,an10,an11, &
                                                as00,as01,as10,as11, &
                                                bb00,bb01,bb10,bb11
        integer::i, j, iter
        integer::i_o, j_o, i_e, j_e, ii, jj, mm, nn
        integer::i0,i1,j0,j1
        !check counter
        integer,dimension(0:md,0:nd):: icheck

        ! ----------------
        !   SOR algorithm
        ! ----------------

        !convert local matrix  p to f00 f01 f10 f11
        i_e= mod(m,2)  ! 0 for m=even / 1 for m=odd
        j_e= mod(n,2)  ! 0 for n=even / 1 for n=odd
        i_o= 1-i_e     ! 0 for m=odd / 1 for m=even  
        j_o= 1-j_e     ! 0 for n=odd / 1 for n=even
        mm =(m+i_e)/2 ! p(1,*)=f*1(*,1) for m=odd / p(1,*)=f*0(1,*) for m=even
        nn =(n+j_e)/2 ! p(*,1)=f1*(*,1) for n=odd / p(1,*)=f0*(*,1) for n=even    
                                ! (0,0) <= f**(jj,ii) <= (nn+1,mm+1)

        !check write #0
        ! write(*,*) "#0 m,n= ",m,n," =>   mm,nn= ",mm,nn
        ! write(*,*) "i_e, i_o, j_e, j_o =",i_e, i_o, j_e, j_o
        ! write(*,*) ((ii-i_e)*2-i_o, (ii-i_e)*2-i_o+1,ii=1,mm)
        ! write(*,*) ((jj-j_e)*2-j_o, (jj-j_e)*2-j_o+1,jj=1,nn)
        !check write #0

        ! check counter set #1
        do i=1,m
        do j=1,n
        icheck(i,j)=0
        end do
        end do

        ! convert data  p => f00,f01,f10,f11
        do ii = 1, mm
        do jj = 1, nn
        i0= (ii-1)*2 +i_o
        i1= (ii-1)*2 +i_o+1
        j0= (jj-1)*2 +j_o
        j1= (jj-1)*2 +j_o+1

        f00(jj,ii)=p(i0,j0)
        f01(jj,ii)=p(i1,j0)
        f10(jj,ii)=p(i0,j1)
        f11(jj,ii)=p(i1,j1)
        end do
        end do

        ! convert matrix components  ap,ae,aw,an,an,bb
        do ii = 1, mm
        do jj = 1, nn
        i0= (ii-1)*2 +i_o
        i1= (ii-1)*2 +i_o+1
        j0= (jj-1)*2 +j_o
        j1= (jj-1)*2 +j_o+1
        ! check counter
        icheck(i0,j0)=icheck(i0,j0)+1
        icheck(i1,j0)=icheck(i1,j0)+1
        icheck(i0,j1)=icheck(i0,j1)+1
        icheck(i1,j1)=icheck(i1,j1)+1

        ap00(jj,ii)=ap(i0,j0)
        ap01(jj,ii)=ap(i1,j0)
        ap10(jj,ii)=ap(i0,j1)
        ap11(jj,ii)=ap(i1,j1)

        ae00(jj,ii)=ae(i0,j0)
        ae01(jj,ii)=ae(i1,j0)
        ae10(jj,ii)=ae(i0,j1)
        ae11(jj,ii)=ae(i1,j1)

        aw00(jj,ii)=aw(i0,j0)
        aw01(jj,ii)=aw(i1,j0)
        aw10(jj,ii)=aw(i0,j1)
        aw11(jj,ii)=aw(i1,j1)

        an00(jj,ii)=an(i0,j0)
        an01(jj,ii)=an(i1,j0)
        an10(jj,ii)=an(i0,j1)
        an11(jj,ii)=an(i1,j1)

        as00(jj,ii)=as(i0,j0)
        as01(jj,ii)=as(i1,j0)
        as10(jj,ii)=as(i0,j1)
        as11(jj,ii)=as(i1,j1)

        bb00(jj,ii)=bb(i0,j0)
        bb01(jj,ii)=bb(i1,j0)
        bb10(jj,ii)=bb(i0,j1)
        bb11(jj,ii)=bb(i1,j1)
        end do
        end do

        !check write #2
        ! write(*,*) "#2"
        ! do i=1,m
        ! do j=1,n
        ! if(icheck(i,j).ne.1) write(*,*) "#2 i,j, icheck=", i,j, icheck(i,j)
        ! end do 
        ! end do


        do iter = 1, iter_max   ! -- SOR iteration start

        error=0.0

        ! default periodic condition in y-direction
        if(j_e.eq.0)then  ! n=even j_e=0 & j_o=1
        do ii = 1, mm
        f10(   0,ii)=f10(nn,ii)
        f11(   0,ii)=f11(nn,ii)
        f00(nn+1,ii)=f00(1 ,ii)
        f01(nn+1,ii)=f01(1 ,ii)
        end do
                        else  ! n=odd j_e=1 & j_o=0  
        do ii = 1, mm
        f00(   1,ii)=f10(nn,ii)
        f01(   1,ii)=f11(nn,ii)
        f00(nn+1,ii)=f10(1 ,ii)
        f01(nn+1,ii)=f11(1 ,ii)
        end do
                        end if

        ! default Dirichlet(i=m+1) / Neumann(i=0) condition in x-direction
        if(i_e.eq.0)then  ! m=even i_e=0 & i_o=1
        do jj = 1, nn
        f01(jj,0   )=f00(jj,1 )  ! Neumann
        f11(jj,0   )=f10(jj,1 )  ! Neumann
        f00(jj,mm+1)=0.          ! Dirichlet
        f10(jj,mm+1)=0.          ! Dirichlet
        end do
                        else  ! m=odd i_e=1 & i_o=0  
        do jj = 1, nn
        f00(jj,1   )=f01(jj,1 )  ! Neumann
        f10(jj,1   )=f11(jj,1 )  ! Neumann
        f00(jj,mm+1)=0.          ! Dirichlet
        f10(jj,mm+1)=0.          ! Dirichlet
        end do
                        end if

        ! NO need old data replace

        ! f00 corretion

        do ii = 1, mm
        do jj = 1, nn
        correct =( bb00(jj,ii)		              &
                        -ae00(jj,ii)*f01(jj  ,ii  )     &
                                -aw00(jj,ii)*f01(jj  ,ii-1)	  &
                        -an00(jj,ii)*f10(jj  ,ii  )     &
                                -as00(jj,ii)*f10(jj-1,ii  )	  &
                )/ap00(jj,ii)

        !error = max(error, abs(f01(jj,ii)-correct))

        f00(jj,ii)= f00(jj,ii)*(1.-relux_factor) &
                        + correct*relux_factor
        end do
        end do

        ! f11 corretion
        do ii = 1, mm
        do jj = 1, nn
        correct =( bb11(jj,ii)					  &
                        -ae11(jj,ii)*f10(jj  ,ii+1)     &
                                -aw11(jj,ii)*f10(jj  ,ii  )	  &
                        -an11(jj,ii)*f01(jj+1,ii  )     &
                                -as11(jj,ii)*f01(jj  ,ii  )	  &
                )/ap11(jj,ii)

        ! error = max(error, abs(f11(jj,ii)-correct))

        f11(jj,ii)= f11(jj,ii)*(1.-relux_factor) &
                        + correct*relux_factor
        end do
        end do

        ! f01 corretion
        do ii = 1, mm
        do jj = 1, nn
        correct =( bb01(jj,ii)					  &
                        -ae01(jj,ii)*f00(jj  ,ii+1)     &
                                -aw01(jj,ii)*f00(jj  ,ii  )	  &
                        -an01(jj,ii)*f11(jj  ,ii  )     &
                                -as01(jj,ii)*f11(jj-1,ii  )	  &
                )/ap01(jj,ii)

        !error = max(error, abs(f01(jj,ii)-correct))

        f01(jj,ii)= f01(jj,ii)*(1.-relux_factor) &
                        + correct*relux_factor
        end do
        end do

        ! f10 corretion
        do ii = 1, mm
        do jj = 1, nn
        correct =( bb10(jj,ii)					                    &
                        -ae10(jj,ii)*f11(jj  ,ii  )     &
                                -aw10(jj,ii)*f11(jj  ,ii-1)	  &
                        -an10(jj,ii)*f00(jj+1,ii  )     &
                                -as10(jj,ii)*f00(jj  ,ii  )	  &
                )/ap10(jj,ii)

        !error = max(error, abs(f10(jj,ii)-correct))

        f10(jj,ii)= f10(jj,ii)*(1.-relux_factor) &
                        + correct*relux_factor
        end do
        end do

        ! write(*,*)'CHECK iteration no.', iter,'  -- error=', error

        end do   ! SOR intereation finish  


        ! reconvert data  p <= f00,f01,f10,f11
        do ii = 1, mm
        do jj = 1, nn
        i0= (ii-1)*2 +i_o
        i1= (ii-1)*2 +i_o+1
        j0= (jj-1)*2 +j_o
        j1= (jj-1)*2 +j_o+1
        p(i0,j0)=f00(jj,ii)
        p(i1,j0)=f01(jj,ii)
        p(i0,j1)=f10(jj,ii)
        p(i1,j1)=f11(jj,ii)
        end do
        end do

        ! default Dirichlet(i=m+1) / Neumann(i=0) condition in x-direction
        do j=1,n
        p(0,j)  = p(1,j)
        p(m+1,j)= 0.
        end do
        ! default periodic condition in y-direction
        do i=1,m+1
        p(i,0)  = p(i,n)
        p(i,n+1)= p(i,1)
        end do

        ! write(*,*)'SOR_X iteration no.', iter-1,'  -- error=', error

        ! ----------------

        return
end subroutine solve_matrix_x

subroutine  solve_matrix_x_vec_omp (p, ap, ae, aw, an, as, bb, m, n, relux_factor, iter_max)
        use global_2d
        implicit none
        real,intent(inout),dimension(0:md,0:nd)::p
        real,intent(in),dimension(0:md,0:nd)::ap, ae, aw, an, as, bb
        integer,intent(in)::m, n
        real,intent(in)::relux_factor
        integer,intent(in)::iter_max 
        ! local variables
        real::error, correct
        real,dimension(0:md/2+1,0:nd/2+1)::f00,f01,f10,f11
        real,dimension(0:md/2+1,0:nd/2+1)::ap00,ap01,ap10,ap11, &
                                                ae00,ae01,ae10,ae11, &
                                                aw00,aw01,aw10,aw11, &
                                                an00,an01,an10,an11, &
                                                as00,as01,as10,as11, &
                                                bb00,bb01,bb10,bb11
        integer::i, j, iter
        integer::i_o, j_o, i_e, j_e, ii, jj, mm, nn
        integer::i0,i1,j0,j1
        !check counter
        integer,dimension(0:md,0:nd):: icheck

        !$omp parallel private(iter, i, j, k, ii, jj) &
        !$omp & shared(iter_max, relux_factor, m, n, mm, nn) &
        !!$omp & shared(error, p_old, p, ap, ae, aw, an, as, bb) &
        !$omp & shared(p_old, p, ap, ae, aw, an, as, bb) &
        !$omp & shared(i_o, j_o, i_e, j_e) &
        !$omp & shared(ap00, ap01, ap10, ap11) &
        !$omp & shared(ae00, ae01, ae10, ae11) &
        !$omp & shared(aw00, aw01, aw10, aw11) &
        !$omp & shared(an00, an01, an10, an11) &
        !$omp & shared(as00, as01, as10, as11) &
        !$omp & shared(bb00, bb01, bb10, bb11) &
        !$omp & default(none)

        ! ----------------
        !   SOR algorithm
        ! ----------------

        !convert local matrix  p to f00 f01 f10 f11
        i_e= mod(m,2)  ! 0 for m=even / 1 for m=odd
        j_e= mod(n,2)  ! 0 for n=even / 1 for n=odd
        i_o= 1-i_e     ! 0 for m=odd / 1 for m=even  
        j_o= 1-j_e     ! 0 for n=odd / 1 for n=even
        mm =(m+i_e)/2 ! p(1,*)=f*1(*,1) for m=odd / p(1,*)=f*0(1,*) for m=even
        nn =(n+j_e)/2 ! p(*,1)=f1*(*,1) for n=odd / p(1,*)=f0*(*,1) for n=even    
                                ! (0,0) <= f**(jj,ii) <= (nn+1,mm+1)

        !check write #0
        ! write(*,*) "#0 m,n= ",m,n," =>   mm,nn= ",mm,nn
        ! write(*,*) "i_e, i_o, j_e, j_o =",i_e, i_o, j_e, j_o
        ! write(*,*) ((ii-i_e)*2-i_o, (ii-i_e)*2-i_o+1,ii=1,mm)
        ! write(*,*) ((jj-j_e)*2-j_o, (jj-j_e)*2-j_o+1,jj=1,nn)
        !check write #0

        ! check counter set #1
        ! do i=1,m
        ! do j=1,n
        ! icheck(i,j)=0
        ! end do
        ! end do

        ! convert data  p => f00,f01,f10,f11
        !$omp do
        do ii = 1, mm
        do jj = 1, nn
        i0= (ii-1)*2 +i_o
        i1= (ii-1)*2 +i_o+1
        j0= (jj-1)*2 +j_o
        j1= (jj-1)*2 +j_o+1

        f00(jj,ii)=p(i0,j0)
        f01(jj,ii)=p(i1,j0)
        f10(jj,ii)=p(i0,j1)
        f11(jj,ii)=p(i1,j1)
        end do
        end do
        !$omp end do

        ! convert matrix components  ap,ae,aw,an,an,bb
        !$omp do
        do ii = 1, mm
        do jj = 1, nn
        i0= (ii-1)*2 +i_o
        i1= (ii-1)*2 +i_o+1
        j0= (jj-1)*2 +j_o
        j1= (jj-1)*2 +j_o+1
        ! check counter
        ! icheck(i0,j0)=icheck(i0,j0)+1
        ! icheck(i1,j0)=icheck(i1,j0)+1
        ! icheck(i0,j1)=icheck(i0,j1)+1
        ! icheck(i1,j1)=icheck(i1,j1)+1

        ap00(jj,ii)=ap(i0,j0)
        ap01(jj,ii)=ap(i1,j0)
        ap10(jj,ii)=ap(i0,j1)
        ap11(jj,ii)=ap(i1,j1)

        ae00(jj,ii)=ae(i0,j0)
        ae01(jj,ii)=ae(i1,j0)
        ae10(jj,ii)=ae(i0,j1)
        ae11(jj,ii)=ae(i1,j1)

        aw00(jj,ii)=aw(i0,j0)
        aw01(jj,ii)=aw(i1,j0)
        aw10(jj,ii)=aw(i0,j1)
        aw11(jj,ii)=aw(i1,j1)

        an00(jj,ii)=an(i0,j0)
        an01(jj,ii)=an(i1,j0)
        an10(jj,ii)=an(i0,j1)
        an11(jj,ii)=an(i1,j1)

        as00(jj,ii)=as(i0,j0)
        as01(jj,ii)=as(i1,j0)
        as10(jj,ii)=as(i0,j1)
        as11(jj,ii)=as(i1,j1)

        bb00(jj,ii)=bb(i0,j0)
        bb01(jj,ii)=bb(i1,j0)
        bb10(jj,ii)=bb(i0,j1)
        bb11(jj,ii)=bb(i1,j1)
        end do
        end do
        !$omp end do

        !check write #2
        ! write(*,*) "#2"
        ! do i=1,m
        ! do j=1,n
        ! if(icheck(i,j).ne.1) write(*,*) "#2 i,j, icheck=", i,j, icheck(i,j)
        ! end do 
        ! end do


        do iter = 1, iter_max   ! -- SOR iteration start

        error=0.0

        ! default periodic condition in y-direction
        if(j_e.eq.0)then  ! n=even j_e=0 & j_o=1
        !$omp do
        do ii = 1, mm
        f10(   0,ii)=f10(nn,ii)
        f11(   0,ii)=f11(nn,ii)
        f00(nn+1,ii)=f00(1 ,ii)
        f01(nn+1,ii)=f01(1 ,ii)
        end do
        !$omp end do
        else  ! n=odd j_e=1 & j_o=0  
        !$omp do
        do ii = 1, mm
        f00(   1,ii)=f10(nn,ii)
        f01(   1,ii)=f11(nn,ii)
        f00(nn+1,ii)=f10(1 ,ii)
        f01(nn+1,ii)=f11(1 ,ii)
        end do
        end if
        !$omp end do

        ! default Dirichlet(i=m+1) / Neumann(i=0) condition in x-direction
        if(i_e.eq.0)then  ! m=even i_e=0 & i_o=1
        !$omp do
        do jj = 1, nn
        f01(jj,0   )=f00(jj,1 )  ! Neumann
        f11(jj,0   )=f10(jj,1 )  ! Neumann
        f00(jj,mm+1)=0.          ! Dirichlet
        f10(jj,mm+1)=0.          ! Dirichlet
        end do
        !$omp end do
        else  ! m=odd i_e=1 & i_o=0  
        !$omp do
        do jj = 1, nn
        f00(jj,1   )=f01(jj,1 )  ! Neumann
        f10(jj,1   )=f11(jj,1 )  ! Neumann
        f00(jj,mm+1)=0.          ! Dirichlet
        f10(jj,mm+1)=0.          ! Dirichlet
        end do
        !$omp end do
        end if

        ! NO need old data replace

        ! f00 corretion
        !$omp do
        do ii = 1, mm
        do jj = 1, nn
        correct =( bb00(jj,ii)		              &
                        -ae00(jj,ii)*f01(jj  ,ii  )     &
                                -aw00(jj,ii)*f01(jj  ,ii-1)	  &
                        -an00(jj,ii)*f10(jj  ,ii  )     &
                                -as00(jj,ii)*f10(jj-1,ii  )	  &
                )/ap00(jj,ii)

        !error = max(error, abs(f01(jj,ii)-correct))

        f00(jj,ii)= f00(jj,ii)*(1.-relux_factor) &
                        + correct*relux_factor
        end do
        end do
        !$omp end do

        ! f11 corretion
        !$omp do
        do ii = 1, mm
        do jj = 1, nn
        correct =( bb11(jj,ii)					  &
                        -ae11(jj,ii)*f10(jj  ,ii+1)     &
                                -aw11(jj,ii)*f10(jj  ,ii  )	  &
                        -an11(jj,ii)*f01(jj+1,ii  )     &
                                -as11(jj,ii)*f01(jj  ,ii  )	  &
                )/ap11(jj,ii)

        ! error = max(error, abs(f11(jj,ii)-correct))

        f11(jj,ii)= f11(jj,ii)*(1.-relux_factor) &
                        + correct*relux_factor
        end do
        end do
        !$omp end do

        ! f01 corretion
        !$omp do
        do ii = 1, mm
        do jj = 1, nn
        correct =( bb01(jj,ii)					  &
                        -ae01(jj,ii)*f00(jj  ,ii+1)     &
                                -aw01(jj,ii)*f00(jj  ,ii  )	  &
                        -an01(jj,ii)*f11(jj  ,ii  )     &
                                -as01(jj,ii)*f11(jj-1,ii  )	  &
                )/ap01(jj,ii)

        !error = max(error, abs(f01(jj,ii)-correct))

        f01(jj,ii)= f01(jj,ii)*(1.-relux_factor) &
                        + correct*relux_factor
        end do
        end do
        !$omp end do

        ! f10 corretion
        !$omp do
        do ii = 1, mm
        do jj = 1, nn
        correct =( bb10(jj,ii)					                    &
                        -ae10(jj,ii)*f11(jj  ,ii  )     &
                                -aw10(jj,ii)*f11(jj  ,ii-1)	  &
                        -an10(jj,ii)*f00(jj+1,ii  )     &
                                -as10(jj,ii)*f00(jj  ,ii  )	  &
                )/ap10(jj,ii)

        !error = max(error, abs(f10(jj,ii)-correct))

        f10(jj,ii)= f10(jj,ii)*(1.-relux_factor) &
                        + correct*relux_factor
        end do
        end do
        !$omp end do

        ! write(*,*)'CHECK iteration no.', iter,'  -- error=', error

        end do   ! SOR intereation finish  


        ! reconvert data  p <= f00,f01,f10,f11
        !$omp do
        do ii = 1, mm
        do jj = 1, nn
        i0= (ii-1)*2 +i_o
        i1= (ii-1)*2 +i_o+1
        j0= (jj-1)*2 +j_o
        j1= (jj-1)*2 +j_o+1
        p(i0,j0)=f00(jj,ii)
        p(i1,j0)=f01(jj,ii)
        p(i0,j1)=f10(jj,ii)
        p(i1,j1)=f11(jj,ii)
        end do
        end do
        !$omp end do

        ! default Dirichlet(i=m+1) / Neumann(i=0) condition in x-direction
        !$omp do
        do j=1,n
        p(0,j)  = p(1,j)
        p(m+1,j)= 0.
        end do
        !$omp end do
        ! default periodic condition in y-direction
        !$omp do
        do i=1,m+1
        p(i,0)  = p(i,n)
        p(i,n+1)= p(i,1)
        end do
        !$omp end do

        ! write(*,*)'SOR_X iteration no.', iter-1,'  -- error=', error

        ! ----------------

        return
end subroutine solve_matrix_x_vec_omp

end module solver
