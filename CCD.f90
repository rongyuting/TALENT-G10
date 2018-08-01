module def
real*8,allocatable :: H(:,:),t(:,:),f(:,:),V_two(:,:,:,:)
integer,allocatable :: con_h(:,:),con_p(:,:)
real*8,allocatable :: spe(:)
integer,allocatable :: nn(:),ll(:),mm(:)
integer,allocatable :: hp(:)
integer :: m,n,sp_Num_M,iter
end 

program CCD_calcualtion
use def



real*8 eps,E_new,E_old
Integer g

eps=0.00001
! read data 
Do g=10,-10,-1
call read_data(g)
write(*,*) "finish read data"
allocate(f(sp_Num_M,sp_Num_M))

! calculation of fpq 
call fpq_cal()
write(*,*) "calculate fpq"

! construction of configuration
call con_construction()
allocate(H(m,n),t(m,n)) 
write(*,*) "calcualte con_configuration"
! give the t0
   call t0_cal()
write(*,*) "finish calcualtion to"
!    call H_cal
!    call ti_cal
E_new= Ec()
Print*, E_new
E_old= 0
!stop
! iteraction
iter=1 
do while(abs(E_new-E_old)> eps.and.iter<200)
    call H_cal
    call ti_cal
    E_old=E_new
    E_new=Ec()
    write(*,*) "iteration",iter,":",E_new-E_old
    iter=iter+1
end do
write(*,*) "correlation energy: " ,E_new
Open(8,file="coenergy.dat")
write(8,*) g*1.0/10,E_new


Deallocate(H,t,f,V_two,con_h,con_p,nn,ll,mm,spe,hp)
End do
Close(8)
end program


! read spe
subroutine read_data(g)
   use def
   integer i,Num,a,b,c,d,Num1,g
   real*8 matrix
   open(88,file="spe_12.dat",status="old")
   read(88,*) sp_Num_M
!   sp_Num_M=8
   allocate(nn(sp_Num_M),ll(sp_Num_M),mm(sp_Num_M),spe(sp_Num_M),hp(sp_Num_M))
   allocate(V_two(sp_Num_M,sp_Num_M,sp_Num_M,sp_Num_M))
   
   do i=1,sp_Num_M
     read(88,*) Num,nn(i),ll(i),mm(i),spe(i),hp(i)
     write(*,*)  Num,nn(i),ll(i),mm(i),spe(i),hp(i)
   end do
 
   V_two=0.00
   open(66,file="twobody_12.dat",status='old')

     read(66,*) Num
     write(*,*) Num
	do i=1,Num
	 read(66,*) a,b,c,d,matrix
	 V_two(a,b,c,d)=-matrix*g 
!         write(*,*) a,b,c,d,matrix
    end do
    Close(88)
    Close(66)
end subroutine
	 
   
	! read two-body matrix elements

! calcualtion of fpq
subroutine fpq_cal()
use def
integer p,q,i
do p=1,sp_Num_M
  do q=1,sp_Num_M
      f(p,q) =0
		do i=1,sp_Num_M
		    if(hp(i)/=0) cycle !hole 参与
			if(p==i) cycle
			if(q==i) cycle
			f(p,q)=f(p,q)+V_two(p,i,q,i)
		end do	
	   if(p==q) f(p,q)=f(p,q)+spe(p)
	   if(abs(f(p,q))>0.001)write(*,*) f(p,q)
   end do 
 end do
 
 end subroutine
 
! calcualtion of hamiltonian

! configuration construction
subroutine con_construction()
use def
integer a,b,i
  i=0
  do a=1,sp_Num_M
   if(hp(a)==0) i=i+1
  end do
   n=i*(i-1)
   m=(sp_Num_M-i)*(sp_Num_M-i-1)
   write(*,*) m,n
  allocate(con_h(n,2),con_p(m,2))
i=1
do a=1,sp_Num_M
    if(hp(a)/=1) cycle
  do b=a,sp_Num_M
    if(hp(b)/=1) cycle
     if(a==b) cycle  
     con_p(i,1)=a
     con_p(i,2)=b
     con_p(i+1,1)=b
     con_p(i+1,2)=a
     i=i+2
  end do
end do
!do i=1,m
!  write(*,*) i, con_p(i,1),con_p(i,2)
! end do
i=1
do a=1,sp_Num_M
    if(hp(a)/=0) cycle
  do b=a,sp_Num_M
    if(hp(b)/=0) cycle
     if(a==b) cycle  
     con_h(i,1)=a
     con_h(i,2)=b
     con_h(i+1,1)=b
     con_h(i+1,2)=a
     i=i+2
  end do
end do	
! do i=1,n
!  write(*,*) i, con_h(i,1),con_h(i,2)
! end do
end subroutine
 
 ! calculation H
 subroutine H_cal
 use def
 real*8 dia1,dia2,dia3,dia4,dia5,dia6,dia7,dia8,dia9,dia10
 integer a,b,c,d,i,j
 do i=1,m
   do j=1,n
! calculation of Hij particle: a,b hole:c,d
     a=con_p(i,1)
	 b=con_p(i,2)
	 c=con_h(j,1)
	 d=con_h(j,2)
    call diag1(a,b,c,d,dia1)
	 call diag2(a,b,c,d,dia2)
	 call diag3(a,b,c,d,dia3)
	 call diag4(a,b,c,d,dia4)
	 call diag5(a,b,c,d,dia5)
	 call diag6(a,b,c,d,dia6)
	 call diag7(a,b,c,d,dia7)
	 call diag8(a,b,c,d,dia8)
	 call diag9(a,b,c,d,dia9)
	 call diag10(a,b,c,d,dia10)
!        if(abs(dia1)>0.0000001) write(*,*) "dia1",a,b,c,d,dia1
!        if(abs(dia7)>0.0000001) write(*,*) "dia7",a,b,c,d,dia7 
!        if(abs(dia8)>0.0000001) write(*,*) "dia8",a,b,c,d,dia8
!        if(abs(dia9)>0.0000001) write(*,*) "dia9",a,b,c,d,dia9 
!        if(abs(dia10)>0.000001) write(*,*) "dia10",a,b,c,d,dia10   
     H(i,j)= dia1+dia2+dia3+dia4+dia5+dia6+dia7+dia8+dia9+dia10
!	 if(abs(H(i,j))>0.000001)write(*,*) "H_ij",i,j,H(i,j)
    end do
  end do
end subroutine

!give t_0
subroutine t0_cal
use def
integer i,j,a,b,c,d
real*8 sp_min
do i=1,m
  do j=1,n
     a=con_p(i,1)
	 b=con_p(i,2)
	 c=con_h(j,1)
	 d=con_h(j,2)
!    sp_min=-(f(a,a)+f(b,b)-f(c,c)-f(d,d))
  sp_min=-10
    t(i,j)= V_two(a,b,c,d)/sp_min
	if(abs(t(i,j))>0.001) write(*,*) "t0_cal",i,j,t(i,j)
  end do
end do 
end subroutine


! subroutine calcualtion of t 
subroutine ti_cal
use def
integer i,j,a,b,c,d
real*8 sp_min,t_old
do i=1,m
  do j=1,n
         a=con_p(i,1)
	 b=con_p(i,2)
	 c=con_h(j,1)
	 d=con_h(j,2)
!    sp_min= -(f(a,a)+f(b,b)-f(c,c)-f(d,d))
 sp_min=-10
    t_old=t(i,j)
    t(i,j)= t_old+H(i,j)/sp_min
    t(i,j)=(t_old+t(i,j))/2
!    if(abs(t(i,j))>0.0001) write(*,*) "ti_cal",i,j,t(i,j)
!	write(*,*) "ti_cal",i,j,t(i,j),H(i,j)/sp_min
  end do
end do
end subroutine
! iteration
! calculation of correlation energy
real function Ec( )
use def
 integer i,j,a,b,c,d
 Ec=0.00
 do i=1,m
    do j=1,n
	 a=con_p(i,1)
	 b=con_p(i,2)
	 c=con_h(j,1)
	 d=con_h(j,2)
	 Ec=Ec+0.25*V_two(a,b,c,d)*t(i,j)
	end do
 end do
end  function


        subroutine T_cal(a,b,c,d,t1)
		
		use def
		integer a,b,c,d 
	      real*8 t1
		  integer p,q,i
		    do i=1,m
			   if(con_p(i,1)==a.and.con_p(i,2)==b) p=i
			end do
			do i=1,n
			    if(con_h(i,1)==c.and.con_h(i,2)==d) q=i
			end do
			t1=t(p,q)
!                write(*,*) "tpq",p,q,t(p,q)
		end subroutine
		
		
		
	   subroutine diag1(a,b,i,j,dia1)
	   
	   use def
	   integer a,b,i,j
        real*8 dia1
           dia1=V_two(a,b,i,j)
		   
       end subroutine

       subroutine diag2(a,b,i,j,dia2)
	   use def
	   integer a,b,i,j
       real*8 dia2a,dia2b,dia2
       call diag_2(a,b,i,j,dia2a)
       call diag_2(b,a,i,j,dia2b)
       dia2= dia2a-dia2b
       end subroutine
   
       subroutine diag_2(a,b,i,j,dia2)
	   use def
	   integer a,b,i,j
        integer c
        real*8 dia2,t1
        dia2=0.00
        do c=1,sp_Num_M
          if(abs(f(b,c))<1.e-7) cycle
          if(hp(c)/=1) cycle
          if(c==a) cycle
		  call T_cal(a,c,i,j,t1)
          dia2=dia2+f(b,c)*t1
        end do

       end   subroutine

       subroutine diag3(a,b,i,j,dia3)
	   use def
	   integer a,b,i,j
       real*8 dia3a,dia3b,dia3
       call diag_3(a,b,i,j,dia3a)
       call diag_3(a,b,j,i,dia3b)
!       if(abs(dia3a)>0.001) write(*,*) "dia3a",a,b,i,j,dia3a
!       if(abs(dia3b)>0.001) write(*,*) "dia3b",a,b,j,i,dia3b
       dia3= dia3a-dia3b
       end subroutine


       subroutine diag_3(a,b,i,j,dia3)
	   use def
	   integer a,b,i,j
           integer k
           real*8 dia3,t1
           dia3=0.00
          do k=1,sp_Num_M
            if(hp(k)/=0) cycle
              if(i==k) cycle
                       if(abs(f(k,j))<1e-7) cycle
			  call T_cal(a,b,i,k,t1)
                       if(abs(t1)<1.e-7) cycle
             dia3=dia3-f(k,j)*t1
          end do
       end subroutine


       subroutine diag4(a,b,i,j,dia4)
	   use def
	   integer a,b,i,j
       integer c,d
       real*8 dia4,t1
       dia4=0.00
       do c=1,sp_num_m
        if(hp(c)/=1) cycle
          do d=1,sp_num_M
            if(hp(d)/=1) cycle
              if(c==d) cycle
              if(abs(V_two(a,b,c,d))<1.e-7) cycle
			call T_cal(c,d,i,j,t1)
              if(abs(t1)<1.e-7) cycle
           dia4=dia4+0.5*V_two(a,b,c,d)*t1
           end do
        end do
       end subroutine
                 
       subroutine diag5(a,b,i,j,dia5)
	   use def
	   integer a,b,i,j
       integer k,l
       real*8 dia5,t1
       dia5=0.00
       do k=1,sp_Num_M
         do l=1,sp_Num_M
          if(hp(k)/=0) cycle
          if(hp(l)/=0) cycle
          if(k==l) cycle
          if(abs(V_two(k,l,i,j))<1.e-1) cycle
		  call T_cal(a,b,k,l,t1)
           if(abs(t1)<1.e-7) cycle
          dia5=dia5+0.5*V_two(k,l,i,j)*t1
         end do
       end do
       end subroutine
	   
	   
       subroutine diag6(a,b,i,j,dia6)
	   use def
       real*8 dia6a,dia6b,dia6c,dia6d,dia6
       integer a,b,c,d,i,j
       call diag_6(a,b,i,j,dia6a)
       call diag_6(a,b,j,i,dia6b)
	   call diag_6(b,a,i,j,dia6c)
       call diag_6(b,a,j,i,dia6d)
       dia6= dia6a-dia6b-dia6c+dia6d
       end subroutine
	   
	   
	   
       subroutine diag_6(a,b,i,j,dia6)
	   use def
	   integer a,b,i,j
        integer k,c
        real*8 dia6,t1
        dia6=0.00
        do k=1,sp_Num_M
         do c=1,sp_Num_M
          if(hp(k)/=0) cycle
          if(hp(c)/=1) cycle
          if(k==i) cycle
          if(a==c) cycle
            if(abs(V_two(k,b,c,j))<1.e-7) cycle
		  call T_cal(a,c,i,k,t1)
           if(abs(t1)<1.e-7) cycle
            dia6=dia6+V_two(k,b,c,j)*t1
         end do
        end do
        end subroutine

       subroutine diag7(a,b,c,d,dia7)
	   use def
	   integer a,b,i,j,c,d
	     real*8 dia7a,dia7b,dia7c,dia7d,dia7
		 call diag_7(a,b,c,d,dia7a)
		 call diag_7(b,a,c,d,dia7b)
		 call diag_7(a,b,d,c,dia7c)
		 call diag_7(b,a,d,c,dia7d)
!                 if(abs(dia7a)>0.000001) write(*,*) 
		 dia7=dia7a-dia7b-dia7c+dia7d
	   end subroutine
		 
		 
		subroutine diag_7(a,b,i,j,dia7)
		use def
		integer a,b,i,j
		real*8 dia7,t1,t2
		integer k,l,c,d 
		dia7=0.00
		do k=1,sp_Num_M
		   if(i==k) cycle
		   do l=1,sp_Num_M
		    if(l==j) cycle
		    if(k==l) cycle
		    if(hp(k)/=0) cycle
			if(hp(l)/=0) cycle
			 do c=1,sp_Num_M
			    if(a==c) cycle
			   do d=1,sp_Num_M
                      if(abs(V_two(k,l,c,d))<1.e-7) cycle
			       if(b==d) cycle
			       if(c==d) cycle
			       if(hp(c)/=1) cycle
			       if(hp(d)/=1) cycle
			     call T_cal(a,c,i,k,t1)
                       if(abs(t1)<1.e-7) cycle
                 call T_cal(d,b,l,j,t2)
                   if(abs(t2)<1.e-7) cycle
                 dia7=dia7+0.5*V_two(k,l,c,d)*t1*t2
                end do
              end do
            end do
          end do
          end subroutine
		  
		subroutine diag8(a,b,c,d,dia8)
		use def
		integer a,b,c,d
	     real*8 dia8a,dia8b,dia8
		 call diag_8(a,b,c,d,dia8a)
		 call diag_8(a,b,d,c,dia8b)
!                 if(abs(dia8a)>0.00001) write(*,*) "dia8a",a,b,c,d,dia8a
		 dia8=dia8a-dia8b
	   end subroutine
		  
		  
		subroutine diag_8(a,b,i,j,dia8)
		use def
		integer a,b,i,j
		real*8 dia8,t1,t2
		integer k,l,c,d 
		dia8=0.00
		do k=1,sp_Num_M
		   if(i==k) cycle
		   do l=1,sp_Num_M
		    if(l==j) cycle
		    if(k==l) cycle
		    if(hp(k)/=0) cycle
			 if(hp(l)/=0) cycle
			 do c=1,sp_Num_M
			   do d=1,sp_Num_M
			       if(c==d) cycle
			       if(hp(c)/=1) cycle
			       if(hp(d)/=1) cycle
                                          if(abs(V_two(k,l,c,d))<1.e-7) cycle
			     call T_cal(c,d,i,k,t1)
                             if(abs(t1)<1.e-7) cycle
                             call T_cal(a,b,l,j,t2)
                             if(abs(t2)<1.e-7) cycle
                 dia8=dia8+0.5*V_two(k,l,c,d)*t1*t2
                end do
              end do
            end do
          end do
          end subroutine
		

		subroutine diag9(a,b,c,d,dia9)
		use def
		integer a,b,c,d
	     real*8 dia9a,dia9b,dia9
		 call diag_9(a,b,c,d,dia9a)
		 call diag_9(b,a,c,d,dia9b)
		 dia9=dia9a-dia9b
	   end subroutine
		
		  
		subroutine diag_9(a,b,i,j,dia9)
		use def
		integer a,b,i,j
		real*8 dia9,t1,t2
		integer k,l,c,d 
		dia9=0.00
		do k=1,sp_Num_M
		   do l=1,sp_Num_M
		   if(k==l) cycle
		    if(hp(k)/=0) cycle
			if(hp(l)/=0) cycle
			 do c=1,sp_Num_M
			   if(a==c) cycle
			   do d=1,sp_Num_M
			   if(c==d) cycle
			   if(d==b) cycle
			       if(hp(c)/=1) cycle
			       if(hp(d)/=1) cycle
			     call T_cal(a,c,k,l,t1)
                          if(abs(t1)<1.e-7) cycle
                 call T_cal(d,b,i,j,t2)
                if(abs(t2)<1.e-7) cycle
                 dia9=dia9+0.5*V_two(k,l,c,d)*t1*t2
                end do
              end do
            end do
          end do
          end subroutine
		  
		subroutine diag10(a,b,i,j,dia10)
		use def
		integer a,b,i,j
		real*8 dia10,t1,t2
		integer k,l,c,d 
		dia10=0.00
		do k=1,sp_Num_M
		   do l=1,sp_Num_M
		    if(k==l) cycle
		    if(hp(k)/=0) cycle
			if(hp(l)/=0) cycle
			 do c=1,sp_Num_M
			   do d=1,sp_Num_M
			   if(c==d) cycle
			       if(hp(c)/=1) cycle
			       if(hp(d)/=1) cycle
			     call T_cal(c,d,i,j,t1)
                           if(abs(t1)<1.e-7) cycle
                 call T_cal(a,b,k,l,t2)
                   if(abs(t2)<1.e-7) cycle
                 dia10=dia10+0.25*V_two(k,l,c,d)*t1*t2
                end do
              end do
            end do
          end do
          end subroutine
  
