program saint_venant_mac_cormack
  ! ================================================================
  ! PARTE 1. Inicio y precision numerica
  ! Este programa calcula flujo en un canal TRAPEZOIDAL.
  ! Esquema: MacCormack Explicito (Forma Conservativa)
  ! ================================================================
  use iso_fortran_env, only: real64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  implicit none
  integer, parameter :: dp=real64

  ! Variables base
  real(dp) :: bw,z,L,dx,dt,tmax,g,nman,S0,A0,Q0,Qup,ratio,t,cfl,h0,v0,c
  ! Vectores de estado
  real(dp), allocatable :: A(:),Q(:),Ap(:),Qp(:),Ac(:),Qc(:),An(:),Qn(:)
  ! Vectores de flujos y fuerzas
  real(dp), allocatable :: press(:),press_p(:),F2(:),F2p(:),Sf_vec(:),Sfp_vec(:)
  integer :: nnodes,nsteps,step,i,ios,iu
  real(dp) :: visc = 0.05_dp ! Viscosidad artificial para amortiguar oscilaciones del resalto

  print '(A)','Saint-Venant Canal Trapezoidal: MacCormack Explicito Conservativo.'
  call pedir('Ancho base bw (m): ',bw,0.0_dp,1000.0_dp)
  call pedir('Talud z (H:1V) (0 rectangular, >0 trapezoidal): ',z,0.0_dp,10.0_dp,.true.)
  call pedir('Longitud L (m): ',L,0.0_dp,1.0e6_dp)
  call pedir('Paso espacial dx (m): ',dx,0.0_dp,L)
  ratio=L/dx
  nnodes=nint(ratio)+1

  call pedir('Paso temporal dt (s): ',dt,0.0_dp,1.e9_dp)
  call pedir('Tiempo maximo (s): ',tmax,0.0_dp,1.e9_dp)
  call pedir('Gravedad g (m/s2): ',g,0.0_dp,40.0_dp)
  call pedir('Coeficiente Manning n: ',nman,0.0_dp,1.0_dp,.true.)
  call pedir('Pendiente S0 (m/m): ',S0,-1.0_dp,1.0_dp,.true.)
  call pedir('Area inicial A0 (m2): ',A0,0.0_dp,1.e6_dp)
  call pedir('Caudal inicial Q0 (m3/s): ',Q0,-1.e6_dp,1.e6_dp,.true.)
  call pedir('Caudal contorno aguas arriba Qup (m3/s): ',Qup,-1.e6_dp,1.e6_dp,.true.)

  ! Calculo del tirante inicial h0 a partir de A0
  if (z > 0.0_dp) then
      h0 = (-bw + sqrt(bw**2 + 4.0_dp*z*A0)) / (2.0_dp*z)
  else
      h0 = A0 / bw
  end if

  v0 = Q0/A0
  c = sqrt(g * (A0 / (bw + 2.0_dp*z*h0)))
  cfl = dt*(abs(v0)+c)/dx
  if(cfl>1.0_dp) then
    write(*,'(A,ES12.4)') 'CUIDADO: El esquema explícito será inestable. CFL > 1.0 (CFL = ', cfl
  end if

  nsteps=nint(tmax/dt)
  allocate(A(nnodes),Q(nnodes),Ap(nnodes),Qp(nnodes),Ac(nnodes),Qc(nnodes),An(nnodes),Qn(nnodes))
  allocate(press(nnodes),press_p(nnodes),F2(nnodes),F2p(nnodes),Sf_vec(nnodes),Sfp_vec(nnodes))
  A=A0; Q=Q0; Q(1)=Qup

  open(newunit=iu,file='resultados_mac_cormack_trapecio.csv',status='replace',action='write',iostat=ios)
  write(iu,'(A)') 'paso;tiempo;nodo;x;A;Q;h;v;Rh;Sf'
  call guardar(0,0.0_dp)

  do step=1,nsteps
    t = step * dt

    ! --- PREDICTOR (Diferencias hacia ADELANTE) ---
    call presion_trapecio(A, press)
    do i=1, nnodes
       Sf_vec(i) = get_sf(A(i), Q(i))
       F2(i) = (Q(i)**2)/A(i) + g * press(i)
    end do

    do i=1, nnodes-1
      Ap(i) = A(i) - (dt/dx)*(Q(i+1) - Q(i))
      Qp(i) = Q(i) - (dt/dx)*(F2(i+1) - F2(i)) + g*A(i)*dt*(S0 - Sf_vec(i))
    end do
    ! Extrapolación de borde explícito aguas abajo
    Ap(nnodes) = A(nnodes) - (dt/dx)*(Q(nnodes) - Q(nnodes-1))
    Qp(nnodes) = Q(nnodes) - (dt/dx)*(F2(nnodes) - F2(nnodes-1)) + g*A(nnodes)*dt*(S0 - Sf_vec(nnodes))

    call validar(Ap,Qp,'predictor MacCormack')

    ! --- CORRECTOR (Diferencias hacia ATRAS) ---
    call presion_trapecio(Ap, press_p)
    do i=1, nnodes
       Sfp_vec(i) = get_sf(Ap(i), Qp(i))
       F2p(i) = (Qp(i)**2)/Ap(i) + g * press_p(i)
    end do

    do i=2, nnodes
      Ac(i) = Ap(i) - (dt/dx)*(Qp(i) - Qp(i-1))
      Qc(i) = Qp(i) - (dt/dx)*(F2p(i) - F2p(i-1)) + g*Ap(i)*dt*(S0 - Sfp_vec(i))
    end do
    ! Extrapolación de borde explícito aguas arriba
    Ac(1) = Ap(1) - (dt/dx)*(Qp(2) - Qp(1))
    Qc(1) = Qp(1) - (dt/dx)*(F2p(2) - F2p(1)) + g*Ap(1)*dt*(S0 - Sfp_vec(1))

    ! --- PROMEDIO Y ACTUALIZACION ---
    do i=1, nnodes
      An(i) = 0.5_dp * (A(i) + Ac(i))
      Qn(i) = 0.5_dp * (Q(i) + Qc(i))
    end do

    ! --- VISCOSIDAD ARTIFICIAL (Suaviza el frente de choque del resalto) ---
    if (visc > 0.0_dp) then
      do i = 2, nnodes-1
         An(i) = An(i) + visc * (A(i+1) - 2.0_dp * A(i) + A(i-1))
         Qn(i) = Qn(i) + visc * (Q(i+1) - 2.0_dp * Q(i) + Q(i-1))
      end do
    end if

    ! --- CONDICIONES DE FRONTERA REALES ---
    An(1) = A0  ! Condición de compuerta

    ! Cierre súbito de compuerta aguas abajo (genera el aumento del tirante aguas arriba).
    Qn(nnodes) = 0.0_dp

    call validar(An,Qn,'promedio predictor-corrector')
    A=An; Q=Qn
    call guardar(step,t)
  end do

  close(iu)
  write(*,'(A,I0,A,I0,A,ES12.4)') 'Listo. Nodos: ',nnodes,'; pasos: ',nsteps,'; CFL inicial: ',cfl

contains

  real(dp) function get_h(area)
    real(dp), intent(in) :: area
    if (z > 0.0_dp) then
      get_h = (-bw + sqrt(max(0.0_dp, bw**2 + 4.0_dp*z*area))) / (2.0_dp*z)
    else
      get_h = area / bw
    end if
  end function

  real(dp) function get_Rh(area)
    real(dp), intent(in) :: area
    real(dp) :: h_val, perim
    h_val = get_h(area)
    perim = bw + 2.0_dp*h_val*sqrt(1.0_dp + z**2)
    get_Rh = area / perim
  end function

  subroutine pedir(msg,x,lo,hi,inclusive)
    character(*),intent(in)::msg
    real(dp),intent(out)::x
    real(dp),intent(in)::lo,hi
    logical,intent(in),optional::inclusive
    integer::st
    logical::inc
    inc=.false.; if(present(inclusive)) inc=inclusive
    do
      write(*,'(A)',advance='no') msg
      read(*,*,iostat=st)x
      if(st==0.and.ieee_is_finite(x)) then
        if((x>lo.or.(inc.and.x>=lo)).and.x<=hi) return
      end if
      print '(A)','Valor invalido.'
    end do
  end subroutine

  subroutine presion_trapecio(aa,pp)
    real(dp),intent(in)::aa(:)
    real(dp),intent(out)::pp(:)
    integer :: k
    real(dp) :: h_val
    do k=1, size(aa)
       h_val = get_h(aa(k))
       pp(k) = (h_val**2)*(bw/2.0_dp + z*h_val/3.0_dp)
    end do
  end subroutine

  real(dp) function get_sf(area, caudal)
    real(dp),intent(in)::area, caudal
    real(dp)::vv, rh
    vv=caudal/area
    rh=get_Rh(area)
    ! Se asegura que el radio hidráulico no sea cero para evitar divisiones
    get_sf=nman**2*vv*abs(vv)/(max(1.e-12_dp, rh)**(4.0_dp/3.0_dp))
  end function

  subroutine validar(aa,qq,label)
    real(dp),intent(in)::aa(:),qq(:)
    character(*),intent(in)::label
    integer::k
    do k=1,size(aa)
      if(.not.ieee_is_finite(aa(k)).or..not.ieee_is_finite(qq(k)).or.aa(k)<=0) then
        write(*,'(A,A,A,I0)') 'ERROR: estado invalido en ',label,', nodo ',k
        stop 1
      end if
    end do
  end subroutine

  subroutine guardar(paso,tiempo)
    integer,intent(in)::paso
    real(dp),intent(in)::tiempo
    integer::k
    real(dp)::hh,vv,rh,fr
    do k=1,nnodes
      hh=get_h(A(k)); vv=Q(k)/A(k); rh=get_Rh(A(k))
      fr=nman**2*vv*abs(vv)/(max(1.e-12_dp, rh)**(4.0_dp/3.0_dp))
      write(iu,'(I0,A,ES12.2E3,A,I0,A,ES12.2E3,6(A,ES12.2E3))') paso,';',tiempo,';',k,';',real(k-1,dp)*dx,';', &
        A(k),';',Q(k),';',hh,';',vv,';',rh,';',fr
    end do
  end subroutine
end program
