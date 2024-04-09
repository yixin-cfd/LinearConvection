! 线性对流方程 u_t + u_x = 0; u(x, 0) = sin(x); 周期性条件; x in [0, 2pi]; 精确解：u(x,t) = sin(x -t)
! u_x: 采用五阶迎风格式
! u_t: 采用三阶龙格库塔法
program main
real*8:: t,dx,dt,pi, CFL 
integer:: m,i  
real*8, allocatable:: u0(:), x(:),  res(:), exact(:)
print*, "solve u_t + u_x = 0; period; u(x, 0) = sin(x)."
print*, "Exact: sin(x-t)"

print*, "Please enter the number of discrete points m and the time t:"
read(*, *) m, t
allocate(u0(m), x(m), res(m), exact(m))
! 计算初值
pi = acos(-1.0)
do i=1, m
x(i) = 2*pi*(i-1.d0)/(m-1.d0); ! 计算 x
u0(i) = sin(x(i));              ! 计算初值
enddo

dx = 2.d0*pi/(m-1.d0)              ! 计算 dx
CFL = 0.001d0                       ! 设置 CFL 数
dt = CFL /dx                    ! 计算dt

print*, "dx=",dx, "dt=", dt
 
! solve 数值解
CALL solve(res, u0, x, m, t, dx, dt)

! 计算精确解
do i=1, m
    exact(i) = sin(x(i) - t)
enddo

!---------------------- 后处理

open(99, file="linearConvResult.plt")
write(99, *) "variables=x, cal, exact"

do i=1,m
    write(99, "(4E20.12)") x(i), res(i), exact(i)
enddo

print*, "output: linearConvResult.plt"


deallocate(u0, x, res, exact)
end


function u_x(u, dx, m, j)
!---- 计算uj_x ------
! u: 场变量
! x: 
! m: 离散点数量
! j: 第几个点的导数
implicit none
integer m, j, idx
real*8:: u(m)
real*8:: u_x, dx

!=== 中心差分
! if(j > 1 .and. j < m) then 
!     u_x = (u(j+1)-u(j-1))/(2.d0*dx)
! else if(j == 1) then
!     u_x = (u(2) - u(1))/dx
! else 
!     u_x = (u(m)-u(m-1))/dx
! endif

!=== 五阶迎风格式
if(j>=4 .and. j <= (m-2)) then
u_x = (-2.d0*u(j-3) + 15.0*u(j-2) -60.d0 *u(j-1) + 20.d0*u(j) + 30.d0*u(j+1) -3.d0*u(j+2))/(60.d0*dx)

else if (j>1 .and. j<4) then ! 低阶格式：中心差分
u_x = (u(j+1)-u(j-1))/(2.d0*dx)

else if (j == 1) then           ! 第一个点的导数
u_x = 2.d0*(u(2) - u(1))/dx - (u(3)-u(1))/(2.d0*dx) 

else if(j == (m-1)) then        ! 低阶格式：中心差分
u_x = (u(j+1)-u(j-1))/(2.d0*dx)

else if(j == m) then            ! 最后一个点的导数       
    u_x = 2.d0*(u(m)-u(m-1))/dx - (u(m) -u(m-2))/(2.d0*dx)

endif

end


subroutine L(NL, u, dx, m)
!--- 计算非线性项(加上负号)
! NL: 返回的非线性项
! u: 场变量
! x: 坐标
! m: 离散点总数
implicit none
real*8:: u_x, dx
real*8:: NL(m), u(m)
integer:: i, m

do i=1, m
    NL(i) = -u_x(u, dx, m, i)
enddo

end


subroutine solve(res, u0, x, m, t, dx, dt)
!---求解方程
! res: t 时刻的解
! u0: 初始值
! x: 
! m:点总数
! t:真解
! dt: 时间离散
implicit none
integer m, i
real*8:: res(m), u0(m), x(m), U(m), U_temp(m), U1(m), U2(m)
real*8:: t,  dt, dx, tSum

tSum = 0.d0     ! 记录时间

U(1:m) = u0(1:m)    ! 赋初值


do while (tSum < t)
    ! 计算 U1
    CALL L(U_temp, U, dx, m)
    do i = 1, m
        U1(i) = U(i) + dt*U_temp(i)
    enddo
    ! 计算 U2
    CALL L(U_temp, U1, dx, m)
    do i=1, m
        U2(i) = 3.d0*U(i)/4.d0 + (1.d0/4.d0)*(U1(i) + dt*U_temp(i))
    enddo 
    ! 下一时刻的 U
    CALL L(U_temp, U2, dx, m)
    do i=1, m
        U(i) = U(i)/3.d0 + (2.d0/3.d0)*(U2(i) + dt*U_temp(i))
    enddo

    U(1) = (U(1) + U(m))/2.d0   ! 周期条件
    U(m) = U(1)
    tSum = tSum + dt
enddo

! 赋值给res
res(1:m) = U(1:m)
end
