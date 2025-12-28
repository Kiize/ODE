program simulation
    real*8, dimension(6) :: b4, b5

    ! Coefficients for runge_kutta order 4th and 5th

    b4(1) = 37.0d0/378.0d0
    b4(2) = 0.0d0
    b4(3) = 250.0d0/621.0d0
    b4(4) = 125.0d0/594.0d0
    b4(5) = 0.0d0
    b4(6) = 512.0d0/1771.0d0

    b5(1) = 2825.0d0/27648.0d0
    b5(2) = 0.0d0
    b5(3) = 18575.0d0/48384.0d0
    b5(4) = 13525.0d0/55296.0d0
    b5(5) = 277.0d0/14336.0d0
    b5(6) = 1.0d0/4.0d0
    
end program simulation

subroutine runge_kutta(y, m, t, h, b)
    real*8, dimension(m) :: y
    real*8 :: t, h, tmp
    real*8, dimension(6) :: rk, b, c
    real*8, dimension(6, 6) :: a

    ! Coefficients for time part evolution

    c(1) = 0.0d0
    c(2) = 1.0d0/5.0d0
    c(3) = 3.0d0/10.0d0
    c(4) = 3.0d0/5.0d0
    c(5) = 1.0d0
    c(6) = 7.0d0/8.0d0

    ! Coefficients for spatial part evolution

    ! Initialize a to 0

    do i = 1, 6
        do j = 1, 6
            a(i, j) = 0.0d0
        enddo   
    enddo    

    a(2, 1) = 1.0d0/5.0d0
    a(3, 1) = 3.0d0/40.0d0
    a(4, 1) = 3.0d0/10.0d0
    a(5, 1) = -11.0d0/54.0d0
    a(6, 1) = 1631.0d0/55296.0d0

    a(3, 2) = 9.0d0/40.0d0
    a(4, 2) = -9.0d0/10.0d0
    a(5, 2) = 5.0d0/2.0d0
    a(6, 2) = 175.0d0/512.0d0

    a(4, 3) = 6.0d0/5.0d0
    a(5, 3) = -70.0d0/27.0d0
    a(6, 3) = 575.0d0/13824.0d0

    a(5, 4) = 35.0d0/27.0d0
    a(6, 4) = 44275.0d0/110592.0d0

    a(6, 5) = 253.0d0/4096.0d0




    do i = 1, 6
        tmp = 0.0d0

        do j = 1, i - 1
            tmp = tmp + a(i, j) * k(j)
        enddo

        rk(i) = f(t + c(i) * h, y + h * tmp)
    enddo

    y = y + h * (b * rk)

end subroutine

subroutine rk5()

end subroutine
