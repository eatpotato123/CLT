PROGRAM parallel_01
 
    USE omp_lib
    IMPLICIT NONE
 
    INTEGER :: i,j
    INTEGER(4) :: time_begin, time_end, time_rate
    REAL, DIMENSION(1:50,1:50) :: f, g
    REAL :: k
    WRITE(*,*) '开始进行串行计算'
!>@ 1、通过串行计算获得两个矩阵的初始化计算
    CALL system_clock(time_begin,time_rate)
    DO i = 1, 50
        DO j = 1, 50
            f(i,j) = i*j
            k = k + 1
        END DO
    END DO
 
    DO i = 1, 50
        DO j = 1, 50
            g(i,j) = i*j + 1
            k = k + 1
        END DO
    END DO
 
    CALL system_clock(time_end,time_rate)
    WRITE(*,*) 'The value of k after serial computing is: ', k
    WRITE(*,*) 'The time wasted on serial computing is: ',(time_end - time_begin)/time_rate
 
    WRITE(*,*)
    WRITE(*,*)
    WRITE(*,*) '开始进行第一类串行计算—SECTIONS'
 
!>@ 2、通过块并行计算获得两个矩阵的初始化计算
    k = 0 ! 重新初始化k的值 
    CALL system_clock(time_begin,time_rate)
    CALL omp_set_num_threads(2)
    !$omp parallel
    !$omp sections private(i,j,k)
        !$omp section
             DO i = 1, 50
                DO j = 1, 50
                    f(i,j) = i*j
                    k = k + 1
                END DO
            END DO
            WRITE(*,*) 'The value of k after parallel computing is: ', k,', and it comes from the thread of ',omp_get_thread_num()
  
        !$omp section
            DO i = 1, 50
                DO j = 1, 50
                    g(i,j) = i*j + 1
                    k = k + 1
                END DO
            END DO
            WRITE(*,*) 'The value of k after parallel computing is: ', k,', and it comes from the thread of ',omp_get_thread_num()
    !$omp end sections
    !$omp end parallel
 
    CALL system_clock(time_end,time_rate)
    WRITE(*,*) 'The time wasted on the first class parallel computing is: ',(time_end - time_begin)/time_rate
    WRITE(*,*)
    WRITE(*,*)
    WRITE(*,*) '开始进行第二类并行计算—DO'
 
!>@ 3、通过DO循环实现两个矩阵的初始化计算
    k = 0 ! 重新初始化k的值
    CALL system_clock(time_begin,time_rate) 
    !$omp parallel private(k,j,i)
    !$omp do
        DO i = 1, 50
            DO j = 1, 50
                f(i,j) = i*j
                k = k + 1
                ! 去掉注释后，可现实每一次循环所在的线程ID
 !               WRITE(*,*) 'The value of k after parallel computing is: ', k,', and it comes from the thread of ',omp_get_thread_num()
            END DO
        END DO
    !$omp end do
    !$omp end parallel
 
    !$omp parallel private(k,j,i)
    !$omp do
        DO i = 1, 50
            DO j = 1, 50
                g(i,j) = i*j
                k = k + 1
!                WRITE(*,*) 'The value of k after parallel computing is: ', k,', and it comes from the thread of ',omp_get_thread_num()
            END DO
        END DO
    !$omp end do
    !$omp end parallel
    CALL system_clock(time_end,time_rate)
    WRITE(*,*) 'The time wasted on the first class parallel computing is: ',(time_end - time_begin)/time_rate
END PROGRAM