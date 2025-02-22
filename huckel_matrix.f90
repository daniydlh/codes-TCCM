program read_matrix
    implicit none

    integer :: dim, i, j
    character(len=20) :: type
    real*8, allocatable :: matrix(:,:), eigen_values(:), eigen_vectors(:,:)
    real*8 :: beta1, beta2
    integer :: len, wide, chains, chain_len
    ! Open the input file
    open(3, file='input.inp')

    ! Read matrix type and dimension
    read(3, *) type
    read(3, *) beta1
    read(3, *) beta2
    read(3, *) len
    read(3, *) wide
    close(3)

    !dimension = total number of atoms

    chains = len + 1
    chain_len = wide * 2
    dim = chains * chain_len

    ! Allocate the matrix based on the dimension
    allocate(matrix(dim, dim),eigen_values(dim),eigen_vectors(dim,dim))
    matrix = 0.0
       
    ! Print the type and dimension
    print*, 'Matrix type:', type
    print*, 'Matrix dimension:', dim
    print*, 'Beta1:', beta1
    print*, 'Beta2:', beta2


    ! Construct the matrix
    if (type .eq. 'open') then
        write(*, '(A)') 'Open matrix'
        ! Fill the open matrix (no loop closure)
        do i = 1, dim
            do j = 1, dim
                if (j .eq. i+1) then
                    if (mod(i,2) .eq. 0) then
                                matrix(i, j) = beta2
                                matrix(j, i) = beta2
                    end if

                    if (mod(i,2) .eq. 1) then
                                matrix(i, j) = beta1
                                matrix(j, i) = beta1
                    end if
                    
                end if
            end do
        end do

        elseif (type .eq. 'closed') then
        write(*, '(A)') 'Closed matrix'
        ! Fill the closed matrix (with loop closure)
        do i = 1, dim
            do j = 1, dim
                if (j .eq. i+1) then
                        if (mod(i,2) .eq. 0) then
                                matrix(i, j) = beta2
                                matrix(j, i) = beta2
                        end if
                        if (mod(i,2) .eq. 1) then
                                matrix(i, j) = beta1
                                matrix(j, i) = beta1
                        end if

                end if
            end do
        end do

        ! Add connections to close the loop
        
        if (mod(dim,2) .eq. 0) then
                matrix(1, dim) = beta2
                matrix(dim, 1) = beta2
        else
                matrix(1, dim) = beta2
                matrix(dim, 1) = beta1
        end if 
        
        


    elseif (type .eq. 'tube') then
        write(*, '(A)') 'Closed matrix'
        ! Fill the closed matrix
        !do k = 0, chains
                do i = 1, dim
                        do j = 1, dim
                                if (i-j .eq. chain_len .and. mod(i,2) /= 0 .and. mod(j,2) /= 0) then
                                       matrix(i,j) = beta1
                                       matrix(j,i) = beta1
                                end if

                                if (j-i .eq. -1) then
                                                if (mod(j, chain_len) .eq. 0 .and. mod(i, chain_len) .eq. 1) then

                                                        matrix(i, j) = 0
                                                        matrix(j, i) = 0
                                                        
                                                else 
                                                        matrix(i, j) = beta1
                                                        matrix(j, i) = beta1
                                                end if
                                end if

                        end do
                end do
         !end do


        

    else
        write(*, '(A)') 'You did not choose the correct type'
        stop
    end if

    write(*, '(A)') 'Matrix:'
    do i = 1, dim
        write(*, '(F6.2)', advance="no") matrix(i, 1)
        do j = 2, dim
            write(*, '(F6.2)', advance="no") matrix(i, j)
        end do
        print * 
    end do

    eigen_vectors(:,:) = matrix(:,:)
    call diagonalize_matrix(dim,eigen_vectors,eigen_values)
        open (2,file='eigen_values')
        
        do i = 1,dim
                write (2,'(I4,2x,f16.8)') i, eigen_values(i)
        end do
        close(2)

        open(3,file='eigen_vectors.out')
        do i =1,dim
                write(3,'(I4,2x,1000f16.8)')i, eigen_vectors(i,:)
        end do
        close(3)
end program read_matrix

subroutine diagonalize_matrix(N,A,e)

      ! Diagonalize a square matrix
    
      implicit none
    
      ! Input variables
    
      integer,intent(in)            :: N
      double precision,intent(inout):: A(N,N)
      double precision,intent(out)  :: e(N)
    
      ! Local variables
    
      integer                       :: lwork,info
      integer                       :: i
      double precision,allocatable  :: work(:)
      
    
      ! Memory allocation
    
      allocate(work(3*N))
      lwork = size(work)
    
      call dsyev('V','U',N,A,N,e,work,lwork,info)
     
      if(info /= 0) then 
        write(*,'(a)') 'Problem in diagonalize_matrix (dsyev)!!'
        stop
      endif
      
      do i = 1 , N
        if (abs(e(i)) < 1e-10) e(i) = 0
      end do  

end subroutine diagonalize_matrix
