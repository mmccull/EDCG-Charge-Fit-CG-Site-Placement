


subroutine compute_pca_matrix(coord,nAtoms,nSteps,refAvg,pcaMat,nCg,edcgVar,edcgAvg)
        implicit none
        integer nAtoms
        integer nSteps
        integer nCg
        real coord(nAtoms,3,nSteps)
        real refAvg(nAtoms,3)
        real pcaMat(nAtoms,nAtoms)
        real edcgVar
        real edcgAvg

        ! Align the coordinates to their average (done iteratively)
        call align_to_avg(coord,nAtoms,nSteps,refAvg)

        ! Compute the covariance matrix
        call compute_covar(coord,nAtoms,nSteps,refAvg,pcaMat,nCg,edcgVar,edcgAvg)

endsubroutine compute_pca_matrix


subroutine compute_covar(coord,nAtoms,nSteps,avgCoord,pcaMat,nCg,edcgVar,edcgAvg)
        implicit none
        integer nAtoms
        integer nSteps
        integer nCg
        real eigenvalues(3*nAtoms)
        real coord(nAtoms,3,nSteps)
        real avgCoord(nAtoms,3)
        real pcaMat(nAtoms,nAtoms)
        real tempPcaMat(nAtoms,nAtoms)
        real temp
        integer i, j, k, l
        integer atom1, atom2
        integer step
        integer index1, index2
!        real, allocatable :: eigenvalues(:)     ! will contain eigenvalues of covariance matrix after DSYEVD subroutine
        real, allocatable :: work(:)            ! work array for DSYEVD subroutine
        integer info                                     ! output integer from DSYEVD which describes if routine completed successfully
        integer lwork                                    ! size of double work array for DSYEVD to be computed below (depends on nAtoms)
        real edcgAvg
        real edcgAvg2
        real edcgVar
        real val

        ! zero pca matrices and summed averages
        pcaMat = 0
        tempPcaMat = 0
        edcgAvg = 0
        edcgAvg2 = 0

        !compute covariance matrix
        do step=1,nSteps
                ! zero the summed variable
                temp = 0
                ! first compute the diagonal elements
                do i=1,nAtoms
                        do j=1,3
                                val  = (coord(i,j,step)-avgCoord(i,j))*(coord(i,j,step)-avgCoord(i,j))
                                tempPcaMat(i,i) = tempPcaMat(i,i) + val*val
                        enddo
                enddo
                do i=1,nAtoms-1
                        do k=i+1,nAtoms
                                do j=1,3
                                        val = (coord(i,j,step)-avgCoord(i,j))*(coord(k,l,step)-avgCoord(k,j))
                                        tempPcaMat(i,k) = tempPcaMat(i,k)+ val*val
                                enddo
                                temp = temp + tempPcaMat(i,k)
                        enddo
                enddo
                temp = temp/real(3*nCg)
!                print*, "Step", step, "PCA sum:", temp
                edcgAvg = edcgAvg+temp 
                edcgAvg2 = edcgAvg2+temp*temp
                pcaMat = pcaMat + tempPcaMat
        enddo

        ! Finalize averages and variance
        edcgAvg = edcgAvg / real(nSteps)
        edcgAvg2 = edcgAvg2 / real(nSteps)
        edcgVar = edcgAvg2-edcgAvg*edcgAvg 

        !time average and finalize covariance
        do i=1,nAtoms
                do k=i,nAtoms
                        pcaMat(i,k) = pcaMat(i,k)/real(nSteps)
                enddo
        enddo

        !diagonalize covariance                  
        !run lapack routine to compute the eigenvalues and eigenvectors of a double precision real symmetric matrix.  Entered only the upper half (hence 'U') of
        !the symmetric matrix.
        !set work array sizes (equations from lapack website http://www.netlib.org/lapack/double/ssyevd.f)
!        lwork = 1 + 6*3*nAtoms + 2*(3*nAtoms)**2
        !allocate arrays for lapack subroutine
!        allocate(work(lwork))
!        call ssyev( 'V', 'U', 3*nAtoms, tempPcaMat, 3*nAtoms, eigenvalues, work, lwork, info)

        !now compile the pcaMat from 3*nCg-6 modes from tempPcaMat
!        do i=3*nAtoms,3*nAtoms-(3*nCg-5),-1
!                do atom1 = 1, nAtoms
!                        do atom2 = atom1, nAtoms
!                                do j=1,3
!                                        index1 = (atom1-1)*3+j
!                                        index2 = (atom2-1)*3+j
!                                        pcaMat(atom1,atom2) = pcaMat(atom1,atom2)+tempPcaMat(index1,i)*eigenvalues(i)*tempPcaMat(index2,i)
!                                enddo
!                        enddo
!                enddo
!        enddo

endsubroutine compute_covar

subroutine align_to_avg(coord,nAtoms,nSteps,refAvg)
        implicit none
        integer, parameter :: maxIter = 100
        real, parameter :: thresh = 1E-5
        integer nAtoms
        integer nSteps
        real coord(nAtoms,3,nSteps)
        real refAvg(nAtoms,3)
        real newAvg(nAtoms,3)
        real residual
        integer iter
        integer atom
        integer step
        integer j
        
        do iter = 1, maxIter
                newAvg = 0
                do step=1,nsteps

                        !orient to reference using Kabsch algorithm
                        call orient_to_reference(coord(:,:,step),nAtoms,refAvg)

                        !average aligned coordinates
                        newAvg = newAvg + coord(:,:,step)

                enddo

                !Time average the positions of each atom
                newAvg = newAvg/real(nSteps) 

                !residual is just the RMSD of the average from this step and the last
                call compute_residual(newAvg,nAtoms,refAvg,residual)

                !update reference
                refAvg = newAvg

                !check convergence
                if (residual < thresh) then
                        exit
                endif
        
        enddo

endsubroutine align_to_avg

! This subroutine calculates the optimal rotation matrix to rotate the new set of coordinates to the reference set using the Kabsch algorithm
! This subroutine assumes that refCoord and alpha_pos have their COM at the origin 
subroutine orient_to_reference(coord,nAtoms,refCoord)
   implicit none
   integer, parameter :: lwork=99
   integer nAtoms
   real coord(nAtoms,3)
   real refCoord(nAtoms,3)
   integer i, j
   real  rot(3,3)
   real  covar(3,3)
   real  work(lwork)
   real  singular_values(3)
   real  right_eigv_T(3,3)
   real  left_eigv(3,3)
   real  right_handed(3,3)
   real  det, sign_det
   integer info


   !The first step in the Kabsch algorithm is to compute a 3x3 covariance matrix between new set of coordinates and reference set of coordinates
   !compute the covariance matrix between this time point and the reference
   covar = matmul(transpose(coord),refCoord)

   !compute the determinant of covar
   call determinant_three_three(covar,det)

   !We now perform a singular value decomposition of the covariance matrix using the lapack subroutine
   call sgesvd('a','a',3,3,covar,3,singular_values,left_eigv,3,right_eigv_T,3,work,lwork,info)

   !This is to make sure we are using a right handed coordinate system
   if (det<0) then
      sign_det=-1
   else
      sign_det=1
   endif

   do i=1, 3
      do j=1, 3
         if (i.eq.3 .and. j.eq.3 ) then
            right_handed(i,j)=sign_det
         elseif (i.eq.j) then
            right_handed(i,j)=1
         else
            right_handed(i,j)=0
         endif
      enddo
   enddo

   !Can now compute the rotation matrix
   left_eigv = matmul(left_eigv,right_handed)
   rot = matmul(left_eigv,right_eigv_T)

   !Rotate the coordinates:
   coord = matmul(coord,rot)

endsubroutine orient_to_reference

!calculate the determinant of a three by three matrix
subroutine determinant_three_three(mat,det)
   implicit none
   real  mat(3,3)
   real  det

   det = mat(1,1)*mat(2,2)*mat(3,3) - mat(1,1)*mat(2,3)*mat(3,2) - mat(1,2)*mat(2,1)*mat(3,3) + mat(1,2)*mat(2,3)*mat(3,1)+&
   mat(1,3)*mat(2,1)*mat(3,2) - mat(1,3)*mat(2,2)*mat(3,1)

endsubroutine determinant_three_three


subroutine compute_residual(coord,nAtoms,refCoord,residual)
   implicit none
   integer nAtoms
   real coord(nAtoms,3)
   real refCoord(nAtoms,3)
   real residual
   real temp
   integer i
   integer j

   residual = 0
   do i=1,nAtoms

      do j=1,3
         temp = coord(i,j)-refCoord(i,j)
         residual = residual + temp*temp
      enddo

   enddo

   residual = sqrt(residual/dble(nAtoms))

endsubroutine compute_residual   

