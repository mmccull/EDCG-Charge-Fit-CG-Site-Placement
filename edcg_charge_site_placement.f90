!!requires lapack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!  Modules !!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module inputData

        integer deltaStep

endmodule inputData

module openmp

        integer np !number of threads

endmodule openmp

module minData

        integer saSteps ! number of simulated annealing steps
        integer sdSteps ! number of steepest descent steps
        integer nSets   ! number of copies

        real, allocatable :: minChi2(:,:)
        integer, allocatable :: minBoundaryRes(:,:)
        real, allocatable :: minCgCharges(:,:)

endmodule minData

module atomData

        integer nAtoms
        real, allocatable :: atomPos(:,:,:)
        real, allocatable :: atomCharges(:)
        real, allocatable :: atomMasses(:)
        real, allocatable :: caPos(:,:,:)
        real, allocatable :: caAvgPos(:,:)
        real atomEspVar
        real edcgVar
        real edcgAvg

        integer nRes                               ! number of residues in psf file
        integer, allocatable :: atomResNumber(:)   ! residue numbers directly from psf file
        integer, allocatable :: resNumber(:)       ! residue numbers renumbered to be sequential
        integer, allocatable :: resAtomStart(:)    ! first atom number in each residue.  To be used for boundary atoms
        character*4, allocatable :: atomName(:)

        character*80 atomPsfFile
        character*80 atomDcdFile

        integer nSteps

endmodule atomData


module cgData

        integer nCg
        real, allocatable :: cgPos(:,:,:)
        real, allocatable :: cgCharges(:)
        real, allocatable :: pcaMat(:,:)
        character*80 cgDcdFile

endmodule cgData

module output

        character*80 outputFile

endmodule output

module timing

        real (kind=8) ti,t1,t2
        real (kind=8) readTime,computeTime,accumTime, fitTime

endmodule timing

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!  Main Program !!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program cg_charge_fit
        use atomData
        use cgData
        use output
        use timing
        use openmp
        use minData
        implicit none
        real, parameter :: thresh = 1E-4
        real, allocatable :: A(:,:)
        real, allocatable :: B(:,:)
        real, allocatable :: C(:,:)
        integer, allocatable :: boundaryRes(:)
        real chargeChi2
        real edcgChi2
        real totalChi2
        real lambda
        real diff
        real (kind=8) omp_get_wtime
        integer iter
        integer set
        real check
        real tempFactor
        real deltaTemp
        real overallMin
        integer minSet
        integer i, j, cg, step
        integer res1, res2
        integer accept

        ti = omp_get_wtime()

        ! obtain user input from command line
        call parse_command_line(atomPsfFile,atomDcdFile,cgDcdFile,outputFile,nCg,lambda)
        open(25,file=outputFile)
        nSets=1 
        saSteps=1   
        sdSteps=1  

        ! read psf file to obtain nAtoms etc
        call read_psf_file

        ! allocate arrays for charge fitting
        allocate(A(nCg,nCg),B(nCg,nAtoms),C(nAtoms,nAtoms),cgCharges(nCg))
        ! allocate PCA matrix for EDCG
        allocate(pcaMat(nRes,nRes))
        ! allocate boundary residue arrays.  These will be updated and minimized wrt
        allocate(boundaryRes(nCg-1),minBoundaryRes(nCg-1,nSets))

        ! read the trajectory and compute the atomPos*atomPos matrix, C, for charge fitting
        call read_atom_trajectory(C)

        !Compute PCA matrix for EDCG procedure
        call compute_pca_matrix(caPos,nRes,nSteps,caAvgPos,pcaMat,nCg,edcgVar,edcgAvg)
        print*, "EDCG Variance:", edcgVar, "EDCG Average:", edcgAvg

        !allocate CG arrays
        allocate(cgPos(nCg,3,nSteps),minCgCharges(nCg,nSets),minChi2(nSets,3))

        write(25,'("Number of Sets:",i6)') nSets
        write(25,'("Number of Simulated Annealing Steps:",i6)') saSteps
        write(25,'("Number of Steepest Descent Steps:",i6)') sdSteps

        do set=1,nSets

                tempFactor = 0.001
                deltaTemp = (tempFactor)/real(saSteps+1)
                accept = 0
                write(25,'("Starting simulated annealing for set ",i5)') set
                do iter=1,saSteps

                        !generate boundary atoms
                        call gen_rand_ordered_seq(nCg-1,1,nRes-1,boundaryRes)

                        !make CG trajectory
                        call create_CG_traj(atomPos,atomMasses,nAtoms,resAtomStart,nRes,cgPos,nCg,nSteps,boundaryRes)

                        !Accumulate A and B
                        call compute_A_B_matrices(atomPos,nAtoms,cgPos,nCg,A,B,nSteps)

                        !Fit charges and compute charge residual
                        call fit_charges(A, B, C, atomCharges, cgCharges, nAtoms, nCg,chargeChi2)
                        chargeChi2 = chargeChi2/real(nSteps)/atomEspVar
                        print*, "Charge residual:", chargeChi2

                        !compute EDCG residual
                        call edcg_residual(pcaMat,nRes,boundaryRes,nCg,edcgChi2)
                        print*, "EDCG reisdual:", edcgChi2
                        edcgChi2=edcgChi2/edcgVar
                        print*, "EDCG reisdual:", edcgChi2
                        !compute total residual
                        totalChi2 = lambda*chargeChi2 + (1-lambda)*edcgChi2

                        call random_number(check)
                        if (iter==1 .or. totalChi2 < minChi2(set,1) .or. exp(-(chargeChi2-minChi2(set,1))/tempFactor)>=check) then   ! need to change this to be simulated annealing
                                minBoundaryRes(:,set) = boundaryRes
                                minCgCharges(:,set) = cgCharges
                                minChi2(set,1)=totalChi2
                                minChi2(set,2)=chargeChi2
                                minChi2(set,3)=edcgChi2
                                accept = accept+1
                        endif
                        tempFactor = tempFactor-deltaTemp

                enddo

                write(25,'("    Min Total Residual:",f20.5," Charge residual:", f20.5, " EDCG residual:",f20.5 )') (minChi2(set,j),j=1,3)
                write(25,'("    Acceptance percentage:",f8.3)') accept/real(saSteps)*100
                write(25,'("    Starting steepest descent for set ",i5)') set
                flush(25)
                flush(6)
                do iter=1,sdSteps
                        diff = 0
                        do cg=1,nCg-1
                                !move one boundary atom
                                boundaryRes = minBoundaryRes(:,set)
                                if ((boundaryRes(cg)+1) < nRes) then ! .and. (cg==(nCg-1) .or. (boundaryRes(cg)+1)<boundaryRes(cg+1))) then
                                        boundaryRes(cg) = boundaryRes(cg)+1

                                        !update the CG trajectory to represent the change in boundary atoms
                                        call create_CG_traj(atomPos,atomMasses,nAtoms,resAtomStart,nRes,cgPos,nCg,nSteps,boundaryRes)

                                        !Accumulate A and B
                                        call compute_A_B_matrices(atomPos,nAtoms,cgPos,nCg,A,B,nSteps)

                                        !Fit charges
                                        call fit_charges(A, B, C, atomCharges, cgCharges, nAtoms, nCg,chargeChi2)
                                        chargeChi2 = chargeChi2/real(nSteps)/atomEspVar
                                        print*, "Charge residual:", chargeChi2

                                        !compute EDCG residual
                                        call edcg_residual(pcaMat,nRes,boundaryRes,nCg,edcgChi2)
                                        edcgChi2=edcgChi2/edcgVar
                                        print*, "EDCG residual:", edcgChi2

                                        !compute total residual
                                        totalChi2 = lambda*chargeChi2 + (1-lambda)*edcgChi2

                                        if (totalChi2 < minChi2(set,1)) then
                                                diff = diff + minChi2(set,1)-totalChi2
                                                minBoundaryRes(:,set) = boundaryRes
                                                minCgCharges(:,set) = cgCharges
                                                minChi2(set,1)=totalChi2
                                                minChi2(set,2)=chargeChi2
                                                minChi2(set,3)=EdcgChi2
                                        endif
                                        !reset boundaryRes
                                        boundaryRes(cg) = boundaryRes(cg)-1
                                endif
                                !move one boundary atom backward
                                if (boundaryRes(cg) > 1 .and. (cg==1 .or. (boundaryRes(cg)-1)>boundaryRes(cg-1))) then
                                        boundaryRes(cg) = boundaryRes(cg)-1

                                        !make CG trajectory
                                        call create_CG_traj(atomPos,atomMasses,nAtoms,resAtomStart,nRes,cgPos,nCg,nSteps,boundaryRes)

                                        !Accumulate A and B
                                        call compute_A_B_matrices(atomPos,nAtoms,cgPos,nCg,A,B,nSteps)

                                        !Fit charges
                                        call fit_charges(A, B, C, atomCharges, cgCharges, nAtoms, nCg,chargeChi2)
                                        chargeChi2 = chargeChi2/real(nSteps)/atomEspVar
                                        print*, "Charge residual:", chargeChi2

                                        !compute EDCG residual
                                        call edcg_residual(pcaMat,nRes,boundaryRes,nCg,edcgChi2)
                                        edcgChi2=(edcgChi2-edcgAvg)**2/edcgVar
                                        print*, "EDCG residual:", edcgChi2

                                        !compute total residual
                                        totalChi2 = lambda*chargeChi2 + (1-lambda)*edcgChi2

                                        if (totalChi2 < minChi2(set,1)) then 
                                                diff = diff + minChi2(set,1)-totalChi2
                                                minBoundaryRes(:,set) = boundaryRes
                                                minCgCharges(:,set) = cgCharges
                                                minChi2(set,1)=totalChi2
                                                minChi2(set,2)=chargeChi2
                                                minChi2(set,3)=EdcgChi2
                                        endif
                                endif !forward step

                        enddo  !boundary atoms

                        if (diff<thresh) exit ! converged

                enddo   ! steepest descent

                write(*,'("Min Residual:",f20.5," Charge Residual:", f20.5," EDCG Residual:", f20.5)') (minChi2(set,j),j=1,3)
                write(25,'("     Min Residual:",f20.5," Charge Residual:", f20.5," EDCG Residual:", f20.5)') (minChi2(set,j),j=1,3)
                write(25,'("    Total CG charge:",f12.3," Total atomic charge:",f12.3)') sum(minCgCharges(:,set)),sum(atomCharges)
                write(25,'("    Min CG Charges:")')
                do i=1,nCg
                        write(25,'(f10.5)') minCgCharges(i,set)
                enddo
                write(25,'("    Min CG Boundary Residues:")')
                do i=1,nCg-1
                        write(25,'(i10)') minBoundaryRes(i,set)
                enddo

                if (set==1 .or. minChi2(set,1)<overallMin) then
                        overallMin = minChi2(set,1)
                        minSet = set
                endif
                flush(25)
                flush(6)
        enddo       ! sets
        write(25,'("Overall Minimum Residual:",f20.5, "Charge:", f20.5, "EDCG:", f20.5)') (minChi2(minSet,j),j=1,3)
        close(25)


        !Now to print out min info from Sets
        open(31,file="set_min_bas.dat")
                do i=1,nCg-1
                        write(31,'(i10)') minBoundaryRes(i,minSet)
                enddo
        close(31)
        open(32,file="set_min_charges.dat")
                do i=1,nCg
                        write(32,'(f10.5)') minCgCharges(i,minSet)
                enddo
        close(32)
        open(35,file="minCgtraj.xyz")
                boundaryRes = minBoundaryRes(:,minSet)
                call create_CG_traj(atomPos,atomMasses,nAtoms,resAtomStart,nRes,cgPos,nCg,nSteps,boundaryRes)
                do step=1,nSteps
                        write(35,'(i5)') nCg
                        write(35,'("CG coordinates from minimized charge/EDCG procedure")')
                        do cg=1,nCg
                                write(35,'("C  ", 3f8.3)') (cgPos(cg,j,step),j=1,3)
                        enddo
                enddo

        close(35)



        t2 = omp_get_wtime()
        write(*,'("Total time elapsed:",f8.3)') t2-ti

endprogram cg_charge_fit



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! Subroutines  !!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!subroutine to create CG traj from a list of boundary residues.  It requires the atomic coordinates, residue start atoms, boundary residues and number
!of CG sites as input.  It spits out the center of mass CG coordinates into cgPos.
subroutine create_CG_traj(atomPos,atomMasses,nAtoms,resAtomStart,nRes,cgPos,nCg,nSteps,boundaryRes)
        implicit none
        integer nAtoms
        integer nSteps
        integer nCg
        integer nRes
        real atomPos(nAtoms,3,nSteps)
        real atomMasses(nAtoms)
        integer resAtomStart(nRes)
        real cgPos(nCg,3,nSteps)
        integer boundaryRes(nCg-1)
        integer startAtom, stopAtom
        integer cg, atom, step
        real temp(3)
        real mass
        integer j


        do step=1,nSteps
                startAtom=1
                do cg=1,nCg
                
                        temp=0
                        mass=0
                        if (cg<nCg) then
                                stopAtom = resAtomStart(boundaryRes(cg)+1)-1
                        else
                                stopAtom = nAtoms
                        endif

                        do atom=startAtom,stopAtom
                                temp = temp + atomMasses(atom)*atomPos(atom,:,step)
                                mass = mass + atomMasses(atom)
                        enddo
                        cgPos(cg,:,step) = temp/mass
                        if (cg<nCg) then
                                startAtom = resAtomStart(boundaryRes(cg)+1)
                        endif
                enddo
        enddo

endsubroutine create_CG_traj

!subroutine to update CG traj due to a change of one boundary residue
!of CG sites as input.  It spits out the center of mass CG coordinates into cgPos.
subroutine update_CG_traj(atomPos,atomMasses,nAtoms,resAtomStart,nRes,cgPos,nCg,nSteps,boundaryRes,updatedBR)
        implicit none
        integer nAtoms
        integer nSteps
        integer nCg
        integer nRes
        integer updatedBR
        real atomPos(nAtoms,3,nSteps)
        real atomMasses(nAtoms)
        integer resAtomStart(nRes)
        real cgPos(nCg,3,nSteps)
        integer boundaryRes(nCg-1)
        integer startAtom, stopAtom
        integer cg, atom, step
        real temp(3)
        real mass
        integer j


        do step=1,nSteps
                startAtom=1
                do cg=updatedBR,updatedBR+1
                
                        temp=0
                        mass=0
                        if (cg<nCg) then
                                stopAtom = resAtomStart(boundaryRes(cg)+1)-1
                        else
                                stopAtom = nAtoms
                        endif

                        do atom=startAtom,stopAtom
                                temp = temp + atomMasses(atom)*atomPos(atom,:,step)
                                mass = mass + atomMasses(atom)
                        enddo
                        cgPos(cg,:,step) = temp/mass
                        if (cg<nCg) then
                                startAtom = resAtomStart(boundaryRes(cg)+1)
                        endif
                enddo
        enddo

endsubroutine update_CG_traj



!generate an ordered random sequence of integers numSeq long between integers minSeq and maxSeq
subroutine gen_rand_ordered_seq(numSeq,minSeq,maxSeq,seq)
        implicit none
        integer numSeq
        integer minSeq
        integer maxSeq
        integer seq(numSeq)
        real temp
        integer tempInt
        integer i, j, k
        logical diff

        do i=1,numSeq
                diff = .false.
                ! genererate a random integer between minSeq and maxSeq and make sure we don't have it already
                do while (diff.eqv..false.) 
                        call random_number(temp)
                        tempInt = int(temp*(maxSeq-minSeq))+minSeq
                        diff = .true.
                        if (i.gt.1) then
                                do j=1,i-1
                                        if (tempInt==seq(j)) then
                                                diff = .false.
                                                exit 
                                        endif
                                enddo
                        endif
                enddo
                seq(i) = tempInt
                ! if we have more than one element we need to sort them in ascending order
                if (i.gt.1) then
                        do j=1,i-1
                                !check to see if the current element is less than element j 
                                if (seq(i)<seq(j)) then
                                        !shift the entire array from element j to i circularly 1 to the right
                                        seq(j:i)=cshift(seq(j:i),-1)
                                        exit
                                endif
                        enddo
                endif
        enddo
       
endsubroutine gen_rand_ordered_seq

subroutine read_atom_trajectory(C)
        use inputData
        use cgData
        use atomData
        use openmp
        use timing
        implicit none
        integer step
        integer atom1, atom2, j
        real C(nAtoms,nAtoms)
        real CStep(nAtoms,nAtoms)
        real atomChargesM(nAtoms,1)
        real, allocatable :: atomEsp(:)
        real atomEspAvg
        real atomEspAvg2
        real tempMat(1,1)
        real dist, temp
        real caAvg(3)
        integer caCount


        C=0
        atomEspAvg = 0
        atomEspAvg2 = 0
        atomChargesM(:,1) = atomCharges
        ! Read atom dcd header and grab nAtoms, nSteps
        call read_dcd_header(atomDcdFile,nAtoms,nSteps,20)
        nSteps=3 !! MM temp
        allocate(atomPos(nAtoms,3,nSteps),atomEsp(nSteps),caPos(nRes,3,nSteps),caAvgPos(nRes,3))

        caAvgPos = 0

        do step=1,nSteps
                CStep=0
                call read_dcd_step(atomPos(:,:,step),nAtoms,20)
                !Atom-atom distance for later
!!                !$omp parallel private(atom1,atom2,j,dist,temp) SHARED(atomPos,C) num_threads(np)
!!                !$omp do schedule(dynamic,10)
                caCount = 1
                caAvg = 0
                do atom1 = 1, nAtoms-1
                        if (atomName(atom1) .eq. "CA  ") then
                                caPos(caCount,:,step) = atomPos(atom1,:,step)
                                caAvg = caAvg + atomPos(atom1,:,step)
                                caCount = caCount + 1
                        endif
                        do atom2 = atom1+1,nAtoms
                                dist = 0
                                do j=1,3
                                        temp = atomPos(atom1,j,step)-atomPos(atom2,j,step)
                                        dist = dist + temp*temp
                                enddo
                                dist = sqrt(dist)
                                CStep(atom1,atom2) = CStep(atom1,atom2)-dist
                                ! symmetrize the matrix
                                CStep(atom2,atom1) = CStep(atom1,atom2)
                        enddo
                enddo
                !compute atom ESP for this step
                tempMat = matmul(transpose(atomChargesM),matmul(CStep,atomChargesM))
                atomEsp(step) = tempMat(1,1)
                print*, "Step:", step, "Atom ESP:", atomEsp(step)
                atomEspAvg = atomEspAvg + atomEsp(step)
                atomEspAvg2 = atomEspAvg2 + atomEsp(step)*atomEsp(step)
                !add to C matrix 
                C = C + CStep
                !average CA position and center values
                caAvg = caAvg/real(nRes)
                do j=1,3
                        caPos(:,j,step) = caPos(:,j,step) - caAvg(j)
                enddo
                caAvgPos = caAvgPos(:,:) + caPos(:,:,step)

        enddo

        caAvgPos = caAvgPos/real(nSteps)

        close(20)
        !finish the average 
        atomEspAvg = atomEspAvg/real(nSteps)
        atomEspAvg2 = atomEspAvg2/real(nSteps)
        ! compute the variance
        atomEspVar = atomEspAvg2-atomEspAvg*atomEspAvg
        print*, "Charge Normalization factor:", atomEspVar

endsubroutine read_atom_trajectory


!read command line information
subroutine parse_command_line(atomPsfFile,atomDcdFile,cgDcdFile,outputFile,nCg,lambda)
        use inputData
        use openmp
        implicit none
        integer i
        character*30 arg
        character*80 atomPsfFile
        character*80 atomDcdFile
        character*80 cgDcdFile
        character*80 outputFile
        logical deltaStepFlag
        logical atomPsfFlag
        logical atomDcdFlag
        logical cgDcdFlag
        logical nCgFlag
        logical npFlag
        logical outputFlag
        logical lambdaFlag
        integer nCg
        real lambda

        atomPsfFlag = .false.
        atomDcdFlag = .false.
        cgDcdFlag = .false.
        nCgFlag = .false.
        npFlag = .false.
        deltaStepFlag = .false.
        outputFlag=.false.
        lambdaFlag=.false.
        i=1
        do 
   
                call get_command_argument(i, arg) 
   
                select case (arg) 
   
                case ('-ncg')
                        i = i+1
                        call get_command_argument(i,arg)
                        read(arg,'(i10)') nCg
                        nCgFlag=.true.
                case ('-np')
                        i = i+1
                        call get_command_argument(i,arg)
                        read(arg,'(i10)') np
                        npFlag=.true.
                case ('-l')
                        i = i+1
                        call get_command_argument(i,arg)
                        read(arg,'(f8.3)') lambda
                        lambdaFlag=.true.
                case ('-psf')
                        i = i+1
                        call get_command_argument(i,atomPsfFile)
                        atomPsfFlag=.true.
                case ('-o')
                        i = i+1
                        call get_command_argument(i,outputFile)
                        outputFlag=.true.
                case ('-dcd')
                        i = i+1
                        call get_command_argument(i,atomDcdFile)
                        atomDcdFlag=.true.
                case ('-step') 
                        i = i+1
                        call get_command_argument(i,arg)
                        read(arg,'(i10)') deltaStep
                        deltaStepFlag = .true.
                case default 
                        print '(a,a,/)', 'Unrecognized command-line option: ', arg 
                        print*, 'Usage: fit_charges.x -ncg [# cg sites] -psf [atom psf file] -dcd [atom dcd file] -step [delta step size] -o [outputfile name] -np [# of threads]'
                        stop 
                end select 
                i = i+1
                if (i.ge.command_argument_count()) exit
        enddo

        if (nCgFlag.eqv..false.) then
                write(*,'("Must provide number of CG sites using command line argument -ncg [number of CG sites]")')
                stop
        endif
        if (atomPsfFlag.eqv..false.) then
                write(*,'("Must provide a psf file using command line argument -psf [psf file name]")')
                stop
        endif
        if (atomDcdFlag.eqv..false.) then
                write(*,'("Must provide a atom dcd file using command line argument -adcd [atom dcd file name]")')
                stop
        endif
        if (outputFlag.eqv..false.) then
                write(*,'("No output file name provided.  Please provide an output file name with the following command line argument: -o [output file name]")')
                stop
        endif
        if (deltaStepFlag.eqv..false.) then
                write(*,'("Using default step size of 1.  Change this with command line argument -step [delta step size]")')
                deltaStep=1
        endif
        if (npFlag.eqv..false.) then
                write(*,'("Using a default charge residual weight of 0.5 (equial charge and EDCG weightt).  Change with command line: -l [charge residual weight]")')
                lambda = 0.5
        endif
        if (npFlag.eqv..false.) then
                write(*,'("Running on 1 thread.  Change this with command line argument -np [number of threads]")')
                np = 1
        else
                write(*,'("Running on",i3," thread(s)")') np
        endif

endsubroutine parse_command_line


!Read atomic charges from some file
subroutine read_psf_file
        use atomData
        implicit none
        integer atom, j
        character*6 check           !character to check if NATOM is in the line
        character*8 numChar         !character to read number of atoms.  must be converted to an integer
        character*4 atomCheck
        character*24 posChar
        integer ios
        integer resCount

        !open the psf file
        open(10,file=atomPsfFile)

        !run through the header of the file and look for the number of atoms
        do 

                read(10,'(a8,2x,a6)') numChar, check

                !if we read the number of atoms exit this do loop
                if (check.eq.'NATOM ') then
                        !this converts the character natoms_char to the integer natoms
                        read(numChar,*) nAtoms
                        write(*,*) "Number of atoms=", nAtoms
                        !Now that we know the number of atoms we must allocate the arrays
                        allocate(atomCharges(nAtoms),atomResNumber(nAtoms),resNumber(nAtoms),atomMasses(nAtoms),atomName(nAtoms))
                        !Now we loop through the number of atoms and read the pertinent information
                        resCount=0
                        do atom=1,nAtoms

                                read(10,'(14x,i4,6x,a4,6x,f10.6,4x,f10.4)') atomResNumber(atom),atomName(atom),atomCharges(atom),atomMasses(atom)
                                ! make residue numbers sequential
                                if (atomResNumber(atom) .ne. atomResNumber(atom-1)) then
                                        resCount = resCount + 1
                                endif
                                resNumber(atom) = resCount

                        enddo
                elseif (check.eq.'NBOND:') then
                        exit
                endif

        enddo

        close(10)

        nRes = resCount
        print*, "Number of residues:", nRes
        !Create array of residue start atom numbers
        allocate(resAtomStart(nRes))
        resCount=1
        do atom=1,nAtoms

                if (atom==1 .or. resNumber(atom) .ne. resNumber(atom-1)) then
                        resAtomStart(resCount) = atom
                        resCount = resCount +1
                endif

        enddo

        write(*,*) "Total Charge of the system = ", sum(atomCharges)

endsubroutine read_psf_file


!requires lapack
subroutine compute_A_B_matrices(atomPos,nAtoms,cgPos,nCg,A,B,nSteps)
        use openmp
        implicit none
        integer nAtoms
        integer nCg
        integer nSteps
        real atomPos(nAtoms,3,nSteps)
        real cgPos(nCg,3,nSteps)
        real A(nCg,nCg)
        real B(nCg,nAtoms)
        real dist, temp
        !loop indeces
        integer cgSite1, cgSite2
        integer j
        integer atom1, atom2
        integer step

        A=0
        B=0

        !Compute the distance between the CG sites
        do step=1,nSteps
                do cgSite1 = 1, nCg-1
                        do cgSite2 = cgSite1+1,nCg
                                dist = 0
                                do j=1,3
                                        temp = cgPos(cgSite1,j,step)-cgPos(cgSite2,j,step)
                                        dist = dist + temp*temp
                                enddo
                                dist = sqrt(dist)
                                A(cgSite1,cgSite2) = A(cgSite1,cgSite2)-dist
                                !symmetrize the matrix
                                A(cgSite2,cgSite1) = A(cgSite1,cgSite2)
                        enddo
                enddo

                !Compute the distance between cg sites and atoms
                do cgSite1 = 1, nCg
                        do atom1 = 1,nAtoms
                                dist = 0
                                do j=1,3
                                        temp = cgPos(cgSite1,j,step)-atomPos(atom1,j,step)
                                        dist = dist + temp*temp
                                enddo
                                B(cgSite1,atom1) = B(cgSite1,atom1)-sqrt(dist)
                        enddo
                enddo
        enddo

endsubroutine compute_A_B_matrices

!requires lapack
!updates A and B matrices due to change of one boundary residue
subroutine update_A_B_matrices(atomPos,nAtoms,cgPos,nCg,A,B,nSteps,updatedBR)
        use openmp
        implicit none
        integer nAtoms
        integer nCg
        integer nSteps
        integer updatedBR
        real atomPos(nAtoms,3,nSteps)
        real cgPos(nCg,3,nSteps)
        real A(nCg,nCg)
        real B(nCg,nAtoms)
        real dist, temp
        !loop indeces
        integer cgSite1, cgSite2
        integer j
        integer atom1, atom2
        integer step

        !zero the appropriate rows and columns
        A(updatedBR,:)=0
        A(updatedBR+1,:)=0
        A(:,updatedBR)=0
        A(:,updatedBR+1)=0
        B(updatedBR,:)=0
        B(updatedBR+1,:)=0
        !Compute the distance between the CG sites
        do step=1,nSteps
                do cgSite1 = updatedBR, updatedBR+1
                        do cgSite2 = 1,nCg
                                if (cgSite1 .ne. cgSite2) then
                                        dist = 0
                                        do j=1,3
                                                temp = cgPos(cgSite1,j,step)-cgPos(cgSite2,j,step)
                                                dist = dist + temp*temp
                                        enddo
                                        dist = sqrt(dist)
                                        A(cgSite1,cgSite2) = A(cgSite1,cgSite2)-dist
                                        !symmetrize the matrix
                                        A(cgSite2,cgSite1) = A(cgSite1,cgSite2)
                                endif
                        enddo
                enddo

                do cgSite1 = updatedBR, updatedBR+1
                        do atom1 = 1,nAtoms
                                dist = 0
                                do j=1,3
                                        temp = cgPos(cgSite1,j,step)-atomPos(atom1,j,step)
                                        dist = dist + temp*temp
                                enddo
                                B(cgSite1,atom1) = B(cgSite1,atom1)-sqrt(dist)
                        enddo
                enddo
        enddo

endsubroutine update_A_B_matrices


!subroutine to compute CG charges from input matrices A and B.  The residual sum of squares is calculated with aid of all-atom matrix C
subroutine fit_charges(A, B, C, atomCharges, cgCharges, nAtoms, nCg, rss)
        implicit none
        integer nAtoms
        integer nCg
        real A(nCg,nCg)
        real (kind=8) ATemp(nCg,nCg)
        real D(nCg-1,nCg)
        real B(nCg,nAtoms)
        real C(nAtoms,nAtoms)
        real BTemp(nCg,nAtoms)
        real cgCharges(nCg,1)
        real atomCharges(nAtoms)
        real atomChargesM(nAtoms,1)
        real (kind=8) newB(nCg)
        real temp(1,1)
        real rss
        !lapack routine variables
        real (kind=8) workQuery(1)
        real (kind=8), allocatable :: work(:)
        integer info
        integer lwork
        integer j, i, k

        !First we need to modify A and B to have the correct matrix properties 
        !create D matrices
        D=0
        do i=1,nCg-1
                D(i,i)=1
                D(i,i+1)=-1
        enddo 
        !multiply A by D0 and B by D1 giving new matrices A and B the correct behavior
        ATemp(1:(nCg-1),:) = matmul(dble(D),dble(A))
        BTemp(1:(nCg-1),:) = matmul(D,B)
        !generate new matrices with last line having 1.0s forcing charge conservation
        ATemp(nCg,:) = 1.0
        BTemp(nCg,:) = 1.0

        !Now use lapack routine to solve least squares problem A*cgCharges=B*atomCharges
        newB = matmul(dble(BTemp),dble(atomCharges))
        call dgels("N",nCg,nCg,1,ATemp,nCg,newB,nCg,workQuery,-1,info)
        lwork = int(workQuery(1))
        allocate(work(lwork))
        call dgels("N",nCg,nCg,1,ATemp,nCg,newB,nCg,work,lwork,info)
        deallocate(work)
        
        cgCharges(:,1) = real(newB(1:nCg))
        atomChargesM(:,1) = atomCharges

        temp = matmul(transpose(atomChargesM),matmul(C,atomChargesM))+matmul(transpose(cgCharges),matmul(A,cgCharges))-2*matmul(transpose(cgCharges),matmul(B,atomChargesM))
        rss = temp(1,1)

!        write(*,'(a17,f12.3," Total atomic charge:",f12.3)') "Total CG Charge:",sum(cgCharges),sum(atomCharges)
!        write(*,'("Residual Sum of Squares:",f8.3,f8.3)') newB(nCg+1), rss(1,1)


endsubroutine fit_charges

!subroutine to compute edcg residual given the pca matrix and boundary atoms
subroutine edcg_residual(pcaMat,nRes,boundaryRes,nCg,edcgChi2)
        implicit none
        integer nRes
        real pcaMat(nRes*3,nRes*3)
        integer nCg
        integer boundaryRes(nCg-1)
        real edcgChi2

        integer res1,res2
        integer cg
        integer startRes, stopRes
        integer index1, index2
        integer i, j

        edcgChi2=0
        startRes=1
        do cg=1,nCg
                if (cg<nCg) then
                        stopRes = boundaryRes(cg)
                else 
                        stopRes = nRes
                endif
                do res1=startRes,stopRes-1
                        do i=1,3
                                index1 = (res1-1)*3+i
                                do res2 = res1,stopRes
                                        do j=1,3
                                                index2 = (res2-1)*3+j
                                                edcgChi2 = edcgChi2 + pcaMat(index1,index1)-2.0*pcaMat(index1,index2)+pcaMat(index2,index2)
                                        enddo
                                enddo        
                        enddo
                enddo 

                startRes = stopRes+1
        enddo

endsubroutine edcg_residual

