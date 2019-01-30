      Subroutine User_Mod ( IDTask, iMod, IsUndr, iStep, iTer, & ! 5
                              iEl, Int, X, Y, Z, Time0, dTime, & ! 7
                              Props, Sig0, Swp0, StVar0, dEps, & ! 5
                              D, Bulk_W, Sig, Swp, StVar, ipl, & ! 6
                              nStat, NonSym, iStrsDep, iTimeDep, & ! 4
                              iTang, iPrjDir, iPrjLen, iAbort ) ! 4 
                                                                  ! = 31
!
! Purpose: User supplied soil model
!          Example: iModel=1 : Bi-Linear Elastic
!
!  Depending on IDTask, 1 : Initialize state variables
!                       2 : calculate stresses,
!                       3 : calculate material stiffness matrix
!                       4 : return number of state variables
!                       5 : inquire matrix properties
!                           return switch for non-symmetrinoParamsD-matrix
!                           stress/time dependent matrix
!                       6 : calculate elastinoParamsmaterial stiffness matrix
!
! Arguments:
!          I/O  Type
!  IDTask   I   I    : see above
!  iModel   I   I    : model number (1..10)
!  IsUndr   I   I    : = 1 for undrained, 0 otherwise
!  iStep    I   I    : Global step number
!  iter     I   I    : Global iteration number
!  iel      I   I    : Global element number
!  Int      I   I    : Global integration point number
!  X        I   R    : X-Position of integration point
!  Y        I   R    : Y-Position of integration point
!  Z        I   R    : Z-Position of integration point
!  Time0    I   R    : Time at start of step
!  dTime    I   R    : Time increment
!  Props    I   R()  : List with model parameters
!  Sig0     I   R()  : Stresses at start of step
!  Swp0     I   R    : Excess pore pressure start of step
!  StVar0   I   R()  : State variable at start of step
!  dEps     I   R()  : Strain increment
!  D       I/O  R(,) : Material stiffness matrix
!  BulkW   I/O  R    : Bulkmodulus for water (undrained only)
!  Sig      O   R()  : Resulting stresses
!  Swp      O   R    : Resulting excess pore pressure
!  StVar    O   R()  : Resulting values state variables
!  ipl      O   I    : Plasticity indicator
!  nStat    O   I    : Number of state variables
!  NonSym   O   I    : Non-SymmetrinoParamsD-matrix ?
!  iStrsDep O   I    : =1 for stress dependent D-matrix
!  iTimeDep O   I    : =1 for time dependent D-matrix
!  iTang    O   I    : =1 for tangent matrix
!  iAbort   O   I    : =1 to force stopping of calculation

      Implicit Double Precision (A-H, O-Z)

      Dimension Props(*), Sig0(*), StVar0(*), dEps(*), D(6,6), Sig(*), StVar(*), iPrjDir(*)
 
      Data iounit / 0 /
      Save iounit
      integer i, j
      Character*100 BaseName
      Character(8)  :: date
      Character(10) :: time
      !DEC$ ATTRIBUTES DLLExport, StdCall, reference :: User_Mod
      
      BaseName = 'example'
      ! Possibly open a file for debugging purposes
      If (iounit.Eq.0) Then

        Call Open_Dbg_File( iPrjDir, iPrjLen, BaseName )

        Write(1,*)'File 1 opened: ', Trim(baseName)
        ! maybe write some more info on version to debug file ?
        !Write(1,*)'Compiled : ',__DATE__,' ', __TIME__  - Probabaly wrong for compiler

        ! using keyword arguments
        call date_and_time(date,time)
        call date_and_time(DATE=date)
        call date_and_time(TIME=time)

        Write(1,*)'Compiled : ', DATE, TIME

        1050   format ( 1X,A,1x,I0 )
        iounit = 1
        Call WriVec(1,'Props',Props,50)
        Call Flush(1)
      End If

      Call WriIvl( -1, 'iounit',iounit )
      Call WriIvl( -1, 'IDTask',IDTask )

      nStat = 0

      If (IDTask .Eq. 1) Then ! Initialize state variables StVar0
        ! Initialise state variables based on K0 conditions
        ! We could set p = grout pressure if suitable?
        p = (Sig0(1) + Sig0(2) + Sig0(3)) / 3
        StVar0(1) = Min(StVar0(1), p)
      End If  ! IDTask = 1

      If (IDTask .Eq. 2) Then ! Calculate the constitutive stresses Sig (and Swp)

        Do i=1,6
          Sig(i) = Sig0(i)
          Do j=1,6
            Sig(i) = Sig(i) + D(i,j) * dEps(j)
          End Do
        End Do  

      End If  ! IDTask = 2

      If (IDTask .Eq. 3 .Or. IDTask .Eq. 6) Then ! Create effective material stiffness/elastinoParamsstiffness matrix D[]
      ! This is where the magic happens!
      ! 1. Use stress boundary Props[1] to select E and v from user inputs Props[2-5]
        SigB = Props(1)
        E1 = Props(2)
        v1 = Props(3)
        E2 = Props(4)
        v2 = Props(5) 
      ! SigD is the deviator stress
        SigD = Sig(1) - (Sig(2)+ Sig(3)) / 2
        If (SigD .Gt. SigB) Then ! Material is stiff
          E = E2
          v = v2
        Else ! Material is less stiff
          E = E1
          v = v1
        End If

      ! Calculate Matrix Terms
        G = 0.5 * E / (1.0 + v)
        Fac= 2 * G / (1.0 - 2 * v)  ! Make sure that v < 0.5 !!!!
        Term1 = Fac* (1 - v)
        Term2 = Fac* v
        D(1,1) = Term1
        D(1,2) = Term2        
        D(1,3) = Term2
        D(2,1) = Term2
        D(2,2) = Term1
        D(2,3) = Term2
        D(3,1) = Term2
        D(3,2) = Term2
        D(3,3) = Term1
        D(4,4) = G
        D(5,5) = G
        D(6,6) = G

      End If  ! IDTask = 3 & 6

      If (IDTask .Eq. 4) Then ! Number of state parameters
        nStat    = 1
      End If  ! IDTask = 4

      If (IDTask .Eq. 5) Then ! matrix type
        NonSym   = 0  ! 1 for non-symmetrinoParamsD-matrix
        iStrsDep = 1  ! 1 for stress dependent D-matrix
        iTang    = 1  ! 1 for tangent D-matrix
        iTimeDep = 0  ! 1 for time dependent D-matrix
      End If  ! IDTask = 5

      Return

!      Case Default
!          Write(1,*) 'invalid model number in UsrMod', iMod
!          Write(1,*) 'IDTask: ',IDTask
!          Stop 'invalid model number in UsrMod'
!          iAbort=1
!          Return
!      End Select ! iMod

      If (IDTask .Eq. 5.And.iel+int.Eq.2) Then
        Write(1,*)'nStat   : ',nStat
        Write(1,*)'NonSym  : ',NonSym
        Write(1,*)'StrsDep : ',iStrsDep
        Write(1,*)'TimeDep : ',iTimeDep
        Write(1,*)'Tangent : ',iTang
      End If
      If (IDTask == -333 .And. iel+int == -1234) Then
!        Write(1,*)'IDTask: ',IDTask,' iStep,iTer',iStep,iTer
!        Call Flush(1)
      End If
      Call WriIvl( -1, 'IDTask end',IDTask )
      Return
     End ! End of model
! OTHER Subroutines in this file:

!  Subroutine Open_Dbg_File( iPrjDir, iPrjLen, BaseName )   
!  Subroutine Add_Str_Length( aString )
!  Subroutine GetModelCount( nMod )
!  Subroutine GetModelName ( iModel , ModelName )
!  Subroutine GetParamCount( iModel , nParam )
!  Subroutine GetParamName ( iModel , iParam, ParamName )
!  Subroutine GetParamUnit ( iModel , iParam, Units )
!  Subroutine GetStateVarCount( iModel , nVar )
!  Subroutine GetStateVarName ( iModel , iVar, Name )
!  Subroutine GetStateVarUnit ( iModel , iVar, Unit )
!

     Subroutine Open_Dbg_File( iPrjDir, iPrjLen, BaseName )
     Implicit None

     Integer, intent(in) :: iPrjLen, iPrjDir(*)
     Character*(*), intent(in):: BaseName

     Character*255 PrjDir, Dbg_Name

     Integer i, nErr, ios

     PrjDir=' '
     Do i=1,iPrjLen
       PrjDir(i:i) = Char( iPrjDir(i) )
     End Do

     Dbg_Name=PrjDir(:iPrjLen)//'data.'//trim(BaseName)//'.rr0'
     nErr=0
   1 Continue
       Open( Unit= 1, File= Dbg_Name,iostat=ios)
       If (ios.Eq.0) Close(Unit=1,Status='delete',iostat=ios)

       If (ios.Ne.0) Then
         !
         ! in case of error try ...udsmex1 or udsmex2 or ..
         !
         nErr=nErr+1
         Dbg_Name=PrjDir(:iPrjLen)//'data.'//     trim(BaseName)//char(48+nErr)//'.rr0'

         If (nErr.Lt.10) Goto 1
       End If

     Open( Unit= 1, File= Dbg_Name, Status= 'OLD')

     End Subroutine Open_Dbg_File

     !******************************************************************************

     Subroutine Add_Str_Length( aString )
     Implicit None
     Character*(*) aString
     Character *255 tString
     Integer Lt              ! length of incoming string

     ! routine should add the length of the string as the first character

     tString = aString
     Lt      = Len_Trim(tString)
     aString = Char(Lt) // tString(1:Lt)

     End Subroutine Add_Str_Length


      Subroutine GetModelCount(nMod)
      !
      ! Return the maximum model number (iModel) in this DLL
      !
      Implicit None
      Integer (Kind=4) nMod
      !DEC$ ATTRIBUTES DLLExport, StdCall, reference :: GetModelCount

      nMod = 1 ! Maximum model number (iModel) in current DLL

      Return
      End ! GetModelCount

      Subroutine GetModelName( iModel , ModelName )
      !
      ! Return the name of the different models
      !
      Implicit None
      Integer  iModel
      Character (Len= * ) ModelName
      Character (Len=255) tName
      !DEC$ ATTRIBUTES DLLExport, StdCall, reference :: GetModelName

      Select Case (iModel)
        Case (1)
          ModelName = "Bi-linear elastic"
        Case Default
          ModelName = "not in DLL"
      End Select

      Call Add_Str_Length( ModelName )

      Return
      End ! GetModelName

      Subroutine GetParamCount( iModel , noParams)
      !
      ! Return the number of parameters of the different models
      !
      Implicit None
      Integer iModel, noParams
      !DEC$ ATTRIBUTES DLLExport, StdCall, reference :: GetParamCount
          noParams= 5

      Return
      End ! GetParamCount

      Subroutine GetParamName (iModel, iParam, ParamName)

      Implicit None
      Integer iModel, iParam
      Character (Len=255) ParamName
      !DEC$ ATTRIBUTES DLLExport, StdCall, reference :: GetParamName

      Select Case (iModel)
        Case (1)
          Select Case (iParam)
            Case (1)
              ParamName = "@s#_B#" ! SigmaB
            Case (2)
              ParamName = "E_1#"   ! E1
            Case (3)
              ParamName = "@n#_1#" ! v1
            Case (4)
              ParamName = "E_2#"   ! E2
            Case (5)
              ParamName = "@n#_2#" ! v2
            Case Default
              ParamName = "???"   
          End Select
 
        Case Default
          ! model not in DLL
          ParamName = "N/A"     
      End Select

      Call Add_Str_Length( ParamName )

      Return
      End ! GetParamAndUnit


      Subroutine GetParamUnit( iModel , iParam, Units )
      !
      ! Return the parameters name and units of the different models
      !
      ! Units: use F for force unit
      !            L for length unit
      !            T for time unit
      !
      Implicit None
      Integer iModel, iParam
      Character (Len=255) Units
      !DEC$ ATTRIBUTES DLLExport, StdCall, reference :: GetParamUnit

      Select Case (iModel)
        Case (1)
          Select Case (iParam)
            Case (1)
              Units     = "F/L^2#"      ! SigmaB
            Case (2)
              Units   = "F/L^2#"      ! E1
            Case (3)
              Units   = "-"        ! v1
            Case (4)
              Units    = "F/L^2#"      ! E2
            Case (5)
              Units    = "-"        ! v2
            Case Default
              Units    = "???"
          End Select
 
        Case Default
          ! model not in DLL
          Units     = "N/A"
      End Select

      Call Add_Str_Length( Units )

      Return
      End ! GetParamAndUnit


      Subroutine GetStateVarCount( iModel , nVar )
      !
      ! Return the number of state variables of the different models
      !
      Implicit None
      Integer iModel, nVar
      !DEC$ ATTRIBUTES DLLExport, StdCall, reference :: GetStateVarCount
      Select Case (iModel)
      Case (1)
        nVar = 0
      Case Default
        nVar = 0
      End Select

      Return
      End


      Subroutine GetStateVarName( iModel , iVar, Name)
      !
      ! Return the name and unit of the different state variables of the different models
      !
      Implicit None
      Integer iModel, iVar
      Character (Len=255) Name
      !DEC$ ATTRIBUTES DLLExport, StdCall, reference :: GetStateVarName
      Select Case (iModel)
      Case (1)
        Select Case (iVar)
          Case Default
            Name="N/A"
        End Select
      Case Default
        Name="N/A"
      End Select

      Call Add_Str_Length ( Name )

      Return
      End

      Subroutine GetStateVarUnit( iModel , iVar, Unit )
      !
      ! Return the name and unit of the different state variables of the different models
      !
      Implicit None
      Integer iModel, iVar
      Character (Len=255) Unit
      !DEC$ ATTRIBUTES DLLExport, StdCall, reference :: GetStateVarUnit
      Select Case (iModel)
      Case (1)
        Select Case (iVar)
          Case Default
            Unit = "?"
        End Select
      Case Default
      Unit = "?"
      End Select

      Call Add_Str_Length( Unit )
      
      Return
      End

! Be painfully Lazy and C+P All USRLIB.for subroutines into this file

      Subroutine MZEROR(R,K)
!
!***********************************************************************
!
!     Function: To make a real array R with dimension K to zero
!
!***********************************************************************
!
      Implicit Double Precision (A-H,O-Z)
      Dimension R(*)

      Do J=1,K
        R(J) = 0.0D0
      End Do

      Return
      End


      Subroutine MZEROI(I,K)
!
!***********************************************************************
!
!     Function: To make an integre array I with Dimension K to zero
!
!***********************************************************************
!
      Dimension I(*)

      Do J=1,K
        I(J)=0
      End Do

      Return
      End

      Subroutine SETRVAL(R,K,V)
!
!***********************************************************************
!
!     Function: To fill a real array R with Dimension K with value V
!
!***********************************************************************
!
      Implicit Double Precision (A-H,O-Z)
      Dimension R(*)

      Do J=1,K
        R(J)=V
      End Do

      Return
      End

      Subroutine SETIVAL(I,K,IV)
!
!***********************************************************************
!
!     Function: To fill an integer array I with Dimension K with value IV
!
!***********************************************************************
!
      Implicit Double Precision (A-H,O-Z)
      Dimension I(*)

      Do J=1,K
        I(J)=IV
      End Do

      Return
      End

      Subroutine COPYIVEC(I1,I2,K)
!
!***********************************************************************
!
!     Function: To copy an integer array I1 with Dimension K to I2
!
!***********************************************************************
!
      Implicit Double Precision (A-H,O-Z)
      Dimension I1(*),I2(*)

      Do  J=1,K
        I2(J)=I1(J)
      End Do

      Return
      End

      Subroutine COPYRVEC(R1,R2,K)
!
!***********************************************************************
!
!     Function: To copy a Double array R1 with Dimension K to R2
!
!***********************************************************************
!
      Implicit Double Precision (A-H,O-Z)
      Dimension R1(*),R2(*)

      Do J=1,K
        R2(J)=R1(J)
      End Do

      Return
      End


      Logical Function IS0ARR(A,N)
!
!***********************************************************************
!    Function :  To check whether a real array contains only zero values.
!                When an array contains only zero's is might not need to be
!                written to the XXX file.
!                exit Function when first non-zero value occured or when
!                all elements are checked and are zero.
!
!    Input:  A : array to be checked
!            N : number of elements in array that should be checked
!
!    Output : .TRUE.  when all elements are 0
!             .FALSE. when at least one element is not zero
!
!    Called by :  Subroutine TOBXX
!
!***********************************************************************
!
      Implicit Double Precision (A-H,O-Z)
      Dimension A(*)
      Is0Arr=.False.
      Do I=1,N
        If ( A(I) .Ne. 0 ) Return
      End Do
      Is0Arr=.True.
      Return
      End

      Logical Function IS0IARR(IARR,N)
!
!***********************************************************************
!    Function :  To check whether a integer array contains only zero values.
!                Similar to IS0ARR
!
!    Input:  IARR : array to be checked
!            N    : number of elements in array that should be checked
!
!    Output : .TRUE.  when all elements are 0
!             .FALSE. when at least one element is not zero
!
!    Called by :  Subroutine TOBXX
!
!***********************************************************************
!
      Implicit Double Precision (A-H,O-Z)
      Dimension IARR(*)

      Is0IArr=.False.
      Do I=1,N
        If ( IARR(I) .Ne. 0 ) Return
      End Do
      Is0IArr=.True.
      Return
      End
!***********************************************************************
      Subroutine MulVec(V,F,K)
!***********************************************************************
!
!     Function: To multiply a real vector V with dimension K by F
!
!***********************************************************************
!
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION V(*)

      Do J=1,K
        V(J)=F*V(J)
      End Do

      Return
      End     ! Subroutine Mulvec
!***********************************************************************
      Subroutine MatVec(xMat,IM,Vec,N,VecR)
!***********************************************************************
!
!     Calculate VecR = xMat*Vec
!
! I   xMat  : (Square) Matrix (IM,*)
! I   Vec   : Vector
! I   N     : Number of rows/colums
! O   VecR  : Resulting vector
!
!***********************************************************************
      Implicit Double Precision (A-H,O-Z)
      Dimension xMat(IM,*),Vec(*),VecR(*)
!***********************************************************************
      Do I=1,N
        X=0
        Do J=1,N
          X=X+xMat(I,J)*Vec(J)
        End Do
        VecR(I)=X
      End Do
      Return
      End    ! Subroutine MatVec

!***********************************************************************
      Subroutine AddVec(Vec1,Vec2,R1,R2,N,VecR)
!***********************************************************************
!
!     Calculate VecR() = R1*Vec1()+R2*Vec2()
!
! I   Vec1,
! I   Vec2  : Vectors
! I   R1,R2 : Multipliers
! I   N     : Number of rows
! O   VecR  : Resulting vector
!
!***********************************************************************
      Implicit Double Precision (A-H,O-Z)
      Dimension Vec1(*),Vec2(*),VecR(*)
!***********************************************************************
      Do I=1,N
        X=R1*Vec1(I)+R2*Vec2(I)
        VecR(I)=X
      End Do
      Return
      End    ! Subroutine AddVec
!
!***********************************************************************
      Double Precision Function DInProd(A,B,N)
!***********************************************************************
!
!     Returns the Inproduct of two vectors
!
! I   A,B  : Two vectors
! I   N    : Used length of vectors
!***********************************************************************
      Implicit Double Precision (A-H,O-Z)
      Dimension A(*),B(*)
!***********************************************************************

      X = 0
      Do I=1,N
        X = X + A(I)*B(I)
      End Do
      DInProd = X
      Return
      End     ! Function DInProd
!
!***********************************************************************
      Subroutine MatMat(xMat1,Id1,xMat2,Id2,nR1,nC2,nC1,xMatR,IdR)
!***********************************************************************
!
!     Calculate xMatR = xMat1*xMat2
!
! I   xMat1 : Matrix (Id1,*)
! I   xMat2 : Matrix (Id2,*)
! I   nR1   : Number of rows in resulting matrix    (= No rows in xMat1)
! I   nC2   : Number of columns in resulting matrix (= No cols in xMat2)
! I   nC1   : Number of columns in matrix xMat1
!             = Number  rows    in matrix xMat2
! O   xMatR : Resulting matrix (IdR,*)
!
!***********************************************************************
      Implicit Double Precision (A-H,O-Z)
      Dimension xMat1(Id1,*),xMat2(Id2,*),xMatR(IdR,*)
!**********************************************************************

      Do I=1,nR1
        Do J=1,nC2
          X=0
          Do K=1,nC1
            X=X+xMat1(I,K)*xMat2(K,J)
          End Do
          xMatR(I,J)=X
        End Do
      End Do

      Return
      End     ! Subroutine MatMat

!***********************************************************************
      Subroutine MatMatSq(n, xMat1, xMat2, xMatR)
!***********************************************************************
!
!     Calculate xMatR = xMat1*xMat2 for square matrices, size n
!
! I   n     : Dimension of matrices
! I   xMat1 : Matrix (n,*)
! I   xMat2 : Matrix (n,*)
! O   xMatR : Resulting matrix (n,*)
!
!***********************************************************************
      Implicit Double Precision (A-H,O-Z)
      Dimension xMat1(n,*),xMat2(n,*),xMatR(n,*)
!**********************************************************************

      Do I=1,n
        Do J=1,n
          X=0
          Do K=1,n
            X=X+xMat1(I,K)*xMat2(K,J)
          End Do
          xMatR(I,J)=X
        End Do
      End Do

      Return
      End     ! Subroutine MatMatSq

!***********************************************************************
      Subroutine WriVal ( io, C , V )
!***********************************************************************
!
! Write (Double) value to file unit io (when io>0)
!
!***********************************************************************
!
      Implicit Double Precision (A-H,O-Z)
      Character C*(*)

      If (io.Le.0) Return

      Write(io,*) C,V
    1 Format( A,3x, 1x,1p,e12.5)
      Return
      End
!***********************************************************************
      Subroutine WriIVl ( io, C , I )
!***********************************************************************
!
! Write (integer) value to file unit io (when io>0)
!
!***********************************************************************
      Implicit Double Precision (A-H,O-Z)
      Character C*(*)

      If (io.Le.0) Return

      Write(io,*) C,I
    1 Format( A,3x, 1x,I6)
      Return
      End
!***********************************************************************
      Subroutine WriIVc ( io, C , iV , n )
!***********************************************************************
!
! Write (integer) vector to file unit io (when io>0)
!
!***********************************************************************
      Character C*(*)
      Dimension iV(*)

      If (io.Le.0) Return

      Write(io,*) C
      Write(io,1) (iv(i),i=1,n)
    1 Format( ( 2(3x,5i4) ) )
      Return
      End
!***********************************************************************
      Subroutine WriVec ( io, C , V , n )
!***********************************************************************
!
! Write (Double) vector to file unit io (when io>0)
! 6 values per line
!***********************************************************************
     Implicit Double Precision (A-H,O-Z)
      Character C*(*)
      Dimension V(*)

      If (io.Le.0) Return

      If (Len_Trim(C).Le.6) Then
        Write(io,2) C,( V(i),i=1,n)
      Else
        Write(io,*) C
        Write(io,1) ( V(i),i=1,n)
      End If
    1 Format( ( 2(1x, 3(1x,1p,e10.3) ) ) )
    2 Format( A, ( T7, 2(1x, 3(1x,1p,e10.3) ) ) )
      Return
      End
!***********************************************************************
      Subroutine WriVec5( io, C , V , n )
!***********************************************************************
!
! Write (Double) vector to file unit io (when io>0)
! 5 values per line
!***********************************************************************
      Implicit Double Precision (A-H,O-Z)
      Character C*(*)
      Dimension V(*)

      If (io.Le.0) Return

      Write(io,*) C
      Write(io,1) ( V(i),i=1,n)
    1 Format( 5(1x,1p,e12.5) )
      Return
      End
!***********************************************************************
      Subroutine WriMat ( io, C , V , nd, nr, nc )
!***********************************************************************
!
! Write (Double) matrix to file unit io (when io>0)
! 6 values per line
!***********************************************************************
      Implicit Double Precision (A-H,O-Z)
      Character C*(*)
      Dimension V(nd,*)

      If (io.Le.0) Return

      Write(io,*) C
      Do j=1,nr
        Write(io,1) j,( V(j,i),i=1,nc)
      End Do
    1 Format(i4, (  T7,2(1x, 3(1x,1p,e10.3) ) ) )
      Return
      End
!***********************************************************************

      Subroutine MatInvPiv(Aorig,B,N)
      Implicit Double Precision (A-H,O-Z)
      Dimension Aorig(n,*), B(n,*),A(:,:)
      Allocatable :: A
      Allocate ( A(n,n) ) ! No error checking !!
      Call CopyRVec(AOrig, A, n*n )
      Call MZeroR(B,n*n)
      Do i=1,n
        B(i,i) = 1d0
      End Do
      Do I=1,n
        T=A(I,I)
        iPiv=i
        Do j=i+1,n
          If ( Abs(A(j,i)) .Gt. Abs(A(iPiv,i))  ) iPiv=j
        End Do
        If (iPiv.Ne.i) Then
          Do j=1,n
            x         = A( i  ,j)
            A( i  ,j) = A(iPiv,j)
            A(iPiv,j) = x
            x         = B( i  ,j)
            B( i  ,j) = B(iPiv,j)
            B(iPiv,j) = x
          End Do
          T=A(I,I)
        End If
        Do J=1,n
          A(I,J)=A(I,J)/T
          B(I,J)=B(I,J)/T
        End Do
        Do K=1,n
          If (K.Ne.I) Then
            T=A(K,I)
            Do J=1,n
              A(K,J)=A(K,J)-T*A(I,J)
              B(K,J)=B(K,J)-T*B(I,J)
            End Do
          End If
        End Do
      End Do
      DeAllocate ( A  )
      Return
      End ! MatinvPiv

!***********************************************************************
      Subroutine PrnSig(IOpt,S,xN1,xN2,xN3,S1,S2,S3,P,Q)
      Implicit Double Precision (A-H,O-Z)
      Dimension S(*),xN1(*),xN2(*),xN3(*)

      If (iOpt.Eq.1) Then
        Call Eig_3(0,S,xN1,xN2,xN3,S1,S2,S3,P,Q) ! with Eigenvectors
      Else
        Call Eig_3a(0,S,S1,S2,S3,P,Q) ! no Eigenvectors
      End If
      Return
      End
!***********************************************************************
      Subroutine Eig_3(iOpt,St,xN1,xN2,xN3,S1,S2,S3,P,Q)
      Implicit Double Precision (A-H,O-Z)
      Dimension St(6),A(3,3),V(3,3), xN1(3),xN2(3),xN3(3)
      !
      ! Get Eigenvalues/Eigenvectors for 3*3 matrix
      ! Wim Bomhof 15/11/'01
      ! PGB : adaption to Principal stress calculation
      !
      ! Applied on principal stresses, directions
      ! Stress vector St(): XX, YY, ZZ, XY, YZ, ZX
      !
      A(1,1) = St(1) ! xx
      A(1,2) = St(4) ! xy = yx
      A(1,3) = St(6) ! zx = xz

      A(2,1) = St(4) ! xy = yx
      A(2,2) = St(2) ! yy
      A(2,3) = St(5) ! zy = yz

      A(3,1) = St(6) ! zx = xz
      A(3,2) = St(5) ! zy = yz
      A(3,3) = St(3) ! zz

      ! Set V to unity matrix
      V(1,1) = 1
      V(2,1) = 0
      V(3,1) = 0

      V(1,2) = 0
      V(2,2) = 1
      V(3,2) = 0

      V(1,3) = 0
      V(2,3) = 0
      V(3,3) = 1


      abs_max_s=0.0
      Do i=1,3
        Do j=1,3
          if (abs(a(i,j)) .Gt. abs_max_s) abs_max_s=abs(a(i,j))
        End Do
      End Do
      Tol = 1d-20 * abs_max_s
      it = 0
      itmax = 50
      Do While ( it.Lt.itMax .And.abs(a(1,2))+abs(a(2,3))+abs(a(1,3)) .Gt. Tol )
        it=it+1
        Do k=1,3
          If (k .Eq. 1) Then
            ip=1
            iq=2
          Else If (k .Eq.2) Then
            ip=2
            iq=3
          Else
            ip=1
            iq=3
          End If
!          If (a(ip,iq) .Ne. 0.0) Then
          If (Abs(a(ip,iq)) .Gt. 1d-50) Then         ! ongelijk nul ?
            tau=(a(iq,iq)-a(ip,ip))/(2.0*a(ip,iq))
            If (tau .Ge.0.0) Then
              sign_tau=1.0
            Else
              sign_tau=-1.0
            End If
            t=sign_tau/(abs(tau)+sqrt(1.0+tau*tau))
            c=1.0/sqrt(1.0+t*t)
            s=t*c
            a1p=c*a(1,ip)-s*a(1,iq)
            a2p=c*a(2,ip)-s*a(2,iq)
            a3p=c*a(3,ip)-s*a(3,iq)
            a(1,iq)=s*a(1,ip)+c*a(1,iq)
            a(2,iq)=s*a(2,ip)+c*a(2,iq)
            a(3,iq)=s*a(3,ip)+c*a(3,iq)
            a(1,ip)=a1p
            a(2,ip)=a2p
            a(3,ip)=a3p

            v1p=c*v(1,ip)-s*v(1,iq)
            v2p=c*v(2,ip)-s*v(2,iq)
            v3p=c*v(3,ip)-s*v(3,iq)
            v(1,iq)=s*v(1,ip)+c*v(1,iq)
            v(2,iq)=s*v(2,ip)+c*v(2,iq)
            v(3,iq)=s*v(3,ip)+c*v(3,iq)
            v(1,ip)=v1p
            v(2,ip)=v2p
            v(3,ip)=v3p

            ap1=c*a(ip,1)-s*a(iq,1)
            ap2=c*a(ip,2)-s*a(iq,2)
            ap3=c*a(ip,3)-s*a(iq,3)
            a(iq,1)=s*a(ip,1)+c*a(iq,1)
            a(iq,2)=s*a(ip,2)+c*a(iq,2)
            a(iq,3)=s*a(ip,3)+c*a(iq,3)
            a(ip,1)=ap1
            a(ip,2)=ap2
            a(ip,3)=ap3
          End If ! a(ip,iq)<>0
        End Do ! k
      End Do ! While
      ! principal values on diagonal of a
      S1 = a(1,1)
      S2 = a(2,2)
      S3 = a(3,3)
      ! Derived invariants
      P = (S1+S2+S3)/3
      Q = Sqrt( ( (S1-S2)**2 + (S2-S3)**2 + (S3-S1)**2 ) / 2 )

      ! Sort eigenvalues S1 <= S2 <= S3
      is1 = 1
      is2 = 2
      is3 = 3
      if (s1.Gt.s2) Then
        t   = s2
        s2  = s1
        s1  = t
        it  = is2
        is2 = is1
        is1 = it
      End If
      if (s2.Gt.s3) Then
        t   = s3
        s3  = s2
        s2  = t
        it  = is3
        is3 = is2
        is2 = it
      End If
      if (s1.Gt.s2) Then
        t   = s2
        s2  = s1
        s1  = t
        it  = is2
        is2 = is1
        is1 = it
      End If
      Do i=1,3
        xN1(i) = v(i,is1) ! first  column
        xN2(i) = v(i,is2) ! second column
        xN3(i) = v(i,is3) ! third  column
      End Do
      Return
      End ! Eig_3

      Subroutine Eig_3a(iOpt,St,S1,S2,S3,P,Q) ! xN1,xN2,xN3,
      Implicit Double Precision (A-H,O-Z)
      Dimension St(6),A(3,3)   !  V(3,3),xN1(3),xN2(3),xN3(3)
      !
      ! Get Eigenvalues ( no Eigenvectors) for 3*3 matrix
      ! Wim Bomhof 15/11/'01
      !
      ! Applied on principal stresses, directions
      ! Stress vector XX, YY, ZZ, XY, YZ, ZX
      !
      A(1,1) = St(1) ! xx
      A(1,2) = St(4) ! xy = yx
      A(1,3) = St(6) ! zx = xz

      A(2,1) = St(4) ! xy = yx
      A(2,2) = St(2) ! yy
      A(2,3) = St(5) ! zy = yz

      A(3,1) = St(6) ! zx = xz
      A(3,2) = St(5) ! zy = yz
      A(3,3) = St(3) ! zz

      abs_max_s=0.0
      Do i=1,3
        Do j=1,3
          if (abs(a(i,j)) .Gt. abs_max_s) abs_max_s=abs(a(i,j))
        End Do
      End Do
      Tol = 1d-20 * abs_max_s
      If (iOpt.Eq.1) Tol = 1d-50*abs_max_s
      it=0
      itmax = 50
      Do While ( it.lt.itmax .And.abs(a(1,2))+abs(a(2,3))+abs(a(1,3)) .Gt. Tol )

        it=it+1
        Do k=1,3
          If (k .Eq. 1) Then
            ip=1
            iq=2
          Else If (k .Eq.2) Then
            ip=2
            iq=3
          Else
            ip=1
            iq=3
          End If
!          If (a(ip,iq) .Ne. 0.0) Then         ! ongelijk nul ?
          If (Abs(a(ip,iq)) .Gt. 1d-50) Then         ! ongelijk nul ?
            tau=(a(iq,iq)-a(ip,ip))/(2.0*a(ip,iq))
            If (tau .Ge.0.0) Then
              sign_tau=1.0
            Else
              sign_tau=-1.0
            End If
            t=sign_tau/(abs(tau)+sqrt(1.0+tau*tau))
            c=1.0/sqrt(1.0+t*t)
            s=t*c
            a1p=c*a(1,ip)-s*a(1,iq)
            a2p=c*a(2,ip)-s*a(2,iq)
            a3p=c*a(3,ip)-s*a(3,iq)
            a(1,iq)=s*a(1,ip)+c*a(1,iq)
            a(2,iq)=s*a(2,ip)+c*a(2,iq)
            a(3,iq)=s*a(3,ip)+c*a(3,iq)
            a(1,ip)=a1p
            a(2,ip)=a2p
            a(3,ip)=a3p

            ap1=c*a(ip,1)-s*a(iq,1)
            ap2=c*a(ip,2)-s*a(iq,2)
            ap3=c*a(ip,3)-s*a(iq,3)
            a(iq,1)=s*a(ip,1)+c*a(iq,1)
            a(iq,2)=s*a(ip,2)+c*a(iq,2)
            a(iq,3)=s*a(ip,3)+c*a(iq,3)
            a(ip,1)=ap1
            a(ip,2)=ap2
            a(ip,3)=ap3
          End If ! a(ip,iq)<>0
        End Do ! k
      End Do ! While
      ! principal values on diagonal of a
      S1 = a(1,1)
      S2 = a(2,2)
      S3 = a(3,3)
      ! Derived invariants
      P = (S1+S2+S3)/3
      Q = Sqrt( ( (S1-S2)**2 + (S2-S3)**2 + (S3-S1)**2 ) / 2 )

      if (s1.Gt.s2) Then
        t   = s2
        s2  = s1
        s1  = t
      End If
      if (s2.Gt.s3) Then
        t   = s3
        s3  = s2
        s2  = t
      End If
      if (s1.Gt.s2) Then
        t   = s2
        s2  = s1
        s1  = t
      End If
      Return
      End ! Eig_3a

!
!***********************************************************************
      Logical Function LEqual(A,B,Eps)
!***********************************************************************
!
!     Returns .TRUE.  when two real values are (almost) equal,
!             .FALSE. otherwise
!
! I   A,B  : Two real values to be compared
! I   Eps  : Toleration (Magnitude ~= 1E-5)
!***********************************************************************
      Implicit Double Precision (A-H,O-Z)
!***********************************************************************
      LEqual =.True.
      If (A .Eq. B) Return
      If (DAbs(A-B) .LT. 0.5D0*Eps*( DAbs(A) + DAbs(B) + Eps ) )Return
      LEqual =.False.
      Return
      End     ! function LEqual
!
!***********************************************************************
      Subroutine CrossProd(xN1,xN2,xN3)
!***********************************************************************
!
!     Returns cross product of xN1 and xN2
!
! I   xN1,xN2 : Two basic vectors
! O   xN3     : Resulting vector
!***********************************************************************
      Implicit Double Precision (A-H,O-Z)
      Dimension xN1(*),xN2(*),xN3(*)
!***********************************************************************

      xN3(1) = xN1(2)*xN2(3) - xN1(3)*xN2(2)
      xN3(2) = xN1(3)*xN2(1) - xN1(1)*xN2(3)
      xN3(3) = xN1(1)*xN2(2) - xN1(2)*xN2(1)

      Return
      End     ! Subroutine CrossProd
!
!***********************************************************************
      Double Precision Function ArcSin(X,ie)
!***********************************************************************
!
!     Returns the Arc Sine of X
!
! I   X : Input value
!
!     Note : In stead of using default routine DASIN we use this one
!            because �X� can be slightly beyond 1 and this will give
!            a RTE using DASIN(X)
!
!***********************************************************************
      Implicit Double Precision (A-H,O-Z)
!***********************************************************************
      Ie=0
      S = (1-X*X)
!      If (S .Lt. -1E-10) Ie=1
!      If (S .Lt. -1E-10) Write(*,1) X,S
!      If (S .Lt. -1E-10) Write(2,1) X,S
    1 Format(' ArcSin(',1x,1p,e13.5e3,') , S =',1x,1p,e13.5e3)
      If (S.LT.0) S = 0
      S = DSQRT(S)
      ArcSin = DATan2(X,S)
      Return
      End     ! function ArcSin
!
!***********************************************************************
      Subroutine CarSig(S1,S2,S3,xN1,xN2,xN3,SNew)
!***********************************************************************
!
!     Returns the Cartesian stresses using the principal stresses S1..S3
!     and the principal directions
!
! I   S1..S3   : Principal stresses
! I   xN1..xN3 : Principal directions (xNi for Si)
!
!***********************************************************************
      Implicit Double Precision (A-H,O-Z)
      Dimension xN1(*),xN2(*),xN3(*),SNew(*)
      Dimension SM(3,3),T(3,3),TT(3,3),STT(3,3)
!***********************************************************************
!
!**** Fill transformation (rotation) matrix
!
      Do I=1,3
        T(I,1) = xN1(I)
        T(I,2) = xN2(I)
        T(I,3) = xN3(I)
        TT(1,I) = T(I,1)
        TT(2,I) = T(I,2)
        TT(3,I) = T(I,3)
      End Do
!      Call MatTranspose(T,3,TT,3,3,3)

      Call MZeroR(SM,9)
      SM(1,1) = S1
      SM(2,2) = S2
      SM(3,3) = S3
!
!**** SMnew = T*SM*TT
!
      Call MatMat(SM ,3,  TT,3 , 3,3,3 ,STT,3)
      Call MatMat( T ,3, STT,3 , 3,3,3 ,SM ,3)
!     Call MatMatSq(3, SM,  TT, STT )   ! STT = SM*TT
!     Call MatMatSq(3,  T, STT, SM  )   ! SM  =  T*STT
!
!**** Extract cartesian stress vector from stress matrix
!
      Do I=1,3
        SNew(I) = SM(I,I)
      End Do
      SNew(4) = SM(2,1)
      SNew(5) = SM(3,2)
      SNew(6) = SM(3,1)

      Return
      End     ! Subroutine CarSig
!**********************************************************************
      subroutine setveclen(xn,n,xl)
!**********************************************************************
      Implicit Double Precision (A-H,O-Z)
      Dimension xN(*)
      x=0
      do i=1,n
        x=x+xn(i)**2
      end do
      if (x.Ne.0) Then
        f=xl/sqrt(x)
        do i=1,3
          xn(i)=xn(i)*f
        end do
      end if
      return
      end ! setveclen

!**********************************************************************
! End Of file
!**********************************************************************
