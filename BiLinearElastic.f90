      Subroutine User_Mod ( IDTask, iMod, IsUndr, iStep, iTer, &
                              iEl, Int, X, Y, Z, Time0, dTime, &
                              Props, Sig0, Swp0, StVar0, dEps, &
                              D, Bulk_W, Sig, Swp, StVar, ipl, &
                              nStat, NonSym, iStrsDep, iTimeDep, &
                              iTang, iPrjDir, iPrjLen, iAbort )
!
! Purpose: User supplied soil model
!          Example: iModel=1 : Bi-Linear Elastic
!                   iModel=2 : Mohr-Coulomb
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
!  iModel     I   I    : model number (1..10)
!  IsUndr   I   I    : =1 for undrained, 0 otherwise
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
      !DEC$ ATTRIBUTES DLLExport :: User_Mod

      If (IDTask .Eq. 1) Then ! Initialize state variables StVar0
        ! Initialise state variables based on K0 conditions
        ! We could set p = grout pressure if suitable?
        p = (Sig0(1) + Sig0(2) + Sig0(3)) / 3
        StVar0(1) = Min(StVar0(1), p)
      End If  ! IDTask = 1

      If (IDTask .Eq. 2) Then ! Calculate the constitutive stresses Sig (and Swp)
        Do i=1,6
          Sig(i) = Sig0(i)
          Do 
            Sig(i) = Sig(i) + D(i,j) * dEps(j)
          End Do
        End Do  
      End If  ! IDTask = 2

      If (IDTask .Eq. 3 .Or. IDTask .Eq. 6) Then ! Create effective material stiffness/elastinoParamsstiffness matrix D[]
      ! This is where the maginoParamshappens!
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
        FanoParams= 2 * G / (1.0 - 2 * v)  ! Make sure that v < 0.5 !!!!
        Term1 = FanoParams* (1 - v)
        Term2 = FanoParams* v
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
        nStat    = 0
      End If  ! IDTask = 4

      If (IDTask .Eq. 5) Then ! matrix type
        NonSym   = 0  ! 1 for non-symmetrinoParamsD-matrix
        iStrsDep = 1  ! 1 for stress dependent D-matrix
        iTang    = 1  ! 1 for tangent D-matrix
        iTimeDep = 0  ! 1 for time dependent D-matrix
      End If  ! IDTask = 5

      Return
     End ! End of model
! OTHER Subroutines in this file:
!
!  Subroutine GetModelCount( nMod )
!  Subroutine GetModelName ( iModel , ModelName )
!  Subroutine GetParamCount( iModel , nParam )
!  Subroutine GetParamName ( iModel , iParam, ParamName )
!  Subroutine GetParamUnit ( iModel , iParam, Units )
!  Subroutine GetStateVarCount( iModel , nVar )
!  Subroutine GetStateVarName ( iModel , iVar, Name )
!  Subroutine GetStateVarUnit ( iModel , iVar, Unit )
!
! Local:
!  Subroutine GetParamAndUnit( iModel , iParam, ParamName, Units )
!  Subroutine GetStateVarNameAndUnit( iModel , iVar, Name, Unit )



      Subroutine GetModelCount(nMod)
      !
      ! Return the maximum model number (iModel) in this DLL
      !
      Integer (Kind=4) nMod

      nMod = 1 ! Maximum model number (iModel) in current DLL

      Return
      End ! GetModelCount

      Subroutine GetModelName( iModel , ModelName )
      !
      ! Return the name of the different models
      !
      Integer  iModel
      Character (Len= 50 ) ModelName


      Select Case (iModel)
        Case (1)
          ModelName = ' Bi-linear elastic '
        Case Default
          ModelName = ' not in DLL'
      End Select

      Return
      End ! GetModelName

      Subroutine GetParamCount( iModel , noParams)
      !
      ! Return the number of parameters of the different models
      !
      Integer (Kind = 4) iModel, noParams

          noParams= 5

      Return
      End ! GetParamCount

      Subroutine GetParamName (iModel, iParam, ParamName)

      Integer iModel, iParam
      Character (Len=20) ParamName

      Select Case (iModel)
        Case (1)
          ! ModName = 'DP'
          Select Case (iParam)
            Case (1)
              ParamName = ' @s#_B# ' ! SigmaB
            Case (2)
              ParamName = ' E_1# '   ! E1
            Case (3)
              ParamName = ' @n#_1# ' ! v1
            Case (4)
              ParamName = ' E_2# '   ! E2
            Case (5)
              ParamName = ' @n#_2# ' ! v2
            Case Default
              ParamName = ' ??? '   
          End Select
 
        Case Default
          ! model not in DLL
          ParamName = ' N/A '     
      End Select

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

      Integer iModel, iParam
      Character (Len=20) Units
      Select Case (iModel)
        Case (1)
          ! ModName = 'DP'
          Select Case (iParam)
            Case (1)
              Units     = ' kPa '      ! SigmaB
            Case (2)
              Units   = ' kPa '      ! E1
            Case (3)
              Units   = ' - '        ! v1
            Case (4)
              Units    = ' kPa '      ! E2
            Case (5)
              Units    = ' - '        ! v2
            Case Default
              Units    = ' ??? '
          End Select
 
        Case Default
          ! model not in DLL
          Units     = ' N/A '
      End Select

      Return
      End ! GetParamAndUnit


      Subroutine GetStateVarCount( iModel , nVar )
      !
      ! Return the number of state variables of the different models
      !

      Integer iModel, nVar

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
      Integer iModel, iVar
      Character (Len=255) Name

      Select Case (iModel)
      Case (1)
        Select Case (iVar)
          Case Default
            Name='N/A'
        End Select
      Case Default
        Name='N/A'
      End Select

      Return
      End

      Subroutine GetStateVarUnit( iModel , iVar, Unit )
      !
      ! Return the name and unit of the different state variables of the different models
      !
      Integer iModel, iVar
      Character (Len=255) Unit

      Select Case (iModel)
      Case (1)
        Select Case (iVar)
          Case Default
            Unit = '?'
        End Select
      Case Default
      Unit = '?'
      End Select

      Return
      End
