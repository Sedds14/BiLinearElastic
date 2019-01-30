      Subroutine User_Mod ( IDTask, iMod, IsUndr,
     *                      iStep, iTer, iEl, Int,
     *                      X, Y, Z,
     *                      Time0, dTime,
     *                      Props, Sig0, Swp0, StVar0,
     *                      dEps, D, BulkW,
     *                      Sig, Swp, StVar, ipl,
     *                      nStat, NonSym, iStrsDep, iTimeDep,iTang,
     *                      iPrjDir, iPrjLen, iAbort )
!
! Purpose: User supplied soil model
!          Example: iMod=1 : Drucker-Prager
!                   iMod=2 : Mohr-Coulomb
!
!  Depending on IDTask, 1 : Initialize state variables
!                       2 : calculate stresses,
!                       3 : calculate material stiffness matrix
!                       4 : return number of state variables
!                       5 : inquire matrix properties
!                           return switch for non-symmetric D-matrix
!                           stress/time dependent matrix
!                       6 : calculate elastic material stiffness matrix
! Arguments:
!          I/O  Type
!  IDTask   I   I    : see above
!  iMod     I   I    : model number (1..10)
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
!  NonSym   O   I    : Non-Symmetric D-matrix ?
!  iStrsDep O   I    : =1 for stress dependent D-matrix
!  iTimeDep O   I    : =1 for time dependent D-matrix
!  iTang    O   I    : =1 for tangent matrix
!  iAbort   O   I    : =1 to force stopping of calculation
!
      Implicit Double Precision (A-H, O-Z)
!
      Dimension Props(*), Sig0(*), StVar0(*), dEps(*), D(6,6),
     *          Sig(*),   StVar(*), iPrjDir(*)

      Data iounit / 0 /
      Save iounit
!
!---  Local variables
!
      Character*100 BaseName

      !DEC$ ATTRIBUTES DLLExport, StdCall, reference :: User_Mod

      BaseName = 'example'
      ! Possibly open a file for debugging purposes
      If (iounit.Eq.0) Then

        Call Open_Dbg_File( iPrjDir, iPrjLen, BaseName )

        Write(1,*)'File 1 opened: ', Trim(baseName)
        ! maybe write some more info on version to debug file ?
        write(1,*)'Compiled : ',__DATE__,' ', __TIME__
!DEC$ IF DEFINED(_X86_)
        ! this 32-bit ??
        Write(1,1050)'IF32',__INTEL_COMPILER
!DEC$ ELSE
        ! this 64-bit ??
        Write(1,1050)'IF64', __INTEL_COMPILER
!DEC$ ENDIF
 1050   format ( 1X,A,1x,I0 )
        iounit = 1
        Call WriVec(1,'Props',Props,50)
        Call Flush(1)
      End If

      Call WriIvl( -1, 'iounit',iounit )
      Call WriIvl( -1, 'IDTask',IDTask )
      Select Case (iMod)
        Case (1)   ! DP
          Call MyModel1( IDTask, iMod, IsUndr, iStep, iTer, iEl, Int,
     *                   X, Y, Z, Time0, dTime,
     *                   Props, Sig0, Swp0, StVar0,
     *                   dEps, D, BulkW, Sig, Swp, StVar, ipl,
     *                   nStat, NonSym, iStrsDep, iTimeDep, iTang,
     *                   iAbort )
        Case (2)   ! MC+Tension cut-off
          Call MyMod_MC( IDTask, iMod, IsUndr, iStep, iTer, iEl, Int,
     *                   X, Y, Z, Time0, dTime,
     *                   Props, Sig0, Swp0, StVar0,
     *                   dEps, D, BulkW, Sig, Swp, StVar, ipl,
     *                   nStat, NonSym, iStrsDep, iTimeDep, iTang,
     *                   iAbort )
!        Case (3)
!          Call MyModel3( IDTask, ....
        Case Default
          Write(1,*) 'invalid model number in UsrMod', iMod
          Write(1,*) 'IDTask: ',IDTask
          Stop 'invalid model number in UsrMod'
          iAbort=1
          Return
      End Select ! iMod
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
      End ! User_Mod

! **********************************************************************

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
          Dbg_Name=PrjDir(:iPrjLen)//'data.'//
     *      trim(BaseName)//char(48+nErr)//'.rr0'

          If (nErr.Lt.10) Goto 1
        End If

      Open( Unit= 1, File= Dbg_Name,blocksize=4096)

      End Subroutine Open_Dbg_File

! **********************************************************************

!
! Interfaces to routines to report parameter names, counts etc.
!

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

      Subroutine GetModelCount( nMod )
      !
      ! Return the maximum model number (nMod) in this DLL
      !
      Implicit None
      Integer (Kind=4) nMod

      !DEC$ ATTRIBUTES DLLExport, StdCall, reference  :: GetModelCount

      Call Get_Model_Count( nMod )

      Return
      End ! GetModelCount

      Subroutine GetModelName( iMod , ModelName )
      !
      ! Return the name of the different models
      !
      Implicit None
      Integer  iMod
      Character (Len= * ) ModelName
      Character (Len=255) tName
      !DEC$ ATTRIBUTES DLLExport, StdCall, reference  :: GetModelName

      Call Get_Model_Name( iMod , ModelName )
      Call Add_Str_Length( ModelName )

      End ! GetModelName

      Subroutine GetParamCount( iMod , nParam )
      !
      ! Return the number of parameters of the different models
      !
      Implicit None
      Integer  iMod, nParam

      !DEC$ ATTRIBUTES DLLExport, StdCall, reference  :: GetParamCount

      Call Get_Param_Count( iMod , nParam )

      End ! GetParamCount

      Subroutine GetParamName( iMod , iParam, ParamName )
      !
      ! Return the parameters name of the different models
      !
      Implicit None
      Integer  iMod, iParam
      Character (Len=255) ParamName, Units

      !DEC$ ATTRIBUTES DLLExport, StdCall, reference  :: GetParamName

      Call GetParamAndUnit(iMod,iParam,ParamName,Units)
      Call Add_Str_Length( ParamName )

      End ! GetParamName

      Subroutine GetParamUnit( iMod , iParam, Units )
      !
      ! Return the units of the different parameters of the different models
      !
      Implicit None
      Integer  iMod, iParam
      Character (Len=255) ParamName, Units

      !DEC$ ATTRIBUTES DLLExport, StdCall, reference  :: GetParamUnit

      Call GetParamAndUnit(iMod,iParam,ParamName,Units)
      Call Add_Str_Length( Units )

      End ! GetParamUnit

      Subroutine GetStateVarCount( iMod , nVar )
      !
      ! Return the number of state variables of the different models
      !
      Implicit None
      Integer  iMod, nVar

      !DEC$ ATTRIBUTES DLLExport, StdCall, reference :: GetStateVarCount

      Call Get_StateVar_Count( iMod , nVar )

      End ! GetStateVarCount 

      Subroutine GetStateVarName( iMod , iVar, Name )
      !
      ! Return the name of the different state variables
      ! of the different models
      !
      Implicit None
      Integer  iMod, iVar
      Character (Len=255) Name, Unit

      !DEC$ ATTRIBUTES DLLExport, StdCall, reference  :: GetStateVarName

      Call GetStateVarNameAndUnit( iMod , iVar, Name, Unit )
      Call Add_Str_Length( Name )

      End ! GetStateVarName

      Subroutine GetStateVarUnit( iMod , iVar, Unit )
      !
      ! Return the units of the different state variables of the different models
      !
      Implicit None
      Integer  iMod, iVar
      Character (Len=255) Name, Unit

      !DEC$ ATTRIBUTES DLLExport, StdCall, reference  :: GetStateVarUnit

      Call GetStateVarNameAndUnit( iMod , iVar, Name, Unit )
      Call Add_Str_Length( Unit )

      End ! GetStateVarUnit
