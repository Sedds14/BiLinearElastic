! Subroutines in this file:
!
!  Subroutine GetModelCount( nMod )
!  Subroutine GetModelName ( iMod , ModelName )
!  Subroutine GetParamCount( iMod , nParam )
!  Subroutine GetParamName ( iMod , iParam, ParamName )
!  Subroutine GetParamUnit ( iMod , iParam, Units )
!  Subroutine GetStateVarCount( iMod , nVar )
!  Subroutine GetStateVarName ( iMod , iVar, Name )
!  Subroutine GetStateVarUnit ( iMod , iVar, Unit )
!
! Local:
!  Subroutine GetParamAndUnit( iMod , iParam, ParamName, Units )
!  Subroutine GetStateVarNameAndUnit( iMod , iVar, Name, Unit )



      Subroutine Get_Model_Count(nMod)
      !
      ! Return the maximum model number (iMod) in this DLL
      !
      Integer (Kind=4) nMod

      nMod = 2 ! Maximum model number (iMod) in current DLL

      Return
      End ! GetModelCount

      Subroutine Get_Model_Name( iMod , ModelName )
      !
      ! Return the name of the different models
      !
      Integer  iMod
      Character (Len= * ) ModelName


      Select Case (iMod)
        Case (1)
          ModelName = 'DP'
        Case (2)
          ModelName = 'MC'
        Case Default
          ModelName = 'not in DLL'
      End Select

      Return
      End ! Get_Model_Name

      Subroutine Get_Param_Count( iMod , nParam )
      !
      ! Return the number of parameters of the different models
      !

      Select Case (iMod)
        Case ( 1 )
          nParam = 5
        Case ( 2 )
          nParam = 6
        Case Default
          nParam = 0
      End Select
      Return
      End ! Get_Param_Count

      Subroutine GetParamAndUnit( iMod , iParam, ParamName, Units )
      !
      ! Return the parameters name and units of the different models
      !
      ! Units: use F for force unit
      !            L for length unit
      !            T for time unit
      !
      Character (Len=255) ParamName, Units, tName
      Select Case (iMod)
        Case (1)
          ! ModName = 'DP'
          Select Case (iParam)
            Case (1)
              ParamName = 'G'     ; Units     = 'F/L^2#'
            Case (2)
              ParamName = 'nu'    ; Units     = '-'
            Case (3)
              ParamName = 'Alpha' ; Units     = '-'
            Case (4)
              ParamName = 'Beta'  ; Units     = '-'
            Case (5)
              ParamName = 'C'     ; Units     = 'F/L^2#'
            Case Default
              ParamName = '???'   ; Units     = '???'
          End Select
        Case (2)
          ! ModName = 'MC'
          Select Case (iParam)
            Case (1)
              ParamName = 'G'     ; Units     = 'F/L^2#'
            Case (2)
              ParamName = 'nu'    ; Units     = '-'
            Case (3)
              ParamName = 'C'     ; Units     = 'F/L^2#'
            Case (4)
              ParamName = 'Phi'   ; Units     = 'deg'
            Case (5)
              ParamName = 'Psi'   ; Units     = 'deg'
            Case (6)
              ParamName = 'Tens'  ; Units     = 'F/L^2#'
            Case Default
              ParamName = '???'   ; Units     = '???'
          End Select
        Case Default
          ! model not in DLL
          ParamName = ' N/A '     ; Units     = ' N/A '
      End Select

      Return
      End ! GetParamAndUnit


      Subroutine Get_StateVar_Count( iMod , nVar )
      !
      ! Return the number of state variables of the different models
      !

      Select Case (iMod)
      Case (1)
        nVar = 0
      Case (2)
        nVar = 0
      Case Default
        nVar = 0
      End Select

      Return
      End


      Subroutine GetStateVarNameAndUnit( iMod , iVar, Name, Unit )
      !
      ! Return the name and unit of the different state variables of the different models
      !
      Character (Len=255) Name, Unit

      Select Case (iMod)
      Case (1)
        Select Case (iVar)
          Case (1)
            Name = 'var_1#'         ; Unit = '-'
          Case (2)
            Name = 'var_2#'         ; Unit = '-'
          Case Default
            Name='N/A'              ; Unit = '?'
        End Select
      Case Default
        Name='N/A'                  ; Unit = '?'
      End Select

      Return
      End



