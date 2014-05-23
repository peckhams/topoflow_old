


#***************************************************************
def Make_RTS_File_for_Qnet_SW(Qn_SW_RTS_file, RTI_file,
                              start_month, start_day, start_hour,
                              stop_month, stop_day, stop_hour,
                              timestep, time_zone,
                              T_air, T_air_file, T_air_type,
                              RH, RH_file, RH_type,
                              albedo, albedo_file, albedo_type,
                              dust_att, dust_att_file, dust_att_type,
                              factor, factor_file, factor_type,
                              slope, slope_file, slope_type,
                              aspect, aspect_file, aspect_type,
                              lon_deg, lon_file, lon_type,
                              lat_deg, lat_file, lat_type,
                              MSG_BOX_ID=None):

#-------------------------------------------------------
#Notes:  If time is before local sunrise or after local
#        sunset then Qnet should be zero.
#-------------------------------------------------------

#----------------------------
#Read info from the RTI file
#----------------------------
    n_params = 37
    _opt = (MSG_BOX_ID,)
    def _ret():
        _optrv = zip(_opt, [MSG_BOX_ID])
        _rv = [Qn_SW_RTS_file, RTI_file, start_month, start_day, start_hour, stop_month, stop_day, stop_hour, timestep, time_zone, T_air, T_air_file, T_air_type, RH, RH_file, RH_type, albedo, albedo_file, albedo_type, dust_att, dust_att_file, dust_att_type, factor, factor_file, factor_type, slope, slope_file, slope_type, aspect, aspect_file, aspect_type, lon_deg, lon_file, lon_type, lat_deg, lat_file, lat_type]
        _rv += [_o[1] for _o in _optrv if _o[0] is not None]
        return tuple(_rv)
    
    info = Read_RTI_File(RTI_file)
    nx = info.ncols
    ny = info.nrows
    
    #------------------------------
    #Create grids of lats and lons
    #------------------------------
    if (idl_func.n_elements(lon_deg) == 0):    
        lon_deg = Longitude_Grid(info)
        lat_deg = Latitude_Grid(info)
    else:    
        if (lon_type == 2):    
            lon_deg = Read_Grid(lon_file, _TYPE='FLOAT', SILENT=True)
        if (lat_type == 2):    
            lat_deg = Read_Grid(lat_file, _TYPE='FLOAT', SILENT=True)
    
    #------------------------------------
    #If data type is grid, read the grid
    #------------------------------------
    if (T_air_type == 2):    
        T_air = Read_Grid(T_air_file, _TYPE='FLOAT', SILENT=True)
    if (RH_type == 2):    
        RH = Read_Grid(RH_file, _TYPE='FLOAT', SILENT=True)
    if (albedo_type == 2):    
        albedo = Read_Grid(albedo_file, _TYPE='FLOAT', SILENT=True)
    if (dust_att_type == 2):    
        dust_att = Read_Grid(dust_att_file, _TYPE='FLOAT', SILENT=True)
    if (factor_type == 2):    
        factor = Read_Grid(factor_file, _TYPE='FLOAT', SILENT=True)
    
    #----------------------------------------
    #Read slope grid, convert to slope angle
    #NB!  RT slope grids have NaNs on edges.
    #----------------------------------------
    slopes = Read_Grid(slope_file, _TYPE='FLOAT', SILENT=True)
    beta = arctan(slopes)
    twopi = float64(2) * numpy.pi
    beta = (twopi + beta) % twopi
    #---------------------------------------
    w_nan = I2PY_w = where(ravel(isfinite(beta) == 0))[0]
    n_nan = size(I2PY_w)
    if (n_nan != 0):    
        beta[w_nan] = float64(0)
    w_bad = I2PY_w = where(ravel(logical_or((beta < float64(0)), (beta > numpy.pi / float64(2)))))[0]
    n_bad = size(I2PY_w)
    if (n_bad != 0):    
        msg = array(['ERROR:  Some slope angles are out of range.', ' '])
        result = GUI_Message(msg, INFO=True, TITLE='ERROR MESSAGE')
        return _ret()
    
    #----------------------------------------------
    #Read aspect grid; Alpha must be CW from north  ;****************
    #NB!  RT aspect grids have NaNs on edges.
    #----------------------------------------------
    aspects = Read_Grid(aspect_file, _TYPE='FLOAT', SILENT=True)
    alpha = (numpy.pi / float64(2)) - aspects
    twopi = float64(2) * numpy.pi
    alpha = (twopi + alpha) % twopi
    #----------------------------------------
    w_nan = I2PY_w = where(ravel(isfinite(alpha) == 0))[0]
    n_nan = size(I2PY_w)
    if (n_nan != 0):    
        alpha[w_nan] = float64(0)
    
    #---------------------------------
    #Open files for input variables ?
    #---------------------------------
    if logical_or((T_air_type == 1), (T_air_type == 3)):    
        T_air_unit = TF_Get_LUN(T_air_file)
        file_T_air_unit = open(T_air_file, 'rb')
        I2PY_SWAP_ENDIAN = Not_Same_Byte_Order(info.byte_order)
        if (T_air_type == 3):    
            T_air = zeros([ny, nx], dtype='Float32')
    #-----------------------------------------------------------
    if logical_or((RH_type == 1), (RH_type == 3)):    
        RH_unit = TF_Get_LUN(RH_file)
        file_RH_unit = open(RH_file, 'rb')
        I2PY_SWAP_ENDIAN = Not_Same_Byte_Order(info.byte_order)
        if (RH_type == 3):    
            RH = zeros([ny, nx], dtype='Float32')
    #-----------------------------------------------------------
    if logical_or((albedo_type == 1), (albedo_type == 3)):    
        albedo_unit = TF_Get_LUN(albedo_file)
        file_albedo_unit = open(albedo_file, 'rb')
        I2PY_SWAP_ENDIAN = Not_Same_Byte_Order(info.byte_order)
        if (albedo_type == 3):    
            albedo = zeros([ny, nx], dtype='Float32')
    #-----------------------------------------------------------
    if logical_or((dust_att_type == 1), (dust_att_type == 3)):    
        dust_att_unit = TF_Get_LUN(dust_att_file)
        file_dust_att_unit = open(dust_att_file, 'rb')
        I2PY_SWAP_ENDIAN = Not_Same_Byte_Order(info.byte_order)
        if (dust_att_type == 3):    
            dust_att = zeros([ny, nx], dtype='Float32')
    #-----------------------------------------------------------
    if logical_or((factor_type == 1), (factor_type == 3)):    
        factor_unit = TF_Get_LUN(factor_file)
        file_factor_unit = open(factor_file, 'rb')
        I2PY_SWAP_ENDIAN = Not_Same_Byte_Order(info.byte_order)
        if (factor_type == 3):    
            factor = zeros([ny, nx], dtype='Float32')
    
    #-----------------------
    #Open RTS_file to write
    #-----------------------
    Qn_SW_RTS_unit = TF_Get_LUN(Qn_SW_RTS_file)
    file_Qn_SW_RTS_unit = open(Qn_SW_RTS_file, 'wb')
    I2PY_SWAP_ENDIAN = Not_Same_Byte_Order(info.byte_order)
    
    #----------------------------------------------
    #Get start & stop times as decimal Julian days
    #----------------------------------------------
    start_Julian_day = Julian_Day(start_month, start_day)
    stop_Julian_day = Julian_Day(stop_month, stop_day)
    start_time = start_Julian_day + (start_hour / float64(24))
    stop_time = stop_Julian_day + (stop_hour / float64(24))
    
    #------------------------------------
    #Convert timestep from hours to days
    #------------------------------------
    timestep_JD = (timestep / float64(24))
    
    #------------------------------------
    #Create the RTS file, frame by frame
    #------------------------------------
    for Julian_day in arange(start_time, (stop_time)+(timestep_JD), timestep_JD):
        #----------------------------------------
        #Compute the offset from True Solar Noon
        #clock_hour is in 24-hour military time
        #but it can have a decimal part.
        #----------------------------------------
        clock_hour = (Julian_day - int16(Julian_day)) * float64(24)
        solar_noon = True_Solar_Noon(Julian_day, lon_deg, time_zone)
        t_offset = (clock_hour - solar_noon)    #[hours]
        
        #-------------------------
        #Write a progress message
        #-------------------------
        if ((MSG_BOX_ID not in [0,None])):    
            jstr = 'Day = ' + TF_String(int16(Julian_day))
            hstr = 'Hour = ' + TF_String(clock_hour, FORMAT='(F5.2)')
            mstr = jstr + ', ' + hstr
            widget_control(MSG_BOX_ID, set_value=mstr)
        
        #------------------------------
        #Read next values from files ?
        #------------------------------
        if (T_air_type == 1):    
            T_air = idl_func.readf(T_air_unit, T_air)
        if (T_air_type == 3):    
            T_air = fromfile(file_T_air_unit, count=size(T_air), dtype=str(T_air.dtype))
            if (I2PY_SWAP_ENDIAN):
                array(T_air, copy=0).byteswap(True)
        #------------------------------------------------------
        if (RH_type == 1):    
            RH = idl_func.readf(RH_unit, RH)
        if (RH_type == 3):    
            RH = fromfile(file_RH_unit, count=size(RH), dtype=str(RH.dtype))
            if (I2PY_SWAP_ENDIAN):
                array(RH, copy=0).byteswap(True)
        #------------------------------------------------------
        if (albedo_type == 1):    
            albedo = idl_func.readf(albedo_unit, albedo)
        if (albedo_type == 3):    
            albedo = fromfile(file_albedo_unit, count=size(albedo), dtype=str(albedo.dtype))
            if (I2PY_SWAP_ENDIAN):
                array(albedo, copy=0).byteswap(True)
        #-----------------------------------------------------------
        if (dust_att_type == 1):    
            dust_att = idl_func.readf(dust_att_unit, dust_att)
        if (dust_att_type == 3):    
            dust_att = fromfile(file_dust_att_unit, count=size(dust_att), dtype=str(dust_att.dtype))
            if (I2PY_SWAP_ENDIAN):
                array(dust_att, copy=0).byteswap(True)
        #-----------------------------------------------------------
        if (factor_type == 1):    
            factor = idl_func.readf(factor_unit, factor)
        if (factor_type == 3):    
            factor = fromfile(file_factor_unit, count=size(factor), dtype=str(factor.dtype))
            if (I2PY_SWAP_ENDIAN):
                array(factor, copy=0).byteswap(True)
        
        #-----------------------------------------
        #Compute min/max Sunrise and Sunset times
        #-----------------------------------------
        #sunrise_offset = Sunrise_Offset_Slope(lat_deg, Julian_day, alpha, beta)
        #sunrise_time   = solar_noon + sunrise_offset
        #---------------------------------------------
        #sunset_offset  = Sunset_Offset_Slope(lat_deg, Julian_day, alpha, beta)
        #sunset_time    = solar_noon + sunset_offset
        #---------------------------------------------
        #print,'solar_noon     = ', min(solar_noon),   max(solar_noon)
        #print,'sunrise_offset = ', min(sunrise_offset), max(sunrise_offset)
        #print,'sunset_offset  = ', min(sunset_offset), max(sunset_offset)
        #print,'sunrise_time   = ', min(sunrise_time), max(sunrise_time)
        #print,'sunset_time    = ', min(sunset_time),  max(sunset_time)
        #print,'-----------------------------------------------------'
        
        #--------------------------------
        # Compute Qnet_SW for this time
        #--------------------------------
        Qnet_SW = Clear_Sky_Radiation(lat_deg, Julian_day, T_air, RH,
                                      t_offset, alpha, beta, albedo, dust_att)
        
        #-------------------------------------------
        #Multiply by an optional factor that can be
        #used to account for cloud or canopy cover
        #-------------------------------------------
        Qnet_SW = (Qnet_SW * factor)
        
        #---------------------------
        #Make sure result is a grid
        #---------------------------
        if (idl_func.n_elements(Qnet_SW) == 1):    
            Qnet_SW = Qnet_SW + zeros([ny, nx], dtype='Float32')
        
        #------------------------
        #Write frame to RTS file
        #------------------------
        if (I2PY_SWAP_ENDIAN):
            array(float32(Qnet_SW), copy=0).byteswap(True)
        float32(Qnet_SW).tofile(file_Qn_SW_RTS_unit)
    
    
    #-------------------------
    #Write a finished message
    #-------------------------
    if ((MSG_BOX_ID not in [0,None])):    
        widget_control(MSG_BOX_ID, set_value='Finished.')
    
    #----------------
    #Close the files
    #----------------
    file_Qn_SW_RTS_unit.close()
    
    
    return _ret()
#  Make_RTS_File_for_Qnet_SW
#****************************************************************
def Get_Time_Zone_List(time_zones, time_zones2=None):

    n_params = 2 - [time_zones2].count(None)
    time_zones = array(['GMT – 12', 'GMT – 11', 'GMT – 10', 'GMT – 9', 'GMT – 8', 'GMT – 7', 'GMT – 6', 'GMT – 5', 'GMT – 4', 'GMT – 3', 'GMT – 2', 'GMT – 1', 'GMT    ', 'GMT + 1', 'GMT + 2', 'GMT + 3', 'GMT + 4', 'GMT + 5', 'GMT + 6', 'GMT + 7', 'GMT + 8', 'GMT + 9', 'GMT + 10', 'GMT + 11', 'GMT + 12'])
    time_zones = time_zones + ' hours'
    time_zones[11] = 'GMT - 1 hour'
    time_zones[12] = 'GMT    '
    time_zones[13] = 'GMT + 1 hour'
    
    
    time_zones2 = array(['GMT-12: International Data Line West', 'GMT-11: Midway Island, Samoa', 'GMT-10: Hawaii', 'GMT-9:  Alaska', 'GMT-8:  Pacific Time (US & Canada), Tijuana ', 'GMT-7:  Mountain Time (US & Canada)', 'GMT-6:  Central Time (US & Canada), Central America', 'GMT-5:  Eastern Time (US & Canada), Bogota, Lima', 'GMT-4:  Atlantic Time (Canada), Santiago', 'GMT-3:  Brasilia, Buenos Aires, Greenland ', 'GMT-2:  Mid-Atlantic ', 'GMT-1:  Azores, Cape Verde Island ', 'GMT:    Greenwich Mean Time, London, Casablanca ', 'GMT+1:  Amsterdam, Berlin, Madrid, Paris, Rome ', 'GMT+2:  Athens, Beirut, Cairo, Istanbul, Minsk ', 'GMT+3:  Baghdad, Kuwait, Moscow, Nairobi, Riyadh ', 'GMT+4:  Abu Dhabi, Baku, Muscat, Tbilisi, Yerevan ', 'GMT+5:  Ekaterinburg, Islamabad, Tashkent ', 'GMT+6:  Almaty, Astana, Dhaka, Novosibirsk ', 'GMT+7:  Bangkok, Hanoi, Jakarta, Krasnoyarsk ', 'GMT+8:  Beijing, Hong Kong, Taipei ', 'GMT+9:  Osaka, Tokyo, Seoul, Yakutsk ', 'GMT+10: Brisbane, Canberra, Guam, Vladivostok ', 'GMT+11: Magadan, New Caledonia, Solomon Is., ', 'GMT+12: Auckland, Fiji, Kamchatka, Wellington '])
    
#  Get_Time_Zone_List
#***************************************************************
def GUI_Make_Qnet_SW_File_event(event):

#-----------
#Error trap
#-----------
    n_params = 1
    def _ret():  return event
    
    # CATCH, status
    OK = Trace_Error(status, event)
    if logical_not(OK):    
        return _ret()
    
    Get_Event_Uvalue(event, uvalue, state)
    
    if uvalue == 'START_MONTH':    
        state.start_month = (event.index + 1)
        
        #*****************
    elif uvalue == 'STOP_MONTH':    
        state.stop_month = (event.index + 1)
        
        #****************
    elif uvalue == 'TIME_ZONE':    
        state.time_zone = (event.index - 12)
        
        #***********************************************************
    elif uvalue == 'T_AIR_TYPE':    
        state.T_air_type = event.index
    elif uvalue == 'RH_TYPE':    
        state.RH_type = event.index
    elif uvalue == 'ALBEDO_TYPE':    
        state.albedo_type = event.index
    elif uvalue == 'DUST_ATT_TYPE':    
        state.dust_att_type = event.index
    elif uvalue == 'FACTOR_TYPE':    
        state.factor_type = event.index
    elif uvalue == 'SLOPE_TYPE':    
        state.slope_type = event.index
    elif uvalue == 'ASPECT_TYPE':    
        state.aspect_type = event.index
    elif uvalue == 'LON_TYPE':    
        state.lon_type = (event.index * uint8(2))
    elif uvalue == 'LAT_TYPE':    
        state.lat_type = (event.index * uint8(2))
        #***********************************************************
        
        #************
    elif uvalue == 'START':    
        #---------------------------------
        #Get months selected via droplist
        #---------------------------------
        start_month = state.start_month
        stop_month = state.stop_month
        time_zone = state.time_zone
        #------------------------
        #Read start day and hour
        #------------------------
        Read_Text_Box(state.start_day_ID, start_day, OK, INTEGER=True)
        if logical_not(OK):    
            return _ret()
        Read_Text_Box(state.start_hour_ID, start_hour, OK, FLOAT=True)
        if logical_not(OK):    
            return _ret()
        #-----------------------
        #Read stop day and hour
        #-----------------------
        Read_Text_Box(state.stop_day_ID, stop_day, OK, INTEGER=True)
        if logical_not(OK):    
            return _ret()
        Read_Text_Box(state.stop_hour_ID, stop_hour, OK, FLOAT=True)
        if logical_not(OK):    
            return _ret()
        #---------------------------
        #Read the timestep in hours
        #---------------------------
        Read_Text_Box(state.timestep_ID, timestep, OK, FLOAT=True)
        if logical_not(OK):    
            return _ret()
        
        #--------------------------------------
        #Read name of new RTS file for Qnet_SW
        #--------------------------------------
        Read_Text_Box(state.Qn_SW_RTS_file_ID, Qn_SW_RTS_file, OK, TEXT=True)
        if logical_not(OK):    
            return _ret()
        Check_Overwrite(Qn_SW_RTS_file, OK)
        if logical_not(OK):    
            return _ret()
        #----------------------
        #Read name of RTI file
        #----------------------
        Read_Text_Box(state.RTI_file_ID, RTI_file, OK, TEXT=True)
        if logical_not(OK):    
            return _ret()
        else:    
            OK = File_Found(RTI_file)
        if logical_not(OK):    
            return _ret()
        
        #--------------------------
        #Read additional variables
        #--------------------------
        Read_Input_Type(state.T_air_type, state.T_air_ID, T_air, OK, filename=T_air_file)
        if logical_not(OK):    
            return _ret()
        #-----------------------------------------------------------------
        Read_Input_Type(state.RH_type, state.RH_ID, RH, OK, filename=RH_file)
        if logical_not(OK):    
            return _ret()
        #-----------------------------------------------------------------
        Read_Input_Type(state.albedo_type, state.albedo_ID, albedo, OK, filename=albedo_file)
        if logical_not(OK):    
            return _ret()
        #-----------------------------------------------------------------
        Read_Input_Type(state.dust_att_type, state.dust_att_ID, dust_att, OK, filename=dust_att_file)
        if logical_not(OK):    
            return _ret()
        #-----------------------------------------------------------------
        Read_Input_Type(state.factor_type, state.factor_ID, factor, OK, filename=factor_file)
        if logical_not(OK):    
            return _ret()
        #-----------------------------------------------------------------
        Read_Input_Type(state.slope_type, state.slope_ID, slope, OK, filename=slope_file)
        if logical_not(OK):    
            return _ret()
        #-----------------------------------------------------------------
        Read_Input_Type(state.aspect_type, state.aspect_ID, aspect, OK, filename=aspect_file)
        if logical_not(OK):    
            return _ret()
        
        #-------------------
        #Read lon and lat ?
        #-------------------
        if (state.lon_ID != int32(0)):    
            Read_Input_Type(state.lon_type, state.lon_ID, lon, OK, filename=lon_file)
            if logical_not(OK):    
                return _ret()
            #-------------------------------------------------------
            Read_Input_Type(state.lat_type, state.lat_ID, lat, OK, filename=lat_file)
            if logical_not(OK):    
                return _ret()
        
        #--------------------------------------------------------
        #Get number of frames in new RTS file from time settings
        #--------------------------------------------------------
        #Get start & stop times as decimal Julian days
        #----------------------------------------------
        Julian_start_day = Julian_Day(start_month, start_day, start_hour)
        Julian_stop_day = Julian_Day(stop_month, stop_day, stop_hour)
        timestep_JD = (timestep / float64(24))   #([hours] -> [days])
        nf_match = int32((Julian_stop_day - Julian_start_day) / timestep_JD)
        #-----------------------------------------------
        #Check if number of frames in RTS files matches
        #the number obtained from the time settings
        #-----------------------------------------------
        if (state.T_air_type == 3):    
            Check_Number_of_Frames(T_air_file, RTI_file, nf_match, OK)
        if logical_not(OK):    
            return _ret()
        #-----------------------------------------------------------------
        if (state.RH_type == 3):    
            Check_Number_of_Frames(RH_file, RTI_file, nf_match, OK)
        if logical_not(OK):    
            return _ret()
        #-----------------------------------------------------------------
        if (state.albedo_type == 3):    
            Check_Number_of_Frames(albedo_file, RTI_file, nf_match, OK)
        if logical_not(OK):    
            return _ret()
        #-----------------------------------------------------------------
        if (state.dust_att_type == 3):    
            Check_Number_of_Frames(dust_att_file, RTI_file, nf_match, OK)
        if logical_not(OK):    
            return _ret()
        #-----------------------------------------------------------------
        if (state.factor_type == 3):    
            Check_Number_of_Frames(factor_file, RTI_file, nf_match, OK)
        if logical_not(OK):    
            return _ret()
        
        #---------------------------------------
        #Call routines that create the RTS file
        #---------------------------------------
        widget_control(event.ID, sensitive=0)    #(disable button)
        Qn_SW_RTS_file, RTI_file, start_month, start_day, start_hour, stop_month, stop_day, stop_hour, timestep, time_zone, T_air, T_air_file, state.T_air_type, RH, RH_file, state.RH_type, albedo, albedo_file, state.albedo_type, dust_att, dust_att_file, state.dust_att_type, factor, factor_file, state.factor_type, slope, slope_file, state.slope_type, aspect, aspect_file, state.aspect_type, lon, lon_file, state.lon_type, lat, lat_file, state.lat_type, state.msg_box_ID = Make_RTS_File_for_Qnet_SW(Qn_SW_RTS_file, RTI_file, start_month, start_day, start_hour, stop_month, stop_day, stop_hour, timestep, time_zone, T_air, T_air_file, state.T_air_type, RH, RH_file, state.RH_type, albedo, albedo_file, state.albedo_type, dust_att, dust_att_file, state.dust_att_type, factor, factor_file, state.factor_type, slope, slope_file, state.slope_type, aspect, aspect_file, state.aspect_type, lon, lon_file, state.lon_type, lat, lat_file, state.lat_type, MSG_BOX_ID=state.msg_box_ID)
        widget_control(event.ID, sensitive=1)    #(enable button)
        
        #-------------------------
        #Show a "finished" dialog
        #-------------------------
        msg = array(['Finished creating Qnet-SW file.', ' '])
        result = GUI_Message(msg, INFO=True, TITLE="Finished")
        #*** Close_Dialog, event.top
        
        #***********
    elif uvalue == 'HELP':    
        Show_HTML_Help('shortwave_calc.htm')
        
        #************
    elif uvalue == 'CLOSE':    
        Close_Dialog(event.top)
        
    else:    
        dum = int16(0)
    
    
    if logical_and((uvalue != 'CLOSE'), (uvalue != 'START')):    
        widget_control(event.top, set_uvalue=state)
    
    
    return _ret()
#  GUI_Make_Qnet_SW_File_event
#****************************************************************
def GUI_Make_Qnet_SW_File(leader):

#-----------
#Error trap
#-----------
    n_params = 1
    def _ret():  return leader
    
    status = No_Catch()
    OK = Check_Error_Status(status)
    if logical_not(OK):    
        return _ret()
    
    if (idl_func.n_elements(leader) == 0):    
        leader = int32(0)
    
    #-----------------------------------
    #Get current values from main state
    #-----------------------------------
    Get_TLB_State(leader, mstate, ALIVE)
    if logical_not(ALIVE):    
        return _ret()
    #---------------------------------------
    prefix = mstate.run_vars.prefix
    Qn_SW_RTS_file = prefix + '_Qn-SW.rts'
    RTI_file = prefix + '.rti'
    info, OK = Read_RTI_File(RTI_file, True)    #(to get pixel_geom)
    if logical_not(OK):    
        return _ret()
    
    #--------------------------
    #Get default string values
    #--------------------------
    start_day_str = ' 1 '
    start_hour_str = ' 0.0 '
    stop_day_str = ' 2 '
    stop_hour_str = ' 0.0 '
    timestep_str = ' 1.0 '
    
    #------------------------------------
    #Structure to store selected options
    #------------------------------------
    state = idl_func.bunch(leader_ID=leader, msg_box_ID=int32(0), Qn_SW_RTS_file_ID=int32(0), RTI_file_ID=int32(0), start_month=1, start_day_ID=int32(0), start_hour_ID=int32(0), stop_month=1, stop_day_ID=int32(0), stop_hour_ID=int32(0), timestep_ID=int32(0), time_zone=0, T_air_ID=int32(0), T_air_type=uint8(3), RH_ID=int32(0), RH_type=uint8(0), albedo_ID=int32(0), albedo_type=uint8(0), dust_att_ID=int32(0), dust_att_type=uint8(0), factor_ID=int32(0), factor_type=uint8(0), slope_ID=int32(0), slope_type=uint8(2), aspect_ID=int32(0), aspect_type=uint8(2), lon_ID=int32(0), lon_type=uint8(0), lat_ID=int32(0), lat_type=uint8(0))
    
    ngap = int16(6)
    XS = int16(24)
    months = array([' January ', ' February ', ' March ', ' April ', ' May ', ' June ', ' July ', ' August ', ' September ', ' October ', ' November ', ' December '])
    Get_Time_Zone_List(time_zones)
    
    #-----------------
    #Main base widget
    #-----------------
    Create_TLB(MB, TITLE='Make RTS File for Qnet_SW', COLUMN=True, LEADER=leader)
    B1 = widget_base(MB, COLUMN=True, SPACE=1, FRAME=True)
    B2 = widget_base(MB, COLUMN=True, SPACE=1, FRAME=True)
    B3 = widget_base(MB, COLUMN=True, SPACE=1, FRAME=True)
    BOT = widget_base(MB, ROW=True, SPACE=1)
    
    #---------------------------
    #Get the start day and hour
    #---------------------------
    ST = widget_base(B1, ROW=True, SPACE=1)
    SM = widget_base(ST, ROW=True, SPACE=ngap)
    SM1 = widget_label(SM, VALUE='Start Month: ')
    SM2 = widget_droplist(SM, VALUE=months, UVALUE='START_MONTH')
    widget_control(SM2, set_droplist_select=0)
    #SM2 = widget_text(SM, VALUE=start_month_str, UVALUE='NONE', $
    #                  /EDITABLE, XSIZE=5)
    #--------------------------------------------------------------
    SD = widget_base(ST, ROW=True, SPACE=ngap)
    SD1 = widget_label(SD, VALUE='Day:')
    SD2 = widget_text(SD, VALUE=start_day_str, UVALUE='NONE', EDITABLE=True, XSIZE=5)
    state.start_day_ID = SD2
    #--------------------------------------------------------------
    SH = widget_base(ST, ROW=True, SPACE=ngap)
    SH1 = widget_label(SH, VALUE='Hour:')
    SH2 = widget_text(SH, VALUE=start_hour_str, UVALUE='NONE', EDITABLE=True, XSIZE=5)
    state.start_hour_ID = SH2
    
    #--------------------------
    #Get the stop day and hour
    #--------------------------
    ET = widget_base(B1, ROW=True, SPACE=1)
    EM = widget_base(ET, ROW=True, SPACE=ngap)
    EM1 = widget_label(EM, VALUE='Stop Month: ')
    EM2 = widget_droplist(EM, VALUE=months, UVALUE='STOP_MONTH')
    widget_control(EM2, set_droplist_select=0)
    #EM2 = widget_text(EM, VALUE=stop_month_str, UVALUE='NONE', $
    #                  /EDITABLE, XSIZE=5)
    #--------------------------------------------------------------
    ED = widget_base(ET, ROW=True, SPACE=ngap)
    ED1 = widget_label(ED, VALUE='Day:')
    ED2 = widget_text(ED, VALUE=stop_day_str, UVALUE='NONE', EDITABLE=True, XSIZE=5)
    state.stop_day_ID = ED2
    #--------------------------------------------------------------
    EH = widget_base(ET, ROW=True, SPACE=ngap)
    EH1 = widget_label(EH, VALUE='Hour:')
    EH2 = widget_text(EH, VALUE=stop_hour_str, UVALUE='NONE', EDITABLE=True, XSIZE=5)
    state.stop_hour_ID = EH2
    
    #-----------------
    #Get the timestep
    #-----------------
    TS0 = widget_base(B1, ROW=True, SPACE=1)
    TS = widget_base(TS0, ROW=True, SPACE=ngap)
    TS1 = widget_label(TS, VALUE='Timestep: ')
    TS2 = widget_text(TS, VALUE=timestep_str, UVALUE='NONE', EDITABLE=True, XSIZE=6)
    TS3 = widget_label(TS, VALUE='[hours] ')
    state.timestep_ID = TS2
    
    #------------------
    #Get the time zone
    #------------------
    TZ = widget_base(TS0, ROW=True, SPACE=ngap)
    TZ1 = widget_label(TZ, VALUE='  Time zone: ')
    TZ2 = widget_droplist(TZ, VALUE=time_zones, UVALUE='TIME_ZONE')
    widget_control(TZ2, set_droplist_select=12)
    
    #------------------
    #Get the time zone
    #------------------
    #TZ0 = widget_base(B1, /ROW, SPACE=1)
    #TZ = widget_base(TZ0, /ROW, SPACE=ngap)
    #  TZ1 = widget_label(TZ, VALUE='Time zone: ')
    #  TZ2 = widget_droplist(TZ, VALUE=time_zones2, UVALUE='TIME_ZONE')
    #  widget_control, TZ2, set_droplist_select=12
    
    #------------------
    #Align the widgets
    #------------------
    Align_Text_Boxes(array([SM1, EM1, TS1]))
    Align_Text_Boxes(array([SM2, EM2]))
    
    #---------------------
    #Add some blank space
    #----------------------
    #** P0 = widget_label(B1, VALUE=' ')
    
    
    #------------------------
    #Get the input variables
    #------------------------
    ngap = int16(6)
    XS = int16(22)
    types = Model_Input_Types()
    gtype = array([' Grid '])
    sg_type = array([' Scalar ', ' Grid '])
    #------------------------------------
    T_air_str = prefix + '_Tair.rts'
    RH_str = '0.30'
    albedo_str = '0.8'
    dust_att_str = '0.08'
    factor_str = '1.0'
    slope_str = prefix + '_slope.rtg'
    aspect_str = prefix + '_aspect.rtg'
    lon_str = '0.0'
    lat_str = '0.0'
    
    #-------------------
    #Get the parameters
    #-------------------
    A1 = widget_base(B2, ROW=True, SPACE=ngap)
    A11 = widget_label(A1, VALUE='Variable: ', UVALUE='NONE')
    A12 = widget_label(A1, VALUE='Type: ', UVALUE='NONE')
    A13 = widget_label(A1, VALUE='Scalar or Grid Filename: ', UVALUE='NONE')
    A14 = widget_label(A1, VALUE='Units: ', UVALUE='NONE')
    #--------------------------------------------------------------------------
    TA = widget_base(B2, ROW=True, SPACE=ngap)
    TA1 = widget_label(TA, VALUE='T_air: ', UVALUE='NONE')
    TA2 = widget_droplist(TA, VALUE=types, UVALUE='T_AIR_TYPE')
    TA3 = widget_text(TA, VALUE=T_air_str, UVALUE='NONE', EDITABLE=True, XSIZE=XS)
    TA4 = widget_label(TA, VALUE='[deg C]', UVALUE='NONE')
    state.T_air_ID = TA3
    widget_control(TA2, set_droplist_select=3)
    #--------------------------------------------------------------------------
    RH = widget_base(B2, ROW=True, SPACE=ngap)
    RH1 = widget_label(RH, VALUE='RH: ', UVALUE='NONE')
    RH2 = widget_droplist(RH, VALUE=types, UVALUE='RH_TYPE')
    RH3 = widget_text(RH, VALUE=RH_str, UVALUE='NONE', EDITABLE=True, XSIZE=XS)
    RH4 = widget_label(RH, VALUE='[none]  in [0,1]', UVALUE='NONE')
    state.RH_ID = RH3
    #--------------------------------------------------------------------------
    #CC = widget_base(B2, /ROW, SPACE=ngap)
    #  CC1 = widget_label(CC, VALUE='cloud frac.: ', UVALUE='NONE')
    #  CC2 = widget_droplist(CC, VALUE=types, UVALUE='C_TYPE')
    #  CC3 = widget_text(CC, VALUE=C_str, UVALUE='NONE', /EDITABLE, XSIZE=XS)
    #  CC4 = widget_label(CC, VALUE='[none]  in [0,1]', UVALUE='NONE')
    #  state.C_ID = CC3
    #--------------------------------------------------------------------------
    #FF = widget_base(B2, /ROW, SPACE=ngap)
    #  FF1 = widget_label(FF, VALUE='canopy frac.: ', UVALUE='NONE')
    #  FF2 = widget_droplist(FF, VALUE=types, UVALUE='F_TYPE')
    #  FF3 = widget_text(FF, VALUE=F_str, UVALUE='NONE', /EDITABLE, XSIZE=XS)
    #  FF4 = widget_label(FF, VALUE='[none]  in [0,1]', UVALUE='NONE')
    #  state.F_ID = FF3
    #--------------------------------------------------------------------------
    AL = widget_base(B2, ROW=True, SPACE=ngap)
    AL1 = widget_label(AL, VALUE='albedo: ', UVALUE='NONE')
    AL2 = widget_droplist(AL, VALUE=types, UVALUE='ALBEDO_TYPE')
    AL3 = widget_text(AL, VALUE=albedo_str, UVALUE='NONE', EDITABLE=True, XSIZE=XS)
    AL4 = widget_label(AL, VALUE='[none]  in [0,1]', UVALUE='NONE')
    state.albedo_ID = AL3
    #*** widget_control, AL2, set_droplist_select=3
    #--------------------------------------------------------------------------
    DA = widget_base(B2, ROW=True, SPACE=ngap)
    DA1 = widget_label(DA, VALUE='dust att.: ', UVALUE='NONE')
    DA2 = widget_droplist(DA, VALUE=types, UVALUE='DUST_ATT_TYPE')
    DA3 = widget_text(DA, VALUE=dust_att_str, UVALUE='NONE', EDITABLE=True, XSIZE=XS)
    DA4 = widget_label(DA, VALUE='[none]  in [0,1]', UVALUE='NONE')
    state.dust_att_ID = DA3
    #--------------------------------------------------------------------------
    FA = widget_base(B2, ROW=True, SPACE=ngap)
    FA1 = widget_label(FA, VALUE='factor: ', UVALUE='NONE')
    FA2 = widget_droplist(FA, VALUE=types, UVALUE='FACTOR_TYPE')
    FA3 = widget_text(FA, VALUE=factor_str, UVALUE='NONE', EDITABLE=True, XSIZE=XS)
    FA4 = widget_label(FA, VALUE='[none]  in [0,1]', UVALUE='NONE')
    state.factor_ID = FA3
    #--------------------------------------------------------------------------
    SL = widget_base(B2, ROW=True, SPACE=ngap)
    SL1 = widget_label(SL, VALUE='slope: ', UVALUE='NONE')
    SL2 = widget_droplist(SL, VALUE=gtype, UVALUE='SLOPE_TYPE')
    SL3 = widget_text(SL, VALUE=slope_str, UVALUE='NONE', EDITABLE=True, XSIZE=XS)
    SL4 = widget_label(SL, VALUE='[m/m]', UVALUE='NONE')
    state.slope_ID = SL3
    #--------------------------------------------------------------------------
    _as = widget_base(B2, ROW=True, SPACE=ngap)
    AS1 = widget_label(_as, VALUE='aspect: ', UVALUE='NONE')
    AS2 = widget_droplist(_as, VALUE=gtype, UVALUE='ASPECT_TYPE')
    AS3 = widget_text(_as, VALUE=aspect_str, UVALUE='NONE', EDITABLE=True, XSIZE=XS)
    AS4 = widget_label(_as, VALUE='[radians]', UVALUE='NONE')
    state.aspect_ID = AS3
    
    #-------------------------------------------------
    #Read longitude and latitude (if not in RTI file)
    #-------------------------------------------------
    if (info.pixel_geom != uint8(0)):    
        LO = widget_base(B2, ROW=True, SPACE=ngap)
        LO1 = widget_label(LO, VALUE='longitude: ', UVALUE='NONE')
        LO2 = widget_droplist(LO, VALUE=sg_type, UVALUE='LON_TYPE')
        LO3 = widget_text(LO, VALUE=lon_str, EDITABLE=True, XSIZE=XS)
        LO4 = widget_label(LO, VALUE='[deg]', UVALUE='NONE')
        state.lon_ID = LO3
        #--------------------------------------------------------------------------
        LA = widget_base(B2, ROW=True, SPACE=ngap)
        LA1 = widget_label(LA, VALUE='latitude: ', UVALUE='NONE')
        LA2 = widget_droplist(LA, VALUE=sg_type, UVALUE='LAT_TYPE')
        LA3 = widget_text(LA, VALUE=lat_str, EDITABLE=True, XSIZE=XS)
        LA4 = widget_label(LA, VALUE='[deg]', UVALUE='NONE')
        state.lat_ID = LA3
    
    #---------------------
    #Add some blank space
    #---------------------
    P1 = widget_label(B2, VALUE=' ')
    
    #------------------
    #Align the widgets
    #------------------
    if (info.pixel_geom != uint8(0)):    
        Align_Text_Boxes(array([A11, TA1, RH1, AL1, DA1, FA1, SL1, AS1, LO1, LA1]))
        Align_Text_Boxes(array([A12, TA2, RH2, AL2, DA2, FA2, SL2, AS2, LO2, LA2]))
        Align_Text_Boxes(array([A13, TA3, RH3, AL3, DA3, FA3, SL3, AS3, LO3, LA3]))
    else:    
        Align_Text_Boxes(array([A11, TA1, RH1, AL1, DA1, FA1, SL1, AS1]))
        Align_Text_Boxes(array([A12, TA2, RH2, AL2, DA2, FA2, SL2, AS2]))
        Align_Text_Boxes(array([A13, TA3, RH3, AL3, DA3, FA3, SL3, AS3]))
    
    #-------------------------------------
    #Get name of new RTS file for Qnet_SW
    #-------------------------------------
    NF = widget_base(B3, ROW=True, SPACE=ngap)
    NF1 = widget_label(NF, VALUE='New RTS file for Qnet_SW: ')
    NF2 = widget_text(NF, VALUE=Qn_SW_RTS_file, UVALUE='NONE', EDITABLE=True, XSIZE=XS)
    NF3 = widget_label(NF, VALUE=' [W / m^2] ')
    state.Qn_SW_RTS_file_ID = NF2
    
    #-----------------------------------
    #Get name of new RTS file for e_air
    #-----------------------------------
    #EA = widget_base(B3, /ROW, SPACE=ngap)
    #  EA1 = widget_label(EA, VALUE='New RTS file for e_air: ')
    #  EA2 = widget_text(EA, VALUE=e_air_RTS_file, UVALUE='NONE', $
    #                    /EDITABLE, XSIZE=XS)
    #  EA3 = widget_label(EA, VALUE=' [mbar] ')
    #  state.e_air_RTS_file_ID = EA2
    
    #-------------------------------------
    #Get name of new RTS file for em_air
    #-------------------------------------
    #EM = widget_base(B3, /ROW, SPACE=ngap)
    #  EM1 = widget_label(EM, VALUE='New RTS file for em_air: ')
    #  EM2 = widget_text(EM, VALUE=em_air_RTS_file, UVALUE='NONE', $
    #                    /EDITABLE, XSIZE=XS)
    #  EM3 = widget_label(EM, VALUE=' [none] ')
    #  state.em_air_RTS_file_ID = EM2
    
    #-----------------------
    #Get name of RTI file
    #-----------------------
    RF = widget_base(B3, ROW=True, SPACE=ngap)
    RF1 = widget_label(RF, VALUE='Existing RTI file: ')
    RF2 = widget_text(RF, VALUE=RTI_file, UVALUE='NONE', EDITABLE=True, XSIZE=XS)
    state.RTI_file_ID = RF2
    
    #---------------------
    #Add some blank space
    #---------------------
    P1 = widget_label(B3, VALUE=' ')
    
    #------------------
    #Align the widgets
    #------------------
    Align_Text_Boxes(array([NF1, RF1]))
    Align_Text_Boxes(array([NF2, RF2]))
    
    #------------------
    #Bottom button bar
    #------------------
    CW_Button_Bar(BOT, START=True, HELP=True, CLOSE=True)
    
    #---------------------
    #A status message box
    #---------------------
    MS = widget_base(BOT, ROW=True, SPACE=ngap)
    MS1 = widget_label(MS, VALUE='Status: ')
    MS2 = widget_text(MS, VALUE='Ready.', XSIZE=20)
    state.msg_box_ID = MS2
    
    #------------------------------------
    #Realize widgets and wait for events
    #------------------------------------
    XOFF = int16(480)
    Realize_TLB(MB, state, 'GUI_Make_Qnet_SW_File', XOFF=XOFF)
    
    
    return _ret()
#  GUI_Make_Qnet_SW_File
#***************************************************************
def Compare_em_air_Methods():

#------------------------------------------------------------
#Notes:  There seems to be two different methods that are
#        commonly used to compute the vapor pressure of air,
#        e_air, and then the emissivity of air, em_air, for
#        use in longwave radiation calculations.
#        This routine compares them graphically.
#------------------------------------------------------------
    n_params = 0
    T_air = arange(80, dtype='Float32') - float32(40.0)   #[Celsius]  (-40 to 40)
    RH = float32(1.0)
    C2K = float32(273.15)
    
    #------------------------
    #Brutsaert (1975) method
    #------------------------
    e_air1 = RH * float32(0.611) * exp((float32(17.3) * T_air) / (T_air + float32(237.3)))  #[kPa]
    em_air1 = float32(1.72) * (e_air1 / (T_air + C2K)) ** (float32(1.0) / 7)
    
    #--------------------------
    #Satterlund (1979) method
    #--------------------------
    #NB! e_air has units of Pa
    #--------------------------
    e_air2 = RH * float64(10) ** (float32(11.40) - (float64(2353) / (T_air + C2K)))   #[Pa]
    eterm = exp(-float64(1) * (e_air2 / float64(100)) ** ((T_air + C2K) / float64(2016)))
    em_air2 = float64(1.08) * (float64(1) - eterm)
    
    #--------------------------
    #Plot the two e_air curves
    #------------------------------
    #These two agree quite closely
    #------------------------------
    matplotlib.pyplot.figure(figsize=(8, 6), dpi=80)
    matplotlib.pyplot.show()
    matplotlib.pyplot.plot(T_air, e_air1)
    matplotlib.pyplot.show()
    oplot(T_air, (e_air2 / float32(1000.0)), psym=-3)   #[Pa -> kPa]
    
    #---------------------------
    #Plot the two em_air curves
    #------------------------------------------------
    #These two don't agree very well for some reason
    #------------------------------------------------
    matplotlib.pyplot.figure(figsize=(8, 6), dpi=80)
    matplotlib.pyplot.show()
    matplotlib.pyplot.plot(T_air, em_air1)
    matplotlib.pyplot.show()
    oplot(T_air, em_air2, psym=-3)
    
#  Compare_em_air_Methods
#***************************************************************
def Make_RTS_File_for_Qnet_LW(Qn_LW_RTS_file, e_air_RTS_file,
                              em_air_RTS_file, RTI_file,
                              start_month, start_day, start_hour,
                              stop_month, stop_day, stop_hour,
                              timestep,
                              T_air, T_air_file, T_air_type,
                              RH, RH_file, RH_type,
                              C, C_file, C_type,
                              F, F_file, F_type,
                              T_surf, T_surf_file, T_surf_type,
                              em_surf, em_surf_file, em_surf_type,
                              MSG_BOX_ID=None):

#----------------------------------------------------------------
#Notes:  Net longwave radiation is computed using the
#        Stefan-Boltzman law.  All four data types
#        should be allowed (scalar, time series, grid or
#        grid stack).

#        Qnet_LW = (LW_in - LW_out)
#        LW_in   = em_air  * sigma * (T_air  + 273.15)^4
#        LW_out  = em_surf * sigma * (T_surf + 273.15)^4
#
#        Temperatures in [deg_C] must be converted to
#        [deg_K].  Recall that absolute zero occurs at
#        0 [deg_K] or -273.15 [deg_C].

#----------------------------------------------------------------
# First, e_air is computed as:
#   e_air = RH * 0.611 * exp[(17.3 * T_air) / (T_air + 237.3)]
# Then, em_air is computed as:
#   em_air = (1 - F) * 1.72 * [e_air / (T_air + 273.15)]^(1/7) *
#             (1 + 0.22 * C^2) + F
#----------------------------------------------------------------
    n_params = 29
    _opt = (MSG_BOX_ID,)
    def _ret():
        _optrv = zip(_opt, [MSG_BOX_ID])
        _rv = [Qn_LW_RTS_file, e_air_RTS_file, em_air_RTS_file, RTI_file, start_month, start_day, start_hour, stop_month, stop_day, stop_hour, timestep, T_air, T_air_file, T_air_type, RH, RH_file, RH_type, C, C_file, C_type, F, F_file, F_type, T_surf, T_surf_file, T_surf_type, em_surf, em_surf_file, em_surf_type]
        _rv += [_o[1] for _o in _optrv if _o[0] is not None]
        return tuple(_rv)
    
    # FORWARD_FUNCTION Vapor_Pressure
    
    #** sigma = 5.6697e-8   ;[W/(m^2 K^4)]  (Stefan-Boltzman constant)
    sigma = float32(5.67E-8)   #[W/(m^2 K^4)]  (Stefan-Boltzman constant)
    C2K = float32(273.15)    #(add to convert [deg_C] to [deg_K])
    
    #----------------------------
    #Read info from the RTI file
    #----------------------------
    info = Read_RTI_File(RTI_file)
    nx = info.ncols
    ny = info.nrows
    
    #------------------------------------
    #If data type is grid, read the grid
    #------------------------------------
    if (T_air_type == 2):    
        T_air = Read_Grid(T_air_file, _TYPE='FLOAT', SILENT=True)
    if (RH_type == 2):    
        RH = Read_Grid(RH_file, _TYPE='FLOAT', SILENT=True)
    if (C_type == 2):    
        C = Read_Grid(C_file, _TYPE='FLOAT', SILENT=True)
    if (F_type == 2):    
        F = Read_Grid(F_file, _TYPE='FLOAT', SILENT=True)
    #-----------------------------------------------------------
    if (T_surf_type == 2):    
        T_surf = Read_Grid(T_surf_file, _TYPE='FLOAT', SILENT=True)
    if (em_surf_type == 2):    
        em_surf = Read_Grid(em_surf_file, _TYPE='FLOAT', SILENT=True)
    
    #---------------------------------
    #Open files for input variables ?
    #---------------------------------
    if logical_or((T_air_type == 1), (T_air_type == 3)):    
        T_air_unit = TF_Get_LUN(T_air_file)
        file_T_air_unit = open(T_air_file, 'rb')
        I2PY_SWAP_ENDIAN = Not_Same_Byte_Order(info.byte_order)
        if (T_air_type == 3):    
            T_air = zeros([ny, nx], dtype='Float32')
    #-----------------------------------------------------------
    if logical_or((RH_type == 1), (RH_type == 3)):    
        RH_unit = TF_Get_LUN(RH_file)
        file_RH_unit = open(RH_file, 'rb')
        I2PY_SWAP_ENDIAN = Not_Same_Byte_Order(info.byte_order)
        if (RH_type == 3):    
            RH = zeros([ny, nx], dtype='Float32')
    #-----------------------------------------------------------
    if logical_or((C_type == 1), (C_type == 3)):    
        C_unit = TF_Get_LUN(C_file)
        file_C_unit = open(C_file, 'rb')
        I2PY_SWAP_ENDIAN = Not_Same_Byte_Order(info.byte_order)
        if (C_type == 3):    
            C = zeros([ny, nx], dtype='Float32')
    #-----------------------------------------------------------
    if logical_or((F_type == 1), (F_type == 3)):    
        F_unit = TF_Get_LUN(F_file)
        file_F_unit = open(F_file, 'rb')
        I2PY_SWAP_ENDIAN = Not_Same_Byte_Order(info.byte_order)
        if (F_type == 3):    
            F = zeros([ny, nx], dtype='Float32')
    #-----------------------------------------------------------
    if logical_or((T_surf_type == 1), (T_surf_type == 3)):    
        T_surf_unit = TF_Get_LUN(T_surf_file)
        file_T_surf_unit = open(T_surf_file, 'rb')
        I2PY_SWAP_ENDIAN = Not_Same_Byte_Order(info.byte_order)
        if (T_surf_type == 3):    
            T_surf = zeros([ny, nx], dtype='Float32')
    #-----------------------------------------------------------
    if logical_or((em_surf_type == 1), (em_surf_type == 3)):    
        em_surf_unit = TF_Get_LUN(em_surf_file)
        file_em_surf_unit = open(em_surf_file, 'rb')
        I2PY_SWAP_ENDIAN = Not_Same_Byte_Order(info.byte_order)
        if (em_surf_type == 3):    
            em_surf = zeros([ny, nx], dtype='Float32')
    
    #------------------------
    #Open RTS_files to write
    #------------------------
    Qn_LW_RTS_unit = TF_Get_LUN(Qn_LW_RTS_file)
    file_Qn_LW_RTS_unit = open(Qn_LW_RTS_file, 'wb')
    I2PY_SWAP_ENDIAN = Not_Same_Byte_Order(info.byte_order)
    #--------------------------------------------------------
    e_air_RTS_unit = TF_Get_LUN(e_air_RTS_file)
    file_e_air_RTS_unit = open(e_air_RTS_file, 'wb')
    I2PY_SWAP_ENDIAN = Not_Same_Byte_Order(info.byte_order)
    #--------------------------------------------------------
    em_air_RTS_unit = TF_Get_LUN(em_air_RTS_file)
    file_em_air_RTS_unit = open(em_air_RTS_file, 'wb')
    I2PY_SWAP_ENDIAN = Not_Same_Byte_Order(info.byte_order)
    
    #----------------------------------------------
    #Get start & stop times as decimal Julian days
    #----------------------------------------------
    start_Julian_day = Julian_Day(start_month, start_day)
    stop_Julian_day = Julian_Day(stop_month, stop_day)
    start_time = start_Julian_day + (start_hour / float64(24))
    stop_time = stop_Julian_day + (stop_hour / float64(24))
    
    #------------------------------------
    #Convert timestep from hours to days
    #------------------------------------
    timestep_JD = (timestep / float64(24))
    
    #------------------------------------
    #Create the RTS file, frame by frame
    #------------------------------------
    for Julian_day in arange(start_time, (stop_time)+(timestep_JD), timestep_JD):
    #----------------------------------------
    #Compute the offset from True Solar Noon
    #clock_hour is in 24-hour military time
    #but it can have a decimal part.
    #----------------------------------------
        clock_hour = (Julian_day - int16(Julian_day)) * float64(24)
        
        #-------------------------
        #Write a progress message
        #-------------------------
        if ((MSG_BOX_ID not in [0,None])):    
            jstr = 'Day = ' + TF_String(int16(Julian_day))
            hstr = 'Hour = ' + TF_String(clock_hour, FORMAT='(F5.2)')
            mstr = jstr + ', ' + hstr
            widget_control(MSG_BOX_ID, set_value=mstr)
        
        #------------------------------
        #Read next values from files ?
        #------------------------------
        if (T_air_type == 1):    
            T_air = idl_func.readf(T_air_unit, T_air)
        if (T_air_type == 3):    
            T_air = fromfile(file_T_air_unit, count=size(T_air), dtype=str(T_air.dtype))
            if (I2PY_SWAP_ENDIAN):
                array(T_air, copy=0).byteswap(True)
        #--------------------------------------------------------
        if (RH_type == 1):    
            RH = idl_func.readf(RH_unit, RH)
        if (RH_type == 3):    
            RH = fromfile(file_RH_unit, count=size(RH), dtype=str(RH.dtype))
            if (I2PY_SWAP_ENDIAN):
                array(RH, copy=0).byteswap(True)
        #--------------------------------------------------------
        if (C_type == 1):    
            C = idl_func.readf(C_unit, C)
        if (C_type == 3):    
            C = fromfile(file_C_unit, count=size(C), dtype=str(C.dtype))
            if (I2PY_SWAP_ENDIAN):
                array(C, copy=0).byteswap(True)
        #--------------------------------------------------------
        if (F_type == 1):    
            F = idl_func.readf(F_unit, F)
        if (F_type == 3):    
            F = fromfile(file_F_unit, count=size(F), dtype=str(F.dtype))
            if (I2PY_SWAP_ENDIAN):
                array(F, copy=0).byteswap(True)
        #--------------------------------------------------------
        if (T_surf_type == 1):    
            T_surf = idl_func.readf(T_surf_unit, T_surf)
        if (T_surf_type == 3):    
            T_surf = fromfile(file_T_surf_unit, count=size(T_surf), dtype=str(T_surf.dtype))
            if (I2PY_SWAP_ENDIAN):
                array(T_surf, copy=0).byteswap(True)
        #--------------------------------------------------------
        if (em_surf_type == 1):    
            em_surf = idl_func.readf(em_surf_unit, em_surf)
        if (em_surf_type == 3):    
            em_surf = fromfile(file_em_surf_unit, count=size(em_surf), dtype=str(em_surf.dtype))
            if (I2PY_SWAP_ENDIAN):
                array(em_surf, copy=0).byteswap(True)
        
        #-------------------------------------------------------
        #Brutsaert (1975) method for computing emissivity
        #of the air, em-air.  (From Dingman (2002, p. 196))
        #See notes for Vapor_Pressure function in formulas.pro.
        #RH = relative humidity [unitless]
        #-------------------------------------------------------
        #NB!  Temperatures are assumed to be given with units
        #     of degrees Celsius and are converted to Kelvin
        #     wherever necessary by adding C2K = 273.15.
        #-------------------------------------------------------
        #NB!  I'm not sure about how F is added at end because
        #     of how the equation is printed in Dingman (2002)
        #-------------------------------------------------------
        e_air = Vapor_Pressure(T_air, RH)   #[kPa]
        term1 = (float32(1.0) - F) * float32(1.72) * (e_air / (T_air + C2K)) ** (float32(1.0) / 7)
        term2 = (float32(1.0) + (float32(0.22) * C ** float32(2.0)))
        em_air = (term1 * term2) + F   #***  DOUBLE CHECK  ***
        
        #-------------------------------------------------------
        #Convert e_air from kPa to mbar before saving, since
        #those are the units used in the snowmelt energy balance
        #equations obtained from Zhang et al. (2000).
        #-------------------------------------------------------
        #NB!  100 kPa = 1 bar = 1000 mbars
        #      => 1 kPa = 10 mbars
        #----------------------------------
        e_air = (e_air * float32(10.0))   #[mbars]
        
        #------------------------------------------------------
        #NB!  This formula for vapor pressure as a function of
        #     temperature compares well with the Brutsaert
        #     formula above, see Compare_em_air_Methods.
        #     However, must pay close attention to whether
        #     equations require units of kPa, Pa, or mbar.
        #     Note also that we must add C2K to T_air below.
        #------------------------------------------------------
        #Satterlund (1979) method for computing the emissivity
        #of the air, em_air, that is intended to "correct
        #apparent deficiencies in this formulation at air
        #temperatures below 0 degrees C" (see G. Liston)
        #Liston cites Aase and Idso(1978), Satterlund (1979)
        #------------------------------------------------------
        #e_air  = Vapor_Pressure(T_air, RH, /SATT, /MBAR)   ;[mbar]
        #eterm  = exp(-1d * (e_air_mb)^((T_air + C2K) / 2016d))
        #em_air = 1.08d * (1d - eterm)
        
        #-----------------------------------
        #Compute Qnet_LW grid for this time
        #-----------------------------------
        LW_in = em_air * sigma * (T_air + C2K) ** float64(4)
        LW_out = em_surf * sigma * (T_surf + C2K) ** float64(4)
        LW_out = LW_out + ((float32(1.0) - em_surf) * LW_in)
        Qnet_LW = (LW_in - LW_out)
        
        #---------------------------
        #Make sure result is a grid
        #---------------------------
        if (idl_func.n_elements(Qnet_LW) == 1):    
            Qnet_LW = Qnet_LW + zeros([ny, nx], dtype='Float32')
        if (idl_func.n_elements(e_air) == 1):    
            e_air = e_air + zeros([ny, nx], dtype='Float32')
        if (idl_func.n_elements(em_air) == 1):    
            em_air = em_air + zeros([ny, nx], dtype='Float32')
        
        #--------------------------
        #Write frames to RTS files
        #--------------------------
        if (I2PY_SWAP_ENDIAN):
            array(float32(Qnet_LW), copy=0).byteswap(True)
        float32(Qnet_LW).tofile(file_Qn_LW_RTS_unit)
        if (I2PY_SWAP_ENDIAN):
            array(float32(e_air), copy=0).byteswap(True)
        float32(e_air).tofile(file_e_air_RTS_unit)
        if (I2PY_SWAP_ENDIAN):
            array(float32(em_air), copy=0).byteswap(True)
        float32(em_air).tofile(file_em_air_RTS_unit)
    
    #-------------------------
    #Write a finished message
    #-------------------------
    if ((MSG_BOX_ID not in [0,None])):    
        widget_control(MSG_BOX_ID, set_value='Finished.')
    
    #--------------------
    #Close input files ?
    #--------------------
    if logical_or((T_air_type == 1), (T_air_type == 3)):    
        file_T_air_unit.close()
    if logical_or((RH_type == 1), (RH_type == 3)):    
        file_RH_unit.close()
    if logical_or((C_type == 1), (C_type == 3)):    
        file_C_unit.close()
    if logical_or((F_type == 1), (F_type == 3)):    
        file_F_unit.close()
    if logical_or((T_surf_type == 1), (T_surf_type == 3)):    
        file_T_surf_unit.close()
    if logical_or((em_surf_type == 1), (em_surf_type == 3)):    
        file_em_surf_unit.close()
    
    #-----------------------
    #Close the output files
    #-----------------------
    file_Qn_LW_RTS_unit.close()
    file_e_air_RTS_unit.close()
    file_em_air_RTS_unit.close()
    
    
    return _ret()
#  Make_RTS_File_for_Qnet_LW
#***************************************************************
def Check_Number_of_Frames(RTS_file, RTI_file, nf_match, OK):

#-----------
#Error trap
#-----------
    n_params = 4
    def _ret():  return (RTS_file, RTI_file, nf_match, OK)
    
    # CATCH, status
    OK = Trace_Error(status, event)
    if logical_not(OK):    
        return _ret()
    
    # FORWARD_FUNCTION Number_of_Frames
    OK = uint8(1)
    
    #-----------------------------------
    #Count number of frames in RTS file
    #-----------------------------------
    nf = Number_Of_Frames(RTS_file, RTI_file)
    
    #--------------------------
    #Issue a warning message ?
    #--------------------------
    if (nf != nf_match):    
        msg = array(['ERROR: ', ' ', 'Number of frames in the RTS file: ', '   ' + RTS_file, ' ', 'does not match the number of frames in the RTS', 'file to be created, based on the time settings. ', ' ', 'Number of frames in RTS file = ' + TF_String(nf), 'Number of frames in new file = ' + TF_String(nf_match), ' '])
        GUI_Error_Message(msg)
        OK = uint8(0)
    
    
    return _ret()
#  Check_Number_of_Frames
#***************************************************************
def GUI_Make_Qnet_LW_File_event(event):

#-----------
#Error trap
#-----------
    n_params = 1
    def _ret():  return event
    
    # CATCH, status
    OK = Trace_Error(status, event)
    if logical_not(OK):    
        return _ret()
    
    Get_Event_Uvalue(event, uvalue, state)
    
    if uvalue == 'START_MONTH':    
        state.start_month = (event.index + 1)
        
        #*****************
    elif uvalue == 'STOP_MONTH':    
        state.stop_month = (event.index + 1)
        
        #**************************************************
    elif uvalue == 'T_AIR_TYPE':    
        state.T_air_type = event.index
    elif uvalue == 'RH_TYPE':    
        state.RH_type = event.index
    elif uvalue == 'C_TYPE':    
        state.C_type = event.index
    elif uvalue == 'F_TYPE':    
        state.F_type = event.index
    elif uvalue == 'T_SURF_TYPE':    
        state.T_surf_type = event.index
    elif uvalue == 'EM_SURF_TYPE':    
        state.em_surf_type = event.index
        #**************************************************
        
        #************
    elif uvalue == 'START':    
        #---------------------------------
        #Get months selected via droplist
        #---------------------------------
        start_month = state.start_month
        stop_month = state.stop_month
        #------------------------
        #Read start day and hour
        #------------------------
        Read_Text_Box(state.start_day_ID, start_day, OK, INTEGER=True)
        if logical_not(OK):    
            return _ret()
        Read_Text_Box(state.start_hour_ID, start_hour, OK, FLOAT=True)
        if logical_not(OK):    
            return _ret()
        #-----------------------
        #Read stop day and hour
        #-----------------------
        Read_Text_Box(state.stop_day_ID, stop_day, OK, INTEGER=True)
        if logical_not(OK):    
            return _ret()
        Read_Text_Box(state.stop_hour_ID, stop_hour, OK, FLOAT=True)
        if logical_not(OK):    
            return _ret()
        #---------------------------
        #Read the timestep in hours
        #---------------------------
        Read_Text_Box(state.timestep_ID, timestep, OK, FLOAT=True)
        if logical_not(OK):    
            return _ret()
        
        #--------------------------------------
        #Read name of new RTS file for Qnet_LW
        #--------------------------------------
        Read_Text_Box(state.Qn_LW_RTS_file_ID, Qn_LW_RTS_file, OK, TEXT=True)
        if logical_not(OK):    
            return _ret()
        Check_Overwrite(Qn_LW_RTS_file, OK)
        if logical_not(OK):    
            return _ret()
        #------------------------------------
        #Read name of new RTS file for e_air
        #------------------------------------
        Read_Text_Box(state.e_air_RTS_file_ID, e_air_RTS_file, OK, TEXT=True)
        if logical_not(OK):    
            return _ret()
        Check_Overwrite(e_air_RTS_file, OK)
        if logical_not(OK):    
            return _ret()
        #-------------------------------------
        #Read name of new RTS file for em_air
        #-------------------------------------
        Read_Text_Box(state.em_air_RTS_file_ID, em_air_RTS_file, OK, TEXT=True)
        if logical_not(OK):    
            return _ret()
        Check_Overwrite(em_air_RTS_file, OK)
        if logical_not(OK):    
            return _ret()
        #----------------------
        #Read name of RTI file
        #----------------------
        Read_Text_Box(state.RTI_file_ID, RTI_file, OK, TEXT=True)
        if logical_not(OK):    
            return _ret()
        else:    
            OK = File_Found(RTI_file)
        if logical_not(OK):    
            return _ret()
        
        #--------------------------
        #Read additional variables
        #--------------------------
        Read_Input_Type(state.T_air_type, state.T_air_ID, T_air, OK, filename=T_air_file)
        if logical_not(OK):    
            return _ret()
        #-----------------------------------------------------------------
        Read_Input_Type(state.RH_type, state.RH_ID, RH, OK, filename=RH_file)
        if logical_not(OK):    
            return _ret()
        #-----------------------------------------------------------------
        Read_Input_Type(state.C_type, state.C_ID, C, OK, filename=C_file)
        if logical_not(OK):    
            return _ret()
        #-----------------------------------------------------------------
        Read_Input_Type(state.F_type, state.F_ID, F, OK, filename=F_file)
        if logical_not(OK):    
            return _ret()
        #-----------------------------------------------------------------
        Read_Input_Type(state.T_surf_type, state.T_surf_ID, T_surf, OK, filename=T_surf_file)
        if logical_not(OK):    
            return _ret()
        #-----------------------------------------------------------------
        Read_Input_Type(state.em_surf_type, state.em_surf_ID, em_surf, OK, filename=em_surf_file)
        if logical_not(OK):    
            return _ret()
        
        #--------------------------------------------------------
        #Get number of frames in new RTS file from time settings
        #--------------------------------------------------------
        #Get start & stop times as decimal Julian days
        #----------------------------------------------
        Julian_start_day = Julian_Day(start_month, start_day, start_hour)
        Julian_stop_day = Julian_Day(stop_month, stop_day, stop_hour)
        timestep_JD = (timestep / float64(24))   #([hours] -> [days])
        nf_match = int32((Julian_stop_day - Julian_start_day) / timestep_JD)
        #-----------------------------------------------
        #Check if number of frames in RTS files matches
        #the number obtained from the time settings
        #-----------------------------------------------
        if (state.T_air_type == 3):    
            T_air_file, RTI_file, nf_match, OK = Check_Number_of_Frames(T_air_file, RTI_file, nf_match, OK)
        if logical_not(OK):    
            return _ret()
        #----------------------------------------------------------------
        if (state.RH_type == 3):    
            RH_file, RTI_file, nf_match, OK = Check_Number_of_Frames(RH_file, RTI_file, nf_match, OK)
        if logical_not(OK):    
            return _ret()
        #----------------------------------------------------------------
        if (state.C_type == 3):    
            C_file, RTI_file, nf_match, OK = Check_Number_of_Frames(C_file, RTI_file, nf_match, OK)
        if logical_not(OK):    
            return _ret()
        #----------------------------------------------------------------
        if (state.F_type == 3):    
            F_file, RTI_file, nf_match, OK = Check_Number_of_Frames(F_file, RTI_file, nf_match, OK)
        if logical_not(OK):    
            return _ret()
        #----------------------------------------------------------------
        if (state.T_surf_type == 3):    
            T_surf_file, RTI_file, nf_match, OK = Check_Number_of_Frames(T_surf_file, RTI_file, nf_match, OK)
        if logical_not(OK):    
            return _ret()
        #----------------------------------------------------------------
        if (state.em_surf_type == 3):    
            em_surf_file, RTI_file, nf_match, OK = Check_Number_of_Frames(em_surf_file, RTI_file, nf_match, OK)
        if logical_not(OK):    
            return _ret()
        
        #---------------------------------------
        #Call routines that create the RTS file
        #---------------------------------------
        widget_control(event.ID, sensitive=0)    #(disable button)
        Qn_LW_RTS_file, e_air_RTS_file, em_air_RTS_file, RTI_file, start_month, start_day, start_hour, stop_month, stop_day, stop_hour, timestep, T_air, T_air_file, state.T_air_type, RH, RH_file, state.RH_type, C, C_file, state.C_type, F, F_file, state.F_type, T_surf, T_surf_file, state.T_surf_type, em_surf, em_surf_file, state.em_surf_type, state.msg_box_ID = Make_RTS_File_for_Qnet_LW(Qn_LW_RTS_file, e_air_RTS_file, em_air_RTS_file, RTI_file, start_month, start_day, start_hour, stop_month, stop_day, stop_hour, timestep, T_air, T_air_file, state.T_air_type, RH, RH_file, state.RH_type, C, C_file, state.C_type, F, F_file, state.F_type, T_surf, T_surf_file, state.T_surf_type, em_surf, em_surf_file, state.em_surf_type, MSG_BOX_ID=state.msg_box_ID)
        widget_control(event.ID, sensitive=1)    #(enable button)
        
        #-------------------------
        #Show a "finished" dialog
        #-------------------------
        msg = array(['Finished creating RTS files for', 'Qnet_LW, e_air and em_air.', ' '])
        result = GUI_Message(msg, INFO=True, TITLE="Finished")
        #*** Close_Dialog, event.top
        
        #***********
    elif uvalue == 'HELP':    
        Show_HTML_Help('longwave_calc.htm')
        
        #************
    elif uvalue == 'CLOSE':    
        Close_Dialog(event.top)
        
    else:    
        dum = int16(0)
    
    
    if logical_and((uvalue != 'CLOSE'), (uvalue != 'START')):    
        widget_control(event.top, set_uvalue=state)
    
    
    return _ret()
#  GUI_Make_Qnet_LW_File_event
#****************************************************************
def GUI_Make_Qnet_LW_File(leader):

#-----------
#Error trap
#-----------
    n_params = 1
    def _ret():  return leader
    
    status = No_Catch()
    OK = Check_Error_Status(status)
    if logical_not(OK):    
        return _ret()
    
    if (idl_func.n_elements(leader) == 0):    
        leader = int32(0)
    
    #-----------------------------------
    #Get current values from main state
    #-----------------------------------
    Get_TLB_State(leader, mstate, ALIVE)
    if logical_not(ALIVE):    
        return _ret()
    #----------------------------------------
    prefix = mstate.run_vars.prefix
    Qn_LW_RTS_file = prefix + '_Qn-LW.rts'
    e_air_RTS_file = prefix + '_e-air.rts'
    em_air_RTS_file = prefix + '_em-air.rts'
    RTI_file = prefix + '.rti'
    #-----------------------------------------
    start_day_str = ' 1 '
    start_hour_str = ' 0.0 '
    stop_day_str = ' 2 '
    stop_hour_str = ' 0.0 '
    timestep_str = ' 1.0 '
    
    #------------------------------------
    #Structure to store selected options
    #------------------------------------
    state = idl_func.bunch(leader_ID=leader, msg_box_ID=int32(0), RTI_file_ID=int32(0), Qn_LW_RTS_file_ID=int32(0), e_air_RTS_file_ID=int32(0), em_air_RTS_file_ID=int32(0), start_month=1, start_day_ID=int32(0), start_hour_ID=int32(0), stop_month=1, stop_day_ID=int32(0), stop_hour_ID=int32(0), timestep_ID=int32(0), T_air_ID=int32(0), T_air_type=uint8(3), RH_ID=int32(0), RH_type=uint8(0), C_ID=int32(0), C_type=uint8(0), F_ID=int32(0), F_type=uint8(0), T_surf_ID=int32(0), T_surf_type=uint8(3), em_surf_ID=int32(0), em_surf_type=uint8(0))
    
    ngap = int16(6)
    XS = int16(24)
    months = array([' January ', ' February ', ' March ', ' April ', ' May ', ' June ', ' July ', ' August ', ' September ', ' October ', ' November ', ' December '])
    #** Get_Time_Zone_List, time_zones
    
    #-----------------
    #Main base widget
    #-----------------
    Create_TLB(MB, TITLE='Make RTS File for Qnet_LW', COLUMN=True, LEADER=leader)
    B1 = widget_base(MB, COLUMN=True, SPACE=1, FRAME=True)
    B2 = widget_base(MB, COLUMN=True, SPACE=1, FRAME=True)
    B3 = widget_base(MB, COLUMN=True, SPACE=1, FRAME=True)
    BOT = widget_base(MB, ROW=True, SPACE=1)
    
    #---------------------------
    #Get the start day and hour
    #---------------------------
    ST = widget_base(B1, ROW=True, SPACE=1)
    SM = widget_base(ST, ROW=True, SPACE=ngap)
    SM1 = widget_label(SM, VALUE='Start Month: ')
    SM2 = widget_droplist(SM, VALUE=months, UVALUE='START_MONTH')
    widget_control(SM2, set_droplist_select=0)
    #SM2 = widget_text(SM, VALUE=start_month_str, UVALUE='NONE', $
    #                  /EDITABLE, XSIZE=5)
    #----------------------------------------------------------------
    SD = widget_base(ST, ROW=True, SPACE=ngap)
    SD1 = widget_label(SD, VALUE='Day:')
    SD2 = widget_text(SD, VALUE=start_day_str, UVALUE='NONE', EDITABLE=True, XSIZE=5)
    state.start_day_ID = SD2
    #----------------------------------------------------------------
    SH = widget_base(ST, ROW=True, SPACE=ngap)
    SH1 = widget_label(SH, VALUE='Hour:')
    SH2 = widget_text(SH, VALUE=start_hour_str, UVALUE='NONE', EDITABLE=True, XSIZE=5)
    state.start_hour_ID = SH2
    
    #--------------------------
    #Get the stop day and hour
    #--------------------------
    ET = widget_base(B1, ROW=True, SPACE=1)
    EM = widget_base(ET, ROW=True, SPACE=ngap)
    EM1 = widget_label(EM, VALUE='Stop Month: ')
    EM2 = widget_droplist(EM, VALUE=months, UVALUE='STOP_MONTH')
    widget_control(EM2, set_droplist_select=0)
    #EM2 = widget_text(EM, VALUE=stop_month_str, UVALUE='NONE', $
    #                  /EDITABLE, XSIZE=5)
    #---------------------------------------------------------------
    ED = widget_base(ET, ROW=True, SPACE=ngap)
    ED1 = widget_label(ED, VALUE='Day:')
    ED2 = widget_text(ED, VALUE=stop_day_str, UVALUE='NONE', EDITABLE=True, XSIZE=5)
    state.stop_day_ID = ED2
    #---------------------------------------------------------------
    EH = widget_base(ET, ROW=True, SPACE=ngap)
    EH1 = widget_label(EH, VALUE='Hour:')
    EH2 = widget_text(EH, VALUE=stop_hour_str, UVALUE='NONE', EDITABLE=True, XSIZE=5)
    state.stop_hour_ID = EH2
    
    #-----------------
    #Get the timestep
    #-----------------
    TS0 = widget_base(B1, ROW=True, SPACE=1)
    TS = widget_base(TS0, ROW=True, SPACE=ngap)
    TS1 = widget_label(TS, VALUE='Timestep: ')
    TS2 = widget_text(TS, VALUE=timestep_str, UVALUE='NONE', EDITABLE=True, XSIZE=6)
    TS3 = widget_label(TS, VALUE='[hours] ')
    state.timestep_ID = TS2
    
    #------------------
    #Get the time zone
    #------------------
    #TZ = widget_base(TS0, /ROW, SPACE=ngap)
    #  TZ1 = widget_label(TZ, VALUE='  Time zone: ')
    #  TZ2 = widget_droplist(TZ, VALUE=time_zones, UVALUE='TIME_ZONE')
    #  widget_control, TZ2, set_droplist_select=12
    
    #---------------------
    #Add some blank space
    #----------------------
    #** P0 = widget_label(B1, VALUE=' ')
    
    #------------------
    #Align the widgets
    #------------------
    Align_Text_Boxes(array([SM1, EM1, TS1]))
    Align_Text_Boxes(array([SM2, EM2]))
    
    
    #------------------------
    #Get the input variables
    #------------------------
    ngap = int16(6)
    XS = int16(22)
    types = Model_Input_Types()
    #---------------------------
    RH_str = '0.30'
    C_str = '0.0'
    F_str = '0.0'
    #em_air_str = '0.70'
    em_surf_str = '0.98'
    #------------------------------------
    T_air_str = prefix + '_Tair.rts'
    T_surf_str = prefix + '_Tsurf.rts'
    
    #-------------------
    #Get the parameters
    #-------------------
    A1 = widget_base(B2, ROW=True, SPACE=ngap)
    A11 = widget_label(A1, VALUE='Variable: ', UVALUE='NONE')
    A12 = widget_label(A1, VALUE='Type: ', UVALUE='NONE')
    A13 = widget_label(A1, VALUE='Scalar or Grid Filename: ', UVALUE='NONE')
    A14 = widget_label(A1, VALUE='Units: ', UVALUE='NONE')
    #--------------------------------------------------------------------------
    TA = widget_base(B2, ROW=True, SPACE=ngap)
    TA1 = widget_label(TA, VALUE='T_air: ', UVALUE='NONE')
    TA2 = widget_droplist(TA, VALUE=types, UVALUE='T_AIR_TYPE')
    TA3 = widget_text(TA, VALUE=T_air_str, UVALUE='NONE', EDITABLE=True, XSIZE=XS)
    TA4 = widget_label(TA, VALUE='[deg C]', UVALUE='NONE')
    state.T_air_ID = TA3
    widget_control(TA2, set_droplist_select=3)
    #--------------------------------------------------------------------------
    #EA = widget_base(B2, /ROW, SPACE=ngap)
    #  EA1 = widget_label(EA, VALUE='em_air: ', UVALUE='NONE')
    #  EA2 = widget_droplist(EA, VALUE=types, UVALUE='EM_AIR_TYPE')
    #  EA3 = widget_text(EA, VALUE=em_air_str, UVALUE='NONE', /EDITABLE, XSIZE=XS)
    #  EA4 = widget_label(EA, VALUE='[none]', UVALUE='NONE')
    #  state.em_air_ID = EA3
    #--------------------------------------------------------------------------
    RH = widget_base(B2, ROW=True, SPACE=ngap)
    RH1 = widget_label(RH, VALUE='RH: ', UVALUE='NONE')
    RH2 = widget_droplist(RH, VALUE=types, UVALUE='RH_TYPE')
    RH3 = widget_text(RH, VALUE=RH_str, UVALUE='NONE', EDITABLE=True, XSIZE=XS)
    RH4 = widget_label(RH, VALUE='[none]  in [0,1]', UVALUE='NONE')
    state.RH_ID = RH3
    #--------------------------------------------------------------------------
    CC = widget_base(B2, ROW=True, SPACE=ngap)
    CC1 = widget_label(CC, VALUE='cloud frac.: ', UVALUE='NONE')
    CC2 = widget_droplist(CC, VALUE=types, UVALUE='C_TYPE')
    CC3 = widget_text(CC, VALUE=C_str, UVALUE='NONE', EDITABLE=True, XSIZE=XS)
    CC4 = widget_label(CC, VALUE='[none]  in [0,1]', UVALUE='NONE')
    state.C_ID = CC3
    #--------------------------------------------------------------------------
    FF = widget_base(B2, ROW=True, SPACE=ngap)
    FF1 = widget_label(FF, VALUE='canopy frac.: ', UVALUE='NONE')
    FF2 = widget_droplist(FF, VALUE=types, UVALUE='F_TYPE')
    FF3 = widget_text(FF, VALUE=F_str, UVALUE='NONE', EDITABLE=True, XSIZE=XS)
    FF4 = widget_label(FF, VALUE='[none]  in [0,1]', UVALUE='NONE')
    state.F_ID = FF3
    #--------------------------------------------------------------------------
    TS = widget_base(B2, ROW=True, SPACE=ngap)
    TS1 = widget_label(TS, VALUE='T_surf: ', UVALUE='NONE')
    TS2 = widget_droplist(TS, VALUE=types, UVALUE='T_SURF_TYPE')
    TS3 = widget_text(TS, VALUE=T_surf_str, UVALUE='NONE', EDITABLE=True, XSIZE=XS)
    TS4 = widget_label(TS, VALUE='[deg C]', UVALUE='NONE')
    state.T_surf_ID = TS3
    widget_control(TS2, set_droplist_select=3)
    #--------------------------------------------------------------------------
    ES = widget_base(B2, ROW=True, SPACE=ngap)
    ES1 = widget_label(ES, VALUE='em_surf: ', UVALUE='NONE')
    ES2 = widget_droplist(ES, VALUE=types, UVALUE='EM_SURF_TYPE')
    ES3 = widget_text(ES, VALUE=em_surf_str, UVALUE='NONE', EDITABLE=True, XSIZE=XS)
    ES4 = widget_label(ES, VALUE='[none]', UVALUE='NONE')
    state.em_surf_ID = ES3
    #---------------------
    #Add some blank space
    #---------------------
    P1 = widget_label(B2, VALUE=' ')
    
    #------------------
    #Align the widgets
    #------------------
    Align_Text_Boxes(array([A11, TA1, RH1, CC1, FF1, TS1, ES1]))
    Align_Text_Boxes(array([A12, TA2, RH2, CC2, FF2, TS2, ES2]))
    Align_Text_Boxes(array([A13, TA3, RH3, CC3, FF3, TS3, ES3]))
    
    
    #-------------------------------------
    #Get name of new RTS file for Qnet_LW
    #-------------------------------------
    NF = widget_base(B3, ROW=True, SPACE=ngap)
    NF1 = widget_label(NF, VALUE='New RTS file for Qnet_LW: ')
    NF2 = widget_text(NF, VALUE=Qn_LW_RTS_file, UVALUE='NONE', EDITABLE=True, XSIZE=XS)
    NF3 = widget_label(NF, VALUE=' [W / m^2] ')
    state.Qn_LW_RTS_file_ID = NF2
    
    #-----------------------------------
    #Get name of new RTS file for e_air
    #-----------------------------------
    EA = widget_base(B3, ROW=True, SPACE=ngap)
    EA1 = widget_label(EA, VALUE='New RTS file for e_air: ')
    EA2 = widget_text(EA, VALUE=e_air_RTS_file, UVALUE='NONE', EDITABLE=True, XSIZE=XS)
    EA3 = widget_label(EA, VALUE=' [mbar] ')
    state.e_air_RTS_file_ID = EA2
    
    #-------------------------------------
    #Get name of new RTS file for em_air
    #-------------------------------------
    EM = widget_base(B3, ROW=True, SPACE=ngap)
    EM1 = widget_label(EM, VALUE='New RTS file for em_air: ')
    EM2 = widget_text(EM, VALUE=em_air_RTS_file, UVALUE='NONE', EDITABLE=True, XSIZE=XS)
    EM3 = widget_label(EM, VALUE=' [none] ')
    state.em_air_RTS_file_ID = EM2
    
    #-----------------------
    #Get name of RTI file
    #-----------------------
    RF = widget_base(B3, ROW=True, SPACE=ngap)
    RF1 = widget_label(RF, VALUE='Existing RTI file: ')
    RF2 = widget_text(RF, VALUE=RTI_file, UVALUE='NONE', EDITABLE=True, XSIZE=XS)
    state.RTI_file_ID = RF2
    
    #---------------------
    #Add some blank space
    #---------------------
    P1 = widget_label(B3, VALUE=' ')
    
    #------------------
    #Align the widgets
    #------------------
    Align_Text_Boxes(array([NF1, EA1, EM1, RF1]))
    Align_Text_Boxes(array([NF2, EA2, EM2, RF2]))
    
    #------------------
    #Bottom button bar
    #------------------
    CW_Button_Bar(BOT, START=True, HELP=True, CLOSE=True)
    
    #---------------------
    #A status message box
    #---------------------
    MS = widget_base(BOT, ROW=True, SPACE=ngap)
    MS1 = widget_label(MS, VALUE='Status: ')
    MS2 = widget_text(MS, VALUE='Ready.', XSIZE=20)
    state.msg_box_ID = MS2
    
    #------------------------------------
    #Realize widgets and wait for events
    #------------------------------------
    XOFF = int16(480)
    Realize_TLB(MB, state, 'GUI_Make_Qnet_LW_File', XOFF=XOFF)
    
    
    return _ret()
#  GUI_Make_Qnet_LW_File 
#***************************************************************


