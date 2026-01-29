/******************************************************************************
 * @section DESCRIPTION
 *
 * Print library.
 *****************************************************************************/

#include <vic_driver_shared_all.h>

/******************************************************************************
 * @brief    Print cell data structure.
 *****************************************************************************/
void
print_cell_data(cell_data_struct *cell,
                size_t            nlayers,
                size_t            nfrost)
{
    size_t i;

    // Print state variables
    fprintf(LOG_DEST, "cell_data - states:\n");
    fprintf(LOG_DEST, "\tRa_over :");
    for (i = 0; i < 3; i++) {
        fprintf(LOG_DEST, "\t%f", cell->Ra_over[i]);
    }
    fprintf(LOG_DEST, "\tRa_sub :");
    for (i = 0; i < 3; i++) {
        fprintf(LOG_DEST, "\t%f", cell->Ra_sub[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tasat        : %f\n", cell->asat);
    fprintf(LOG_DEST, "\tCLitter     : %f\n", cell->CLitter);
    fprintf(LOG_DEST, "\tCInter      : %f\n", cell->CInter);
    fprintf(LOG_DEST, "\tCSlow       : %f\n", cell->CSlow);

    fprintf(LOG_DEST, "\trootmoist   : %f\n", cell->rootmoist);
    fprintf(LOG_DEST, "\twetness     : %f\n", cell->wetness);
    fprintf(LOG_DEST, "\tzwt         : %f\n", cell->zwt);
    fprintf(LOG_DEST, "\tzwt_lumped  : %f\n", cell->zwt_lumped);

    // Print fluxes
    fprintf(LOG_DEST, "cell_data - fluxes:\n");
    fprintf(LOG_DEST, "\tbaseflow    : %f\n", cell->baseflow);
    fprintf(LOG_DEST, "\tevap        : %f\n", cell->evap);
    fprintf(LOG_DEST, "\tesoil       : %f\n", cell->esoil);
    fprintf(LOG_DEST, "\ttransp      : %f\n", cell->transp);
    fprintf(LOG_DEST, "\tdewsoil     : %f\n", cell->dewsoil);
    fprintf(LOG_DEST, "\tcanopyevap  : %f\n", cell->canopyevap);
    fprintf(LOG_DEST, "\tvapor_grnd  : %f\n", cell->vapor_grnd);
    fprintf(LOG_DEST, "\tconden_grnd : %f\n", cell->conden_grnd);
    fprintf(LOG_DEST, "\tsnowfrost   : %f\n", cell->snowfrost);
    fprintf(LOG_DEST, "\tsnow_sublim : %f\n", cell->snow_sublim);
    fprintf(LOG_DEST, "\tsoil_inflow : %f\n", cell->soil_inflow);
    fprintf(LOG_DEST, "\trunoff      : %f\n", cell->runoff);
    fprintf(LOG_DEST, "\tRhLitter    : %f\n", cell->RhLitter);
    fprintf(LOG_DEST, "\tRhLitter2Atm: %f\n", cell->RhLitter2Atm);
    fprintf(LOG_DEST, "\tRhInter     : %f\n", cell->RhInter);
    fprintf(LOG_DEST, "\tRhSlow      : %f\n", cell->RhSlow);
    fprintf(LOG_DEST, "\tRhTot       : %f\n", cell->RhTot);
}

/******************************************************************************
 * @brief    Print day-month-year structure.
 *****************************************************************************/
void
print_dmy(dmy_struct *dmy)
{
    fprintf(LOG_DEST, "dmy:\n");
    fprintf(LOG_DEST, "\tday        : %hu\n", dmy->day);
    fprintf(LOG_DEST, "\tday_in_year: %hu\n", dmy->day_in_year);
    fprintf(LOG_DEST, "\tseconds    : %u\n", dmy->dayseconds);
    fprintf(LOG_DEST, "\tmonth      : %hu\n", dmy->month);
    fprintf(LOG_DEST, "\tyear       : %u\n", dmy->year);
}

/******************************************************************************
 * @brief    Print dmy structure as one string.
 *****************************************************************************/
void
sprint_dmy(char       *str,
           dmy_struct *dmy)
{
    sprintf(str,
            "dmy:\n"
            "\tday         : %hu\n"
            "\tday_in_year : %hu\n"
            "\tseconds     : %u\n"
            "\tmonth       : %hu\n"
            "\tyear        : %u\n",
            dmy->day, dmy->day_in_year, dmy->dayseconds, dmy->month, dmy->year);
}

/******************************************************************************
 * @brief    Print energy balance structure.
 *****************************************************************************/
void
print_energy_bal(energy_bal_struct *eb,
                 size_t             nnodes,
                 size_t             nfronts)
{
    size_t i;

    // Print energy_bal - state variables
    fprintf(LOG_DEST, "energy_bal - states:\n");
    fprintf(LOG_DEST, "\tTcanopy    : %f\n", eb->Tcanopy);
    fprintf(LOG_DEST, "\tTair       : %f\n", eb->Tair);
    fprintf(LOG_DEST, "\tTsurf      : %f\n", eb->Tsurf);
    fprintf(LOG_DEST, "\tTgrnd      : %f\n", eb->Tgrnd);
    fprintf(LOG_DEST, "\tTupper     : %f\n", eb->Tupper);
    fprintf(LOG_DEST, "\tTlower     : %f\n", eb->Tlower);
    fprintf(LOG_DEST, "\tFrozenGrnd : %d\n", eb->FrozenGrnd);
    fprintf(LOG_DEST, "\tFrozenOver : %d\n", eb->FrozenOver);
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tCs_node          :");
    for (i = 0; i < nnodes; i++) {
        fprintf(LOG_DEST, "\t%f", eb->Cs_node[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tkappa_node       :");
    for (i = 0; i < nnodes; i++) {
        fprintf(LOG_DEST, "\t%f", eb->kappa_node[i]);
    }
    fprintf(LOG_DEST, "\tT                :");
    for (i = 0; i < nnodes; i++) {
        fprintf(LOG_DEST, "\t%f", eb->T[i]);
    }
    // Print energy_bal - fluxes
    fprintf(LOG_DEST, "energy_bal - fluxes:\n");
    fprintf(LOG_DEST, "\tAPAR_sunlit         : %f\n", eb->APAR_sunlit);
    fprintf(LOG_DEST, "\tAPAR_shade          : %f\n", eb->APAR_shade);
    fprintf(LOG_DEST, "\tadvection           : %f\n", eb->advection);
    fprintf(LOG_DEST, "\tAdvectSub           : %f\n", eb->AdvectSub);
    fprintf(LOG_DEST, "\tAdvectGrnd          : %f\n", eb->AdvectGrnd);
    fprintf(LOG_DEST, "\tAdvectOver          : %f\n", eb->AdvectOver);
    fprintf(LOG_DEST, "\tgrnd_flux           : %f\n", eb->grnd_flux);
    fprintf(LOG_DEST, "\tGroundSub           : %f\n", eb->GroundSub);
    fprintf(LOG_DEST, "\tGroundGrnd          : %f\n", eb->GroundGrnd);
    fprintf(LOG_DEST, "\tlatent              : %f\n", eb->latent);
    fprintf(LOG_DEST, "\tLatentSub           : %f\n", eb->LatentSub);
    fprintf(LOG_DEST, "\tLatentGrnd          : %f\n", eb->LatentGrnd);
    fprintf(LOG_DEST, "\tsensible            : %f\n", eb->sensible);
    fprintf(LOG_DEST, "\tSensibleSub         : %f\n", eb->SensibleSub);
    fprintf(LOG_DEST, "\tSensibleGrnd        : %f\n", eb->SensibleGrnd);
    fprintf(LOG_DEST, "\tSensibleOver        : %f\n", eb->SensibleOver);
    fprintf(LOG_DEST, "\tLatentVapOver       : %f\n", eb->LatentVapOver);
    fprintf(LOG_DEST, "\tLatentVapGrnd       : %f\n", eb->LatentVapGrnd);
    fprintf(LOG_DEST, "\tReflShortSurf       : %f\n", eb->ReflShortSurf);
    fprintf(LOG_DEST, "\tReflShortGrnd       : %f\n", eb->ReflShortGrnd);
    fprintf(LOG_DEST, "\tReflShortSub        : %f\n", eb->ReflShortSub);
    fprintf(LOG_DEST, "\tEmissLongSub        : %f\n", eb->EmissLongSub);
    fprintf(LOG_DEST, "\tEmissLongGrnd       : %f\n", eb->EmissLongGrnd);
    fprintf(LOG_DEST, "\tEmissLongSurf       : %f\n", eb->EmissLongSurf);
    fprintf(LOG_DEST, "\tNetLongSurf         : %f\n", eb->NetLongSurf);
    fprintf(LOG_DEST, "\tNetLongGrnd         : %f\n", eb->NetLongGrnd);
    fprintf(LOG_DEST, "\tNetLongOver         : %f\n", eb->NetLongOver);
    fprintf(LOG_DEST, "\tNetLongSub          : %f\n", eb->NetLongSub);
    fprintf(LOG_DEST, "\tNetShortGrnd        : %f\n", eb->NetShortGrnd);
    fprintf(LOG_DEST, "\tNetShortSurf        : %f\n", eb->NetShortSurf);
    fprintf(LOG_DEST, "\tNetShortSub         : %f\n", eb->NetShortSub);
}

/******************************************************************************
 * @brief    Print forcing type structure.
 *****************************************************************************/
void
print_force_type(force_type_struct *force_type)
{
    fprintf(LOG_DEST, "force_type:\n");
    fprintf(LOG_DEST, "\tSIGNED    : %d\n", force_type->SIGNED);
    fprintf(LOG_DEST, "\tSUPPLIED  : %d\n", force_type->SUPPLIED);
    fprintf(LOG_DEST, "\tmultiplier: %f\n", force_type->multiplier);
}

/******************************************************************************
 * @brief    Print global parameters structure.
 *****************************************************************************/
void
print_global_param(global_param_struct *gp)
{
    size_t i;

    fprintf(LOG_DEST, "global_param:\n");
    fprintf(LOG_DEST, "\twind_h              : %.4f\n", gp->wind_h);
    fprintf(LOG_DEST, "\tresolution          : %.4f\n", gp->resolution);
    fprintf(LOG_DEST, "\tdt                  : %.4f\n", gp->dt);
    fprintf(LOG_DEST, "\tsnow_dt             : %.4f\n", gp->snow_dt);
    fprintf(LOG_DEST, "\trunoff_dt           : %.4f\n", gp->runoff_dt);
    fprintf(LOG_DEST, "\tmodel_steps_per_day : %zu\n", gp->model_steps_per_day);
    fprintf(LOG_DEST, "\tsnow_steps_per_day  : %zu\n", gp->snow_steps_per_day);
    fprintf(LOG_DEST, "\trunoff_steps_per_day: %zu\n",
            gp->runoff_steps_per_day);
    fprintf(LOG_DEST, "\tendday              : %hu\n", gp->endday);
    fprintf(LOG_DEST, "\tendmonth            : %hu\n", gp->endmonth);
    fprintf(LOG_DEST, "\tendyear             : %hu\n", gp->endyear);
    for (i = 0; i < 2; i++) {
        fprintf(LOG_DEST, "\tforceday[%zd]        : %hu\n", i, gp->forceday[i]);
        fprintf(LOG_DEST, "\tforcesec[%zd]        : %u\n", i, gp->forcesec[i]);
        fprintf(LOG_DEST, "\tforcemonth[%zd]      : %hu\n", i,
                gp->forcemonth[i]);
        fprintf(LOG_DEST, "\tforceoffset[%zd]     : %hu\n", i,
                gp->forceoffset[i]);
        fprintf(LOG_DEST, "\tforceskip[%zd]       : %u\n", i, gp->forceskip[i]);
        fprintf(LOG_DEST, "\tforceyear[%zd]       : %hu\n", i,
                gp->forceyear[i]);
    }
    fprintf(LOG_DEST, "\tnrecs               : %zu\n", gp->nrecs);
    fprintf(LOG_DEST, "\tstartday            : %hu\n", gp->startday);
    fprintf(LOG_DEST, "\tstartsec            : %u\n", gp->startsec);
    fprintf(LOG_DEST, "\tstartmonth          : %hu\n", gp->startmonth);
    fprintf(LOG_DEST, "\tstartyear           : %hu\n", gp->startyear);
    fprintf(LOG_DEST, "\tstateday            : %hu\n", gp->stateday);
    fprintf(LOG_DEST, "\tstatemonth          : %hu\n", gp->statemonth);
    fprintf(LOG_DEST, "\tstateyear           : %hu\n", gp->stateyear);
    fprintf(LOG_DEST, "\tstatesec            : %u\n", gp->statesec);
}

/******************************************************************************
 * @brief    Print options structure.
 *****************************************************************************/
void
print_option(option_struct *option)
{
    fprintf(LOG_DEST, "option:\n");
    fprintf(LOG_DEST, "\tBLOWING              : %s\n",
            option->BLOWING ? "true" : "false");
    fprintf(LOG_DEST, "\tBLOWING_VAR_THRESHOLD: %s\n",
            option->BLOWING_VAR_THRESHOLD ? "true" : "false");
    fprintf(LOG_DEST, "\tBLOWING_CALC_PROB    : %s\n",
            option->BLOWING_CALC_PROB ? "true" : "false");
    fprintf(LOG_DEST, "\tBLOWING_SIMPLE       : %s\n",
            option->BLOWING_SIMPLE ? "true" : "false");
    fprintf(LOG_DEST, "\tBLOWING_FETCH        : %s\n",
            option->BLOWING_FETCH ? "true" : "false");
    fprintf(LOG_DEST, "\tBLOWING_SPATIAL_WIND : %s\n",
            option->BLOWING_SPATIAL_WIND ? "true" : "false");
    fprintf(LOG_DEST, "\tCARBON               : %s\n",
            option->CARBON ? "true" : "false");
    fprintf(LOG_DEST, "\tCONTINUEONERROR      : %s\n",
            option->CONTINUEONERROR ? "true" : "false");
    fprintf(LOG_DEST, "\tCORRPREC             : %s\n",
            option->CORRPREC ? "true" : "false");
    fprintf(LOG_DEST, "\tEQUAL_AREA           : %s\n",
            option->EQUAL_AREA ? "true" : "false");
    fprintf(LOG_DEST, "\tEXP_TRANS            : %s\n",
            option->EXP_TRANS ? "true" : "false");
    fprintf(LOG_DEST, "\tFROZEN_SOIL          : %s\n",
            option->FROZEN_SOIL ? "true" : "false");
    fprintf(LOG_DEST, "\tNcanopy              : %zu\n", option->Ncanopy);
    fprintf(LOG_DEST, "\tNlayer               : %zu\n", option->Nlayer);
    fprintf(LOG_DEST, "\tNnode                : %zu\n", option->Nnode);
    fprintf(LOG_DEST, "\tNOFLUX               : %s\n",
            option->NOFLUX ? "true" : "false");
    fprintf(LOG_DEST, "\tNVEGTYPES            : %zu\n", option->NVEGTYPES);
    fprintf(LOG_DEST, "\tRC_MODE              : %d\n", option->RC_MODE);
    fprintf(LOG_DEST, "\tSNOW_DENSITY         : %d\n", option->SNOW_DENSITY);
    fprintf(LOG_DEST, "\tSNOW_BAND            : %zu\n", option->SNOW_BAND);
    fprintf(LOG_DEST, "\tTFALLBACK            : %s\n",
            option->TFALLBACK ? "true" : "false");
    fprintf(LOG_DEST, "\tBASEFLOW             : %d\n", option->BASEFLOW);
    fprintf(LOG_DEST, "\tGRID_DECIMAL         : %d\n", option->GRID_DECIMAL);
    fprintf(LOG_DEST, "\tVEGLIB_PHOTO         : %s\n",
            option->VEGLIB_PHOTO ? "true" : "false");
    fprintf(LOG_DEST, "\tVEGLIB_FCAN          : %s\n",
            option->VEGLIB_FCAN ? "true" : "false");
    fprintf(LOG_DEST, "\tVEGPARAM_ALB         : %s\n",
            option->VEGPARAM_ALB ? "true" : "false");
    fprintf(LOG_DEST, "\tVEGPARAM_LAI         : %s\n",
            option->VEGPARAM_LAI ? "true" : "false");
    fprintf(LOG_DEST, "\tVEGPARAM_FCAN        : %s\n",
            option->VEGPARAM_FCAN ? "true" : "false");
    fprintf(LOG_DEST, "\tALB_SRC              : %d\n", option->ALB_SRC);
    fprintf(LOG_DEST, "\tLAI_SRC              : %d\n", option->LAI_SRC);
    fprintf(LOG_DEST, "\tFCAN_SRC             : %d\n", option->FCAN_SRC);
    fprintf(LOG_DEST, "\tORGANIC_FRACT        : %s\n",
            option->ORGANIC_FRACT ? "true" : "false");
    fprintf(LOG_DEST, "\tBULK_DENSITY_COMB        : %s\n",
            option->BULK_DENSITY_COMB ? "true" : "false");
    fprintf(LOG_DEST, "\tSTATE_FORMAT         : %d\n", option->STATE_FORMAT);
    fprintf(LOG_DEST, "\tINIT_STATE           : %s\n",
            option->INIT_STATE ? "true" : "false");
    fprintf(LOG_DEST, "\tSAVE_STATE           : %s\n",
            option->SAVE_STATE ? "true" : "false");
    fprintf(LOG_DEST, "\tSTATENAME_CESM       : %s\n",
            option->STATENAME_CESM ? "true" : "false");
    fprintf(LOG_DEST, "\tNoutstreams          : %zu\n", option->Noutstreams);
}

/******************************************************************************
 * @brief    Print out data structure.
 *****************************************************************************/
void
print_out_data(double         **out_data,
               metadata_struct *metadata)
{
    size_t i;
    size_t j;

    fprintf(LOG_DEST, "out_data:\n");

    for (i = 0; i < N_OUTVAR_TYPES; i++) {
        fprintf(LOG_DEST, "\tvarname: %s\n", metadata[i].varname);
        fprintf(LOG_DEST, "\t\tnelem: %zu\n", metadata[i].nelem);
        fprintf(LOG_DEST, "\t\tdata:");
        for (j = 0; j < metadata[i].nelem; j++) {
            fprintf(LOG_DEST, "\t%.4f", out_data[i][j]);
        }
        fprintf(LOG_DEST, "\n");
    }
    fprintf(LOG_DEST, "\n");
}

/******************************************************************************
 * @brief    Print stream_file_struct.
 *****************************************************************************/
void
print_stream(stream_struct   *stream,
             metadata_struct *metadata)
{
    size_t       i;
    unsigned int varid;

    fprintf(LOG_DEST, "stream_file_struct:\n");

    fprintf(LOG_DEST, "\tprefix: %s\n", stream->prefix);
    fprintf(LOG_DEST, "\tfilename: %s\n", stream->filename);
    fprintf(LOG_DEST, "\tfh: %p\n", stream->fh);
    fprintf(LOG_DEST, "\tfile_format: %hu\n", stream->file_format);
    fprintf(LOG_DEST, "\tnvars: %zu\n", stream->nvars);
    fprintf(LOG_DEST, "\tngridcells: %zu\n", stream->ngridcells);
    fprintf(LOG_DEST, "\tagg_alarm:\n    ");
    print_alarm(&(stream->agg_alarm));
    fprintf(
        LOG_DEST,
        "\t# \tVARID        \tVARNAME \tTYPE \tMULT \tFORMAT        \tAGGTYPE\n");
    for (i = 0; i < stream->nvars; i++) {
        varid = stream->varid[i];
        fprintf(LOG_DEST, "\t%zu \t%u \t%20s \t%hu \t%f \t%10s \t%hu\n",
                i, varid, metadata[varid].varname,
                stream->type[i], stream->mult[i], stream->format[i],
                stream->aggtype[i]);
    }
    fprintf(LOG_DEST, "\taggdata shape: (%zu, %zu, nelem, 1)\n",
            stream->ngridcells, stream->nvars);

    fprintf(LOG_DEST, "\n");
}

/******************************************************************************
 * @brief    Print stream_file_struct.
 *****************************************************************************/
void
print_alarm(alarm_struct *alarm)
{
    fprintf(LOG_DEST, "alarm_struct:\n");
    fprintf(LOG_DEST, "\tcount: %u\n", alarm->count);
    fprintf(LOG_DEST, "\tfreq: %u\n", alarm->freq);
    fprintf(LOG_DEST, "\tnext_count: %d\n", alarm->next_count);
    fprintf(LOG_DEST, "\tnext_dmy: \n    ");
    print_dmy(&(alarm->next_dmy));
    fprintf(LOG_DEST, "\tn: %d\n", alarm->n);
    fprintf(LOG_DEST, "\tis_subdaily: %s\n",
            alarm->is_subdaily ? "true" : "false");

    fprintf(LOG_DEST, "\n");
}

/******************************************************************************
 * @brief    Print stream_file_struct.
 *****************************************************************************/
void
print_out_metadata(metadata_struct *metadata,
                   size_t           nvars)
{
    size_t i;

    fprintf(LOG_DEST, "metadata_struct: \n");

    for (i = 0; i < nvars; i++) {
        fprintf(LOG_DEST, "\t%s (%zu)\n", metadata[i].varname, i);
        fprintf(LOG_DEST, "\t\tlong_name: %s\n", metadata[i].long_name);
        fprintf(LOG_DEST, "\t\tunits: %s\n", metadata[i].units);
        fprintf(LOG_DEST, "\t\tdescription: %s\n", metadata[i].description);
        fprintf(LOG_DEST, "\t\tnelem: %zu\n", metadata[i].nelem);
    }
    fprintf(LOG_DEST, "\n");
}

/******************************************************************************
 * @brief    print param set structure.
 *****************************************************************************/
void
print_param_set(param_set_struct *param_set)
{
    size_t i;

    fprintf(LOG_DEST, "param_set:\n");
    for (i = 0; i < N_FORCING_TYPES; i++) {
        print_force_type(&(param_set->TYPE[i]));
    }
    fprintf(LOG_DEST, "\tFORCE_DT    : %.4f %.4f\n", param_set->FORCE_DT[0],
            param_set->FORCE_DT[1]);
    fprintf(LOG_DEST, "\tFORCE_ENDIAN: %d %d\n", param_set->FORCE_ENDIAN[0],
            param_set->FORCE_ENDIAN[1]);
    fprintf(LOG_DEST, "\tFORCE_FORMAT: %d %d\n", param_set->FORCE_FORMAT[0],
            param_set->FORCE_FORMAT[1]);
    fprintf(LOG_DEST, "\tFORCE_INDEX :\n");
    for (i = 0; i < N_FORCING_TYPES; i++) {
        fprintf(LOG_DEST, "\t\t%zd: %d %d\n", i, param_set->FORCE_INDEX[0][i],
                param_set->FORCE_INDEX[1][i]);
    }
    fprintf(LOG_DEST, "\tN_TYPES     : %zu %zu\n", param_set->N_TYPES[0],
            param_set->N_TYPES[1]);
}

/******************************************************************************
 * @brief    Print model parameters.
 *****************************************************************************/
void
print_parameters(parameters_struct *param)
{
    fprintf(LOG_DEST, "parameters:\n");
    fprintf(LOG_DEST, "\tLAPSE_RATE: %.4f\n", param->LAPSE_RATE);
    fprintf(LOG_DEST, "\tGAUGE_HEIGHT: %.4f\n", param->GAUGE_HEIGHT);
    fprintf(LOG_DEST, "\tHUGE_RESIST: %.4f\n", param->HUGE_RESIST);
    fprintf(LOG_DEST, "\tALBEDO_BARE_SOIL: %.4f\n", param->ALBEDO_BARE_SOIL);
    fprintf(LOG_DEST, "\tEMISS_GRND: %.4f\n", param->EMISS_GRND);
    fprintf(LOG_DEST, "\tEMISS_VEG: %.4f\n", param->EMISS_VEG);
    fprintf(LOG_DEST, "\tEMISS_ICE: %.4f\n", param->EMISS_ICE);
    fprintf(LOG_DEST, "\tEMISS_SNOW: %.4f\n", param->EMISS_SNOW);
    fprintf(LOG_DEST, "\tEMISS_H2O: %.4f\n", param->EMISS_H2O);
    fprintf(LOG_DEST, "\tSOIL_RESID_MOIST: %.4f\n", param->SOIL_RESID_MOIST);
    fprintf(LOG_DEST, "\tSOIL_SLAB_MOIST_FRACT: %.4f\n",
            param->SOIL_SLAB_MOIST_FRACT);
    fprintf(LOG_DEST, "\tVEG_LAI_SNOW_MULTIPLIER: %.4f\n",
            param->VEG_LAI_SNOW_MULTIPLIER);
    fprintf(LOG_DEST, "\tVEG_MIN_INTERCEPTION_STORAGE: %.4f\n",
            param->VEG_MIN_INTERCEPTION_STORAGE);
    fprintf(LOG_DEST, "\tVEG_LAI_WATER_FACTOR: %.4f\n",
            param->VEG_LAI_WATER_FACTOR);
    fprintf(LOG_DEST, "\tCANOPY_CLOSURE: %.4f\n", param->CANOPY_CLOSURE);
    fprintf(LOG_DEST, "\tCANOPY_RSMAX: %.4f\n", param->CANOPY_RSMAX);
    fprintf(LOG_DEST, "\tCANOPY_VPDMINFACTOR: %.4f\n",
            param->CANOPY_VPDMINFACTOR);
    fprintf(LOG_DEST, "\tSVP1: %.4f\n", param->SVP1);
    fprintf(LOG_DEST, "\tSVP2: %.4f\n", param->SVP2);
    fprintf(LOG_DEST, "\tSVP3: %.4f\n", param->SVP3);
    fprintf(LOG_DEST, "\tPHOTO_OMEGA: %.4f\n", param->PHOTO_OMEGA);
    fprintf(LOG_DEST, "\tPHOTO_LAIMAX: %.4f\n", param->PHOTO_LAIMAX);
    fprintf(LOG_DEST, "\tPHOTO_LAILIMIT: %.4f\n", param->PHOTO_LAILIMIT);
    fprintf(LOG_DEST, "\tPHOTO_LAIMIN: %.4f\n", param->PHOTO_LAIMIN);
    fprintf(LOG_DEST, "\tPHOTO_EPAR: %.4f\n", param->PHOTO_EPAR);
    fprintf(LOG_DEST, "\tPHOTO_FCMAX: %.4f\n", param->PHOTO_FCMAX);
    fprintf(LOG_DEST, "\tPHOTO_FCMIN: %.4f\n", param->PHOTO_FCMIN);
    fprintf(LOG_DEST, "\tPHOTO_ZENITHMIN: %.4f\n", param->PHOTO_ZENITHMIN);
    fprintf(LOG_DEST, "\tPHOTO_ZENITHMINPAR: %.4f\n",
            param->PHOTO_ZENITHMINPAR);
    fprintf(LOG_DEST, "\tPHOTO_ALBSOIPARMIN: %.4f\n",
            param->PHOTO_ALBSOIPARMIN);
    fprintf(LOG_DEST, "\tPHOTO_MINMAXETRANS: %.4f\n",
            param->PHOTO_MINMAXETRANS);
    fprintf(LOG_DEST, "\tPHOTO_MINSTOMCOND: %.4f\n", param->PHOTO_MINSTOMCOND);
    fprintf(LOG_DEST, "\tPHOTO_FCI1C3: %.4f\n", param->PHOTO_FCI1C3);
    fprintf(LOG_DEST, "\tPHOTO_FCI1C4: %.4f\n", param->PHOTO_FCI1C4);
    fprintf(LOG_DEST, "\tPHOTO_OX: %.4f\n", param->PHOTO_OX);
    fprintf(LOG_DEST, "\tPHOTO_KC: %.4f\n", param->PHOTO_KC);
    fprintf(LOG_DEST, "\tPHOTO_KO: %.4f\n", param->PHOTO_KO);
    fprintf(LOG_DEST, "\tPHOTO_EC: %.4f\n", param->PHOTO_EC);
    fprintf(LOG_DEST, "\tPHOTO_EO: %.4f\n", param->PHOTO_EO);
    fprintf(LOG_DEST, "\tPHOTO_EV: %.4f\n", param->PHOTO_EV);
    fprintf(LOG_DEST, "\tPHOTO_ER: %.4f\n", param->PHOTO_ER);
    fprintf(LOG_DEST, "\tPHOTO_ALC3: %.4f\n", param->PHOTO_ALC3);
    fprintf(LOG_DEST, "\tPHOTO_FRDC3: %.4f\n", param->PHOTO_FRDC3);
    fprintf(LOG_DEST, "\tPHOTO_EK: %.4f\n", param->PHOTO_EK);
    fprintf(LOG_DEST, "\tPHOTO_ALC4: %.4f\n", param->PHOTO_ALC4);
    fprintf(LOG_DEST, "\tPHOTO_FRDC4: %.4f\n", param->PHOTO_FRDC4);
    fprintf(LOG_DEST, "\tPHOTO_THETA: %.4f\n", param->PHOTO_THETA);
    fprintf(LOG_DEST, "\tPHOTO_FRLEAF: %.4f\n", param->PHOTO_FRLEAF);
    fprintf(LOG_DEST, "\tPHOTO_FRGROWTH: %.4f\n", param->PHOTO_FRGROWTH);
    fprintf(LOG_DEST, "\tSRESP_E0_LT: %.4f\n", param->SRESP_E0_LT);
    fprintf(LOG_DEST, "\tSRESP_T0_LT: %.4f\n", param->SRESP_T0_LT);
    fprintf(LOG_DEST, "\tSRESP_WMINFM: %.4f\n", param->SRESP_WMINFM);
    fprintf(LOG_DEST, "\tSRESP_WMAXFM: %.4f\n", param->SRESP_WMAXFM);
    fprintf(LOG_DEST, "\tSRESP_WOPTFM: %.4f\n", param->SRESP_WOPTFM);
    fprintf(LOG_DEST, "\tSRESP_RHSAT: %.4f\n", param->SRESP_RHSAT);
    fprintf(LOG_DEST, "\tSRESP_RFACTOR: %.4f\n", param->SRESP_RFACTOR);
    fprintf(LOG_DEST, "\tSRESP_TAULITTER: %.4f\n", param->SRESP_TAULITTER);
    fprintf(LOG_DEST, "\tSRESP_TAUINTER: %.4f\n", param->SRESP_TAUINTER);
    fprintf(LOG_DEST, "\tSRESP_TAUSLOW: %.4f\n", param->SRESP_TAUSLOW);
    fprintf(LOG_DEST, "\tSRESP_FAIR: %.4f\n", param->SRESP_FAIR);
    fprintf(LOG_DEST, "\tSRESP_FINTER: %.4f\n", param->SRESP_FINTER);
    fprintf(LOG_DEST, "\tSNOW_MAX_SURFACE_SWE: %.4f\n",
            param->SNOW_MAX_SURFACE_SWE);
    fprintf(LOG_DEST, "\tSNOW_LIQUID_WATER_CAPACITY: %.4f\n",
            param->SNOW_LIQUID_WATER_CAPACITY);
    fprintf(LOG_DEST, "\tSNOW_NEW_SNOW_DENSITY: %.4f\n",
            param->SNOW_NEW_SNOW_DENSITY);
    fprintf(LOG_DEST, "\tSNOW_NEW_SNOW_DENS_MAX: %.4f\n",
            param->SNOW_NEW_SNOW_DENS_MAX);
    fprintf(LOG_DEST, "\tSNOW_DEPTH_THRES: %.12f\n",
            param->SNOW_DEPTH_THRES);
    fprintf(LOG_DEST, "\tSNOW_DENS_DMLIMIT: %.4f\n", param->SNOW_DENS_DMLIMIT);
    fprintf(LOG_DEST, "\tSNOW_DENS_MAX_CHANGE: %.4f\n",
            param->SNOW_DENS_MAX_CHANGE);
    fprintf(LOG_DEST, "\tSNOW_DENS_ETA0: %.4f\n", param->SNOW_DENS_ETA0);
    fprintf(LOG_DEST, "\tSNOW_DENS_C1: %.4f\n", param->SNOW_DENS_C1);
    fprintf(LOG_DEST, "\tSNOW_DENS_C2: %.4f\n", param->SNOW_DENS_C2);
    fprintf(LOG_DEST, "\tSNOW_DENS_C5: %.4f\n", param->SNOW_DENS_C5);
    fprintf(LOG_DEST, "\tSNOW_DENS_C6: %.4f\n", param->SNOW_DENS_C6);
    fprintf(LOG_DEST, "\tSNOW_DENS_F: %.4f\n", param->SNOW_DENS_F);
    fprintf(LOG_DEST, "\tSNOW_MIN_SWQ_EB_THRES: %.4f\n",
            param->SNOW_MIN_SWQ_EB_THRES);
    fprintf(LOG_DEST, "\tSNOW_A1: %.4f\n", param->SNOW_A1);
    fprintf(LOG_DEST, "\tSNOW_A2: %.4f\n", param->SNOW_A2);
    fprintf(LOG_DEST, "\tSNOW_L1: %.4f\n", param->SNOW_L1);
    fprintf(LOG_DEST, "\tSNOW_L2: %.4f\n", param->SNOW_L2);

    fprintf(LOG_DEST, "\tSNOW_TRACESNOW: %.4f\n", param->SNOW_TRACESNOW);
    fprintf(LOG_DEST, "\tSNOW_CONDUCT: %.4f\n", param->SNOW_CONDUCT);

    fprintf(LOG_DEST, "\tBLOWING_KA: %.4f\n", param->BLOWING_KA);
    fprintf(LOG_DEST, "\tBLOWING_CSALT: %.4f\n", param->BLOWING_CSALT);
    fprintf(LOG_DEST, "\tBLOWING_UTHRESH: %.4f\n", param->BLOWING_UTHRESH);
    fprintf(LOG_DEST, "\tBLOWING_KIN_VIS: %.4f\n", param->BLOWING_KIN_VIS);
    fprintf(LOG_DEST, "\tBLOWING_MAX_ITER: %d\n", param->BLOWING_MAX_ITER);
    fprintf(LOG_DEST, "\tBLOWING_K: %d\n", param->BLOWING_K);
    fprintf(LOG_DEST, "\tBLOWING_SETTLING: %.4f\n", param->BLOWING_SETTLING);
    fprintf(LOG_DEST, "\tBLOWING_NUMINCS: %d\n", param->BLOWING_NUMINCS);
    fprintf(LOG_DEST, "\tTREELINE_TEMPERATURE: %.4f\n",
            param->TREELINE_TEMPERATURE);
    fprintf(LOG_DEST, "\tSNOW_DT: %.4f\n", param->SNOW_DT);
    fprintf(LOG_DEST, "\tSURF_DT: %.4f\n", param->SURF_DT);
    fprintf(LOG_DEST, "\tSOIL_DT: %.4f\n", param->SOIL_DT);
    fprintf(LOG_DEST, "\tCANOPY_DT: %.4f\n", param->CANOPY_DT);
    fprintf(LOG_DEST, "\tCANOPY_VP: %.4f\n", param->CANOPY_VP);
    fprintf(LOG_DEST, "\tFROZEN_MAXITER: %d\n", param->FROZEN_MAXITER);
    fprintf(LOG_DEST, "\tMAX_ITER_GRND_CANOPY: %d\n",
            param->MAX_ITER_GRND_CANOPY);
    fprintf(LOG_DEST, "\tNEWT_RAPH_MAXTRIAL: %d\n", param->NEWT_RAPH_MAXTRIAL);
    fprintf(LOG_DEST, "\tNEWT_RAPH_TOLX: %.4f\n", param->NEWT_RAPH_TOLX);
    fprintf(LOG_DEST, "\tNEWT_RAPH_TOLF: %.4f\n", param->NEWT_RAPH_TOLF);
    fprintf(LOG_DEST, "\tNEWT_RAPH_R_MAX: %.4f\n", param->NEWT_RAPH_R_MAX);
    fprintf(LOG_DEST, "\tNEWT_RAPH_R_MIN: %.4f\n", param->NEWT_RAPH_R_MIN);
    fprintf(LOG_DEST, "\tNEWT_RAPH_RELAX1: %.4f\n", param->NEWT_RAPH_RELAX1);
    fprintf(LOG_DEST, "\tNEWT_RAPH_RELAX2: %.4f\n", param->NEWT_RAPH_RELAX2);
    fprintf(LOG_DEST, "\tNEWT_RAPH_RELAX3: %.4f\n", param->NEWT_RAPH_RELAX3);
    fprintf(LOG_DEST, "\tNEWT_RAPH_EPS2: %.4f\n", param->NEWT_RAPH_EPS2);
    fprintf(LOG_DEST, "\tROOT_BRENT_MAXTRIES: %d\n",
            param->ROOT_BRENT_MAXTRIES);
    fprintf(LOG_DEST, "\tROOT_BRENT_MAXITER: %d\n", param->ROOT_BRENT_MAXITER);
    fprintf(LOG_DEST, "\tROOT_BRENT_TSTEP: %.4f\n", param->ROOT_BRENT_TSTEP);
    fprintf(LOG_DEST, "\tROOT_BRENT_T: %.4f\n", param->ROOT_BRENT_T);
    fprintf(LOG_DEST, "\tFROZEN_MAXITER: %d\n", param->FROZEN_MAXITER);
    fprintf(LOG_DEST, "\tSNOW_NEW_SNOW_ALB: %.4f\n", param->SNOW_NEW_SNOW_ALB);
    fprintf(LOG_DEST, "\tSNOW_ALB_ACCUM_A: %.4f\n", param->SNOW_ALB_ACCUM_A);
    fprintf(LOG_DEST, "\tSNOW_ALB_ACCUM_B: %.4f\n", param->SNOW_ALB_ACCUM_B);
    fprintf(LOG_DEST, "\tSNOW_ALB_THAW_A: %.4f\n", param->SNOW_ALB_THAW_A);
    fprintf(LOG_DEST, "\tSNOW_ALB_THAW_B: %.4f\n", param->SNOW_ALB_THAW_B);
}

/******************************************************************************
 * @brief    Print save data structure.
 *****************************************************************************/
void
print_save_data(save_data_struct *save)
{
    fprintf(LOG_DEST, "save_data:\n");
    fprintf(LOG_DEST, "\ttotal_moist_storage: %.4f\n",
            save->total_moist_storage);
    fprintf(LOG_DEST, "\ttotal_soil_moist: %.4f\n", save->total_soil_moist);
    fprintf(LOG_DEST, "\tsurfstor: %.4f\n", save->surfstor);
    fprintf(LOG_DEST, "\tswe: %.4f\n", save->swe);
    fprintf(LOG_DEST, "\twdew: %.4f\n", save->wdew);
}

/******************************************************************************
 * @brief     Print snow data structure.
 *****************************************************************************/
void
print_snow_data(snow_data_struct *snow)
{
    // Print state variables
    fprintf(LOG_DEST, "snow_data - states:\n");
    fprintf(LOG_DEST, "\talbedo            : %f\n", snow->albedo);
    fprintf(LOG_DEST, "\tcoverage          : %f\n", snow->coverage);
    fprintf(LOG_DEST, "\tdensity           : %f\n", snow->density);
    fprintf(LOG_DEST, "\tsnow_depth        : %f\n", snow->snow_depth);
    fprintf(LOG_DEST, "\tSnowAge           : %f\n", snow->SnowAge);
    fprintf(LOG_DEST, "\told_swq           : %f\n", snow->old_swq);
    fprintf(LOG_DEST, "\tswq               : %f\n", snow->swq);
    fprintf(LOG_DEST, "\tpack_melt         : %f\n", snow->pack_melt);
    fprintf(LOG_DEST, "\tsnow_canopy       : %f\n", snow->snow_canopy);
    fprintf(LOG_DEST, "\tdelta_depth       : %f\n", snow->delta_depth);
    fprintf(LOG_DEST, "\tnew_snow_density  : %f\n", snow->new_snow_density);

}

/******************************************************************************
 * @brief    Print soil_con_struct.
 *****************************************************************************/
void
print_soil_con(soil_con_struct *scon,
               size_t           nlayers,
               size_t           nnodes,
               size_t           nfrost,
               size_t           nbands,
               size_t           nzwt)
{
    size_t i;
    size_t j;

    fprintf(LOG_DEST, "soil_con:\n");
    fprintf(LOG_DEST, "\tFS_ACTIVE             : %s\n",
            scon->FS_ACTIVE ? "true" : "false");
    fprintf(LOG_DEST, "\tKsat                  :");
    for (i = 0; i < nlayers; i++) {
        fprintf(LOG_DEST, "\t%f", scon->Ksat[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tWfc                   :");
    for (i = 0; i < nlayers; i++) {
        fprintf(LOG_DEST, "\t%f", scon->Wfc[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tWpwp                  :");
    for (i = 0; i < nlayers; i++) {
        fprintf(LOG_DEST, "\t%f", scon->Wpwp[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tb_dynamic             : %f\n", scon->b_dynamic);
    fprintf(LOG_DEST, "\tAlbedoPar             : %f\n", scon->AlbedoPar);
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tb_infilt              : %f\n", scon->b_infilt);
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tbubble                :");
    for (i = 0; i < nlayers; i++) {
        fprintf(LOG_DEST, "\t%f", scon->bubble[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tbubble_node           :");
    for (i = 0; i < nnodes; i++) {
        fprintf(LOG_DEST, "\t%f", scon->bubble_node[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tbulk_density          :");
    for (i = 0; i < nlayers; i++) {
        fprintf(LOG_DEST, "\t%f", scon->bulk_density[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tbulk_dens_min         :");
    for (i = 0; i < nlayers; i++) {
        fprintf(LOG_DEST, "\t%f", scon->bulk_dens_min[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tbulk_dens_org       :");
    for (i = 0; i < nlayers; i++) {
        fprintf(LOG_DEST, "\t%f", scon->bulk_dens_org[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tcap_drive             : %f\n", scon->cap_drive);
    fprintf(LOG_DEST, "\tdepth                 :");
    for (i = 0; i < nlayers; i++) {
        fprintf(LOG_DEST, "\t%f", scon->depth[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tdp                    : %f\n", scon->dp);
    fprintf(LOG_DEST, "\tdz_soil               :");
    for (i = 0; i < nnodes; i++) {
        fprintf(LOG_DEST, "\t%f", scon->dz_soil[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tZsum_soil             :");
    for (i = 0; i < nnodes; i++) {
        fprintf(LOG_DEST, "\t%f", scon->Zsum_soil[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\texpt                  :");
    for (i = 0; i < nlayers; i++) {
        fprintf(LOG_DEST, "\t%f", scon->expt[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\texpt_node             :");
    for (i = 0; i < nnodes; i++) {
        fprintf(LOG_DEST, "\t%f", scon->expt_node[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tmax_infil             : %f\n", scon->max_infil);
    fprintf(LOG_DEST, "\tmax_moist             :");
    for (i = 0; i < nlayers; i++) {
        fprintf(LOG_DEST, "\t%f", scon->max_moist[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tDsat                 :");
    for (i = 0; i < nlayers; i++) {
        fprintf(LOG_DEST, "\t%f", scon->Dsat[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tporosity              :");
    for (i = 0; i < nlayers; i++) {
        fprintf(LOG_DEST, "\t%f", scon->porosity[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tporosity_node        :");
    for (i = 0; i < nnodes; i++) {
        fprintf(LOG_DEST, "\t%f", scon->porosity_node[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tquartz              :");
    for (i = 0; i < nlayers; i++) {
        fprintf(LOG_DEST, "\t%f", scon->quartz[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\torganic               :");
    for (i = 0; i < nlayers; i++) {
        fprintf(LOG_DEST, "\t%f", scon->organic[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tresid_moist           :");
    for (i = 0; i < nlayers; i++) {
        fprintf(LOG_DEST, "\t%f", scon->resid_moist[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\trough                 : %f\n", scon->rough);
    fprintf(LOG_DEST, "\tsnow_rough            : %f\n", scon->snow_rough);
    fprintf(LOG_DEST, "\tsoil_density          :");
    for (i = 0; i < nlayers; i++) {
        fprintf(LOG_DEST, "\t%f", scon->soil_density[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tsoil_dens_min         :");
    for (i = 0; i < nlayers; i++) {
        fprintf(LOG_DEST, "\t%f", scon->soil_dens_min[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tsoil_dens_org         :");
    for (i = 0; i < nlayers; i++) {
        fprintf(LOG_DEST, "\t%f", scon->soil_dens_org[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "BandElev                :");
    for (i = 0; i < nbands; i++) {
        fprintf(LOG_DEST, "\t%f", scon->BandElev[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "AreaFract               :");
    for (i = 0; i < nbands; i++) {
        fprintf(LOG_DEST, "\t%f", scon->AreaFract[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "Pfactor               :");
    for (i = 0; i < nbands; i++) {
        fprintf(LOG_DEST, "\t%f", scon->Pfactor[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "Tfactor               :");
    for (i = 0; i < nbands; i++) {
        fprintf(LOG_DEST, "\t%f", scon->Tfactor[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\televation             : %f\n", scon->elevation);
    fprintf(LOG_DEST, "\tlat                   : %f\n", scon->lat);
    fprintf(LOG_DEST, "\tlng                   : %f\n", scon->lng);
    fprintf(LOG_DEST, "\tcell_area             : %f\n", scon->cell_area);
    fprintf(LOG_DEST, "\ttime_zone_lng         : %f\n", scon->time_zone_lng);
    fprintf(LOG_DEST, "\tgridcel               : %d\n", scon->gridcel);
    fprintf(LOG_DEST, "\tzwtvmoist_zwt         :");
    for (i = 0; i < nlayers + 2; i++) {
        for (j = 0; j < nzwt; j++) {
            fprintf(LOG_DEST, "\t%f", scon->zwtvmoist_zwt[i][j]);
        }
        fprintf(LOG_DEST, "\n\t\t\t");
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tzwtvmoist_moist       :");
    for (i = 0; i < nlayers + 2; i++) {
        for (j = 0; j < nzwt; j++) {
            fprintf(LOG_DEST, "\t%f", scon->zwtvmoist_moist[i][j]);
        }
        fprintf(LOG_DEST, "\n\t\t\t");
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tslope                 : %f\n", scon->slope);
    fprintf(LOG_DEST, "\taspect                : %f\n", scon->aspect);
    fprintf(LOG_DEST, "\tehoriz                : %f\n", scon->ehoriz);
    fprintf(LOG_DEST, "\twhoriz                : %f\n", scon->whoriz);
}

/******************************************************************************
 * @brief    Print veg_con structure.
 *****************************************************************************/
void
print_veg_con(veg_con_struct *vcon,
              size_t          nroots,
              char            carbon,
              size_t          ncanopy)
{
    size_t i;

    fprintf(LOG_DEST, "veg_con:\n");
    fprintf(LOG_DEST, "\tCv              : %.4f\n", vcon->Cv);
    fprintf(LOG_DEST, "\troot            :");
    for (i = 0; i < nroots; i++) {
        fprintf(LOG_DEST, "\t%.2f", vcon->root[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tveg_class       : %d\n", vcon->veg_class);
    fprintf(LOG_DEST, "\tvegetat_type_num: %zu\n", vcon->vegetat_type_num);
    if (carbon) {
        fprintf(LOG_DEST, "\tCanopLayerBnd   :");
        for (i = 0; i < ncanopy; i++) {
            fprintf(LOG_DEST, "\t%.2f", vcon->CanopLayerBnd[i]);
        }
    }
}

/******************************************************************************
 * @brief    Print vegetation library variables.
 *****************************************************************************/
void
print_veg_lib(veg_lib_struct *vlib,
              char            carbon)
{
    size_t i;

    fprintf(LOG_DEST, "veg_lib:\n");
    fprintf(LOG_DEST, "\tLAI           :");
    for (i = 0; i < MONTHS_PER_YEAR; i++) {
        fprintf(LOG_DEST, "\t%.2f", vlib->LAI[i]);
    }
    fprintf(LOG_DEST, "\tSAI           :");
    for (i = 0; i < MONTHS_PER_YEAR; i++) {
        fprintf(LOG_DEST, "\t%.2f", vlib->SAI[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\talbedo        :");
    for (i = 0; i < MONTHS_PER_YEAR; i++) {
        fprintf(LOG_DEST, "\t%.2f", vlib->albedo[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tfcanopy        :");
    for (i = 0; i < MONTHS_PER_YEAR; i++) {
        fprintf(LOG_DEST, "\t%.2f", vlib->fcanopy[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tdisplacement  :");
    for (i = 0; i < MONTHS_PER_YEAR; i++) {
        fprintf(LOG_DEST, "\t%.2f", vlib->displacement[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tNVegLibTypes  : %zu\n", vlib->NVegLibTypes);
    fprintf(LOG_DEST, "\trad_atten     : %.4f\n", vlib->rad_atten);
    fprintf(LOG_DEST, "\trmax          : %.4f\n", vlib->rmax);
    fprintf(LOG_DEST, "\trmin          : %.4f\n", vlib->rmin);
    fprintf(LOG_DEST, "\tCanopy_Upper  : %.4f\n", vlib->Canopy_Upper);
    fprintf(LOG_DEST, "\tCanopy_Lower  : %.4f\n", vlib->Canopy_Lower);
    fprintf(LOG_DEST, "\troughness     :");
    for (i = 0; i < MONTHS_PER_YEAR; i++) {
        fprintf(LOG_DEST, "\t%.2f", vlib->roughness[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\ttrunk_ratio   : %.4f\n", vlib->trunk_ratio);
    fprintf(LOG_DEST, "\twind_atten    : %.4f\n", vlib->wind_atten);
    fprintf(LOG_DEST, "\tRGL           : %.4f\n", vlib->RGL);
    fprintf(LOG_DEST, "\tveg_class     : %d\n", vlib->veg_class);
    if (carbon) {
        fprintf(LOG_DEST, "\tCtype         : %d\n", vlib->Ctype);
        fprintf(LOG_DEST, "\tMaxCarboxRate : %.4f\n", vlib->MaxCarboxRate);
        fprintf(LOG_DEST, "\tMaxETransport : %.4f\n", vlib->MaxETransport);
        fprintf(LOG_DEST, "\tCO2Specificity: %.4f\n", vlib->CO2Specificity);
        fprintf(LOG_DEST, "\tLightUseEff   : %.4f\n", vlib->LightUseEff);
        fprintf(LOG_DEST, "\tNscaleFlag    : %s\n",
                vlib->NscaleFlag ? "true" : "false");
        fprintf(LOG_DEST, "\tWnpp_inhib    : %.4f\n", vlib->Wnpp_inhib);
        fprintf(LOG_DEST, "\tNPPfactor_sat : %.4f\n", vlib->NPPfactor_sat);
    }
}

/******************************************************************************
 * @brief    Print vegetation variables.
 *****************************************************************************/
void
print_veg_var(veg_var_struct *vvar,
              size_t          ncanopy,
              size_t          nswband)
{
    extern option_struct options;

    size_t               i;

    // Print state variables
    fprintf(LOG_DEST, "veg_var - states:\n");
    fprintf(LOG_DEST, "\talbedo   : %f\n", vvar->albedo);
    fprintf(LOG_DEST, "\tdisplacement : %f\n", vvar->displacement);
    fprintf(LOG_DEST, "\tfcanopy   : %f\n", vvar->fcanopy);
    fprintf(LOG_DEST, "\tLAI   : %f\n", vvar->LAI);
    fprintf(LOG_DEST, "\tSAI   : %f\n", vvar->SAI);
    fprintf(LOG_DEST, "\tNetLAI   : %f\n", vvar->NetLAI);
    fprintf(LOG_DEST, "\tNetSAI   : %f\n", vvar->NetSAI);
    fprintf(LOG_DEST, "\troughness    : %f\n", vvar->roughness);
    fprintf(LOG_DEST, "\tWdew         : %f\n", vvar->Wdew);
    fprintf(LOG_DEST, "\tMaxSnowInt   : %f\n", vvar->MaxSnowInt);
    fprintf(LOG_DEST, "\tMaxRainInt   : %f\n", vvar->MaxRainInt);
    fprintf(LOG_DEST, "\twetFrac      : %f\n", vvar->wetFrac);

    // Print fluxes
    fprintf(LOG_DEST, "veg_var - fluxes:\n");
    fprintf(LOG_DEST, "\tRainThroughFall  : %f\n", vvar->RainThroughFall);
    fprintf(LOG_DEST, "\tSnowThroughFall  : %f\n", vvar->SnowThroughFall);
    fprintf(LOG_DEST, "\tRainDrip  : %f\n", vvar->RainDrip);
    fprintf(LOG_DEST, "\tSnowDrip  : %f\n", vvar->SnowDrip);
    fprintf(LOG_DEST, "\tint_rain  : %f\n", vvar->int_rain);
    fprintf(LOG_DEST, "\tint_snow  : %f\n", vvar->int_snow);

    if (options.CARBON) {
        // Carbon terms - states
        fprintf(LOG_DEST, "\tAnnualNPP    : %f\n", vvar->AnnualNPP);
        fprintf(LOG_DEST, "\tAnnualNPPPrev: %f\n", vvar->AnnualNPPPrev);
        fprintf(LOG_DEST, "\tCi           : %f\n", vvar->Ci);
        fprintf(LOG_DEST, "\tCiLayer      :");
        for (i = 0; i < ncanopy; i++) {
            fprintf(LOG_DEST, "\t%f", vvar->CiLayer[i]);
        }
        fprintf(LOG_DEST, "\n");
        fprintf(LOG_DEST, "\tNPPfactor    : %f\n", vvar->NPPfactor);
        fprintf(LOG_DEST, "\tNscaleFactor :");
        for (i = 0; i < ncanopy; i++) {
            fprintf(LOG_DEST, "\t%f", vvar->NscaleFactor[i]);
        }
        fprintf(LOG_DEST, "\n");
        fprintf(LOG_DEST, "\trc           : %f\n", vvar->rc);
        fprintf(LOG_DEST, "\trsLayer      :");
        for (i = 0; i < ncanopy; i++) {
            fprintf(LOG_DEST, "\t%f", vvar->rsLayer[i]);
        }
        fprintf(LOG_DEST, "\n");
        // Carbon terms - fluxes
        fprintf(LOG_DEST, "\taPAR         : %f\n", vvar->aPAR);
        fprintf(LOG_DEST, "\taPARLayer    :");
        for (i = 0; i < ncanopy; i++) {
            fprintf(LOG_DEST, "\t%f", vvar->aPARLayer[i]);
        }
        fprintf(LOG_DEST, "\n");
        fprintf(LOG_DEST, "\tGPP          : %f\n", vvar->GPP);
        fprintf(LOG_DEST, "\tLitterfall   : %f\n", vvar->Litterfall);
        for (i = 0; i < ncanopy; i++) {
            fprintf(LOG_DEST, "\t%f", vvar->aPARLayer[i]);
        }
        fprintf(LOG_DEST, "\n");
        fprintf(LOG_DEST, "\tNPP          : %f\n", vvar->NPP);
        fprintf(LOG_DEST, "\tRaut         : %f\n", vvar->Raut);
        fprintf(LOG_DEST, "\tRdark        : %f\n", vvar->Rdark);
        fprintf(LOG_DEST, "\tRgrowth      : %f\n", vvar->Rgrowth);
        fprintf(LOG_DEST, "\tRmaint       : %f\n", vvar->Rmaint);
        fprintf(LOG_DEST, "\tRphoto       : %f\n", vvar->Rphoto);
    }
}
