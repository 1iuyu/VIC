/******************************************************************************
 * @section DESCRIPTION
 *
 * Print library.
 *****************************************************************************/

#include "vic_driver_shared_all.h"

/******************************************************************************
 * @brief    Print cell data structure.
 *****************************************************************************/
void
print_cell_data(cell_data_struct *cell)
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
    fprintf(LOG_DEST, "\trootmoist   : %f\n", cell->rootmoist);
    fprintf(LOG_DEST, "\tzwt         : %f\n", cell->zwt);

    // Print fluxes
    fprintf(LOG_DEST, "cell_data - fluxes:\n");
    fprintf(LOG_DEST, "\tbaseflow    : %f\n", cell->baseflow);
    fprintf(LOG_DEST, "\tevap        : %f\n", cell->evap);
    fprintf(LOG_DEST, "\ttransp      : %f\n", cell->transp);
    fprintf(LOG_DEST, "\tdewsoil     : %f\n", cell->dewsoil);
    fprintf(LOG_DEST, "\tcanopyevap  : %f\n", cell->canopyevap);
    fprintf(LOG_DEST, "\tsnowfrost   : %f\n", cell->snowfrost);
    fprintf(LOG_DEST, "\tsnow_sublim : %f\n", cell->snow_sublim);
    fprintf(LOG_DEST, "\tsoil_inflow : %f\n", cell->soil_inflow);
    fprintf(LOG_DEST, "\trunoff      : %f\n", cell->runoff);
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
                 size_t             nnodes)
{
    size_t i;

    // Print energy_bal - state variables
    fprintf(LOG_DEST, "energy_bal - states:\n");
    fprintf(LOG_DEST, "\tTcanopy    : %f\n", eb->Tcanopy);
    fprintf(LOG_DEST, "\tTfoliage   : %f\n", eb->Tfoliage);
    fprintf(LOG_DEST, "\tTsurf      : %f\n", eb->Tsurf);
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
    fprintf(LOG_DEST, "\tadvection           : %f\n", eb->advection);
    fprintf(LOG_DEST, "\tAdvectSub           : %f\n", eb->AdvectSub);
    fprintf(LOG_DEST, "\tAdvectGrnd          : %f\n", eb->AdvectGrnd);
    fprintf(LOG_DEST, "\tAdvectOver          : %f\n", eb->AdvectOver);
    fprintf(LOG_DEST, "\tgrnd_flux           : %f\n", eb->grnd_flux);
    fprintf(LOG_DEST, "\tlatent              : %f\n", eb->latent);
    fprintf(LOG_DEST, "\tsensible            : %f\n", eb->sensible);
    fprintf(LOG_DEST, "\tSensibleStem        : %f\n", eb->SensibleStem);
    fprintf(LOG_DEST, "\tSensibleLeaf        : %f\n", eb->SensibleLeaf);
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
    fprintf(LOG_DEST, "\tresolution           : %.4f\n", gp->resolution);
    fprintf(LOG_DEST, "\tstep_dt              : %.4f\n", gp->step_dt);
    fprintf(LOG_DEST, "\tmodel_steps_per_day  : %zu\n", gp->model_steps_per_day);
    fprintf(LOG_DEST, "\tendday               : %hu\n", gp->endday);
    fprintf(LOG_DEST, "\tendmonth             : %hu\n", gp->endmonth);
    fprintf(LOG_DEST, "\tendyear              : %hu\n", gp->endyear);
    for (i = 0; i < 2; i++) {
        fprintf(LOG_DEST, "\tforceday[%zd]          : %hu\n", i, gp->forceday[i]);
        fprintf(LOG_DEST, "\tforcesec[%zd]          : %u\n", i, gp->forcesec[i]);
        fprintf(LOG_DEST, "\tforcemonth[%zd]        : %hu\n", i,
                gp->forcemonth[i]);
        fprintf(LOG_DEST, "\tforceoffset[%zd]       : %hu\n", i,
                gp->forceoffset[i]);
        fprintf(LOG_DEST, "\tforceskip[%zd]         : %u\n", i, gp->forceskip[i]);
        fprintf(LOG_DEST, "\tforceyear[%zd]         : %hu\n", i,
                gp->forceyear[i]);
    }
    fprintf(LOG_DEST, "\tnrecs                : %zu\n", gp->nrecs);
    fprintf(LOG_DEST, "\tstartday             : %hu\n", gp->startday);
    fprintf(LOG_DEST, "\tstartsec             : %u\n", gp->startsec);
    fprintf(LOG_DEST, "\tstartmonth           : %hu\n", gp->startmonth);
    fprintf(LOG_DEST, "\tstartyear            : %hu\n", gp->startyear);
    fprintf(LOG_DEST, "\tstateday             : %hu\n", gp->stateday);
    fprintf(LOG_DEST, "\tstatemonth           : %hu\n", gp->statemonth);
    fprintf(LOG_DEST, "\tstateyear            : %hu\n", gp->stateyear);
    fprintf(LOG_DEST, "\tstatesec             : %u\n", gp->statesec);
}

/******************************************************************************
 * @brief    Print options structure.
 *****************************************************************************/
void
print_option(option_struct *option)
{
    fprintf(LOG_DEST, "option:\n");
    fprintf(LOG_DEST, "\tCARBON               : %s\n",
            option->CARBON ? "true" : "false");
    fprintf(LOG_DEST, "\tCONTINUEONERROR      : %s\n",
            option->CONTINUEONERROR ? "true" : "false");
    fprintf(LOG_DEST, "\tCORRPREC             : %s\n",
            option->CORRPREC ? "true" : "false");
    fprintf(LOG_DEST, "\tNlayer               : %zu\n", option->Nlayer);
    fprintf(LOG_DEST, "\tNOFLUX               : %s\n",
            option->NOFLUX ? "true" : "false");
    fprintf(LOG_DEST, "\tNVEGTYPES            : %zu\n", option->NVEGTYPES);
    fprintf(LOG_DEST, "\tSNOW_DENSITY         : %d\n", option->SNOW_DENSITY);
    fprintf(LOG_DEST, "\tSNOW_BAND            : %zu\n", option->SNOW_BAND);
    fprintf(LOG_DEST, "\tTFALLBACK            : %s\n",
            option->TFALLBACK ? "true" : "false");
    fprintf(LOG_DEST, "\tGRID_DECIMAL         : %d\n", option->GRID_DECIMAL);
    fprintf(LOG_DEST, "\tLAI_SRC              : %d\n", option->LAI_SRC);
    fprintf(LOG_DEST, "\tSAI_SRC              : %d\n", option->SAI_SRC);
    fprintf(LOG_DEST, "\tFCAN_SRC             : %d\n", option->FCAN_SRC);
    fprintf(LOG_DEST, "\tPARAM_FROM_SOIL    : %s\n",
            option->PARAM_FROM_SOIL ? "true" : "false");
    fprintf(LOG_DEST, "\tSTATE_FORMAT         : %d\n", option->STATE_FORMAT);
    fprintf(LOG_DEST, "\tINIT_STATE           : %s\n",
            option->INIT_STATE ? "true" : "false");
    fprintf(LOG_DEST, "\tSAVE_STATE           : %s\n",
            option->SAVE_STATE ? "true" : "false");
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
    fprintf(LOG_DEST, "\tVEG_LAI_SNOW_MULTIPLIER: %.4f\n",
            param->VEG_LAI_SNOW_MULTIPLIER);
    fprintf(LOG_DEST, "\tVEG_LAI_WATER_FACTOR: %.4f\n",
            param->VEG_LAI_WATER_FACTOR);
    fprintf(LOG_DEST, "\tCANOPY_CLOSURE: %.4f\n", param->CANOPY_CLOSURE);
    fprintf(LOG_DEST, "\tCANOPY_RSMAX: %.4f\n", param->CANOPY_RSMAX);
    fprintf(LOG_DEST, "\tPHOTO_LRESC3: %.4f\n", param->PHOTO_LRESC3);
    fprintf(LOG_DEST, "\tPHOTO_LRESC4: %.4f\n", param->PHOTO_LRESC4);
    fprintf(LOG_DEST, "\tPHOTO_OX: %.4f\n", param->PHOTO_OX);
    fprintf(LOG_DEST, "\tPHOTO_KC: %.4f\n", param->PHOTO_KC);
    fprintf(LOG_DEST, "\tPHOTO_KO: %.4f\n", param->PHOTO_KO);
    fprintf(LOG_DEST, "\tPHOTO_EC: %.4f\n", param->PHOTO_EC);
    fprintf(LOG_DEST, "\tPHOTO_EO: %.4f\n", param->PHOTO_EO);
    fprintf(LOG_DEST, "\tPHOTO_EV: %.4f\n", param->PHOTO_EV);
    fprintf(LOG_DEST, "\tPHOTO_ER: %.4f\n", param->PHOTO_ER);
    fprintf(LOG_DEST, "\tSNOW_MAX_SURFACE_SWE: %.4f\n",
            param->SNOW_MAX_SURFACE_SWE);
    fprintf(LOG_DEST, "\tSNOW_LIQUID_WATER_CAPACITY: %.4f\n",
            param->SNOW_LIQUID_WATER_CAPACITY);
    fprintf(LOG_DEST, "\tSNOW_NEW_SNOW_DENSITY: %.4f\n",
            param->SNOW_NEW_SNOW_DENSITY);
    fprintf(LOG_DEST, "\tSNOW_NEW_SNOW_DENS_MAX: %.4f\n",
            param->SNOW_NEW_SNOW_DENS_MAX);
    fprintf(LOG_DEST, "\tSNOW_CONDUCT: %.4f\n", param->SNOW_CONDUCT);
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
    fprintf(LOG_DEST, "\tsnow_depth        : %f\n", snow->snow_depth);
    fprintf(LOG_DEST, "\tsnowage           : %f\n", snow->snowage);
    fprintf(LOG_DEST, "\told_swq           : %f\n", snow->last_swq);
    fprintf(LOG_DEST, "\tswq               : %f\n", snow->swq);
    fprintf(LOG_DEST, "\tpack_melt         : %f\n", snow->pack_melt);
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
               size_t           nbands)
{
    size_t i;

    fprintf(LOG_DEST, "soil_con:\n");
    fprintf(LOG_DEST, "\tcapil_drive           : %f\n", scon->capil_drive);
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tAlbedoPar             : %f\n", scon->AlbedoPar);
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tb_dynamic             : %f\n", scon->b_dynamic);
    fprintf(LOG_DEST, "\tb_infilt              : %f\n", scon->b_infilt);
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tbubble_node           :");
    for (i = 0; i < nnodes; i++) {
        fprintf(LOG_DEST, "\t%f", scon->bubble_node[i]);
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
    fprintf(LOG_DEST, "\tdepth                 :");
    for (i = 0; i < nlayers; i++) {
        fprintf(LOG_DEST, "\t%f", scon->depth[i]);
    }
    fprintf(LOG_DEST, "\n");
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
    fprintf(LOG_DEST, "\torganic               :");
    for (i = 0; i < nlayers; i++) {
        fprintf(LOG_DEST, "\t%f", scon->organic_node[i]);
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
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tslope                 : %f\n", scon->slope);
}

/******************************************************************************
 * @brief    Print veg_con structure.
 *****************************************************************************/
void
print_veg_con(veg_con_struct *vcon,
              size_t          nroots)
{
    size_t i;

    fprintf(LOG_DEST, "veg_con:\n");
    fprintf(LOG_DEST, "\tCv              : %.4f\n", vcon->Cv);
    fprintf(LOG_DEST, "\tBandIndex       : %d\n", vcon->BandIndex);
    fprintf(LOG_DEST, "\tNroot           : %zu\n", vcon->Nroot);

    fprintf(LOG_DEST, "\troot            :");
    for (i = 0; i < nroots; i++) {
        fprintf(LOG_DEST, "\t%.2f", vcon->root[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tveg_class       : %d\n", vcon->veg_class);
    fprintf(LOG_DEST, "\tvegetat_type_num: %zu\n", vcon->vegetat_type_num);
}

/******************************************************************************
 * @brief    Print vegetation library variables.
 *****************************************************************************/
void
print_veg_lib(veg_lib_struct *vlib)
{
    size_t i;

    fprintf(LOG_DEST, "veg_lib:\n");
    fprintf(LOG_DEST, "\tLAI           :");
    for (i = 0; i < MONTHS_PER_YEAR; i++) {
        fprintf(LOG_DEST, "\t%.2f", vlib->LAI[i]);
    }
    fprintf(LOG_DEST, "veg_lib:\n");
    fprintf(LOG_DEST, "\tSAI           :");
    for (i = 0; i < MONTHS_PER_YEAR; i++) {
        fprintf(LOG_DEST, "\t%.2f", vlib->SAI[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tfcanopy        :");
    for (i = 0; i < MONTHS_PER_YEAR; i++) {
        fprintf(LOG_DEST, "\t%.2f", vlib->fcanopy[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\treflleaf       :");
    for (i = 0; i < MAX_SWBANDS; i++) {
        fprintf(LOG_DEST, "\t%.2f", vlib->reflleaf[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\treflstem       :");
    for (i = 0; i < MAX_SWBANDS; i++) {
        fprintf(LOG_DEST, "\t%.2f", vlib->reflstem[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\ttransleaf      :");
    for (i = 0; i < MAX_SWBANDS; i++) {
        fprintf(LOG_DEST, "\t%.2f", vlib->transleaf[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\ttransstem      :");
    for (i = 0; i < MAX_SWBANDS; i++) {
        fprintf(LOG_DEST, "\t%.2f", vlib->transstem[i]);
    }

    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tNVegLibTypes  : %zu\n", vlib->NVegLibTypes);
    fprintf(LOG_DEST, "\tCanopy_Radius : %.4f\n", vlib->Canopy_Radius);
    fprintf(LOG_DEST, "\tCanopy_Upper  : %.4f\n", vlib->Canopy_Upper);
    fprintf(LOG_DEST, "\tCanopy_Lower  : %.4f\n", vlib->Canopy_Lower);
    fprintf(LOG_DEST, "\ttrunk_dia     : %.4f\n", vlib->trunk_dia);
    fprintf(LOG_DEST, "\tc_biomass     : %.4f\n", vlib->c_biomass);
    fprintf(LOG_DEST, "\tCOI           : %.4f\n", vlib->COI);
    fprintf(LOG_DEST, "\td_leaf        : %.4f\n", vlib->d_leaf);
    fprintf(LOG_DEST, "\troot_a        : %.4f\n", vlib->root_a);
    fprintf(LOG_DEST, "\troot_b        : %.4f\n", vlib->root_b);
    fprintf(LOG_DEST, "\troot_d        : %.4f\n", vlib->root_d);
    fprintf(LOG_DEST, "\tliq_bioms     : %.4f\n", vlib->liq_bioms);
    fprintf(LOG_DEST, "\tslatop        : %.4f\n", vlib->slatop);
    fprintf(LOG_DEST, "\tstem_num      : %.4f\n", vlib->stem_num);
    fprintf(LOG_DEST, "\tZ0sub_LAImax  : %.4f\n", vlib->Z0sub_LAImax);
    fprintf(LOG_DEST, "\tZ0sub_Cs      : %.4f\n", vlib->Z0sub_Cs);
    fprintf(LOG_DEST, "\tZ0sub_Cr      : %.4f\n", vlib->Z0sub_Cr);
    fprintf(LOG_DEST, "\tZ0sub_c       : %.4f\n", vlib->Z0sub_c);
    fprintf(LOG_DEST, "\tZ0sub_cw      : %.4f\n", vlib->Z0sub_cw);
    fprintf(LOG_DEST, "\tsmpsc         : %.4f\n", vlib->smpsc);
    fprintf(LOG_DEST, "\tsmpso         : %.4f\n", vlib->smpso);

    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tveg_class     : %d\n", vlib->veg_class);
    fprintf(LOG_DEST, "\tCtype         : %d\n", vlib->Ctype);
    fprintf(LOG_DEST, "\tmatric50      : %.4f\n", vlib->matric50);
    fprintf(LOG_DEST, "\tkcano_max     : %.4f\n", vlib->kcano_max);
    fprintf(LOG_DEST, "\tkroot_max     : %.4f\n", vlib->kroot_max);
    fprintf(LOG_DEST, "\ttheta_cj      : %.4f\n", vlib->theta_cj);
    fprintf(LOG_DEST, "\tleaf_CN       : %.4f\n", vlib->leaf_CN);
    fprintf(LOG_DEST, "\tSLA_top       : %.4f\n", vlib->SLA_top);
    fprintf(LOG_DEST, "\tfN_rub        : %.4f\n", vlib->fN_rub);
    fprintf(LOG_DEST, "\tmedlynslope   : %.4f\n", vlib->medlynslope);
    fprintf(LOG_DEST, "\tmedlynint     : %.4f\n", vlib->medlynint);
}

/******************************************************************************
 * @brief    Print vegetation variables.
 *****************************************************************************/
void
print_veg_var(veg_var_struct *vvar,
              size_t          ncanopy)
{
    // Print state variables
    fprintf(LOG_DEST, "veg_var - states:\n");
    fprintf(LOG_DEST, "\tfcanopy   : %f\n", vvar->fcanopy);
    fprintf(LOG_DEST, "\tLAI   : %f\n", vvar->LAI);
    fprintf(LOG_DEST, "\tSAI   : %f\n", vvar->SAI);
    fprintf(LOG_DEST, "\tWdew         : %f\n", vvar->Wdew);
    fprintf(LOG_DEST, "\tMaxSnowInt   : %f\n", vvar->MaxSnowInt);
    fprintf(LOG_DEST, "\tMaxRainInt   : %f\n", vvar->MaxRainInt);
    fprintf(LOG_DEST, "\twetFrac      : %f\n", vvar->wetFrac);
    fprintf(LOG_DEST, "\twetFrac      : %f\n", vvar->dryFrac);
    // Print fluxes
    fprintf(LOG_DEST, "veg_var - fluxes:\n");
    fprintf(LOG_DEST, "\tRainThroughFall  : %f\n", vvar->RainThroughFall);
    fprintf(LOG_DEST, "\tSnowThroughFall  : %f\n", vvar->SnowThroughFall);
    fprintf(LOG_DEST, "\tRainDrip  : %f\n", vvar->RainDrip);
    fprintf(LOG_DEST, "\tSnowDrip  : %f\n", vvar->SnowDrip);
    fprintf(LOG_DEST, "\tint_rain  : %f\n", vvar->int_rain);
    fprintf(LOG_DEST, "\tint_snow  : %f\n", vvar->int_snow);
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tLAI_z      :");
    for (size_t i = 0; i < ncanopy; i++) {
        fprintf(LOG_DEST, "\t%.2f", vvar->LAI_z[i]);
    }
    fprintf(LOG_DEST, "\n");
    fprintf(LOG_DEST, "\tSAI_z      :");
    for (size_t i = 0; i < ncanopy; i++) {
        fprintf(LOG_DEST, "\t%.2f", vvar->SAI_z[i]);
    }    
}
