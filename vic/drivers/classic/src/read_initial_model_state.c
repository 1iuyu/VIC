/******************************************************************************
 * @section DESCRIPTION
 *
 * This subroutine initializes the model state at hour 0 of the date defined in
 * the given state file.
 *
 * Soil moisture, soil thermal, and snowpack variables  are stored for each
 * vegetation type and snow band.  However moisture variables from the
 * distributed precipitation model are averaged so that the model is restarted
 * with mu = 1.
 *****************************************************************************/

#include <vic_driver_classic.h>

/******************************************************************************
 * @brief    This subroutine initializes the model state at hour 0 of the date
 *           defined in the given state file.
 *****************************************************************************/
void
read_initial_model_state(FILE            *init_state,
                         all_vars_struct *all_vars,
                         int              Nveg,
                         int              Nbands,
                         int              cellnum,
                         soil_con_struct *soil_con)
{
    extern option_struct options;

    char                 tmpstr[MAXSTRING];
    int                  veg, iveg;
    int                  band, iband;
    size_t               lidx;
    size_t               nidx;
    int                  tmp_cellnum;
    int                  tmp_Nveg;
    int                  tmp_Nband;
    int                  tmp_char;
    int                  byte, Nbytes;

    cell_data_struct    *cell;
    snow_data_struct    *snow;
    energy_bal_struct   *energy;
    veg_var_struct      *veg_var;

    cell = all_vars->cell;
    veg_var = all_vars->veg_var;
    snow = all_vars->snow;
    energy = all_vars->energy;

    /* read cell information */
    if (options.STATE_FORMAT == BINARY) {
        fread(&tmp_cellnum, sizeof(int), 1, init_state);
        fread(&tmp_Nveg, sizeof(int), 1, init_state);
        fread(&tmp_Nband, sizeof(int), 1, init_state);
        fread(&Nbytes, sizeof(int), 1, init_state);
    }
    else {
        fscanf(init_state, "%d %d %d", &tmp_cellnum, &tmp_Nveg, &tmp_Nband);
    }
    // Skip over unused cell information
    while (tmp_cellnum != cellnum && !feof(init_state)) {
        if (options.STATE_FORMAT == BINARY) {
            // skip rest of current cells info
            for (byte = 0; byte < Nbytes; byte++) {
                fread(&tmp_char, 1, 1, init_state);
            }
            // read info for next cell
            fread(&tmp_cellnum, sizeof(int), 1, init_state);
            fread(&tmp_Nveg, sizeof(int), 1, init_state);
            fread(&tmp_Nband, sizeof(int), 1, init_state);
            fread(&Nbytes, sizeof(int), 1, init_state);
        }
        else {
            // skip rest of current cells info
            fgets(tmpstr, MAXSTRING, init_state); // skip rest of general cell info
            for (veg = 0; veg <= tmp_Nveg; veg++) {
                for (band = 0; band < tmp_Nband; band++) {
                    fgets(tmpstr, MAXSTRING, init_state); // skip snowband info
                }
            }
            // read info for next cell
            fscanf(init_state, "%d %d %d", &tmp_cellnum, &tmp_Nveg, &tmp_Nband);
        } // end if
    } // end while

    if (feof(init_state)) {
        log_err("Requested grid cell (%d) is not in the model state file.",
                cellnum);
    }

    if (tmp_Nveg != Nveg) {
        log_err("The number of vegetation types in cell %d (%d) does not equal"
                "that defined in vegetation parameter file (%d).  Check your "
                "input files.", cellnum, tmp_Nveg, Nveg);
    }
    if (tmp_Nband != Nbands) {
        log_err("The number of snow bands in cell %d (%d) does not equal that"
                "defined in the snow band file (%d).  Check your input files.",
                cellnum, tmp_Nband, Nbands);
    }

    /* Read soil thermal node deltas */
    for (nidx = 0; nidx < options.Nnode; nidx++) {
        if (options.STATE_FORMAT == BINARY) {
            fread(&soil_con->dz_soil[nidx], sizeof(double), 1, init_state);
        }
        else {
            fscanf(init_state, "%lf", &soil_con->dz_soil[nidx]);
        }
    }
    if (options.Nnode == 1) {
        soil_con->dz_soil[0] = 0;
    }

    /* Read soil thermal node depths */
    for (nidx = 0; nidx < options.Nnode; nidx++) {
        if (options.STATE_FORMAT == BINARY) {
            fread(&soil_con->Zsum_soil[nidx], sizeof(double), 1, init_state);
        }
        else {
            fscanf(init_state, "%lf", &soil_con->Zsum_soil[nidx]);
        }
    }
    if (options.Nnode == 1) {
        soil_con->Zsum_soil[0] = 0;
    }
    if (soil_con->Zsum_soil[options.Nnode - 1] - soil_con->dp > DBL_EPSILON) {
        log_warn("Sum of soil nodes (%f) exceeds defined damping depth (%f)."
                 "Resetting damping depth.",
                 soil_con->Zsum_soil[options.Nnode - 1], soil_con->dp);
        soil_con->dp = soil_con->Zsum_soil[options.Nnode - 1];
    }

    /* Input for all vegetation types */
    for (veg = 0; veg <= Nveg; veg++) {
        /* Read cell identification information */
        if (options.STATE_FORMAT == BINARY) {
            if (fread(&iveg, sizeof(int), 1, init_state) != 1) {
                log_err("End of model state file found unexpectedly");
            }
            if (fread(&iband, sizeof(int), 1, init_state) != 1) {
                log_err("End of model state file found unexpectedly");
            }
        }
        else {
            if (fscanf(init_state, "%d %d", &iveg, &iband) == EOF) {
                log_err("End of model state file found unexpectedly");
            }
        }
        if (iveg != veg || iband != band) {
            log_err("The vegetation and snow band indices in the model "
                    "state file (veg = %d, band = %d) do not match those "
                    "currently requested (veg = %d , band = %d).  Model "
                    "state file must be stored with variables for all "
                    "vegetation indexed by variables for all snow bands.",
                    iveg, iband, veg, band);
        }

        /* Read total soil moisture */
        for (nidx = 0; nidx < options.Nnode; nidx++) {
            if (options.STATE_FORMAT == BINARY) {
                if (fread(&cell[veg].moist[nidx],
                          sizeof(double), 1, init_state) != 1) {
                    log_err("End of model state file found unexpectedly");
                }
            }
            else {
                if (fscanf(init_state, " %lf",
                           &cell[veg].moist[nidx]) == EOF) {
                    log_err("End of model state file found unexpectedly");
                }
            }
        }

        /* Read average ice content */
        for (nidx = 0; nidx < options.Nnode; nidx++) {
            if (options.STATE_FORMAT == BINARY) {
                if (fread(&cell[veg].ice[nidx],
                          sizeof(double), 1, init_state) != 1) {
                    log_err("End of model state file found"
                            "unexpectedly");
                }
            }
            else {
                if (fscanf(init_state, " %lf",
                    &cell[veg].ice[nidx]) == EOF) {
                    log_err("End of model state file found"
                            "unexpectedly");
                }
            }

        }

        if (veg < Nveg) {
            /* Read dew storage */
            if (options.STATE_FORMAT == BINARY) {
                if (fread(&veg_var[veg].Wdew, sizeof(double), 1,
                          init_state) != 1) {
                    log_err("End of model state file found unexpectedly");
                }
            }
            else {
                if (fscanf(init_state, " %lf",
                           &veg_var[veg].Wdew) == EOF) {
                    log_err("End of model state file found unexpectedly");
                }
            }

            if (options.CARBON) {
                if (options.STATE_FORMAT == BINARY) {
                    /* Read cumulative annual NPP */
                    if (fread(&(veg_var[veg].AnnualNPP),
                              sizeof(double), 1, init_state) != 1) {
                        log_err("End of model state file found unexpectedly");
                    }
                    if (fread(&(veg_var[veg].AnnualNPPPrev),
                              sizeof(double), 1, init_state) != 1) {
                        log_err("End of model state file found unexpectedly");
                    }
                    /* Read Soil Carbon Storage */
                    if (fread(&(cell[veg].CLitter), sizeof(double), 1,
                              init_state) != 1) {
                        log_err("End of model state file found unexpectedly");
                    }
                    if (fread(&(cell[veg].CInter), sizeof(double), 1,
                              init_state) != 1) {
                        log_err("End of model state file found unexpectedly");
                    }
                    if (fread(&(cell[veg].CSlow), sizeof(double), 1,
                              init_state) != 1) {
                        log_err("End of model state file found unexpectedly");
                    }
                }
                else {
                    /* Read cumulative annual NPP */
                    if (fscanf(init_state, " %lf",
                               &veg_var[veg].AnnualNPP) == EOF) {
                        log_err("End of model state file found unexpectedly");
                    }
                    if (fscanf(init_state, " %lf",
                               &veg_var[veg].AnnualNPPPrev) == EOF) {
                        log_err("End of model state file found unexpectedly");
                    }
                    /* Read Soil Carbon Storage */
                    if (fscanf(init_state, " %lf",
                               &cell[veg].CLitter) == EOF) {
                        log_err("End of model state file found unexpectedly");
                    }
                    if (fscanf(init_state, " %lf",
                               &cell[veg].CInter) == EOF) {
                        log_err("End of model state file found unexpectedly");
                    }
                    if (fscanf(init_state, " %lf",
                               &cell[veg].CSlow) == EOF) {
                        log_err("End of model state file found unexpectedly");
                    }
                }
            }
        }

        /* Read SnowFrac */
        for (nidx = 0; nidx < options.Nnode; nidx++) {
            if (options.STATE_FORMAT == BINARY) {
                if (fread(&cell[veg].SnowFrac[nidx],
                          sizeof(double), 1, init_state) != 1) {
                    log_err("End of model state file found"
                            "unexpectedly");
                }
            }
            else {
                if (fscanf(init_state, " %lf",
                    &cell[veg].SnowFrac[nidx]) == EOF) {
                    log_err("End of model state file found"
                            "unexpectedly");
                }
            }

        }

        /* Read snow data */
        if (options.STATE_FORMAT == BINARY) {
            if (fread(&energy[veg].Tcanopy, sizeof(int), 1,
                      init_state) != 1) {
                log_err("End of model state file found unexpectedly");
            }
            if (fread(&energy[veg].Tair, sizeof(char), 1,
                      init_state) != 1) {
                log_err("End of model state file found unexpectedly");
            }
            if (fread(&energy[veg].Tgrnd, sizeof(double), 1,
                      init_state) != 1) {
                log_err("End of model state file found unexpectedly");
            }
            if (fread(&snow[veg].coverage, sizeof(double), 1,
                      init_state) != 1) {
                log_err("End of model state file found unexpectedly");
            }
            if (fread(&snow[veg].swq, sizeof(double), 1,
                      init_state) != 1) {
                log_err("End of model state file found unexpectedly");
            }
            if (fread(&snow[veg].density, sizeof(double), 1,
                      init_state) != 1) {
                log_err("End of model state file found unexpectedly");
            }
            if (fread(&snow[veg].snow_canopy, sizeof(double), 1,
                      init_state) != 1) {
                log_err("End of model state file found unexpectedly");
            }
            if (fread(&snow[veg].SnowAge, sizeof(double), 1,
                      init_state) != 1) {
                log_err("End of model state file found unexpectedly");
            }
            if (fread(&snow[veg].old_swq, sizeof(double), 1,
                      init_state) != 1) {
                log_err("End of model state file found unexpectedly");
            }
        }
        else {
            if (fscanf(init_state,
                       " %lf %lf %lf %lf %lf %lf %lf %lf %lf",
                       &energy[veg].Tcanopy, &energy[veg].Tair,
                       &energy[veg].Tgrnd,
                       &snow[veg].coverage, &snow[veg].swq,
                       &snow[veg].density,
                       &snow[veg].snow_canopy,
                       &snow[veg].SnowAge,
                       &snow[veg].old_swq) ==
                EOF) {
                log_err("End of model state file found unexpectedly");
            }
        }

        /* Read soil thermal node temperatures */
        for (nidx = 0; nidx < options.Nnode + snow[veg].Nsnow; nidx++) {
            if (options.STATE_FORMAT == BINARY) {
                if (fread(&energy[veg].T[nidx], sizeof(double), 1,
                          init_state) != 1) {
                    log_err("End of model state file found unexpectedly");
                }
            }
            else {
                if (fscanf(init_state, " %lf",
                           &energy[veg].T[nidx]) == EOF) {
                    log_err("End of model state file found unexpectedly");
                }
            }
        }

        /* Read snow layer pack_ice */
        for (lidx = 0; lidx < snow[veg].Nsnow; lidx++) {
            if (options.STATE_FORMAT == BINARY) {
                if (fread(&snow[veg].pack_ice[lidx], sizeof(double), 1,
                          init_state) != 1) {
                    log_err("End of model state file found unexpectedly");
                }
            }
            else {
                if (fscanf(init_state, " %lf",
                           &snow[veg].pack_ice[lidx]) == EOF) {
                    log_err("End of model state file found unexpectedly");
                }
            }
        }

        /* Read Read snow layer pack_liq */
        for (lidx = 0; lidx < snow[veg].Nsnow; lidx++) {
            if (options.STATE_FORMAT == BINARY) {
                if (fread(&snow[veg].pack_liq[lidx], sizeof(double), 1,
                          init_state) != 1) {
                    log_err("End of model state file found unexpectedly");
                }
            }
            else {
                if (fscanf(init_state, " %lf",
                           &snow[veg].pack_liq[lidx]) == EOF) {
                    log_err("End of model state file found unexpectedly");
                }
            }
        }

        /* Read Read snow layer theta_ice */
        for (lidx = 0; lidx < snow[veg].Nsnow; lidx++) {
            if (options.STATE_FORMAT == BINARY) {
                if (fread(&snow[veg].theta_ice[lidx], sizeof(double), 1,
                          init_state) != 1) {
                    log_err("End of model state file found unexpectedly");
                }
            }
            else {
                if (fscanf(init_state, " %lf",
                           &snow[veg].theta_ice[lidx]) == EOF) {
                    log_err("End of model state file found unexpectedly");
                }
            }
        }

        /* Read Read snow layer theta_liq */
        for (lidx = 0; lidx < snow[veg].Nsnow; lidx++) {
            if (options.STATE_FORMAT == BINARY) {
                if (fread(&snow[veg].theta_liq[lidx], sizeof(double), 1,
                          init_state) != 1) {
                    log_err("End of model state file found unexpectedly");
                }
            }
            else {
                if (fscanf(init_state, " %lf",
                           &snow[veg].theta_liq[lidx]) == EOF) {
                    log_err("End of model state file found unexpectedly");
                }
            }
        }

        /* Read Read snow layer porosity */
        for (lidx = 0; lidx < snow[veg].Nsnow; lidx++) {
            if (options.STATE_FORMAT == BINARY) {
                if (fread(&snow[veg].porosity[lidx], sizeof(double), 1,
                          init_state) != 1) {
                    log_err("End of model state file found unexpectedly");
                }
            }
            else {
                if (fscanf(init_state, " %lf",
                           &snow[veg].porosity[lidx]) == EOF) {
                    log_err("End of model state file found unexpectedly");
                }
            }
        }
    }

    // Check that soil moisture does not exceed maximum allowed
    for (veg = 0; veg <= Nveg; veg++) {
        for (nidx = 0; nidx < options.Nnode; nidx++) {
            if (cell[veg].moist[nidx] >
                soil_con->porosity[lidx]) {
                log_warn("Initial soil moisture (%f) exceeds "
                         "maximum (%f) in layer %zu for veg tile "
                         "%d.  Resetting to maximum.",
                         cell[veg].moist[nidx],
                         soil_con->porosity[nidx], nidx, veg);

                cell[veg].ice[nidx] *=
                        soil_con->porosity[nidx] /
                        cell[veg].moist[nidx];

                cell[veg].moist[nidx] =
                    soil_con->porosity[nidx];
            }
            if (cell[veg].ice[nidx] > cell[veg].moist[nidx]) {
                cell[veg].ice[nidx] = cell[veg].moist[nidx];
            }
        }
    }
}
