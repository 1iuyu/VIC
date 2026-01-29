/******************************************************************************
 * @section DESCRIPTION
 *
 * This subroutine saves the model state at second 0 of the date defined in the
 * global control file using STATEDAY, STATEMONTH, and STATEYEAR.  The saved
 * files can then be used to initialize the model to the same state as when the
 * files were created.
 *
 * Soil moisture, soil thermal, and snowpack variables  are stored for each
 * vegetation type and snow band.  However moisture variables from the
 * distributed precipitation model are averaged so that the model is restarted
 * with mu = 1.
 *****************************************************************************/

#include <vic_driver_classic.h>

/******************************************************************************
 * @brief    Saves the model state.
 *****************************************************************************/
void
write_model_state(all_vars_struct *all_vars,
                  veg_con_struct  *veg_con,
                  int              cellnum,
                  filep_struct    *filep,
                  soil_con_struct *soil_con)
{
    extern option_struct options;

    double               tmpval;
    int                  veg;
    int                  band;
    size_t               lidx;
    size_t               nidx;
    int                  Nbands;
    int                  Nbytes;
    int                  Nveg;

    cell_data_struct    *cell;
    snow_data_struct    *snow;
    energy_bal_struct   *energy;
    veg_var_struct      *veg_var;

    Nbands = options.SNOW_BAND;

    cell = all_vars->cell;
    veg_var = all_vars->veg_var;
    snow = all_vars->snow;
    energy = all_vars->energy;
    Nveg = veg_con->vegetat_type_num;
    /* write cell information */
    if (options.STATE_FORMAT == BINARY) {
        fwrite(&cellnum, sizeof(int), 1, filep->statefile);
        fwrite(&Nveg, sizeof(int), 1, filep->statefile);
        fwrite(&Nbands, sizeof(int), 1, filep->statefile);
    }
    else {
        fprintf(filep->statefile, "%i %i %i", cellnum, Nveg, Nbands);
    }
    // This stores the number of bytes from after this value to the end
    // of the line.  DO NOT CHANGE unless you have changed the values
    // written to the state file.
    // IF YOU EDIT THIS FILE: UPDATE THIS VALUE!
    if (options.STATE_FORMAT == BINARY) {
        Nbytes = options.Nnode * sizeof(double) + // dz_node
                 options.Nnode * sizeof(double) + // Zsum_node
                 (Nveg + 1) * 2 * sizeof(int) +   // veg & band
                 (Nveg + 1) * options.Nnode * sizeof(double) + // soil moisture
                 (Nveg + 1) * options.Nnode * sizeof(double) + // soil ice
                 Nveg * sizeof(double); // Wdew
        if (options.CARBON) {
            /* Carbon-specific state vars */
            Nbytes += Nveg * 5 * sizeof(double); // AnnualNPP, AnnualNPPPrev, and 3 soil carbon storages
        }
        Nbytes += (Nveg + 1) * options.Nnode *sizeof(double) + // SnowFrac
                  (Nveg + 1) * sizeof(double) + // Tcanopy
                  (Nveg + 1) * sizeof(double) + // Tair
                  (Nveg + 1) * sizeof(double) + // Tsurf
                  (Nveg + 1) * sizeof(double) + // Tgrnd
                  (Nveg + 1) * sizeof(double) + // Tupper
                  (Nveg + 1) * sizeof(double) + // Tlower
                  (Nveg + 1) * sizeof(double) * 6 + // coverage, swq, density, snow_canopy, SnowAge, old_swq
                  (Nveg + 1) * (snow[0].Nsnow + options.Nnode) * sizeof(double) + // soil temperatures
                  (Nveg + 1) * snow[0].Nsnow * sizeof(double) + // pack_ice
                  (Nveg + 1) * snow[0].Nsnow * sizeof(double) + // pack_liq
                  (Nveg + 1) * snow[0].Nsnow * sizeof(double) + // theta_ice
                  (Nveg + 1) * snow[0].Nsnow * sizeof(double) + // theta_liq
                  (Nveg + 1) * snow[0].Nsnow * sizeof(double) + // porosity
        fwrite(&Nbytes, sizeof(int), 1, filep->statefile);
    }

    /* Write soil thermal node deltas */
    for (nidx = 0; nidx < options.Nnode; nidx++) {
        if (options.STATE_FORMAT == BINARY) {
            fwrite(&soil_con->dz_soil[nidx], sizeof(double), 1,
                   filep->statefile);
        }
        else {
            fprintf(filep->statefile, " "ASCII_STATE_FLOAT_FMT,
                    soil_con->dz_soil[nidx]);
        }
    }
    /* Write soil thermal node depths */
    for (nidx = 0; nidx < options.Nnode; nidx++) {
        if (options.STATE_FORMAT == BINARY) {
            fwrite(&soil_con->Zsum_soil[nidx], sizeof(double), 1,
                   filep->statefile);
        }
        else {
            fprintf(filep->statefile, " "ASCII_STATE_FLOAT_FMT,
                    soil_con->Zsum_soil[nidx]);
        }
    }
    if (options.STATE_FORMAT == ASCII) {
        fprintf(filep->statefile, "\n");
    }

    /* Output for all vegetation types */
    for (veg = 0; veg <= Nveg; veg++) {
        /* Output for all snow bands */
        band = veg_con[veg].BandIndex;
        /* Write cell identification information */
        if (options.STATE_FORMAT == BINARY) {
            fwrite(&veg, sizeof(int), 1, filep->statefile);
            fwrite(&band, sizeof(int), 1, filep->statefile);
        }
        else {
            fprintf(filep->statefile, "%i %i", veg, band);
        }

        /* Write total soil moisture */
        for (lidx = 0; lidx < options.Nnode; lidx++) {
            tmpval = cell[veg].moist[lidx];
            if (options.STATE_FORMAT == BINARY) {
                fwrite(&tmpval, sizeof(double), 1, filep->statefile);
            }
            else {
                fprintf(filep->statefile, " "ASCII_STATE_FLOAT_FMT, tmpval);
            }
        }

        /* Write average ice content */
        for (lidx = 0; lidx < options.Nnode; lidx++) {
            tmpval = cell[veg].ice[lidx];
            if (options.STATE_FORMAT == BINARY) {
                fwrite(&tmpval, sizeof(double), 1, filep->statefile);
            }
            else {
                fprintf(filep->statefile, " "ASCII_STATE_FLOAT_FMT,
                        tmpval);
            }
        }

        if (veg < Nveg) {
            /* Write dew storage */
            tmpval = veg_var[veg].Wdew;
            if (options.STATE_FORMAT == BINARY) {
                fwrite(&tmpval, sizeof(double), 1, filep->statefile);
            }
            else {
                fprintf(filep->statefile, " "ASCII_STATE_FLOAT_FMT, tmpval);
            }
            if (options.CARBON) {
                /* Write cumulative NPP */
                tmpval = veg_var[veg].AnnualNPP;
                if (options.STATE_FORMAT == BINARY) {
                    fwrite(&tmpval, sizeof(double), 1, filep->statefile);
                }
                else {
                    fprintf(filep->statefile, " "ASCII_STATE_FLOAT_FMT,
                            tmpval);
                }
                tmpval = veg_var[veg].AnnualNPPPrev;
                if (options.STATE_FORMAT == BINARY) {
                    fwrite(&tmpval, sizeof(double), 1, filep->statefile);
                }
                else {
                    fprintf(filep->statefile, " "ASCII_STATE_FLOAT_FMT,
                            tmpval);
                }
                /* Write soil carbon storages */
                tmpval = cell[veg].CLitter;
                if (options.STATE_FORMAT == BINARY) {
                    fwrite(&tmpval, sizeof(double), 1, filep->statefile);
                }
                else {
                    fprintf(filep->statefile, " "ASCII_STATE_FLOAT_FMT,
                            tmpval);
                }
                tmpval = cell[veg].CInter;
                if (options.STATE_FORMAT == BINARY) {
                    fwrite(&tmpval, sizeof(double), 1, filep->statefile);
                }
                else {
                    fprintf(filep->statefile, " "ASCII_STATE_FLOAT_FMT,
                            tmpval);
                }
                tmpval = cell[veg].CSlow;
                if (options.STATE_FORMAT == BINARY) {
                    fwrite(&tmpval, sizeof(double), 1, filep->statefile);
                }
                else {
                    fprintf(filep->statefile, " "ASCII_STATE_FLOAT_FMT,
                            tmpval);
                }
            }
        }

        /* Write snow data */
        if (options.STATE_FORMAT == BINARY) {
            fwrite(&snow[veg].coverage, sizeof(int), 1,
                   filep->statefile);
            fwrite(&snow[veg].swq, sizeof(char), 1,
                   filep->statefile);
            fwrite(&snow[veg].density, sizeof(double), 1,
                   filep->statefile);
            fwrite(&snow[veg].snow_canopy, sizeof(double), 1,
                   filep->statefile);
            fwrite(&snow[veg].SnowAge, sizeof(double), 1,
                   filep->statefile);
            fwrite(&snow[veg].old_swq, sizeof(double), 1,
                   filep->statefile);
        }
        else {
            fprintf(filep->statefile,
                    ASCII_STATE_FLOAT_FMT " "ASCII_STATE_FLOAT_FMT " "
                    ASCII_STATE_FLOAT_FMT " "ASCII_STATE_FLOAT_FMT " "
                    ASCII_STATE_FLOAT_FMT " "ASCII_STATE_FLOAT_FMT,
                    snow[veg].coverage, snow[veg].swq,
                    snow[veg].density, snow[veg].snow_canopy,
                    snow[veg].SnowAge, snow[veg].old_swq);
        }

        /* Write thermal node temperatures */
        for (nidx = 0; nidx < options.Nnode + snow[veg].Nsnow; nidx++) {
            if (options.STATE_FORMAT == BINARY) {
                fwrite(&energy[veg].T[nidx], sizeof(double), 1,
                       filep->statefile);
            }
            else {
                fprintf(filep->statefile, " "ASCII_STATE_FLOAT_FMT,
                        energy[veg].T[nidx]);
            }
        }

        /* Write pack_ice */
        for (lidx = 0; lidx < snow[veg].Nsnow; lidx++) {
            if (options.STATE_FORMAT == BINARY) {
                fwrite(&snow[veg].pack_ice[lidx], sizeof(double), 1,
                       filep->statefile);
            }
            else {
                fprintf(filep->statefile, " "ASCII_STATE_FLOAT_FMT,
                        snow[veg].pack_ice[lidx]);
            }
        }

        /* Write pack_liq */
        for (lidx = 0; lidx < snow[veg].Nsnow; lidx++) {
            if (options.STATE_FORMAT == BINARY) {
                fwrite(&snow[veg].pack_liq[lidx], sizeof(double), 1,
                       filep->statefile);
            }
            else {
                fprintf(filep->statefile, " "ASCII_STATE_FLOAT_FMT,
                        snow[veg].pack_liq[lidx]);
            }
        }

        /* Write theta_ice */
        for (lidx = 0; lidx < snow[veg].Nsnow; lidx++) {
            if (options.STATE_FORMAT == BINARY) {
                fwrite(&snow[veg].theta_ice[lidx], sizeof(double), 1,
                       filep->statefile);
            }
            else {
                fprintf(filep->statefile, " "ASCII_STATE_FLOAT_FMT,
                        snow[veg].theta_ice[lidx]);
            }
        }

        /* Write theta_liq */
        for (lidx = 0; lidx < snow[veg].Nsnow; lidx++) {
            if (options.STATE_FORMAT == BINARY) {
                fwrite(&snow[veg].theta_liq[lidx], sizeof(double), 1,
                       filep->statefile);
            }
            else {
                fprintf(filep->statefile, " "ASCII_STATE_FLOAT_FMT,
                        snow[veg].theta_liq[lidx]);
            }
        }

        /* Write snowpack porosity */
        for (lidx = 0; lidx < snow[veg].Nsnow; lidx++) {
            if (options.STATE_FORMAT == BINARY) {
                fwrite(&snow[veg].porosity[lidx], sizeof(double), 1,
                       filep->statefile);
            }
            else {
                fprintf(filep->statefile, " "ASCII_STATE_FLOAT_FMT,
                        snow[veg].porosity[lidx]);
            }
        }

        if (options.STATE_FORMAT == ASCII) {
            fprintf(filep->statefile, "\n");
        }
    }
    /* Force file to be written */
    fflush(filep->statefile);
}
