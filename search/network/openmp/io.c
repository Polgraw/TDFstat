#include <H5Apublic.h>
#include <H5Ipublic.h>
#include <H5Tpublic.h>
#include <stdlib.h>
#include <time.h>
#include <hdf5.h>
#include <hdf5_hl.h>
#include "auxi.h"
#include "struct.h"

#define FORMAT_VERSION  1
#define TRIG_DSET_NAME  "triggers"
#define TRIG_RANK       1


int hdfout_init (char *outname, Command_line_opts *opts, Search_settings *sett, 
     Search_range *s_range, Trigger *sgnlv)
{
     int       i, j;
     hid_t     file, t_dataset, t_space, t_prop, filespace, memspace, space_scalar;
     hid_t     attr;
     hsize_t   dims[2] = {0,0};
     hsize_t   t_dim[TRIG_RANK] = {0};
     hsize_t   t_maxdim[TRIG_RANK] = {H5S_UNLIMITED};
     herr_t    hstat;
     
     hid_t vstr_type_id = H5Tcopy(H5T_C_S1);
     H5Tset_size(vstr_type_id, H5T_VARIABLE);
     
     // ------------------------------------------------------------------------
     // Define data types
     
     // triggers data type (sgnlv)
     hid_t ffstat_type = H5Tvlen_create(H5T_NATIVE_FLOAT);
     hid_t t_tid = H5Tcreate (H5T_COMPOUND, sizeof(Trigger));
     H5Tinsert(t_tid, "m", HOFFSET(Trigger, m), H5T_NATIVE_FLOAT);
     H5Tinsert(t_tid, "n", HOFFSET(Trigger, n), H5T_NATIVE_FLOAT);
     H5Tinsert(t_tid, "s", HOFFSET(Trigger, s), H5T_NATIVE_FLOAT);
     H5Tinsert(t_tid, "ra", HOFFSET(Trigger, ra), H5T_NATIVE_FLOAT);
     H5Tinsert(t_tid, "dec", HOFFSET(Trigger, dec), H5T_NATIVE_FLOAT);
     H5Tinsert(t_tid, "fdot", HOFFSET(Trigger, fdot), H5T_NATIVE_FLOAT);
     H5Tinsert(t_tid, "ffstat", HOFFSET(Trigger, ffstat), ffstat_type);

     // Command line options data type  (opts)
     hid_t cmd_opts_tid = H5Tcreate (H5T_COMPOUND, sizeof(Command_line_opts));
     H5Tinsert(cmd_opts_tid, "checkp_flag", HOFFSET(Command_line_opts, checkp_flag), H5T_NATIVE_INT);
     H5Tinsert(cmd_opts_tid, "veto_flag", HOFFSET(Command_line_opts, veto_flag), H5T_NATIVE_INT);
     //H5Tinsert(cmd_opts_tid, "gen_vlines_flag", HOFFSET(Command_line_opts, gen_vlines_flag), H5T_NATIVE_INT);
     //H5Tinsert(cmd_opts_tid, "help_flag", HOFFSET(Command_line_opts, help_flag), H5T_NATIVE_INT);
     H5Tinsert(cmd_opts_tid, "seg", HOFFSET(Command_line_opts, seg), H5T_NATIVE_INT);
     H5Tinsert(cmd_opts_tid, "band", HOFFSET(Command_line_opts, band), H5T_NATIVE_INT);
     H5Tinsert(cmd_opts_tid, "hemi", HOFFSET(Command_line_opts, hemi), H5T_NATIVE_INT);
     H5Tinsert(cmd_opts_tid, "nod", HOFFSET(Command_line_opts, nod), H5T_NATIVE_INT);
     H5Tinsert(cmd_opts_tid, "thr", HOFFSET(Command_line_opts, thr), H5T_NATIVE_DOUBLE);
     //H5Tinsert(cmd_opts_tid, "fpo_val", HOFFSET(Command_line_opts, fpo_val), H5T_NATIVE_DOUBLE);
     H5Tinsert(cmd_opts_tid, "narrowdown", HOFFSET(Command_line_opts, narrowdown), H5T_NATIVE_DOUBLE);
     H5Tinsert(cmd_opts_tid, "overlap", HOFFSET(Command_line_opts, overlap), H5T_NATIVE_DOUBLE);
     H5Tinsert(cmd_opts_tid, "indir", HOFFSET(Command_line_opts, indir), vstr_type_id);
     H5Tinsert(cmd_opts_tid, "outdir", HOFFSET(Command_line_opts, outdir), vstr_type_id);
     H5Tinsert(cmd_opts_tid, "range_file", HOFFSET(Command_line_opts, range_file), vstr_type_id);
     H5Tinsert(cmd_opts_tid, "grid_file", HOFFSET(Command_line_opts, grid_file), vstr_type_id);
     //H5Tinsert(cmd_opts_tid, "dump_range_file", HOFFSET(Command_line_opts, dump_range_file), H5T_C_S1);
     H5Tinsert(cmd_opts_tid, "usedet", HOFFSET(Command_line_opts, usedet), vstr_type_id);
     H5Tinsert(cmd_opts_tid, "addsig", HOFFSET(Command_line_opts, addsig), vstr_type_id);
     H5Tinsert(cmd_opts_tid, "fstat_norm", HOFFSET(Command_line_opts, fstat_norm), vstr_type_id);
     //H5Tinsert(cmd_opts_tid, "label", HOFFSET(Command_line_opts, label), H5T_C_S1);
     //H5Tinsert(cmd_opts_tid, "state_file", HOFFSET(Command_line_opts, state_file), H5T_C_S1);

     // Search settings data type  (sett)
     hid_t sett_tid = H5Tcreate (H5T_COMPOUND, sizeof(Search_settings));
     H5Tinsert(sett_tid, "fpo", HOFFSET(Search_settings, fpo), H5T_NATIVE_DOUBLE);
     H5Tinsert(sett_tid, "dt", HOFFSET(Search_settings, dt), H5T_NATIVE_DOUBLE);
     H5Tinsert(sett_tid, "B", HOFFSET(Search_settings, B), H5T_NATIVE_DOUBLE);
     H5Tinsert(sett_tid, "oms", HOFFSET(Search_settings, oms), H5T_NATIVE_DOUBLE);
     H5Tinsert(sett_tid, "omr", HOFFSET(Search_settings, omr), H5T_NATIVE_DOUBLE);
     H5Tinsert(sett_tid, "Smin", HOFFSET(Search_settings, Smin), H5T_NATIVE_DOUBLE);
     H5Tinsert(sett_tid, "Smax", HOFFSET(Search_settings, Smax), H5T_NATIVE_DOUBLE);
     H5Tinsert(sett_tid, "sepsm", HOFFSET(Search_settings, sepsm), H5T_NATIVE_DOUBLE);
     H5Tinsert(sett_tid, "cepsm", HOFFSET(Search_settings, cepsm), H5T_NATIVE_DOUBLE);
     H5Tinsert(sett_tid, "nfft", HOFFSET(Search_settings, nfft), H5T_NATIVE_INT);
     H5Tinsert(sett_tid, "nod", HOFFSET(Search_settings, nod), H5T_NATIVE_INT);
     H5Tinsert(sett_tid, "N", HOFFSET(Search_settings, N), H5T_NATIVE_INT);
     H5Tinsert(sett_tid, "nfftf", HOFFSET(Search_settings, nfftf), H5T_NATIVE_INT);
     H5Tinsert(sett_tid, "nmax", HOFFSET(Search_settings, nmax), H5T_NATIVE_INT);
     H5Tinsert(sett_tid, "nmin", HOFFSET(Search_settings, nmin), H5T_NATIVE_INT);
     H5Tinsert(sett_tid, "s", HOFFSET(Search_settings, s), H5T_NATIVE_INT);
     H5Tinsert(sett_tid, "nd", HOFFSET(Search_settings, nd), H5T_NATIVE_INT);
     H5Tinsert(sett_tid, "interpftpad", HOFFSET(Search_settings, interpftpad), H5T_NATIVE_INT);
     H5Tinsert(sett_tid, "fftpad", HOFFSET(Search_settings, fftpad), H5T_NATIVE_INT);
     H5Tinsert(sett_tid, "Ninterp", HOFFSET(Search_settings, Ninterp), H5T_NATIVE_INT);
     H5Tinsert(sett_tid, "nifo", HOFFSET(Search_settings, nifo), H5T_NATIVE_INT);
     H5Tinsert(sett_tid, "numlines_band", HOFFSET(Search_settings, numlines_band), H5T_NATIVE_INT);
     H5Tinsert(sett_tid, "nvlines_all_inband", HOFFSET(Search_settings, nvlines_all_inband), H5T_NATIVE_INT);     
     H5Tinsert(sett_tid, "bufsize", HOFFSET(Search_settings, bufsize), H5T_NATIVE_INT);
     H5Tinsert(sett_tid, "dd", HOFFSET(Search_settings, dd), H5T_NATIVE_INT);
     hid_t M_t = H5Tarray_create2(H5T_NATIVE_DOUBLE, 1, (hsize_t[]){16});
     H5Tinsert(sett_tid, "M", HOFFSET(Search_settings, M), M_t);
     hsize_t dims2d[2] = {sett->nvlines_all_inband, 2};
     hid_t lines_t = H5Tarray_create2(H5T_NATIVE_DOUBLE, 2, dims2d);
     H5Tinsert(sett_tid, "lines", HOFFSET(Search_settings, lines), lines_t);

     // search ranges data type  (s_range)     
     hid_t range_tid = H5Tcreate (H5T_COMPOUND, sizeof(Search_range));
     hid_t t1d_int = H5Tarray_create2(H5T_NATIVE_INT, 1, (hsize_t[]){2});
     hid_t t1d_float = H5Tarray_create2(H5T_NATIVE_FLOAT, 1, (hsize_t[]){2});
     H5Tinsert(range_tid, "pmr", HOFFSET(Search_range, pmr), t1d_int);
     H5Tinsert(range_tid, "fr", HOFFSET(Search_range, fr), t1d_int);
     H5Tinsert(range_tid, "mr", HOFFSET(Search_range, mr), t1d_float);
     H5Tinsert(range_tid, "nr", HOFFSET(Search_range, nr), t1d_float);
     H5Tinsert(range_tid, "spndr", HOFFSET(Search_range, spndr), t1d_float);
     H5Tinsert(range_tid, "mstep", HOFFSET(Search_range, mstep), H5T_NATIVE_FLOAT);
     H5Tinsert(range_tid, "nstep", HOFFSET(Search_range, nstep), H5T_NATIVE_FLOAT);
     H5Tinsert(range_tid, "sstep", HOFFSET(Search_range, sstep), H5T_NATIVE_FLOAT);
     H5Tinsert(range_tid, "mst", HOFFSET(Search_range, mst), H5T_NATIVE_FLOAT);
     H5Tinsert(range_tid, "nst", HOFFSET(Search_range, nst), H5T_NATIVE_FLOAT);
     H5Tinsert(range_tid, "sst", HOFFSET(Search_range, sst), H5T_NATIVE_FLOAT);
     H5Tinsert(range_tid, "pst", HOFFSET(Search_range, pst), H5T_NATIVE_INT);

     // ifo data type  (Detector_settings)
     hid_t ifo_tid = H5Tcreate(H5T_COMPOUND, sizeof(Detector_settings));
     hid_t name_type_id = H5Tcopy(H5T_C_S1);
     H5Tset_size(name_type_id, DETNAME_LENGTH);
     hid_t xdat_type_id = H5Tcopy(H5T_C_S1);
     H5Tset_size(xdat_type_id, FNAME_LENGTH);
     H5Tinsert(ifo_tid, "name", HOFFSET(Detector_settings, name), name_type_id);
     H5Tinsert(ifo_tid, "xdatname", HOFFSET(Detector_settings, xdatname), xdat_type_id);

     // ------------------------------------------------------------------------

     //Create the file.
     file = H5Fcreate(outname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
     
     // Basic attributes
     int fv = FORMAT_VERSION;
     hstat = H5LTset_attribute_int(file, "/", "format_version", &fv, 1);
     hstat = H5LTset_attribute_string(file, "/", "git_commit", CODEVER);

     char datetime_str[80];
     time_t now = time(NULL);
     struct tm *t = localtime(&now);
     strftime(datetime_str, sizeof(datetime_str), "%Y-%m-%d %H:%M:%S", t);
     hstat = H5LTset_attribute_string(file, "/", "t_start", datetime_str);
          
     // write structs as attributes
     hid_t scalar_space_id = H5Screate(H5S_SCALAR);
     
     attr = H5Acreate2(file, "opts", cmd_opts_tid, scalar_space_id, H5P_DEFAULT, H5P_DEFAULT);
     hstat = H5Awrite(attr, cmd_opts_tid, opts);
     H5Aclose(attr);
     
     attr = H5Acreate2(file, "sett", sett_tid, scalar_space_id, H5P_DEFAULT, H5P_DEFAULT);
     hstat = H5Awrite(attr, sett_tid, sett);
     H5Aclose(attr);

     attr = H5Acreate2(file, "s_range", range_tid, scalar_space_id, H5P_DEFAULT, H5P_DEFAULT);
     hstat = H5Awrite(attr, range_tid, s_range);
     H5Aclose(attr);

     hid_t ifo_space_id = H5Screate_simple(1, (hsize_t[]){sett->nifo}, NULL);
     attr = H5Acreate2(file, "ifo", ifo_tid, ifo_space_id, H5P_DEFAULT, H5P_DEFAULT);
     hstat = H5Awrite(attr, ifo_tid, &ifo);
     H5Aclose(attr);

     // Write triggers dataset
     hsize_t chunk_dims[1] = {3};
     t_prop = H5Pcreate(H5P_DATASET_CREATE);
     hstat = H5Pset_chunk(t_prop, TRIG_RANK, chunk_dims);

     t_space = H5Screate_simple(TRIG_RANK, t_dim, t_maxdim);
     t_dataset = H5Dcreate2(file, TRIG_DSET_NAME, t_tid, t_space, H5P_DEFAULT, t_prop, H5P_DEFAULT);
     hstat = H5Dwrite(t_dataset, t_tid, H5S_ALL, H5S_ALL, H5P_DEFAULT, sgnlv);
     
     // Release resources
     H5Pclose(t_prop);
     H5Dclose(t_dataset);
     H5Sclose(t_space);
     H5Fclose(file);
     H5Tclose(ffstat_type);

     return(EXIT_SUCCESS);
}



int hdfout_extend(char *outname, Trigger *sgnlv, int sgnlv_size)
{
     
     // Extend the triggers dataset by sgnlv_size entries from sgnlv bffer
     // Assume the file is already created and contains dataset "triggers"
     // The function closes all HDF5 objects it opens

     int        i, j;
     hid_t      file, t_dataset, t_space, t_prop, filespace, memspace, space_scalar;
     hsize_t    t_dim[TRIG_RANK] = {0};
     hsize_t    t_maxdim[TRIG_RANK] = {H5S_UNLIMITED};
     herr_t     hstat;
     
     // ------------------------------------------------------------------------
     // Define data types
     
     // triggers data type (sgnlv)
     hid_t ffstat_type = H5Tvlen_create(H5T_NATIVE_FLOAT);
     hid_t t_tid = H5Tcreate (H5T_COMPOUND, sizeof(Trigger));
     H5Tinsert(t_tid, "m", HOFFSET(Trigger, m), H5T_NATIVE_FLOAT);
     H5Tinsert(t_tid, "n", HOFFSET(Trigger, n), H5T_NATIVE_FLOAT);
     H5Tinsert(t_tid, "s", HOFFSET(Trigger, s), H5T_NATIVE_FLOAT);
     H5Tinsert(t_tid, "ra", HOFFSET(Trigger, ra), H5T_NATIVE_FLOAT);
     H5Tinsert(t_tid, "dec", HOFFSET(Trigger, dec), H5T_NATIVE_FLOAT);
     H5Tinsert(t_tid, "fdot", HOFFSET(Trigger, fdot), H5T_NATIVE_FLOAT);
     H5Tinsert(t_tid, "ffstat", HOFFSET(Trigger, ffstat), ffstat_type);

     // ------------------------------------------------------------------------

     //Open the file.
     file = H5Fopen(outname, H5F_ACC_RDWR, H5P_DEFAULT);
     
     // Open the dataset.
     t_dataset = H5Dopen2(file, TRIG_DSET_NAME, H5P_DEFAULT);
     
     // Get the dataspace.
     filespace = H5Dget_space(t_dataset);
     
     // Get the current size of the dataset.
     hstat = H5Sget_simple_extent_dims(filespace, t_dim, NULL);
     
     // Extend the dataset.
     t_dim[0] += sgnlv_size;
     hstat = H5Dset_extent(t_dataset, t_dim);
     
     // Close the old dataspace and get the new one.
     hstat = H5Sclose(filespace);
     filespace = H5Dget_space(t_dataset);
     
     // Define hyperslab in the extended portion of the dataset.
     hsize_t start[TRIG_RANK] = {t_dim[0]-sgnlv_size};
     hsize_t count[TRIG_RANK] = {sgnlv_size};
     hstat = H5Sselect_hyperslab(filespace, H5S_SELECT_SET, start, NULL, count, NULL);
     memspace = H5Screate_simple(TRIG_RANK, count, NULL);
     hstat = H5Dwrite(t_dataset, t_tid, memspace, filespace, H5P_DEFAULT, sgnlv);
     
     // Release resources
     H5Sclose(memspace);
     H5Sclose(filespace);
     H5Tclose(ffstat_type);
     H5Dclose(t_dataset);
     H5Fclose(file);
     
     return(EXIT_SUCCESS);
}



int hdfout_finalize(char *outname, int totsgnl, double time_elapsed, int nthreads)
{
     hid_t file;
     herr_t hstat;
     
     //Open the file.
     file = H5Fopen(outname, H5F_ACC_RDWR, H5P_DEFAULT);
     
     char datetime_str[80];
     time_t now = time(NULL);
     struct tm *t = localtime(&now);
     strftime(datetime_str, sizeof(datetime_str), "%Y-%m-%d %H:%M:%S", t);
     hstat = H5LTset_attribute_string(file, "/", "t_end", datetime_str);
     
     hstat = H5LTset_attribute_int(file, "/", "totsgnl", &totsgnl, 1);
     hstat = H5LTset_attribute_double(file, "/", "walltime", &time_elapsed, 1);
     hstat = H5LTset_attribute_int(file, "/", "num_threads", &nthreads, 1);

     H5Fclose(file);

     return(EXIT_SUCCESS);
}
