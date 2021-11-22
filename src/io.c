#include "globals.h"
#include "io.h"
#include <inttypes.h>

/* These two handle the F90 Records required for the file format */
#define WRITE_F90REC  {my_fwrite(&blocksize,sizeof(int),1,fp);}

struct GADGET_Header Header;
int blksize;

void write_header(FILE *fp);

void Write_output()
{
    FILE *fp;

    size_t nBytes = 256 + 0.5 * Param.Ntotal * sizeof(*P)
            + 0.5 * Param.Ntotal * sizeof(*SphP);

    printf("Output : \n"
            "   File Name = %s\n"
            "   File Size ~ %.1f MB\n"
            ,Param.Output_File, nBytes/1e6);

    if (!(fp=fopen(Param.Output_File,"w")))
        fprintf(stderr, "Can't open file %s\n"
                ,Param.Output_File);

    write_header(fp);

    for (int iblock=0; iblock<IO_LASTENTRY; iblock++)
        add_block(fp, (enum iofields) iblock);

    fclose(fp);

    printf("done\n");

    return;
}

void write_header(FILE *fp)
{
    int i, nextblock, blocksize;

    printf("   Header\n");

    /* Set Header */
    for (i=0; i<6; i++) {

        Header.npart[i] = Param.Npart[i];

        Header.mass[i] = Param.Mpart[i];
        Header.npartTotal[i] = Header.npart[i];
    }

    Header.time = 0;
    Header.redshift = 0;
    Header.flag_sfr = 0;
    Header.flag_feedback = 0;
    Header.flag_cooling = 0;
    Header.num_files = 1;
    Header.BoxSize = Param.Boxsize;
    Header.Omega0 = 1;
    Header.OmegaLambda = 0.7;
    Header.HubbleParam = Cosmo.h_100;
    Header.flag_stellarage = 0;
    Header.flag_metals = 0;

    /* Write F90 Unformatted Record */
    blocksize = sizeof(int) + 4 * sizeof(char);
    WRITE_F90REC                    /* <- uses "blocksize" */
    my_fwrite((void *) "HEAD", sizeof(char), 4, fp);
    nextblock = sizeof(Header) + 2 * sizeof(int);
    my_fwrite(&nextblock, sizeof(int), 1, fp);
    WRITE_F90REC

    blocksize = sizeof(Header);
    WRITE_F90REC
    my_fwrite(&Header, blocksize,1,fp);
    WRITE_F90REC

    return;
}

void add_block(FILE *fp, enum iofields iblock)
{
    int  ipart, ibuf, type, offset, imax, nextblock, blocksize;
    void *write_buffer = NULL;
    size_t  nData;

    set_block_info((enum iofields) iblock);

    printf("   Block %d (%s)\n",iblock, Block.Name);

    nData = Block.Ntot*Block.Val_per_element * Block.Bytes_per_element;

    write_buffer = Malloc(nData);

    offset = 0;
    imax = Block.Npart[0];

    for (type=ibuf=0; type<6; type++){

        for (ipart=offset;ipart<imax; ipart++) {

            fill_write_buffer(iblock, write_buffer, ipart, ibuf);

            ibuf += Block.Val_per_element;
        }

        offset += Block.Npart[type];

        if (type < 5)
            imax += Block.Npart[type+1];
    }

    blocksize = sizeof(int) + 4 * sizeof(char);
    WRITE_F90REC
    my_fwrite(&Block.Label, sizeof(char), 4, fp);
    nextblock = nData + 2 * sizeof(int);
    my_fwrite(&nextblock, sizeof(int), 1, fp);
    WRITE_F90REC

    blocksize = nData;
    WRITE_F90REC
    my_fwrite(write_buffer, blocksize, 1, fp);
    blocksize = nData;
    WRITE_F90REC

    free(write_buffer);

    return;
}

void fill_write_buffer(enum iofields blocknr, void *wbuf,
        size_t ipart, size_t ibuf)
{
    int i;

    switch (blocknr) {
        case IO_POS:
            for (i=0; i<3; i++)
                ((float *)wbuf)[ibuf+i] = P[ipart].Pos[i];
        break;
        case IO_VEL:
            for (i=0; i<3; i++)
                ((float *)wbuf)[ibuf+i] = P[ipart].Vel[i];
        break;
        case IO_ID:
            ((unsigned int *)wbuf)[ibuf] = P[ipart].ID;
        break;
        case IO_RHO:
            ((float *)wbuf)[ibuf] = SphP[ipart].Rho;
        break;
        case IO_RHOMODEL:
            ((float *)wbuf)[ibuf] = SphP[ipart].Rho_Model;
        break;
        case IO_HSML:
            ((float *)wbuf)[ibuf] = SphP[ipart].Hsml;
        break;
        case IO_U:
            ((float *)wbuf)[ibuf] = SphP[ipart].U;
        break;
        case IO_BFLD:
            for (i=0; i<3; i++)
                ((float *)wbuf)[ibuf+i] = SphP[ipart].Bfld[i];
        break;
        default:
            Assert(0, "Block not found %d",blocknr);
        break;
    }

    return;
}
void set_block_info(enum iofields blocknr)
{
    int i = 0;

    for (i = Block.Ntot = 0; i < 6; i++)
        Block.Npart[i] = 0;

    strncpy(Block.Name, "                   ",22);

    switch (blocknr) {
        case IO_POS:
        strncpy(Block.Label,"POS ",4);            /* Has to be 4 Letters */
        strncpy(Block.Name, "Coordinates",13);
        for (i=0; i<6; i++)
            Block.Npart[i] = Param.Npart[i];

        Block.Val_per_element = 3;
        Block.Bytes_per_element = sizeof(P[0].Pos[0]);
        break;
        case IO_VEL:
        strncpy(Block.Label,"VEL ",4);
        strncpy(Block.Name, "Velocities",11);
        for (i=0; i<6; i++)
            Block.Npart[i] = Param.Npart[i];
        Block.Val_per_element = 3;
        Block.Bytes_per_element = sizeof(P[0].Vel[0]);
        break;
        case IO_ID:
        strncpy(Block.Label,"ID  ",4);
        strncpy(Block.Name, "ParticleIDs",13);
        for (i=0; i<6; i++)
            Block.Npart[i] = Param.Npart[i];

        Block.Val_per_element = 1;
        Block.Bytes_per_element = sizeof(unsigned int);
        break;
        case IO_RHO:
        strncpy(Block.Label,"RHO ",4);
        strncpy(Block.Name, "Density",16);
        Block.Npart[0] = Param.Npart[0];
        Block.Val_per_element = 1;
        Block.Bytes_per_element = sizeof(SphP[0].Rho);
        break;
        case IO_RHOMODEL:
        strncpy(Block.Label,"RHOM",4);
        strncpy(Block.Name, "Model Density",16);
        Block.Npart[0] = Param.Npart[0];
        Block.Val_per_element = 1;
        Block.Bytes_per_element = sizeof(SphP[0].Rho_Model);
        break;
        case IO_HSML:
        strncpy(Block.Label,"HSML",4);
        strncpy(Block.Name, "SmoothingLength",16);
        Block.Npart[0] = Param.Npart[0];
        Block.Val_per_element = 1;
        Block.Bytes_per_element = sizeof(SphP[0].Hsml);
        break;
        case IO_U:
        strncpy(Block.Label,"U   ",4);
        strncpy(Block.Name, "InternalEnergy",16);
        Block.Npart[0] = Param.Npart[0];
        Block.Val_per_element = 1;
        Block.Bytes_per_element = sizeof(SphP[0].U);
        break;
        case IO_BFLD:
        strncpy(Block.Label,"BFLD ",4);
        strncpy(Block.Name, "MagneticField",16);
        Block.Npart[0] = Param.Npart[0];
        Block.Val_per_element = 3;
        Block.Bytes_per_element = sizeof(SphP[0].Bfld[0]);
        break;
        case IO_LASTENTRY:
        strncpy(Block.Label,"LAST",4);
        strncpy(Block.Name, " ",1);
        for (i=0; i<6; i++)
            Block.Npart[i] = 0;
        Block.Bytes_per_element = 0;
        Block.Val_per_element = 0;
        break;
        default:
            fprintf(stderr, "Block not found %d",blocknr);
            exit(5);
        break;
    }

    for (int i = Block.Ntot = 0; i < 6; i++)
        Block.Ntot += Block.Npart[i];

    return;
}

/* Basic routine to read data from a file
 * */
size_t my_fwrite(void *data, size_t size, size_t nItems, FILE *fp)
{
    size_t nWritten = 0;

    if (nItems*size > 0) {

        if ((nWritten = fwrite(data, size, nItems, fp)) != nItems) {

            fprintf(stderr,"I/O error (fwrite) ");

            fflush(stderr);

            exit(6);
        }

    } else
        nWritten = 0;

    return nWritten;
}

#undef WRITE_F90REC

/* Read a number of tags from an ascii file
 * the comment sign is % */
#define LASTPARAMETERID -1111
#define REAL 1
#define STRING 2
#define INT 3

void Read_param_file(char *filename)
{
    FILE *fd;

    int id[MAXTAGS];
    void *addr[MAXTAGS];
    char tag[MAXTAGS][50];

    char buf[CHARBUFSIZE], buf1[CHARBUFSIZE];
    char buf2[CHARBUFSIZE],buf3[2 * CHARBUFSIZE];
    int tagDone[MAXTAGS]= {0};
    int i, j, nt=0;

    strcpy(tag[nt], "Output_file");
    addr[nt] = &Param.Output_File;
    id[nt++] = STRING;

    strcpy(tag[nt], "Ntotal");
    addr[nt] = &Param.Ntotal;
    id[nt++] = INT;

    strcpy(tag[nt], "Mtotal");
    addr[nt] = &Param.Mtot200;
    id[nt++] = REAL;

    strcpy(tag[nt], "Redshift");
    addr[nt] = &Param.Redshift;
    id[nt++] = REAL;

    strcpy(tag[nt], "Mass_Ratio");
    addr[nt] = &Param.Mass_Ratio;
    id[nt++] = REAL;

    strcpy(tag[nt], "ImpactParam");
    addr[nt] = &Param.Impact_Param;
    id[nt++] = REAL;

    strcpy(tag[nt], "ZeroEOrbitFrac");
    addr[nt] = &Param.Zero_Energy_Orbit_Fraction;
    id[nt++] = REAL;

    strcpy(tag[nt], "Cuspy");
    addr[nt] = &Param.Cuspy;
    id[nt++] = INT;

    strcpy(tag[nt], "beta_0");
    addr[nt] = &Halo[0].Beta;
    id[nt++] = REAL;

    strcpy(tag[nt], "beta_1");
    addr[nt] = &Halo[1].Beta;
    id[nt++] = REAL;

    strcpy(tag[nt], "Bfld_Norm");
    addr[nt] = &Param.Bfld_Norm;
    id[nt++] = REAL;

    strcpy(tag[nt], "Bfld_Eta");
    addr[nt] = &Param.Bfld_Eta;
    id[nt++] = REAL;

    strcpy(tag[nt], "bf_0");
    addr[nt] = &Halo[0].bf;
    id[nt++] = REAL;

    strcpy(tag[nt], "bf_1");
    addr[nt] = &Halo[1].bf;
    id[nt++] = REAL;

    strcpy(tag[nt], "rcut_r200_ratio_0");
    addr[nt] = &Halo[0].Rcut_R200_Ratio;
    id[nt++] = REAL;

    strcpy(tag[nt], "rcut_r200_ratio_1");
    addr[nt] = &Halo[1].Rcut_R200_Ratio;
    id[nt++] = REAL;

    /* System of Units */
    strcpy(tag[nt], "UnitLength_in_cm");
    addr[nt] = &Unit.Length;
    id[nt++] = REAL;

    strcpy(tag[nt], "UnitMass_in_g");
    addr[nt] = &Unit.Mass;
    id[nt++] = REAL;

    strcpy(tag[nt], "UnitVelocity_in_cm_per_s");
    addr[nt] = &Unit.Vel;
    id[nt++] = REAL;

#ifdef GIVEPARAMS
    strcpy(tag[nt], "c_nfw_0");
    addr[nt] = &Halo[0].C_nfw;
    id[nt++] = REAL;

    strcpy(tag[nt], "v_com_0");
    addr[nt] = &Param.VelMerger[0];
    id[nt++] = REAL;

    strcpy(tag[nt], "rc_0");
    addr[nt] = &Halo[0].Rcore;
    id[nt++] = REAL;



    strcpy(tag[nt], "c_nfw_1");
    addr[nt] = &Halo[1].C_nfw;
    id[nt++] = REAL;

    strcpy(tag[nt], "v_com_1");
    addr[nt] = &Param.VelMerger[1];
    id[nt++] = REAL;

    strcpy(tag[nt], "rc_1");
    addr[nt] = &Halo[1].Rcore;
    id[nt++] = REAL;

#endif

#ifdef ADD_THIRD_SUBHALO

    strcpy(tag[nt], "SubFirstMass");
    addr[nt] = &Param.SubFirstMass;
    id[nt++] = REAL;

    strcpy(tag[nt], "SubFirstPos0");
    addr[nt] = &Param.SubFirstPos[0];
    id[nt++] = REAL;

    strcpy(tag[nt], "SubFirstPos1");
    addr[nt] = &Param.SubFirstPos[1];
    id[nt++] = REAL;

    strcpy(tag[nt], "SubFirstPos2");
    addr[nt] = &Param.SubFirstPos[2];
    id[nt++] = REAL;

    strcpy(tag[nt], "SubFirstVel0");
    addr[nt] = &Param.SubFirstVel[0];
    id[nt++] = REAL;

    strcpy(tag[nt], "SubFirstVel1");
    addr[nt] = &Param.SubFirstVel[1];
    id[nt++] = REAL;

    strcpy(tag[nt], "SubFirstVel2");
    addr[nt] = &Param.SubFirstVel[2];
    id[nt++] = REAL;

#endif

#ifdef DOUBLE_BETA_COOL_CORES
    strcpy(tag[nt], "Rho0_Fac");
    addr[nt] = &Param.Rho0_Fac;
    id[nt++] = REAL;

    strcpy(tag[nt], "Rc_Fac");
    addr[nt] = &Param.Rc_Fac;
    id[nt++] = REAL;
#endif

#ifdef TURB_B_FIELD
    strcpy(tag[nt], "Bfld_Scale");
    addr[nt] = &Param.Bfld_Scale;
    id[nt++] = REAL;

    strcpy(tag[nt], "Turb_Spectral_Index");
    addr[nt] = &Param.Spectral_Index;
    id[nt++] = REAL;
#endif

    /* Add above */
    id[nt] = LASTPARAMETERID;

    if ((fd = fopen(filename, "r"))) {
        sprintf(buf, "%s%s", filename, "-usedvalues");
        printf("\nReading Parameter file : %s \n\n",
                filename);
        while (fgets(buf, CHARBUFSIZE, fd)) {
            if (sscanf(buf, "%s%s%s", buf1, buf2,buf3) < 2)
                    continue;

            if (buf1[0] == '%')
                    continue;

                for (i = 0, j = -1; i < nt; i++) {

                    if (!(strcmp(buf1, tag[i])) && (tagDone[i] != 1) ){
                        j = i;
                        tagDone[i] = 1;
                        break;
                    }
                }

                if (j >= 0) {
                    switch (id[j]) {
                    case REAL:
                        *((double *)addr[j]) = atof(buf2);
                        break;
                    case STRING:
                        strcpy((char *)addr[j], buf2);
                        break;
                    case INT:
                        *((int *)addr[j]) = atoi(buf2);
                        break;
                    }
                }
            }
            fclose(fd);
    } else {
        fprintf(stderr, "Parameter file %s not found.\n",
            filename);
        exit(1);
    }

    for (i = 0; i < nt; i++) {
        if (!tagDone[i]) {
            fprintf(stderr,
                "Value for tag '%s' missing in parameter file '%s'.\n",
                tag[i], filename);
            exit(1);
        }
    }

#ifndef GIVEPARAMS // Keep Xm < 0 so we don't screw up Dvel
     //Param.Mass_Ratio = fmin(Param.Mass_Ratio, 1/Param.Mass_Ratio);
    if (Param.Mass_Ratio == 0)
         Sub.First = 1;
    else
        Sub.First = 2;
#endif

     return;
}
#undef REAL
#undef STRING
#undef INT
#undef LASTPARAMETERID
