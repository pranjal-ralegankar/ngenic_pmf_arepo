#include "allvars.h"


struct io_header_1 header1, header;

int WhichSpectrum;


int SphereMode;
int *Local_nx_table;

FILE *FdTmp, *FdTmpInput;

int Nmesh, Nsample;

long long IDStart;



char GlassFile[500];
char FileWithInputSpectrum[500];
char FileWithInputSpectrumB[500];
char FileWithInputSpectrumDM[500];
char FileWithInputSpectrumBpmf[500];
char FileWithInputSpectrumDMpmf[500];
char FileWithTransfer[500];

int GlassTileFac;

double Box;
int Seed;
int Seedb;
double klimitB;
long long TotNumPart;

double B1mpc;
double nB;
double eta;
double lambda_pmf_b, lambda_pmf_dm, lambda_lcdm_b, lambda_lcdm_dm;
double boost;

int NumPart;

int *Slab_to_task;

int NTaskWithN;

struct part_data *P;

int Nglass;

double InitTime;
double Redshift;
double MassTable[6];


char OutputDir[100], FileBase[100];
int NumFilesWrittenInParallel;


int ThisTask, NTask;

int Local_nx, Local_x_start;

int IdStart;

rfftwnd_mpi_plan Inverse_plan;
fftw_real *Disp, *Workspace;
fftw_complex *Cdata;


double UnitTime_in_s, UnitLength_in_cm, UnitMass_in_g, UnitVelocity_in_cm_per_s;
double InputSpectrum_UnitLength_in_cm;
double G, Hubble;
double RhoCrit;

double Omega, OmegaLambda, OmegaDM_2ndSpecies, Sigma8;
double OmegaBaryon, HubbleParam;
double ShapeGamma;
double PrimordialIndex;
double Dplus;			/* growth factor */

#ifdef DIFFERENT_TRANSFER_FUNC
int Type, MinType, MaxType;
#endif

int ReNormalizeInputSpectrum;

int WDM_On;
int WDM_Vtherm_On;
double WDM_PartMass_in_kev;


int NU_On;
int NU_Vtherm_On;
double NU_PartMass_in_ev;

