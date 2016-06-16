// Do forward-backward, value 1
// Do original cnF2freq tree evaluation, value 0
#define DOFB 1
const int SELFBITS = 2;

const int INDCOUNT = 1000000;
const bool DOREMAPDISTANCES = false;
const bool DOINFPROBS = false;
const bool SELFING = false;
const bool RELSKEWS = true;

// F2 with haplotyping
const int NUMGEN = 3;
const int TYPEBITS = (1 << NUMGEN) - 2;
const int TYPESEXES[TYPEBITS + SELFING * SELFBITS + RELSKEWS] = { 0, 0, 1, 1, 0, 1/*, 0, 0*//*, 0*/ };
//const int TYPEGENS[TYPEBITS] = {1, 0, 0, 1, 0, 0};
const int TYPEGENS[TYPEBITS + SELFING * 2] = { 1, 0, 0, 1, 0, 0/*, 2, 2*//*, 3*/ };

const int BITS_W_SELF = TYPEBITS + (SELFING ? SELFBITS : 0);
const int TOTBITS = BITS_W_SELF + (RELSKEWS ? 1 : 0);
const int NUMTYPES = 1 << TOTBITS;
const double EVENGEN = 1.0 / NUMTYPES;
const float MINFACTOR = -1e15;
const unsigned int NUMFLAG2GEN = NUMGEN;
const unsigned int HALFNUMPATHS = 1 << (TYPEBITS / 2);
const unsigned int NUMPATHS = 1 << (TYPEBITS + 1);
const unsigned int NUMSHIFTGEN = NUMGEN - 1;
const unsigned int HALFNUMSHIFTS = 1 << ((1 << (NUMSHIFTGEN - 1)) - 1);
const unsigned int NUMSHIFTS = 1 << ((1 << NUMSHIFTGEN) - 1);
const bool HAPLOTYPING = true;
const int NONSELFNUMTYPES = NUMTYPES >> (SELFING * 2);

// Only NUMTYPES * 3 unless we have RELSKEWS messing things up as well
// TODO: Reorder bits...
const int VALIDSELFNUMTYPES = RELSKEWS ? NUMTYPES : (NUMTYPES - SELFING * (NUMTYPES >> 2));

//#define DOEXTERNFORGCC

#ifdef DOEXTERNFORGCC
#define EXTERNFORGCC extern
#else
#define EXTERNFORGCC
#endif

// F2 with no haplotyping
/*const int NUMGEN = 2;
const int TYPEBITS = (1 << NUMGEN) - 2;
const int TYPESEXES[TYPEBITS] = {0, 1};
const int NUMTYPES = 1 << TYPEBITS;
const double EVENGEN = 1.0 / NUMTYPES;
const float MINFACTOR = -1e15;
const unsigned int NUMFLAG2GEN = 1;
const unsigned int HALFNUMPATHS = 1;
const unsigned int NUMPATHS = 2;
const unsigned int NUMSHIFTGEN = 0;
const unsigned int HALFNUMSHIFTS = 1;
const unsigned int NUMSHIFTS = 1;
const bool HAPLOTYPING = false;*/

// USED IN QTLMAS15
/*
const int NUMGEN = 2;
const int TYPEBITS = (1 << NUMGEN) - 2;
const int TYPESEXES[TYPEBITS] = {0, 1};
const int NUMTYPES = 1 << TYPEBITS;
const double EVENGEN = 1.0 / NUMTYPES;
const float MINFACTOR = -1e15;
const unsigned int NUMFLAG2GEN = NUMGEN;
const unsigned int HALFNUMPATHS = 1 << (TYPEBITS / 2);
const unsigned int NUMPATHS = NUMTYPES << 1;
const unsigned int NUMSHIFTGEN = NUMGEN - 1;
const unsigned int HALFNUMSHIFTS = 1 << ((1 << (NUMSHIFTGEN - 1)) - 1);
const unsigned int NUMSHIFTS = 1 << ((1 << NUMSHIFTGEN) - 1);
/*const unsigned int NUMSHIFTGEN = 0;
const unsigned int HALFNUMSHIFTS = 1;
const unsigned int NUMSHIFTS = 1;*/
//const bool HAPLOTYPING = true;




const int HALFNUMTYPES = 1 << (TYPEBITS / 2);

// Infer corrections for impossible genotypes according to the pedigree
// Also infer values for missing markers from existing information in offspring
// When enabled, results similar to QTL Express without data reduction
// ccoeff does not provide correction inference, so exact result reproduction
// is achieved when this flag is disabled.
bool CORRECTIONINFERENCE = true;
