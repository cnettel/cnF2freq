// cnF2freq, (c) Carl Nettelblad, Department of Cell and Molecular Biology, Uppsala University
// 2008-2014
//
// Modified version of Plantimpute with haplotype skewness values written.
//
// carl.nettelblad@it.uu.se
//
// This code is allowed to be freely used for any commercial or research purpose. If the code is
// used or integrated into another project largely unchanged, attribution to the original author
// should be included. No warranties are given.


#define NDEBUG
// These defines fixed an error in one particular site installation of the Portland compiler.
#define _STLP_EXPOSE_GLOBALS_IMPLEMENTATION 1
#define _REENTRANT 1
#define _SECURE_SCL 0
// For MSCVC
// Note: MSVC OpenMP support is insufficient in current release. Disable OpenMP for compilation
// in MSVC.
//
// Recent releases of g++ on MacOS X and Linux, as well as the Intel C++ compiler on
// Linux (x86 and x64) and Windows have been tested. The Portland compiler collection works only
// with some optimization settings, some race conditions in STL, although some
// workarounds are used.
#define _CRT_SECURE_NO_WARNINGS

#include <vector>
#include <string.h>
#include <stdio.h>

float templgeno[8] = {-1, -0.5,
	0,  0.5,
	0,  0.5,
	1, -0.5};



// _MSC_VER is here to be interpreted as any compiler providing TR1 C++ headers
//#ifdef _MSC_VER
//#include <array>
//#else
// Boost also provides an array implementation, that is largely compatible
#include <boost/array.hpp>
//#endif

#include <boost/math/distributions/binomial.hpp>
#include <boost/units/quantity.hpp>
#include <boost/units/conversion.hpp>
#include <boost/units/make_system.hpp>
#include <boost/units/base_unit.hpp>
#include <boost/units/base_dimension.hpp>
#include <boost/units/static_constant.hpp>
#include <boost/static_assert.hpp>
#include <boost/random/linear_congruential.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>
/*#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi.hpp>*/


#include <errno.h>
#include <assert.h>
#include <stdlib.h>
#include <set>
#include <algorithm>
#include <math.h>
#include <map>
#include <float.h> // use these libraries


using namespace std; // use functions that are part of the standard library
#ifdef _MSC_VER
//using namespace tr1;
#else
using namespace boost;
#endif

using namespace boost;
//using namespace boost::mpi;
using namespace boost::random;
using namespace boost::units;


#define none cnF2freqNONE


#ifndef _MSC_VER
#define _isnan isnan
#define _finite finite
#endif


mt19937 rng;
int myrand(int max)
{	
	uniform_int<> orig(0, max - 1);

	return orig(rng);
}



struct negtrue
{
	const int value;
	negtrue(bool neg, int value) : value(neg ? -1 : (value >= 0 ? value : 0))
	{	
	}

	operator bool() const
	{
		return value < 0;
	}
};

// Declarations used for type safety based on Boost Units
class FactorUnit {};

struct MarkerBaseDimension : base_dimension<MarkerBaseDimension, 931> {};
typedef MarkerBaseDimension::dimension_type MarkerDimension;
struct MarkerValBaseUnit : public base_unit<MarkerValBaseUnit, MarkerDimension, 932> {};
typedef make_system<MarkerValBaseUnit>::type UnitSystem;

typedef unit<MarkerDimension, UnitSystem> MarkerValUnit;

BOOST_UNITS_STATIC_CONSTANT(MarkerValue,MarkerValUnit);  


typedef quantity<MarkerValUnit, int> MarkerVal;
typedef quantity<FactorUnit, float> Factor;

typedef pair<MarkerVal, MarkerVal> MarkerValPair;

const MarkerVal UnknownMarkerVal = (MarkerVal) 0;
const MarkerVal sexmarkerval = 9 * MarkerValue;

const float maxdiff = 0.00005;

#include "settings.h"

bool early = false;

vector<double> markerposes;
vector<double> actrec[2];
//int selfgen = 0;
double genrec[3];
vector<unsigned int> chromstarts;

vector<int> markertranslation;
typedef vector<boost::array<double, 5 > > MWTYPE;

template<class T> class vectorplus;

template<> class vectorplus<MWTYPE >
{
public:
	MWTYPE operator () (const MWTYPE& l, const MWTYPE& r)
	{
		MWTYPE result;
		result.resize(l.size());

		for (int i = 0; i < l.size(); i++)
		{
			for (int j = 0; j < 5; j++)
			{
				result[i][j] = l[i][j] + r[i][j];
			}
		}

		return result;
	}
};

template<class T> class vectorplus<vector<T> >
{
public:
	vector<T> operator () (const vector<T>& l, const vector<T>& r)
	{
		vector<T> result;
		result.resize(l.size());

		for (int i = 0; i < l.size(); i++)
		{
			result[i] = l[i] + r[i];
		}

		return result;
	}
};

/*namespace boost { namespace mpi {

template<class T>
struct is_commutative<vectorplus<vector<T> >, vector<T> > 
: mpl::true_ { };

} }*/

MWTYPE markerweight;
MWTYPE markerweight2;

double discstep = 0.1;
double baserec[2];
int sexc = 1;


// Is a specific marker value a admissible as a match to marker b
// The logic we use currently represents "0" as unknown, anything else as a known
// marker value.
// NOTE, VALUE CHANGED FOR ZEROS!
template<bool zeropropagate> bool markermiss(MarkerVal& a, const MarkerVal b)
{
	// This is the real logic; we do not bind anything at all when zeropropagate is true
	if (zeropropagate) return false;

	if (a == UnknownMarkerVal)
	{
		if (!zeropropagate) a = b;
		return false;
	}
	if (b == UnknownMarkerVal && a != sexmarkerval) return false;

	return a != b;
}


// Move a tree-based binary flag up a generation. The structure of bit flags might look like
// (1)(3)(3), where each group of (3) is in itself (1)(2), where (2) is of course (1)(1).
int upflagit(int flag, int parnum, int genwidth)
{

	if (flag < 0) return flag;
	flag >>= parnum * (genwidth - 1);
	flag &= ((1 << (genwidth - 1)) - 1);

	return flag;
}

struct individ;

int generation = 1;
int shiftflagmode;
int lockpos[NUMSHIFTS];
int quickmark[NUMSHIFTS];
int quickgen[NUMSHIFTS];

template<class T2> class PerStateArray
{
public:
	typedef boost::array<T2, NUMTYPES> T;
};

template<class T2> class StateToStateMatrix
{
public:
	typedef typename PerStateArray<
		typename PerStateArray<T2>::T >
		::T T;
};

// The quick prefixes are caches that retain the last invocation.
// Only used fully when we stop at marker positions exactly, i.e. not a general grid search.
double quickfactor[NUMSHIFTS];
PerStateArray<int>::T quickendmarker[NUMSHIFTS];
PerStateArray<double>::T quickendfactor[NUMSHIFTS];
StateToStateMatrix<double>::T quickendprobs[NUMSHIFTS];
PerStateArray<double>::T quickmem[NUMSHIFTS];

// A hashed store of inheritance pathway branches that are known to be impossible.
// Since we can track the two branches that make up the state in the F_2 individual independently,
// this optimization can reduce part of the cost by sqrt(number of states).
typedef boost::array<boost::array<boost::array<boost::array<boost::array<boost::array<boost::array<int, 4>, HALFNUMSHIFTS>, HALFNUMPATHS + 1>, HALFNUMTYPES>, 2>, 2>, 2> IAT;
IAT impossible;

// A memory structure storing haplo information for later update.
// By keeping essentially thread-independent copies, no critical sections have to
// be acquired during the updates.
boost::array<boost::array<float, 2>, 1000000> haplos;
boost::array<boost::array<map<MarkerVal, float>, 2>, 1000000> infprobs;

// done, factors and cacheprobs all keep track of the same data
// done indicates that a specific index (in the binary tree of blocks of multi-step transitions) is done
// with a "generation id" that's semi-unique, meaning no active clearing of the data structure is performed
vector<int> done[NUMSHIFTS];
// factors contain the mantissas of the extended floating-point representation
vector<PerStateArray<double>::T > factors[NUMSHIFTS];
// cacheprobs contain actual transitions from every possible state to every possible other state
vector<StateToStateMatrix<double>::T > cacheprobs[NUMSHIFTS];
vector<individ*> reltree;
map<individ*, int> relmap; //containing flag2 indices

//#pragma omp threadprivate(realdone, realfactors, realcacheprobs)


#pragma omp threadprivate(generation, done, factors, cacheprobs, shiftflagmode, impossible, haplos, lockpos, quickmark, quickgen, quickmem, quickfactor, quickendfactor, quickendprobs, reltree, relmap, infprobs)

// We put all thread-local structures in our own separate struct. This is because many compilers implementing OpenMP
// use a relatively expensive call for determining thread-local global data, which is not cached over function calls.
// The simple pointer-arithmetic for a lookup within a struct is cheaper.
struct threadblock
{
	int* const generation;
	int* const shiftflagmode;
	int* const quickmark;
	int* const quickgen;
	int* const lockpos;
	double* const quickfactor;
	PerStateArray<double>::T* const quickmem;
	PerStateArray<int>::T* const quickendmarker;
	IAT* const impossible;
	boost::array<boost::array<float, 2>, 1000000>* const haplos;
	vector<int>* const done;
	vector<PerStateArray<double>::T >* const factors;
	vector<StateToStateMatrix<double>::T >* const cacheprobs;	
	PerStateArray<double>::T* const quickendfactor;
	StateToStateMatrix<double>::T* const quickendprobs;
	boost::array<boost::array<map<MarkerVal, float>, 2>, 1000000>* infprobs;

	threadblock() : generation(&::generation), shiftflagmode(&::shiftflagmode), impossible(&::impossible),
		done(::done), factors(::factors), cacheprobs(::cacheprobs), haplos(&::haplos),
		quickmark(::quickmark), quickgen(::quickgen), lockpos(::lockpos), quickmem(::quickmem),
		quickfactor(::quickfactor), quickendfactor(::quickendfactor), quickendprobs(::quickendprobs),
		quickendmarker(::quickendmarker), infprobs(&::infprobs)
	{
	};
};


// Turners are used as mix-ins when probabilities are filtered at the "fixed" marker
// (the one we are actually probing).
// This functionality is used to invert the state in different manner, to test for the
// possibility that the complete haplotype assignment from an arbitrary point until the end
// has been mixed up. This can easily happen as the assignment of haplotype numbers is
// basically arbitrary, so the only thing keeping it in place is the linkage to the original
// defining locus.
class noneturner
{
public:
	void operator() (const PerStateArray<double>::T& probs) const
	{
	};

	bool canquickend() const
	{
		return true;
	}
} none;

class aroundturner
{
	int turn;
	int flagmodeshift;

public:
	aroundturner(int turn)
	{
		if (NUMGEN == 3)
		{
			this->turn = turn & 54;
			flagmodeshift = (turn >> TYPEBITS) | ((turn & 1) ? 2 : 0) |
				((turn & 8) ? 4 : 0);
		}
		else
		{
			this->turn = turn & 3;
			flagmodeshift = (turn >> TYPEBITS);
		}
	}

	bool canquickend() const
	{
		return false;
	}


	void operator() (PerStateArray<double>::T& probs) const
	{
		PerStateArray<double>::T probs2;

		for (int i = 0; i < NUMTYPES; i++)
		{
		  int newval = i ^ turn;
		  if (flagmodeshift)
		    {
		      // Emulating actually switching parental genotypes around
		      //		      newval = ((newval & 1) << 1) | ((newval & 2) >> 1);
		    }
			probs2[newval] = probs[i];
		}

		for (int i = 0; i < NUMTYPES; i++)
		{
			probs[i] = probs2[i];
		}

		// Not *tb out of laziness, the number of calls is limited
		shiftflagmode ^= flagmodeshift;
		// Hacky, tacky
	};
};

// A struct containing some auxiliary arguments to the trackpossible family of functions.
// These parameters are not known at compile-time, hence not template arguments, but they do not
// change over the series of recursive calls.
const struct trackpossibleparams
{
	const float updateval;
	int* const gstr;

	trackpossibleparams() : updateval(0.0f), gstr(0)
	{
	}

	trackpossibleparams(float updateval, int* gstr) : updateval(updateval), gstr(gstr)
	{
	}
} tpdefault;



struct classicstop
{
	int lockpos;
	int genotype;

	classicstop(int lockpos, int genotype) : lockpos(lockpos), genotype(genotype)
	{}

	operator int() const
	{
		return lockpos;
	}


	const int getgenotype(int startmark) const
	{
		return genotype;
	}



	negtrue okstep(const int startmark, const int endmark) const
	{
		bool found = (lockpos <= -1000 - startmark && lockpos > -1000 - endmark) ||
			(lockpos >= markerposes[startmark] && lockpos <= markerposes[endmark]);

		return negtrue(!found, genotype);
	}

	bool fixtofind(int& genotype, double& startpos, double& endpos, int j) const
	{
		bool tofind = (lockpos >= startpos && lockpos <= endpos);
		if (tofind)
		{
			endpos = lockpos;
		}


		if (lockpos == -1000 - j + 1)
		{
			tofind = true;
			endpos = startpos;
		}

		if (tofind)
		{
			genotype = this->genotype;
		}

		return tofind;
	}

};


struct smnonecommon
{
	const int getgenotype(int startmark) const
	{
		return -1;
	}
};

struct nonestop : smnonecommon
{
	operator int() const
	{
		assert(false);

		return -1;
	}

	negtrue okstep(const int startmark, const int endmark) const
	{
		return negtrue(true, 0);
	}

	bool fixtofind(int& genotype, double& startpos, double& endpos, int j) const
	{
		return false;
	}
} NONESTOP;

struct tssmcommon
{
	int lockpos;
	int genotype[2];


	operator int() const
	{
		return lockpos;
	}

	negtrue okstep(const int startmark, const int endmark) const
	{
		// twicestop is only used for markers
		bool found = (lockpos <= -1000 - startmark + 1 && lockpos > -1000 - endmark);

		if (found)
		{
			int index = -lockpos - startmark - 1000;

			if (index == 0 || index == 1) return negtrue(!found, genotype[index]);
		}

		return negtrue(!found, 0);
	}

	bool fixtofind(int& genotype, double& startpos, double& endpos, int j) const
	{
		bool two = lockpos == -1000 - j + 2;
		if (lockpos == -1000 - j + 1 || two)
		{
			endpos = startpos;
			genotype = this->genotype[two];
			return true;
		}

		return false;
	}

protected: tssmcommon(int lockpos) : lockpos(lockpos)
		   {
		   }
};

struct twicestop : tssmcommon
{
	twicestop(int lockpos, int* genotype) : tssmcommon(lockpos)
	{
		this->genotype[0] = genotype[0];
		this->genotype[1] = genotype[1];
	}

	twicestop(int lockpos, int g1, int g2) : tssmcommon(lockpos)
	{
		this->genotype[0] = g1;
		this->genotype[1] = g2;
	}

	const int getgenotype(int startmark) const
	{
		return genotype[-1000 - startmark != lockpos];
	}

};

struct stopmodpair : tssmcommon, smnonecommon
{
	typedef boost::array<boost::array<float, 2>, 2> miniactrecT;
	miniactrecT actrec;

	stopmodpair(int lockpos, miniactrecT actrec) : tssmcommon(lockpos), actrec(actrec)
	{
		genotype[0] = -1;
		genotype[1] = -1;
	}
};

template<class G> const bool canquickend(int startmark, const G &stopdata)
{
	return false;
}

/*template<> const bool canquickend<smnonecommon>(int startmark, const smnonecommon &stopdata)
{
return false;
}*/

template<> const bool canquickend<classicstop>(int startmark, const classicstop &stopdata)
{
	return stopdata.lockpos <= -1000 && stopdata.genotype >= 0;
}

template<> const bool canquickend<twicestop>(int startmark, const twicestop &stopdata)
{
	return stopdata.lockpos == -1000 - startmark + 1;
}


// Should only be called for markers where locking is really done



template<class G> bool firstlockmatch(int lockpos, const G& stopdata)
{
	return stopdata.lockpos == lockpos;
}


template<> bool firstlockmatch<nonestop>(const int lockpos, const nonestop& stopdata)
{
	return false;
}

#ifdef PERMARKERACTREC
template<class G> double getactrec(const G& stopdata, double startpos, double endpos, int k, int j, int gen)
{
	return actrec[k][j];
}
#else
template<class G> double getactrec(const G& stopdata, double startpos, double endpos, int k, int j, int gen)
{
	return genrec[gen];
}
#endif

template<> double getactrec<stopmodpair>(const stopmodpair& stopdata, double startpos, double endpos, int k, int j, int gen)
{
	int index = j - 1 + stopdata.lockpos + 1000;
	//  printf("getactrec: %d %d %d %d\n", k, j, index, stopdata.lockpos);

	if (index == 0 || index == 1) return stopdata.actrec[k][index];

	return actrec[k][j];
}

const int semicount = 3226;


// A structure containing most of the information on an individual
struct individ
{
	// The individual #.
	int n;
	string name;
	// Generation number. Convenient, while not strictly needed.
	int gen;
	// Number of children.
	int children;
	// Am I explicitly created to not contain genotypes?
	bool empty;
	// Parents.
	individ* pars[2];
	//
	vector<individ*> kids;
	// Sex.
	bool sex;
	// Line or strain of origin, should only exist in founders.
	int strain;
	// Marker data as a list of pairs. No specific ordering assumed.
	vector<MarkerValPair > markerdata;
	vector<pair<float, float> > markersure;

	// Temporary storage of all possible marker values, used in fixparents.
  vector<map<MarkerVal, pair<int, double> > > markervals;
	// The haplotype weight, or skewness. Introducing an actual ordering of the value in markerdata.
	vector<float> haploweight;
	// The cost-benefit value of inverting the haplotype assignment from an arbitrary marker point on.
	vector<float> negshift;
	vector<int> lastinved;
	vector<unsigned int> lockstart;
	//vector<boost::array<boost::array<double, 40>, 4 > > semishift;
	vector<boost::array<boost::array<boost::array<boost::array<double, 2>, 2>, 2>, 2> > parinfprobs;
	vector<boost::array<map<pair<MarkerVal, MarkerVal>, double>, 2> > infprobs;
	vector<boost::array<map<pair<MarkerVal, MarkerVal>, double>, 2> > sureinfprobs;
	vector<boost::array<double, 2> > unknowninfprobs;

	vector<int> genotypegrid;

	// Accumulators for haplotype skewness updates.
	vector<double> haplobase;
	vector<double> haplocount;

	individ()
	{
		pars[0] = 0;
		pars[1] = 0;
		strain = 0;
		sex = false;
		empty = false;
	}

	bool arerelated(individ* b, vector<individ*> stack = vector<individ*>(), int gens = 0)
	{
		if (gens > 2) return false;
		if (!this) return false;

		if (b == this) return true;
		if (find(stack.begin(), stack.end(), this) != stack.end())
		{
			return false;
		}

		stack.push_back(this);
		if (stack.size() == 1)
		{
			if (b->arerelated(this, stack, gens + 1)) return true;
		}

		if (pars[0]->arerelated(b, stack, gens + 1)) return true;
		if (pars[1]->arerelated(b, stack, gens + 1)) return true;

		for (int i = 0; i < kids.size(); i++)
		{
			for (int j = 0; j < b->kids.size(); j++)
			{
				if (b->kids[j] == kids[i]) return true;
			}
		}
		//stack.pop_back();

		return false;
	}


	// A wrapper for calling trackpossible. Additional logic introduced to do lookup in the "impossible" tables of branches
	// not allowed at the current point. This will generally mean most branches and reduces the number of states visited considerably
	// in typical cases.
	//
	// The double overload means that the class can be treated as a conventional function call, returning a double, when the pre-lookup
	// is not needed. In the "real" trackpossible method, one instance is called in that way, while the other one is pre-looked up.
	template<bool update, bool zeropropagate> struct recursetrackpossible
	{
		individ* mother;
		int upflagr;
		int upflag2r;
		int upshiftr;
	  double secondval;
		const trackpossibleparams& extparams;
		const int genwidth;
		const MarkerVal markerval;
		const int marker;
		const threadblock& tb;
		int firstpar;
		int* impossibleref;
		int impossibleval;
		bool prelok;

	  recursetrackpossible(individ* mother, const threadblock& tb, MarkerVal markerval, double secondval, int marker, int upflag, int upflag2, int upshift, int genwidth, int f2n, int firstpar, int numrealpar, const trackpossibleparams& extparams) :
	    genwidth(genwidth), extparams(extparams), markerval(markerval), marker(marker), tb(tb), firstpar(firstpar), mother(mother), secondval(secondval)

		{
			upflagr = upflagit(upflag, firstpar, genwidth);
			upflag2r = upflagit(upflag2, firstpar, genwidth >> (NUMGEN - NUMFLAG2GEN));
			upshiftr = upflagit(upshift, firstpar, genwidth >> (NUMGEN - NUMSHIFTGEN));

			prelok = true;
			if (!zeropropagate && genwidth == (1 << (NUMGEN - 1)))
			{
			  if (false)
			    {
				impossibleref = &(*tb.impossible)[*(tb.shiftflagmode) & 1][firstpar][f2n][upflagr][upflag2r + 1][upshiftr][marker & 3];
			    }
				impossibleval = (*tb.generation) * markerposes.size() + marker;

				if (impossibleref && *impossibleref == impossibleval)
				{
					prelok = false;
				}
			}
		}

		operator double()
		{
			if (!prelok)
			{
			  return 0;
			}

			double baseval =
			  mother->pars[firstpar]->trackpossible<update, zeropropagate>(tb, markerval, secondval, marker,
				upflagr,
				upflag2r,
				upshiftr, extparams, genwidth >> 1);

			if (impossibleref && !zeropropagate && !update && genwidth == (1 << (NUMGEN - 1)) && !baseval)
			{
			  *impossibleref = impossibleval;
			}

			return baseval;
		}
	};

	// zeropropagate also implies that there is a gstr value, we propagate zeros to find any possible source strain
	// The main logic of tracking a specific inheritance pathway, computing haplotype weights, and overall feasibility.
	// update: Should haplotype weights be updated?
	// zeropropagate: Should zero marker values be kept, or changed to the actual values they are matched against.
	// threadblock: Reference to all thread-local data.
	// inmarkerval: The marker val we are matching against.
	// marker: The marker number.
	// flag: The actual genotype flag. Note that this is zero-extended to be one bit more than the actual enumeration of flags
	// (lowest bit always 0).
	// flag99: What's mostly called flag2. A specific shift state of marker numbers, or -1 to indicate that all shifts are allowed.
	// localshift: Based on shiftflagmode, the overall mapping *for the complete sequence* of strand numbers to actual parentage for all
	// individuals in the analysis.
	// extparams: external parameters.
	// genwidth: The width of the generation flags.
  template<bool update, bool zeropropagate> double trackpossible(const threadblock& tb, MarkerVal inmarkerval, double secondval, const unsigned int marker,
		const unsigned int flag, const int flag99, int localshift = 0, const trackpossibleparams& extparams = tpdefault,
		const int genwidth = 1 << (NUMGEN - 1)) /*const*/
	{
		if (this == NULL) return 1;

		// TYPEBITS are the ordinary TYPEBITS. Anything set beyond those indicates selfing. g is multiplied by 2 to become flag, hence TYPEBITS + 1
		const int selfval = (flag >> (TYPEBITS + 1));
		const bool selfingNOW = SELFING && (genwidth == (1 << (NUMGEN - 1))) && selfval;

		MarkerVal selfmarker[2];
		float selfsure[2];

		int upflag2 = -1;
		const int upflag = flag >> 1;
		const int upshift = localshift >> 1;
		int f2s = 0;
		const MarkerVal* themarker = selfingNOW ? selfmarker : &markerdata[marker].first;
		const float* themarkersure = selfingNOW ? selfsure : &markersure[marker].first;
		int f2end = 2;

		if (flag99 != -1 && genwidth >> (NUMGEN - NUMFLAG2GEN) > 0)
		{
			upflag2 = flag99 >> 1;
			f2s = flag99;
			f2end = flag99 + 1;			
		}
		//if (this->n == 1633) printf("%d %d %d %d %d\t%d %d %d\n", this->n, localshift, genwidth, NUMGEN - NUMFLAG2GEN, flag99, upflag2, upflag, upshift);

		int firstpar = flag & 1;
		double ok = 0;

		//#pragma ivdep
		// This ivdep is quite complicated, since we actually change markerval, but not in
		// conditions where this alters the results of the other iteration.

		// flag2 determines which value in the tuple to check. flag and localshift determine the haplotype weight value
		// assigned to that test, and which parent to match against that value.
		for (int flag2 = f2s; flag2 < f2end && (HAPLOTYPING || !ok); flag2++)
		{
		  int f2n = (flag2 & 1);
		  if (selfingNOW)
		    {
		      int selfindex = (selfval >> 1) ^ f2n;
		      selfmarker[0] = UnknownMarkerVal;
		      selfmarker[1] = UnknownMarkerVal;
		      selfsure[0] = 0;
		      selfsure[1] = 0;
		      MarkerVal first = markerdata[marker].first;
		      
		      if (!markermiss<false>(first, markerdata[marker].second))
			{
			  selfmarker[selfindex] = first;
			  selfsure[selfindex] = 1 - (1 - markersure[marker].first) * (1 - markersure[marker].second);
			}
		      else
			{
			  selfmarker[selfindex] = markerdata[marker].second;
			  if (markersure[marker].first == 0) { return 0; }
			  selfsure[selfindex] = 1 - (markersure[marker].first) * (1 - markersure[marker].second);
			}
		    }


		  bool allthesame = themarker[0] == themarker[1];
			MarkerVal markerval = inmarkerval;
			double baseval;

			
			int realf2n = f2n;

			double mainsecondval = 0;
			// If this marker value is not compatible, there is no point in trying.
			if (markermiss<zeropropagate>(markerval, themarker[f2n]))
			{
			  baseval = themarkersure[f2n] + (1.0 - themarkersure[f2n]) * secondval;
			}
			else
			{
			  if (themarkersure[f2n]) mainsecondval = (themarkersure[f2n]) / (1.0 - themarkersure[f2n]);
			  baseval = 1.0 - themarkersure[f2n] + themarkersure[f2n] * secondval;
			}

			if (!baseval) continue;			

			// Normalize, in some sense.
			f2n ^= ((firstpar ^ localshift) & 1);

			if (zeropropagate || !genwidth)
			{
				baseval *= 0.5;
			}
			//			else if (/*!empty &&*/ (allthesame && (CORRECTIONINFERENCE) || (themarker[0] == UnknownMarkerVal && themarker[1] == UnknownMarkerVal && themarkersure[0] + themarkersure[1] == 0)))
						else if (/*!empty &&*/ (allthesame && ((CORRECTIONINFERENCE) || (themarkersure[0] == themarkersure[1]))) && !selfingNOW)
			{
				baseval *= ((f2n) ? 1.0 : 0.0);
				//if (baseval == 0.5) printf("%d\n", n);
			}
			else
			{
				if (HAPLOTYPING)
				{
					//if (selfingNOW) baseval *= realf2n == selfindex; else
					baseval *= fabs((f2n ? 1.0 : 0.0) - (selfingNOW  ? 0: haploweight[marker]));
					//if (selfingNOW) baseval *= 0;

					if (themarker[realf2n] == UnknownMarkerVal && selfingNOW && baseval/* && themarker[!realf2n] == markerval*/)
					{						
						//						baseval *= 0.95;
						/*if (baseval > 0.75) baseval = 0.75;*/

						//baseval = 0.5; // This is wicked
					}
				}
				else
				{
					// No haplotype weights, all interpretations allowed.
					baseval *= 0.5;
				}
			}



			if (!baseval)
			{
				continue;
			}

			// If we are at maximum depth, by depth limit or by lack of ancestor
			if (genwidth == HAPLOTYPING || !pars[firstpar])
			{
				if (zeropropagate && extparams.gstr)
				{
					*(extparams.gstr) += (themarker[f2n] == (2 * MarkerValue));
				}
			}

			// There should be some other flag for the actual search depth
			if (genwidth == HAPLOTYPING)
			{
			}
			else
			{
				// Track to the parent generation, creating evaluation objects first
				// These do a lookup in a special hash for combinations known to be 0, to avoid unnecessary calls
				// Both are checked before either branch of the pedigree is actually traced.
				recursetrackpossible<update, zeropropagate> subtrack1 =
				  recursetrackpossible<update, zeropropagate>(this, tb, markerval, mainsecondval, marker,
					upflag,
					upflag2,
					upshift,
					genwidth,
					f2n,
					firstpar,
					0,
					extparams);

				if (subtrack1.prelok && (!zeropropagate || (genwidth == 1 << (NUMGEN - 1))) )
				  {					  
				    double secsecondval = 0;
				    if (themarkersure[!realf2n])
				      {
					baseval *= (1 - themarkersure[!realf2n]);
					secsecondval = themarkersure[!realf2n] / (1 - themarkersure[!realf2n]);
				      }
				    baseval *= recursetrackpossible<update, zeropropagate>(this, tb, themarker[!realf2n], secsecondval,
											   marker,
					upflag,
					upflag2,
					upshift,
					genwidth,
					f2n,
					!firstpar,
					1,
					extparams);

					if (selfingNOW && extparams.gstr && !(selfval & (1 << !firstpar))) *extparams.gstr = 0;
				  }
				if (!baseval) continue;

				if (!(selfingNOW && extparams.gstr && !(selfval & (1 << firstpar))))
					baseval *= subtrack1;
			}

			if (!baseval) continue;

			ok += baseval;

			if (update /*&& !allthesame*/ && !selfingNOW)
			{
				(*tb.haplos)[n][f2n] += extparams.updateval;
				(*tb.infprobs)[n][realf2n][markerval] += extparams.updateval;
			}
		}
		if (selfingNOW && extparams.gstr) *extparams.gstr *= 2;
		return ok;
	}


	// calltrackpossible is a slight wrapper that hides at least some of the interanl parameters needed for the recursion from outside callers
	template<bool update, bool zeropropagate> double calltrackpossible(const threadblock& tb, const MarkerVal* markervals, const unsigned int marker,
		const int genotype, const unsigned int offset, const int flag2, const double updateval = 0.0)
	{		
	  return trackpossible<update, zeropropagate>(tb, UnknownMarkerVal, 0,
			marker, genotype * 2, flag2, *(tb.shiftflagmode), trackpossibleparams(updateval, 0));
	}

	// "Fix" parents, i.e. infer correct marker values for any zero values existing.
	// Depending on some specifics, this can either assume that all existing marker values will be reflected in some offspring individuals,
	// or indeed only infer specifically those things that can be said from certain, i.e. if a zero value exists, and the other branch of
	// the pedigree doesn't allow a marker value found in the offspring, then it can be assumed to have been transmitted from the individual
	// with the zero value. 
	void fixparents(unsigned int marker, bool latephase)
	{
		MarkerValPair& themarker = markerdata[marker];
		individ* parp[3] = {pars[0], pars[1], this};
		int az = 0;

#pragma omp critical (parmarkerval)
		for (int i = 0; i < 3; i++)
		{
			if (parp[i])
			{
				az += parp[i]->markerdata[marker].first == UnknownMarkerVal;
				az += parp[i]->markerdata[marker].second == UnknownMarkerVal;
				while (parp[i]->markervals.size() <= marker)
				{
					parp[i]->markervals.resize(markerdata.size());
				}
			}
		}

		//if (!az) return;		

		double okvals[2] = {0};
		// We only accept an interpretation when it is by exclusion the only possible one. As soon as one intepretation has gained acceptance,
		// no need to retest it.
		for (shiftflagmode = 0; shiftflagmode < 1; shiftflagmode+=1)
		{
			threadblock tb;
			for (int i = 0; i < NUMTYPES; i++)
			{
				for (int flag2 = 0; flag2 < NUMPATHS; flag2++)
				{
					if (okvals[flag2 & 1]) continue;

					double ok = calltrackpossible<false, false>(tb, 0, marker, i, 0, flag2);
					if (!ok) continue;
					okvals[flag2 & 1] += ok;
					if (okvals[0] && okvals[1]) break;
				}
				if (okvals[0] && okvals[1]) break;
			}
		}

		if (!okvals[0] && !okvals[1])
		{
			printf("Clearing %d:%d\n", this->n, marker);
			markerdata[marker] = make_pair(UnknownMarkerVal, UnknownMarkerVal);
			markersure[marker] = make_pair(0.0, 0.0);
		}

		if ((((bool) okvals[0]) ^ ((bool) okvals[1])) || latephase)
		{
			for (int flag2 = 0; flag2 < 2; flag2++)
			{
				if (!okvals[flag2]) continue;

				for (int k = 0; k < 2; k++)
				{
					if (pars[k])
					{
						int u = ((k ^ flag2 /*^ *tb.shiftflagmode*/) & 1);
						if (latephase || (&themarker.first)[u] != UnknownMarkerVal)
#pragma omp critical (parmarkerval)
						{
						  double old1 = 1;
						  int old2 = 0;
						  map<MarkerVal, pair<int, double> >::iterator i = pars[k]->markervals[marker].find((&themarker.first)[u]);
						  if (i != pars[k]->markervals[marker].end())
						    {
						      old1 = i->second.second;
						      old2 = i->second.first;
						    }
						  double probit = (markersure[marker].first + markersure[marker].second);
						  probit /= (1 - probit);
						  pars[k]->markervals[marker][(&themarker.first)[u]] = make_pair(old2 + 1, old1 * probit);
						}
					}
				}
			}
		}		
	}

	// Verifies that an update should take place and performs it. Just a wrapper to calltrackpossible.
	void updatehaplo(const threadblock& tb, const unsigned int marker, const unsigned int i, const int flag2, const double updateval)
	{
		MarkerValPair& themarker = markerdata[marker];		

		double ok = calltrackpossible<false, false>(tb, &themarker.first, marker, i, 0, flag2, updateval);

		if (ok)
		{
			calltrackpossible<true, false>(tb, &themarker.first, marker, i, 0, flag2, updateval);
		}
		else
		{

		}
	}

	// Adjust the probability, i.e. filter all probability values based on the haplotype weights and overall admissibility for the different
	// states.
	void adjustprobs(const threadblock& tb, PerStateArray<double>::T& probs, const unsigned int marker, double& factor, const bool oldruleout, int flag99)
	{
		double sum = 0;
		PerStateArray<double>::T probs2;

		const MarkerValPair& themarker = markerdata[marker];
		const bool ruleout = true; // TODO

		for (int q = 0; q <= (int) !ruleout; q++)
		{
			int f2start = 0;
			int f2end = NUMPATHS;

			// Negative values other than -1 just cause trouble

			// Can we really trust the logic to expand correctly even for zero values?
			//if (flag99 >= 0 || !somezero[marker])
			{
				f2start = flag99;
				f2end = flag99 + 1;
			}

			for (unsigned int i = 0; i < NUMTYPES; i++) // genotype
			{
				probs2[i] = probs[i];

				// We will multiply this already small number with an even smaller number... let's assume it's zero and be done with it.
				if (probs[i] < 1e-200)
				{
					probs[i] = 0;
					continue;
				}

				double realok = 0;

				for (int flag2 = f2start; flag2 < f2end; flag2++)
				{
					realok += calltrackpossible<false, false>(tb, &themarker.first, marker, i, 0, flag2);
				}					

				// TODO UGLY CAPPING
				if (HAPLOTYPING)
				{
					probs[i] *= (double) realok;
				}
				else
				{
					probs[i] *= (bool) realok;
				}
				sum += probs[i];
			}

			if (sum == 0 && !ruleout)
			{
				for (int i = 0; i < NUMTYPES; i++)
				{
					probs[i] = probs2[i];
				}

				// The code sees: This doesn't make sense, ignore this marker!
				// NOTE: If the parent-grandparent genotypes are completely incompatible, the problem
				// is NOT eliminated.

				// Current caching scheme does not execute full elimination. Efficiency gains possible by not allowing
				// for genotyping error in every single marker. A small epsilon should otherwise be introduced whenever
				// any probs[i] is zero.
				markerdata[marker] = make_pair(UnknownMarkerVal, UnknownMarkerVal);
				fprintf(stderr, "Error in %x, marker %d, impossible path\n", this, marker);
				sum = 1;
			}
			else
				break;
		}

		if (sum <= 0)
		  {
		    factor = MINFACTOR;
		  }
		else
		  {
		    // Normalize, update auxiliary exponent
		    for (int i = 0; i < NUMTYPES; i++)
		      {
			probs[i] /= sum;
		      }
		    factor += log(sum);
		  }
	}
		
	// Append a "multi-step" transition. If the cached values (essentially a N * N transition matrix for the steps from startmark to
	// endmark) are missing, calculate them first.
	double fillortake(const threadblock& tb, const int index, const unsigned int startmark, const unsigned int endmark, PerStateArray<double>::T& probs)
	{
		if ((tb.done[*(tb.shiftflagmode)])[index] != (*tb.generation))
		{

			for (unsigned int i = 0; i < NUMTYPES; i++)
			{
				PerStateArray<double>::T probs2 = {{0}};
				probs2[i] = 1;

				// Note ruleout here, could be set to false if we "prime" with an evenly distributed probs first
				(tb.factors[*tb.shiftflagmode])[index][i] = quickanalyze<false, noneturner>(tb, none, startmark,
					endmark,
					NONESTOP,
					-1,
					true,
					probs2);

				double sum = 0;

#pragma ivdep:back
				for (int j = 0; j < NUMTYPES; j++)
				{
					if (!_finite(probs2[j])) probs2[j] = 0.0;
					sum += probs2[j];
				}

				sum *= NUMTYPES;
				if (sum <= 0)
				  {
				    (tb.factors[*tb.shiftflagmode])[index][i] = MINFACTOR;
				    sum = 0;
				  }
				else
				  {
				    (tb.factors[*tb.shiftflagmode])[index][i] += log(sum);
				    sum = 1/sum;
				  }
				double* probdest = &(tb.cacheprobs[*tb.shiftflagmode])[index][i][0];

#pragma ivdep
				for (int j = 0; j < NUMTYPES; j++)
				{
					probdest[j] = probs2[j] * sum;
				}
			}
			(tb.done[*tb.shiftflagmode])[index] = (*tb.generation);
		}

	        double factor = MINFACTOR;
		for (int i = 0; i < NUMTYPES; i++)
		{
		  if (probs[i] == 0.0) continue;

			factor = max(factor, (tb.factors[*tb.shiftflagmode])[index][i]);
		}

		PerStateArray<double>::T probs2 = {{0}};

		for (int i = 0; i < NUMTYPES; i++)
		{
			double step = (tb.factors[*tb.shiftflagmode])[index][i] - factor;
			if (probs[i] == 0.0 || step <= -100.0f) continue;
			double basef = exp((double) step) * probs[i];
			//			if (basef == 0.0 || !_finite(basef)) continue;
			const double* probsource = &(tb.cacheprobs[*tb.shiftflagmode])[index][i][0];
#pragma ivdep
			for (int j = 0; j < NUMTYPES; j++)
			{
				const double term = basef * probsource[j];
				probs2[j] += term;
			}
		}

		double sum = 0;

		for (int i = 0; i < NUMTYPES; i++)
		{
			sum += probs2[i];
		}

		if (sum <= 0)
		  {
		    factor = MINFACTOR;
		    sum = 0;
		  }
		else
		  {
		    factor += log(sum);
		    sum = 1 / sum;
		  }
		    
		for (int i = 0; i < NUMTYPES; i++)
		  {
		    probs[i] = probs2[i] * sum;
		  }


		return factor;
	}

	// Analyze for a specific range, including a possible fixed specific state at some position (determined by stopdata)
	template<bool inclusive, class T, class G> double quickanalyze(const threadblock& tb, const T& turner, unsigned int startmark,
		const unsigned int endmark, const G &stopdata, const int flag2, bool ruleout, PerStateArray<double>::T& probs,
		float minfactor = MINFACTOR)
	{
		unsigned int stepsize;
		double factor = 0;
		bool allowfull = inclusive;
		bool frommem = false;
		if (inclusive && tb.quickgen[*tb.shiftflagmode] == *tb.generation && firstlockmatch(tb.lockpos[*tb.shiftflagmode], stopdata))
		{
			allowfull = true;
			frommem = true;
			factor = tb.quickfactor[*tb.shiftflagmode];
			if (factor <= minfactor)
			{
				return MINFACTOR;
			}

			startmark = tb.quickmark[*tb.shiftflagmode];
//			printf("STARTING FROM MARKER: %d\n", startmark);
			for (int i = 0; i < NUMTYPES; i++)
			{
				probs[i] = tb.quickmem[*tb.shiftflagmode][i];
			}
		}

		// Loop until we have reached the end, with blocks of varying sizes.
		while (startmark < endmark)
		{
			for (stepsize = 1; stepsize < (endmark - startmark + allowfull) &&
				stopdata.okstep(startmark, startmark + stepsize) &&
				!(startmark & (stepsize - 1)); stepsize *= 2);

			// A single step, either due to the fixated genotypee being within this range, or simply because we've gone all the way down
			// the tree.
			if (stepsize <= 2)
			{
				stepsize = 1;

				if (!frommem && !stopdata.okstep(startmark, startmark + 1))
				{
					// If we have a fixated genotype at marker x in one call, it is very likely that the next call will also
					// be a fixated genotype at marker x. Possibly another one, but still fixated. The fixated genotype does not
					// change the probability values leading up to this position, so they are cached.
					tb.quickgen[*tb.shiftflagmode] = *tb.generation;
					tb.quickmark[*tb.shiftflagmode] = startmark;
					tb.lockpos[*tb.shiftflagmode] = stopdata;
					tb.quickfactor[*tb.shiftflagmode] = factor;

					for (int i = 0; i < NUMTYPES; i++)
					{
						tb.quickmem[*tb.shiftflagmode][i] = probs[i];
						tb.quickendmarker[*tb.shiftflagmode][i] = -1;
						tb.quickendfactor[*tb.shiftflagmode][i] = 1.0;
					}
					frommem = true;
				}
				bool willquickend = (frommem && turner.canquickend() && canquickend(startmark, stopdata));
				int genotype = stopdata.getgenotype(startmark);

				if (genotype < 0) willquickend = false;

				if (willquickend)
				{
					if (factor + tb.quickendfactor[*tb.shiftflagmode][genotype] <= minfactor)
					{
						return MINFACTOR;
					}
					// If we are doing a quick end
					factor += realanalyze<4, T>(tb, turner, startmark, startmark + stepsize, stopdata, flag2, ruleout, &probs);
				}
				else
				{
					//printf("Not a quick end %d %d\n", startmark, stopdata.getgenotype(startmark));
					factor += realanalyze<0, T>(tb, turner, startmark, startmark + stepsize, stopdata, flag2, ruleout, &probs);
				}

				// This will work.
				if (!_finite(factor) || factor <= minfactor)
				{
					//			    			    if (startmark) printf("%d;%d;%d : %lf\n", n, shiftflagmode, startmark, factor);
			  if (!_finite(factor)) printf("Non-finiteA %d %d\n", n, startmark);					
return MINFACTOR;
				}
				if (willquickend)
				{
					double wasval = probs[genotype];
					if (tb.quickendfactor[*tb.shiftflagmode][genotype] > 0.0 || tb.quickendmarker[*tb.shiftflagmode][genotype] != startmark)
					{
						// Precalc the full set of transitions from this very marker, for all states, all the way to the end
						// That way, any new analysis starting from this marker can be solved with no recomputation at all.
						// Note that only the latest marker used is stored here, but as we generally loop over all states and possibly
						// multiple flag values, it still helps.
						for (int i = 0; i < NUMTYPES; i++)
						{
							tb.quickendprobs[*tb.shiftflagmode][genotype][i] = 0;				     
						}
						tb.quickendprobs[*tb.shiftflagmode][genotype][genotype] = 1.0;				 
						double superfactor = realanalyze<0 | 2, T>(tb, turner, startmark, startmark + stepsize, NONESTOP, -1, ruleout, &tb.quickendprobs[*tb.shiftflagmode][genotype]);

//						printf("Setting!\n");
						superfactor += quickanalyze<true, T>(tb, turner, startmark + stepsize, endmark, NONESTOP, flag2, ruleout, tb.quickendprobs[*tb.shiftflagmode][genotype],
							// all uses of this precalced data will have the
							// quickfactor component in common, so taking that
							// into account for the limit is not problematic
							// any strong filtering at this specific locus can be
							// handled by the return MINFACTOR line above
							minfactor - tb.quickfactor[*tb.shiftflagmode]);

						tb.quickendmarker[*tb.shiftflagmode][genotype] = startmark;
						tb.quickendfactor[*tb.shiftflagmode][genotype] = superfactor;
					}

					factor += tb.quickendfactor[*tb.shiftflagmode][genotype];
					factor += log(wasval);
					if (factor <= minfactor) return MINFACTOR;

					for (int i = 0; i < NUMTYPES; i++)
					{
						probs[i] = tb.quickendprobs[*tb.shiftflagmode][genotype][i];
					}

					return factor;
				}
			}
			else
			{
				// Use a stored multi-step transition.
				stepsize /= 2;

				int base = 0;
				for (int q = stepsize / 2; q > 1; q /= 2)
				{
					base += markerposes.size() / q;
				}

				base += startmark / stepsize;

				factor += fillortake(tb, base, startmark, startmark + stepsize, probs);
			}
			startmark += stepsize;
			allowfull |= true;

			if (!_finite(factor) || factor <= minfactor)
			{
			  if (!_finite(factor)) printf("Non-finiteB %d %d\n", n, startmark);
				if (!frommem && !stopdata.okstep(startmark, endmark))
				{
					tb.quickgen[*tb.shiftflagmode] = *tb.generation;
					tb.lockpos[*tb.shiftflagmode] = stopdata;
					tb.quickfactor[*tb.shiftflagmode] = MINFACTOR;
				}

				return MINFACTOR;
			}
		}		

		return factor;
	}

	// A wrapper to quickanalyze, preparing the start probability vector.
	template<class T, class G> double doanalyze(const threadblock& tb, const T& turner, const int startmark, const int endmark, const G& stopdata,
		const int flag2, bool ruleout = false, PerStateArray<double>::T* realprobs = 0, float minfactor = MINFACTOR)
	{
		PerStateArray<double>::T fakeprobs;
		PerStateArray<double>::T& probs = realprobs ? *realprobs : fakeprobs;

		if (realprobs)
		{
			//probs = realprobs;
		}
		else
		{
			//probs = fakeprobs;
			double selfingfactors[4];
			if (SELFING)
			{
			  int selfgen = gen - 2;
				selfingfactors[0] = 1.0 / (1 << selfgen);
				selfingfactors[1] = (1 - selfingfactors[0]) * 0.5;
				selfingfactors[2] = (1 - selfingfactors[0]) * 0.5;
				selfingfactors[3] = 0;
			}
			for (int i = 0; i < NUMTYPES; i++)
			{
				fakeprobs[i] = EVENGEN * (SELFING ?
					selfingfactors[i >> TYPEBITS] : 1.0);
			}
		}

		double factor = quickanalyze<true, T>(tb, turner, startmark, endmark, stopdata, flag2, ruleout, probs, minfactor);
		bool small = !_finite(factor) || minfactor >= factor;

		if (!small) adjustprobs(tb, probs, endmark, factor, ruleout, -1); // TODO, the very last marker can be somewhat distorted

		return factor;
	}	

	// This is the actual analyzing code. It works with no caches, and can function independently, but is generally only used to patch in those
	// elements not found in caches by quickanalyze and fillortake.
	//
	// Both transition and emission (through adjustprobs) probabilities are handled here.
	//
	// first bit in updateend signals whether the interval is end-inclusive at endmark
	// the second bit in updateend will be 0 if the interval is end-inclusive at startmark, and 1 IF NOT
	// the third bit will cause the code to quit early, after processing the genotype and turner condition
	template<int updateend, class T, class G> double realanalyze(const threadblock& tb, const T& turner, const int startmark, const int endmark, const G& stopdata,
		const int flag2, const bool ruleout = false, PerStateArray<double>::T* realprobs = 0)
	{
		PerStateArray<double>::T fakeprobs;
		PerStateArray<double>::T& probs = realprobs ? *realprobs : fakeprobs;

		if (realprobs)
		{
			//probs = realprobs;
		}
		else
		{
			//probs = fakeprobs;
			for (int i = 0; i < NUMTYPES; i++)
			{
				fakeprobs[i] = EVENGEN;
			}
		}

		// The *logged* normalization factor
		double factor = 0;

		// Walk over all markers.
		for (int j = startmark + 1; j <= endmark; j++)
		{
			double startpos = markerposes[j - 1];
			double endpos = markerposes[j];
			int genotype = -1;

			bool tofind = stopdata.fixtofind(genotype, startpos, endpos, j);

			int f2use = -1;

			if (tofind)
			{
				f2use = flag2;
			}

			if (genotype != -2)
			{
				// If we are at the very first position, and the specific flag was set, include the emission probabilities for the previous
				// marker. Used to maximize the caching.
				if (!((updateend & 2) && (j == startmark + 1))) adjustprobs(tb, probs, j - 1, factor, ruleout, f2use);
			}
			else
			{
				// We could do some stuff here to avoid excessive adjustprobs calls
				// a -2 genotype does not only mean that all genotypes are allowed, but indeed that the marker data at this marker
				// is ignored!
			}

			// For a specific intra-marker region, we have two cases: the case of a fixated position between the two markers, and the simple case
			// of no fixated position.
			for (int iter = 0; iter <= (int) tofind; iter++)
			{
				// If iter is 1, we have currently handled the transition all the way to the fixated position. Now filter to keep only
				// a single state value positive.
				if (iter)
				{
					turner(probs);
					if (genotype >= 0)
					{
						for (int i = 0; i < NUMTYPES; i++)
						{
							probs[i] = probs[i] * (i == genotype);
						}
					}

					// Were we asked to stop at this very state, in the middle of things?
					if (updateend & 4) return factor;
				}


				double dist = (endpos - startpos);

				// Compute transitions, assuming there is any distance to transfer over.
				if (dist > 0)
				{
					PerStateArray<double>::T probs2 = {{0}};
					double recprob[2 + SELFING][2];

#pragma ivdep
					// Compute recombination probabilities for this specific distance, for the two sexes.
					// (as the sex-dependent marker distance might not be a simple transformation, the actrec
					// data comes into play).
					const int selfgen = gen - 2;
					for (int gen = 0; gen < 2 + SELFING; gen++)
					{
						for (int k = 0; k < 2; k++)
						{
						  recprob[gen][k] = 0.5 * (1.0 - exp((gen == 2 ? selfgen : 1) * getactrec(stopdata, startpos, endpos, k, j, gen) * (dist)));
							//					if (iter == tofind) recprob[k] = max(recprob[k], 1e-5);
						}
					}

					// Precompute values for presence (/lack of) a crossover event for either sex.
					double other[2][2][2];
#pragma ivdep
					for (int gen = 0; gen < 2; gen++)
					{
						for (int m = 0; m < 2; m++)
						{
							for (int k = 0; k < 2; k++)
							{
								double prob = recprob[gen][k];
								if (m) prob = 1.0 - prob;

								other[gen][m][k] = prob;
							}
						}
					}

					boost::array<double, NONSELFNUMTYPES> recombprec;
					
#pragma ivdep
					for (int index = 0; index < NONSELFNUMTYPES; index++)
					{
						recombprec[index] = 1;
					}

					double selfprec[2 * SELFING + 1][2 * SELFING + 1];
					if (SELFING)
					  {
					    selfprec[0][1] = selfprec[0][2] = recprob[2][0];
					    selfprec[0][0] = 1 - 2 * recprob[2][0];
					    selfprec[1][0] = selfgen ? selfprec[0][1] * 2.0 / ((1 << selfgen) - 1) : 1;
					    selfprec[1][2] = selfprec[1][0] * selfprec[0][1]; // TODO, UNDERESTIMATED!
					    selfprec[2][0] = selfprec[1][0];
					    selfprec[2][1] = selfprec[1][2];
					    selfprec[2][2] = selfprec[1][1] = 1 - selfprec[1][0] - selfprec[1][2];
					  }
					else
					  {
					  }


					// Compute probabilities for arbitrary xor values of current and future state
					for (int t = 0; t < TYPEBITS; t++)
					{
						int sex = TYPESEXES[t];
						int gen = TYPEGENS[t];
#pragma ivdep
						for (int index = 0; index < NONSELFNUMTYPES; index++)
						{
							int val = !((index >> t) & 1);
							recombprec[index] *= other[gen][val][sex];
						}
					}


					// Use those xor values
					// For the 4-state model, this is an inefficient way to go about it, but it is quite a bit more efficient for
					// the 64-state model (or beyond).
					for (int from = 0; from < VALIDSELFNUMTYPES; from++) // SELFING condition removes the double-bit set case, which is not allowed
 					{
						if (probs[from] < MINFACTOR || !probs[from]) continue;
						for (int to = 0; to < VALIDSELFNUMTYPES; to++)
						{
							probs2[to] += probs[from] * recombprec[(from ^ to) & (NONSELFNUMTYPES - 1)] * (SELFING ? selfprec[from >> TYPEBITS][to >> TYPEBITS] : 1);
						}
					}
 
					for (int c = 0; c < NUMTYPES; c++)
					{
						probs[c] = probs2[c];
					}
				}
				//				else
				{
					//				if (iter == tofind)
					{
						/*for (int c = 0; c < VALIDSELFNUMTYPES; c++)
						{p
							if (probs[c] < 1e-200) probs[c] = 1e-200;
						}*/
					}
				}

				startpos = endpos;
				endpos = markerposes[j];
			}			
		}

		if (updateend & 1)
		{
			adjustprobs(tb, probs, endmark, factor, ruleout, -1); // TODO
		}

		return factor;
	}
};

// Oh, how we waste memory, in a pseudo-O(1) manner
individ* individer[1000000];
// dous contains those individuals that really should be analyzed
vector<individ*> dous;

// retrieve the individual with a specific number
#ifdef _MSC_VER
individ* const getind(int n, bool real = false) 
#else
individ* const __attribute__((const)) getind(int n, bool real = false) 
#endif
{
	if (!real && !individer[n]) n = 0;
	if (n <= 0) return 0;	

	if (!individer[n])
	{
		printf("Creating %d\n", n);
		individer[n] = new individ();
		individer[n]->n = n;
		individ* ind = individer[n];
		ind->empty = false;
		  ind->markerdata.resize(markerposes.size());
		  ind->haplobase.resize(markerposes.size());
		  ind->haplocount.resize(markerposes.size());
		  ind->haploweight.resize(markerposes.size());
		  ind->negshift.resize(markerposes.size());
      
		  ind->infprobs.resize(markerposes.size());
		  ind->sureinfprobs.resize(markerposes.size());
		  ind->unknowninfprobs.resize(markerposes.size());
		  ind->parinfprobs.resize(markerposes.size());
		  ind->markersure.resize(markerposes.size());
		  for (int x = 0; x < markerposes.size(); x++)
		    {
		      ind->markerdata[x] = make_pair(UnknownMarkerVal, UnknownMarkerVal);
		      ind->haplobase[x] = 0;
		      ind->haplocount[x] = 0;
		      ind->haploweight[x] = 0.5;
		      ind->negshift[x] = 0;
		      ind->markersure[x] = make_pair(0, 0);
		    }
		  
				//		ind->semishift.resize(5000);
		  ind->lastinved.resize(chromstarts.size());
		  ind->lockstart.resize(chromstarts.size());
      
		  for (int i = 0; i < chromstarts.size(); i++)
		{
		  ind->lastinved[i] = -1;
		  ind->lockstart[i] = 0;
		}
      

	}

	return individer[n];
}

// read qtlmas sample data, code used for testing haplotypes, quite a bit of hardcoding present in this function
void readqtlmas()
{
	FILE* inpheno = fopen("phenotype.txt", "r");
	FILE* ingeno = fopen("genotype_cor.txt", "r");
	markertranslation.resize(6000);
	for (int i = 0; i < 6; i++)
	{
		chromstarts.push_back(i * 1000);
		for (int j = 0; j < 1000; j++)
		{
			markertranslation[i * 1000 + j] = i * 1000 + j;
			markerposes.push_back(j * 0.01 * 2);
			for (int t = 0; t < 2; t++)
			{
				actrec[t].push_back(baserec[t]);
			}
		}
	}
	chromstarts.push_back(6000);

	char tlf[16384];
	fgets(tlf, 16384, inpheno);

	for (int indn = 1; indn < 7000; indn++)
	{
		int fa = 0;
		int mo = 0;
		individ* ind = getind(indn);
		ind->pars[0] = getind(fa);
		ind->pars[1] = getind(mo);

		ind->gen = 5;

		ind->sex = 0;
		ind->strain = 1;

		ind->haplobase.resize(markerposes.size());
		ind->haplocount.resize(markerposes.size());
		ind->haploweight.resize(markerposes.size());
		ind->negshift.resize(markerposes.size());
		ind->markerdata.resize(markerposes.size());

		for (int i = 0; i < 6000; i++)
		{
			ind->haploweight[i] = 0.0;
		}
	}

	while (fgets(tlf, 16384, inpheno))
	{
		int indn, fa, mo, sex, gen;
		if (sscanf(tlf, "%d %d %d %d %d", &indn, &fa, &mo, &sex, &gen) < 4) break;


		individ* ind = getind(indn);

		ind->pars[0] = getind(fa);
		ind->pars[1] = getind(mo);


		if (fa && mo && (!ind->pars[0]->haplobase.size() ||
			!ind->pars[1]->haplobase.size()))
		{
			printf("PROBLEM %d %d %d\n", indn, fa, mo);
		}
		ind->gen = gen;

		ind->sex = sex - 1;
		ind->strain = (indn % 2) + 1;

		if (gen >= 2) dous.push_back(ind);

		ind->haplobase.resize(markerposes.size());
		ind->haplocount.resize(markerposes.size());
		ind->haploweight.resize(markerposes.size());
		ind->negshift.resize(markerposes.size());
		ind->markerdata.resize(markerposes.size());
	}

	int indn;
	while (fscanf(ingeno, "%d", &indn) == 1 && indn)
	{
		individ* ind = getind(indn);
		if (!ind->markerdata.size())
		{
			printf("Stopping at %d\n", indn);
			break;
		}

		for (int i = 0; i < 6000; i++)
		{
			int a, b;
			fscanf(ingeno, "%d %d", &a, &b);
			if (i % 10)
			{
				a = 0;
				b = 0;
			}

			ind->markerdata[i] = make_pair(a * MarkerValue, b * MarkerValue);

		}
	}
}

void readqtlmas14()
{
	//	FILE* inpheno = fopen("phenotype.txt", "r");
	// /bubo/home/h24/nettel/
	FILE* ingeno = fopen("genotypes.txt", "r");
	FILE* inminfo = fopen("marker-info.txt", "r");
	FILE* inped = fopen("pedigree.txt", "r");
	markertranslation.resize(10033);
	int n, c, bppos;
	int n2 = 0;
	int oldc = -1;

	fprintf(stderr, "Input: %d %d %d\n", ingeno, inminfo, inped);

	while (fscanf(inminfo, "%d %d %d", &n, &c, &bppos) == 3)
	{
		markertranslation[n - 1] = n - 1;
		if (c != oldc)
		{
			//if (oldc != -1)
			{
				chromstarts.push_back(n - 1);
			}		
			n2 = 0;
			oldc = c;
		}
		markerposes.push_back(bppos / 1000000.0);
		if (markerposes.size() > 1975) break;
		for (int t = 0; t < 2; t++)
		{
			actrec[t].push_back(baserec[t]);
		}
		n2++;
	}
	chromstarts.push_back(n);
	//chromstarts.push_back(n + 1);
	markertranslation[n] = n;
	markertranslation[n + 1] = n + 1;
	markerposes.push_back(0);
	markerposes.push_back(1);
	for (int t = 0; t < 2; t++)
	{
		actrec[t].push_back(baserec[t]);
	}
	for (int t = 0; t < 2; t++)
	{
		actrec[t].push_back(baserec[t]);
	}



	int fa, mo;
	int indn;
	char sex[3];
	while (fscanf(inped, "%d %d %d %s", &indn, &fa, &mo, sex) == 4)
	{
		individ* ind = getind(indn, true);
		ind->pars[0] = getind(fa, true);
		ind->pars[1] = getind(mo, true);

		ind->gen = 0;
		ind->gen += indn > 20;
		ind->gen += indn > 438;
		ind->gen += indn > 1349;
		ind->gen += indn > 2326;

		ind->n = indn;

		//if (ind->gen > 0) dous.push_back(ind);
		if (indn > 2326 && ind->n < 2332) dous.push_back(ind);
		//		if (indn > 3196 && ind->n < 3227) dous.push_back(ind);

		ind->sex = (sex[0] == 'F');
		ind->strain = 1;
		ind->markerdata.resize(markerposes.size());

		//if (ind->n >= 1600)
		{
			ind->haplobase.resize(markerposes.size());
			ind->haplocount.resize(markerposes.size());
			ind->haploweight.resize(markerposes.size());
			ind->negshift.resize(markerposes.size());

			ind->infprobs.resize(markerposes.size());
			ind->sureinfprobs.resize(markerposes.size());
			ind->unknowninfprobs.resize(markerposes.size());
			ind->parinfprobs.resize(markerposes.size());
			ind->markersure.resize(markerposes.size());
			//		ind->semishift.resize(5000);
			ind->lastinved.resize(chromstarts.size());
			ind->lockstart.resize(chromstarts.size());

			for (int i = 0; i < chromstarts.size(); i++)
			{
				ind->lastinved[i] = -1;
				ind->lockstart[i] = 0;
			}

			for (int i = 0; i < markerposes.size(); i++)
			{
				ind->haploweight[i] = 0.5;
			}

		}		

		//		if (ind->n > 2360) break;
	}

	int indread = 0;
	while (fscanf(ingeno, "%d", &indn) == 1 && indn)
	{
		individ* ind = getind(indn, true);
		indread++;
		if (!ind->markerdata.size())
		{
			printf("Stopping at %d\n", indn);
			break;
		}

		for (int i = 0; i < /*markerposes.size() - 2*/ 10031; i++)
		{
			int a, b;
			fscanf(ingeno, "%d %d", &a, &b);
			if (a > 5 || b > 5) fprintf(stderr, "Marker offset mismatch: %d %d %d %d\n", indn, i, a, b);

			if (i < markerposes.size())
			{
				ind->markerdata[i] = make_pair(a * MarkerValue, b * MarkerValue);

				if (ind->gen < 4)
				{
					ind->markerdata[i] = make_pair(UnknownMarkerVal, UnknownMarkerVal);
				}
				else
				  if (i == 1181) fprintf(stderr, "Input check: %d %d %d\n", indn, a, b);
			}

		}
		
		//		if (indn > 2360) break;
	}
	fprintf(stderr, "Individuals read: %d\n", indread);
	getind(1633)->markerdata[4] = make_pair(2 * MarkerValue, 2 * MarkerValue);
	getind(1351)->markerdata[4] = make_pair(2 * MarkerValue, 2 * MarkerValue);
}


// read marker info in a format similar to that accepted by ccoeff
void readmarkerinfo(FILE* in)
{
	int n, m;
	fscanf(in, "%d %d", &n, &m);
	vector<int> count;
	markertranslation.resize(m);

	for (int i = 0, j = 0; i < n; i++)
	{
		fscanf(in, "%d", &m);
		count.push_back(m);

		printf("Number on ch %d: %d \n", i + 1, m);

		for (int k = 0; k < count.end()[-1]; k++)
		{
			fscanf(in, "%d", &m);
			markertranslation[m - 1] = ++j;
		}
	}

	int pos = 0;
	for (int i = 0; i < n; i++)
	{
		chromstarts.push_back(pos);		

		vector<double> part[2];
		for (int t = 0; t < sexc; t++)
		{
			fscanf(in, "%d", &m);

			double j = 0;
			for (int k = 0; k < count[i]; k++)
			{
				double v;
				fscanf(in, "%lf", &v);
				//v = 1;
				j += v;

				for (int p = 0; p < 2 / sexc; p++)
				{
					part[t + p].push_back(j / discstep);
				}
			}
		}

		for (int k = 0; k < count[i]; k++)
		{
			double sum = 0;
			for (int t = 0; t < 2; t++)
			{
				sum += part[t][k];
			}
			sum /= 2.0;

			markerposes.push_back(sum);
			printf("%lf\n", sum);

			for (int t = 0; t < 2; t++)
			{
				double avgdist = 0;
				if (k)
				{
					avgdist = markerposes[pos] - markerposes[pos - 1];
				}

				if (avgdist)
				{					
					actrec[t].push_back(baserec[t] * (part[t][k] - part[t][k - 1]) / avgdist);
				}
				else
				{
					actrec[t].push_back(-1.0);
				}
			}
			pos++;
		}
	}
	chromstarts.push_back(pos);

	printf("%d chromosomes, %d markers, really %d/%d\n", n, m, markerposes.size(), chromstarts.size() - 1);
}

// read a pedigree in a format similar to the one accepted by ccoeff
void readped(FILE* in)
{
	int famsize;

	// The format is a set of full sibships, each started by four founders (possibly some identical), two parents
	while (fscanf(in, "%d", &famsize) == 1)
	{
		for (int i = 0; i < famsize + 6; i++)
		{
			int indn;
			int fa, mo;
			int sex = 0;

			// strain here could just as well have been called line
			// OTOH, in a line-based file format, animal lines can be confusing
			int strain = -1;
			char tlf[255];
			while (fgets(tlf, 255, in))
			{
				if (sscanf(tlf, "%d %d %d %d %d", &indn, &fa, &mo, &sex, &strain) >= 3) break;
			}

			individ* ind = getind(indn, true);
			ind->pars[1] = getind(fa, true);
			ind->pars[0] = getind(mo, true);

			if (ind->pars[0] && ind->pars[1] && ind->pars[0]->sex == 1 && ind->pars[1]->sex == 0) {
				swap(ind->pars[0], ind->pars[1]);
				printf("%d\n", indn);
			}

			ind->sex = sex - 1;
			ind->strain = strain;
			// The generation is defined by the index, the first few are always the first two generations
			ind->gen = 0;
			ind->gen += i >= 4;
			ind->gen += i >= 6;

			// Only analyze the F2s
			if (i >= 6) dous.push_back(ind);
		}
	}
	printf("Pedigree containing %d F2 individuals\n", dous.size());
}

// Read marker data and dimension several data fields dependent on marker data
// These are expected to be in a format similar to the one accepted by ccoeff
void readmarkerdata(FILE* in)
{
	int indn;
	int n = 0;
	while (fscanf(in, "%d", &indn) == 1)
	{
		individ* ind = getind(indn, true);
		n++;
		ind->markerdata.resize(markerposes.size());
		ind->haplobase.resize(markerposes.size());
		ind->haplocount.resize(markerposes.size());
		ind->haploweight.resize(markerposes.size());
		ind->negshift.resize(markerposes.size());

		for (unsigned int i = 0; i < markerposes.size(); i++)
		{
			ind->haploweight[i] = 0.5;
		}

		for (unsigned int i = 0; i < markertranslation.size(); i++)
		{
			int a, b;
			if (fscanf(in, "%d %d", &a, &b) != 2)
			{
				fprintf(stderr, "Marker mismatch: %d\n", indn);
			}
			if (markertranslation[i])
			{
				ind->markerdata[markertranslation[i] - 1] =
					make_pair(a * MarkerValue, b * MarkerValue);
			}
		}
	}
	printf("Marker data parsed for %d individuals\n", n);
}

void readhaploweights(FILE* in)
{
	int indn;
	int end = markerposes[chromstarts[1] - 1];

	while (fscanf(in, "%d", &indn) == 1)
	{
		// TODO chromcount?
		individ* ind = getind(indn, true);
		double num1, num2;
		int a, b;
		for (int i = chromstarts[0]; i < chromstarts[1]; i++)
		{
			fscanf(in, "%lf %d %d %lf", &num1, &a, &b, &num2);
			ind->haploweight[i] = num1;

			if (ind->haploweight[i] <= 0.5)
			{
				swap(ind->markerdata[i].first, ind->markerdata[i].second);
				ind->haploweight[i] = 1 - ind->haploweight[i];
			}
			if (ind->haploweight[i]) ind->haploweight[i] = 0.5;
		}	  
		ind->genotypegrid.resize(end + 1);

		if (!ind->pars[0] && !ind->pars[1])
		{
			for (int i = 0; i <= end; i++)
			{
				ind->genotypegrid[i] = (ind->strain - 1) * 3;
			}
		}
	}

	FILE* indgrid = fopen("indgrid.out", "w");
	FILE* indmarkers = fopen("indmarkers.out", "w");
	for (int gen = 0; gen <= 2; gen++)
	{
		for (int i = 0; i < 1000000; i++)
		{
			individ* ind = getind(i, true);		  

			if ((ind && ind->gen == gen && ind->genotypegrid.size()) /*&& ind->pars[0] && ind->pars[1]*/)
			{
				fprintf(indmarkers, "%d ", i);
				fprintf(indgrid, "%d ", i);		  

				int phase = myrand(4);
				int marker = 0;

				if (gen)
				{
					for (int marker = 0; marker < chromstarts[1]; marker++)
					{
						if (myrand(2))
						{
							swap(ind->markerdata[marker].first, ind->markerdata[marker].second);
						}
					}
				}


				for (int i = 0; i <= end; i++)
				{
					int gval = 0;
					int startmarker = marker;
					if (gen)
					{
						for (int p = 0; p < 2 && gen; p++)
						{
							marker = startmarker;

							int mask = 1 << p;
							bool parhalf = phase & mask;

							gval |= mask * ((bool) (ind->pars[p]->genotypegrid[i] & (1 << parhalf)));

							for (;markerposes[marker] < i + 1 && marker < chromstarts[1]; marker++)
							{
								if ((&ind->markerdata[marker].first)[p] != UnknownMarkerVal)
								{
									(&ind->markerdata[marker].first)[p] = (&ind->pars[p]->markerdata[marker].first)[parhalf];
								}
							}

							phase &= 65535 - mask;

							int high = 100;
							int lo = 0;

							if (marker < chromstarts[1])
							{
								high = 10000;
								if (actrec[p][marker])
									lo = actrec[p][marker] * -5000;
							}

							if (myrand(high) <= lo) parhalf = !parhalf;

							phase |= parhalf * mask;
						}
						ind->genotypegrid[i] = gval;
					}

					fprintf(indgrid, " %d", ind->genotypegrid[i]);			  
				}
				for (int marker = 0; marker < chromstarts[1]; marker++)
				{
				  fprintf(indmarkers, " %d %d", ind->markerdata[marker].first.value(), ind->markerdata[marker].second.value());
				}

				fprintf(indmarkers, "\n");
				fprintf(indgrid, "\n");
			}
			/*else
			{
			for (int i = 0; i <= end; i++)
			{
			fprintf(indgrid, " %d", ind->genotypegrid[i]);
			}
			}*/


		}
	}
}

void lockhaplos(individ* ind, int i)
{
	unsigned int j;
	if (ind->n > 220) return;

	if (ind->lockstart[i] >= chromstarts[i + 1])
	{
		ind->lockstart[i] = 0;
	}
	int bestpos = -1;
	double bestcount = -1;

	for (j = max(chromstarts[i], ind->lockstart[i]); j != chromstarts[i + 1]; j++)
	{
		if (ind->markerdata[j].first == ind->markerdata[j].second) continue;
		if (fabs(ind->haploweight[j] - 0.5) > 0.49999) continue;	

		individ* pars[3] = {ind, ind->pars[0], ind->pars[1]};
		double count = 0;

		for (int p = 0; p < 3; p++)
		{
			if (!pars[p] || pars[p]->markerdata.size() <= j) continue;

			for (int i = 0; i < 2; i++)
			{
			  if ((&(pars[p]->markerdata[j].first))[i] != UnknownMarkerVal)
			    {
			      count += 1.0 - (&(pars[p]->markersure[j].first))[i];
			    }
			}			
			if (p == 0) count *= 2;						
		}

		if (count > bestcount)
		{
			bestcount = count;
			bestpos = j;
		}	
	}

	if (bestpos == -1)
	{
		return;
	}
	j = bestpos;

	printf("Fixing point: %d %d %d\n", ind->n, i + 1, j);
	fflush(stdout);

	if (j != chromstarts[i + 1])
	{
		ind->haploweight[j] = 0;
		ind->lockstart[i] = j + 1;
	}
}


double dosureval(int what, pair<int, double> val)
{
  if (val.second == 0) return 0;

  /*  double toret = log(val.second);
  toret += log(0.5) * (what * 0.5 - val.first);
  toret *= 2;
  toret /= what;*/
  double toret = log(val.second);
  toret /= what;
  toret *= 4;
  toret = exp(toret);
  toret /= (1 + toret);

  return toret;
}

// Some operations performed when marker data has been read, independent of format.
void postmarkerdata()
{
	int any, anyrem;
	bool latephase = false;


	markerweight.resize(markerposes.size());
	// If inference is active, add "new" marker data until all data has been found.
	if (CORRECTIONINFERENCE) do
	{
#pragma omp parallel for schedule(dynamic,32)
		// all haploweights MUST be non-zero at this point, as we do not explore all shiftflagmode values
		for (int i = 1; i < 1000000; i++)
		{
			individ* ind = getind(i);
			ind->children = 0;
		}
		for (int i = 1; i < 1000000; i++)
		{
			individ* ind = getind(i);
			for (int j = 0; j < 2; j++)
			  {
			  if (ind->pars[j])
			    {
			      ind->pars[j]->children++;
			    }
			  }

			if (ind->markerdata.size())
			{
				generation++;

				for (int g = 0; g < ind->markerdata.size(); g++)
				{
					ind->fixparents(g, latephase);
				}
			}
		}

		any = 0;
		anyrem = 0;
		for (int i = 1; i < 1000000; i++)
		{
			individ* ind = getind(i);
			if (ind->sex) continue;

			for (int g = 0; g < (int) ind->markervals.size(); g++)
			{
			  //				if (g == 1) continue;

				int startsize = ind->markervals[g].size();
				const int known =
					(ind->markerdata[g].first != UnknownMarkerVal) +
					(ind->markerdata[g].second != UnknownMarkerVal);

				const int oldany = any;

				if (latephase && known == 2 && !startsize && ind->gen < 1)
				{
					ind->markerdata[g].second = UnknownMarkerVal;
					ind->markerdata[g].first = UnknownMarkerVal;
					ind->markersure[g] = make_pair(0.0, 0.0);
					any++;
					anyrem++;
				}
				else
					if (known == 2) continue;

				ind->markervals[g].erase(UnknownMarkerVal);
				startsize = ind->markervals[g].size();

				// map insert preserves existing value if present
				if (ind->markerdata[g].first  != UnknownMarkerVal) ind->markervals[g].insert(make_pair(ind->markerdata[g].first, make_pair(ind->children, ind->markersure[g].first)));
				if (ind->markerdata[g].second != UnknownMarkerVal) ind->markervals[g].insert(make_pair(ind->markerdata[g].second, make_pair(ind->children, ind->markersure[g].second)));

				if (ind->markervals[g].size() >= 3)
				{
					fprintf(stderr, "Error, too many matches: %d\t%d\n", i, g);
				}
				if (!latephase && ind->markervals[g].size() == 2)
				{
				  int knowncount = ind->markervals[g].begin()->second.first + (++ind->markervals[g].begin())->second.first;
					ind->markerdata[g] = make_pair(ind->markervals[g].begin()->first, (++ind->markervals[g].begin())->first);
					ind->markersure[g] = make_pair(dosureval(knowncount, ind->markervals[g].begin()->second), dosureval(knowncount, (++ind->markervals[g].begin())->second));
					any++;
				}
				if (latephase && ind->markervals[g].size() == 1 && startsize == 1 && known == 1 && !ind->pars[0] && !ind->pars[1])
				{
					if (ind->markerdata[g].first == UnknownMarkerVal || ind->markerdata[g].second == UnknownMarkerVal) any++;
					ind->markerdata[g] = make_pair(ind->markervals[g].begin()->first, ind->markervals[g].begin()->first);		
					ind->markersure[g] = make_pair(dosureval(ind->children, ind->markervals[g].begin()->second),
								       dosureval(ind->children, ind->markervals[g].begin()->second));		       
				} // DANGEROUS ASSUMPTIONS
				else if (!latephase && ind->markervals[g].size() == 1 && known == 0)
				{
					any++;
					ind->markerdata[g] = make_pair(ind->markervals[g].begin()->first, UnknownMarkerVal);
					ind->markersure[g] = make_pair(dosureval(ind->children, ind->markervals[g].begin()->second), 0.0);
				}

				
				if (any != oldany) printf("Correction at %d, marker %d (%d;%d) (%lf;%lf)\n", i, g,
							  ind->markerdata[g].first.value(), ind->markerdata[g].second.value(),
							  ind->markersure[g].first, ind->markersure[g].second);

			}

			for (int g = 0; g < (int) ind->markervals.size(); g++)
			{
				if (ind->markerdata[g].first == sexmarkerval) {
					ind->markerdata[g] = make_pair(ind->markerdata[g].second, ind->markerdata[g].first);
				}
				ind->markervals[g].clear();
			}
		}
		fprintf(stderr, "Number of corrected genotypes: %d\n", any);
		if (latephase)
		{
			latephase = false;
		}
		else
		{
		  /*			if (!any && !latephase)
			{
				any++;
				latephase = true;
				}*/
		}
	}
	while (any > anyrem);

	for (int i = 1; i < 1000000; i++)
	{
		individ* ind = getind(i);

		// Lock the first position
		if (HAPLOTYPING && ind && ind->haploweight.size())
			// These days we do the locking in all generations
		{
			ind->markervals.clear();
			// Lock the haplotype (in an arbitrary manner) for the first marker in each linkage group
			// Propagation would be more efficient if we locked the mid-position (which should really be determined in cM)
			// Doing so would be much more opaque, though...
			for (unsigned int i = 0; i < chromstarts.size() - 1; i++)
			{
			  if (!ind->pars[0] && !ind->pars[1])
				lockhaplos(ind, i);				
			}
		}
	}
}

typedef boost::tuple<individ*, double, int> negshiftcand;
struct inferiorrelated
{
	negshiftcand ourtuple;
	bool anymatch;

	inferiorrelated(negshiftcand init) : ourtuple(init), anymatch(false)
	{
	}

	bool operator() (negshiftcand b)
	{
		if (b.get<0>()->arerelated(ourtuple.get<0>()))
		{
			if (b.get<1>() > ourtuple.get<1>()) return true;
			else
				anymatch = true;
		}

		return false;
	}
};

struct negshifter
{
	int c;
	negshifter(int c) : c(c)
	{
	}

	void operator() (negshiftcand b)
	{
		int minstart = b.get<2>();
		individ* ind = b.get<0>();
		ind->lastinved[c] = minstart;
		for (int p = minstart + 1; p < (int) chromstarts[c + 1]; p++)
		{
			if (p == minstart + 1) fprintf(stdout, "Inv: %d %d\n", ind->n, p);
			if (ind->n < 2327) ind->haploweight[p] = 1.0f - ind->haploweight[p];
		}
	}
};

bool ignoreflag2(int flag2, int g, int q, int flag2ignore, const map<individ*, int>& relmap)
{
  int flag2filter = (1 << 30) - 1;
  // Below lines relied on incorrect assumption of remapping of inheritance
  /*  const int selfval = g >> TYPEBITS;
    if (SELFING && selfval)
      {
	int basefilter = (HALFNUMPATHS - 1) << 1;
        if ((((flag2 ^ (g * 2)) / HALFNUMPATHS) ^ (flag2 ^ (g * 2))) & basefilter)
	  {
	    return true;
	  }
	flag2filter = basefilter * (selfval == 1 ? HALFNUMPATHS : 1);
	flag2filter |= 1;
	//	printf("Selfval is %d, flag2filter is %d, flag2ignore is %d\n", selfval, flag2filter, flag2ignore);
	}*/
    if (flag2 & (flag2ignore & flag2filter)) return true;

	int marker = -q-1000;
	for (map<individ*, int>::const_iterator i = relmap.begin(); i != relmap.end(); i++)
	{
	  int currfilter = (i->second & flag2filter);
		int filtered = ((flag2 ^ (g * 2)) & currfilter);
		// Require ALL bits in the flag to be set, if at least one is set
		if (filtered && filtered != currfilter) return true;
		//if (marker >= 0 && i->first->markerdata[marker].first == UnknownMarkerVal && i->first->markerdata[marker].second == UnknownMarkerVal && (!(flag2 & i->second)))
		if (marker >= 0 && i->first->markerdata[marker].first == i->first->markerdata[marker].second && i->first->markersure[marker].first == i->first->markersure[marker].second && (!(flag2 & currfilter)) && (!SELFING || currfilter != 1 /*|| selfgen == 0*/))
			{
			  				return true;
							}
	}
	return false;
}

// The actual walking over all chromosomes for all individuals in "dous"
// If "full" is set to false, we assume that haplotype inference should be done, over marker positions.
// A full scan is thus not the iteration that takes the most time, but the scan that goes over the full genome grid, not only
// marker positions.
template<bool full> void doit(FILE* out, bool printalot
#ifdef F2MPI
							  , mpi::communicator& world
#endif
							  )
{
	const bool doprint = full;
#ifdef F2MPI
	broadcast(world, actrec, 2, 0);
#endif

	int count = 0;
	vector<vector<boost::array<float, 2> > > realgeno;	

	realgeno.resize(dous.size());

	for (int a = 1; a <= 20; a++)
	{
		individ* ind2 = getind(a);
		/*	    ind2->semishift.resize(5000);
		for (int q = 0; q < 5000; q++)
		{
		for (int z = 0; z < 40; z++)
		{
		for (int p = 0; p < 2; p++)
		{
		for (int r = 0; r < 2; r++)
		{
		ind2->semishift[q][p * 2 + r][z] = (z == (a - 1) * 2 + p);
		}
		}
		}
		*/
	}

	for (unsigned int j = 0; j < dous.size(); j++)
	{
		if (dous[j]->markerdata.size())
		{
			count++;
		}
		else
		{
			dous.erase(dous.begin() + j);
		}
	}

	static int iter = 0;
	iter++;


	map<pair<individ*, individ*>, map<int, boost::array<double, 8> > > nsm;
	if (doprint)
	{
		//fprintf(out, "%d %d\n", count, chromstarts.size() - 1);
	}

	for (int i = 0; i < 1000000; i++)
	{
		individ* ind = getind(i);
		if (!ind) continue;
		ind->children = 0;
		ind->kids.clear();
		if (!ind || !ind->haplocount.size()) continue;

		for (unsigned int j = 0; j < ind->haplocount.size(); j++)
		{
			ind->haplocount[j] = 0.0f;
			ind->haplobase[j] = 0.0f;
		}

		//		fprintf(out, "B%d:%d\n", world.rank(), i);
#ifdef F2MPI
		broadcast(world, ind->haploweight, 0);

		//		fprintf(out, "C%d:%d\n", world.rank(), i);

		world.barrier();
#endif
		//		fflush(out);

	}

	for (int j = 0; j < (int) dous.size(); j++)
	{

		for (int i = 0; i < 2; i++)
		{
			individ* pnow = dous[j]->pars[i];

			if (pnow)
			{
				pnow->children++;
				pnow->kids.push_back(dous[j]);
			}
		}
	}


	for (unsigned int i = 0; i < chromstarts.size() - 1; i++)
	{
		if (doprint)
		{
			//fprintf(out, "%d %d\n", i + 1, (int) markerposes[chromstarts[i + 1] - 1]);
		}
		//printf("Chromosome %d\n", i + 1);

		// The output for all individuals in a specific iteration is stored, as we have parallelized the logic and 
		// want the output to be done in order.
		vector<vector<char> > outqueue;

		outqueue.resize(dous.size());


#pragma omp parallel for schedule(dynamic,1)
		for (int j = 0; j < (int) dous.size(); j++)
		{
#ifdef F2MPI
			if (j % world.size() != world.rank()) continue;			
#endif
			//			realgeno[j].resize(markerposes[chromstarts[1] - 1] + 1);

			generation++;
			threadblock tborig;
			threadblock tb = tborig;


			// Some heaps are not properly synchronized. Putting a critical section here makes the operations not safe,
			// but *safer*.
#pragma omp critical(uglynewhack)
			for (int t = 0; t < NUMSHIFTS; t++)
			{
				factors[t].resize(markerposes.size());
				cacheprobs[t].resize(markerposes.size());
				done[t].resize(markerposes.size());
			}

			if (dous[j]->markerdata.size())
			{				
				int qstart = -1000 - chromstarts[i];
				int qend = -1000 - chromstarts[i + 1];
				int qd = -1;
				int f2s = 0;
				int f2end = NUMPATHS;

				if (!HAPLOTYPING)
				{
					f2s = -1;
					f2end = 0;
				}

				int shifts = 0;
				int shiftend = NUMSHIFTS;

				reltree.clear();
				relmap.clear();
				reltree.push_back(dous[j]);
				relmap[dous[j]] = 1;
				int flag2ignore = 0;

				// Special optimization hardcoded for this population structure, eagerly skipping flags that do not
				// correspond to any inheritance, i.e. if not the full pedigree of 6 individuals back is present.
				int shiftignore = 0;
				if (HAPLOTYPING)
				{
					flag2ignore = 1;
					shiftignore = 1;
					for (int lev1 = 0; lev1 < 2; lev1++)
					{
						individ* lev1i = dous[j]->pars[lev1];
						if (!lev1i) continue;
						int flag2base = 1 << (1 + lev1 * ((1 << (NUMFLAG2GEN - 1)) - 1));
						if (!lev1i->empty)
							{
								flag2ignore |= flag2base;
								relmap[lev1i] |= flag2base;
						}

						bool anypars = false;
						reltree.push_back(lev1i);
						if (NUMGEN > 2)
						{
							for (int lev2 = 0; lev2 < 2; lev2++)
							{
								individ* lev2i = lev1i->pars[lev2];
								if (!lev2i) continue;

								if (!lev2i->empty)
									{
										flag2ignore |= (flag2base << (lev2 + 1));
										relmap[lev2i] |= (flag2base << (lev2 + 1));
								}
								anypars = true;
								reltree.push_back(lev2i);
							}
						}
						if (anypars)
						  {
						    shiftignore |= 2 << lev1;
						  }
					}

					flag2ignore ^= (NUMPATHS - 1);
					shiftignore ^= (NUMSHIFTS - 1);
				}

				sort(reltree.begin(), reltree.end());
				reltree.resize(unique(reltree.begin(), reltree.end()) - reltree.begin());

				bool skipsome = false;
				for (int u = 0; u < reltree.size() && !skipsome; u++)
				{
					for (int p = 0; p < 2; p++)
					{
						if ((&(reltree[u]->markerdata[chromstarts[i]].first))[p] == sexmarkerval)
							skipsome = true;
					}
				}

				if (!HAPLOTYPING)
				{
					reltree.resize(0);
					flag2ignore = 0;
				}

				if (full)
				{
					qstart = (int) markerposes[chromstarts[i]];
					qend = (int) markerposes[chromstarts[i + 1] - 1] + 1;
					qd = 1;
					f2s = -1;
					f2end = 0;
					/*shifts =s 0;
					shiftend = 1;*/
				}

				if (dous[j]->gen < 2) shiftend = min(2, shiftend);

				// BIG SHIFT HACK, removed again as it makes resolving a shared 12 12 ancestry very hard
				bool anygood = false;
				anygood |= dous[j]->pars[0] && !dous[j]->pars[0]->empty;
				anygood |= dous[j]->pars[1] && !dous[j]->pars[1]->empty;
				//if (!anygood) shiftend = 2;
				if (!anygood)
				{
						shiftignore = 7;
						flag2ignore = 0;
				}

				double factor = -1e15;
				double factors[NUMSHIFTS];
				for (shiftflagmode = shifts; shiftflagmode < shiftend; shiftflagmode++)
				{
					if (shiftflagmode & shiftignore) //continue;
							factors[shiftflagmode] = -1e30;
					else
 						factors[shiftflagmode] = dous[j]->doanalyze<noneturner>(tb, none, chromstarts[i], chromstarts[i + 1] - 1, NONESTOP, -1, false, 0, -10000 + factor);
					factor = max(factor, factors[shiftflagmode]);
				}

				// Normalize!
				double realfactor = 0;
				for (int s = shifts; s < shiftend; s++)
				{
					if (s & shiftignore) continue;
					realfactor += exp(factors[s] - factor);
				}
				int countf2i = 0;
				for (int k = 0; k < TYPEBITS + 1; k++)
				{
					if (flag2ignore & (1 << k)) countf2i++;
				}
				factor -= log(1 << countf2i);

								printf("%d,%03d,%03d: %lf\t", dous[j]->n, flag2ignore, shiftignore, factor);
				factor += log(realfactor);
								printf("%lf %d\n", factor, shiftend);
								fflush(stdout);

				//fflush(stdout);
				// Flushing can be useful for debugging, but not for performance!
				// This output can get ugly due to race conditions. One shouldn't rely on it.

				if (_isnan(factor)) continue;

				char lineout[255];

/*				// States are mapped onto values describing the line/strain origin, in the sense of 00, 01, 10 or 11 (old comment)
				PerStateArray<int>::T maptogeno;
				shiftflagmode = 0;
				for (int g = 0; g < NUMTYPES; g++)
				{
					int sum = 0;
					dous[j]->trackpossible<false, true>(tb, UnknownMarkerVal, 0
						, 0, g * 2, 0, 0, trackpossibleparams(0, &sum));u

					maptogeno[g] = sum;
					printf("Maptogeno %d\n", sum);
				}*/

				// Walk over all chromosome positions, whether it be markers (negative q values <= -1000) or grid positions


				for (int q = qstart; q != qend; q+=qd)
				{
					double probs[4] = {0};
					//double mwvals[NUMTYPES][NUMTYPES] = {0};
					//double mwfvals[NUMTYPES] = {0};
					double mwvals[1][1];
					double mwfvals[1];

					double mwval[4] = {0};

					for (int g = 0; g < NUMTYPES; g++)
					{						
						for (shiftflagmode = shifts; shiftflagmode < shiftend; shiftflagmode++)
						{
							if (shiftflagmode & shiftignore) continue;
							if (factor - factors[shiftflagmode] > 10) continue;
							if (q <= -1000 && false)
							{
								double val;

								for (int g2 = 0; g2 < NUMTYPES; g2++)
								{
									val = dous[j]->doanalyze<noneturner>(tb, none, chromstarts[i], chromstarts[i + 1] - 1, twicestop(q, g, g2),
										-1, true, 0, -5000.0 + factor) - factor;


									val = exp(val);	
									mwvals[g][g2] += val;
									if (q == -1010 && false)
									{
										printf("%d: %d -> %d: %lf\n", dous[j]->n, g, g2, val);
									}
									//									fprintf(stderr, "%d:%d %d %d %lf\n", dous[j]->n, q, g, g2, val);

								}

								val = 0;
								/*								if (g == 0)
								{
								stopmodpair::miniactrecT newactrec;
								for (int gs = 0; gs < 2; gs++)
								{
								for (int k = 0; k < 2; k++)
								{
								newactrec[k][gs] = actrec[k][-q - 1000 + 1 + gs];
								}
								}

								double sum = newactrec[0][0] + newactrec[0][1];
								sum /= 2;
								newactrec[0][0] = 0;
								newactrec[0][1] = sum * 2;

								val = dous[j]->doanalyze<noneturner>(tb, none, chromstarts[i], chromstarts[i + 1] - 1, stopmodpair(q, newactrec),
								-1, true, 0, -20.0 + factor) - factor;
								}*/


								mwfvals[g] += val;
							}

							for (int flag2 = f2s; flag2 < f2end; flag2++)
							{
							  if (ignoreflag2(flag2, g, q, flag2ignore, relmap)) continue;
								//if (flag2 & (flag2ignore)) continue;

								int firstpar = 0;
								double val;

								if (q <= -1000 && flag2 != -1)
								{
									// Do a lookup in the impossible structure. If we are at a marker position, those branches that are impossible will have an exact
									// score of 0.
									// If we are NOT at a marker, this information can still be used in trackpossible, but not in here, as there is some recombination
									// added between the marker and fixated position we're checking.
									for (int a = 0; a < 2; a++)
									{
										if (a) firstpar = !firstpar;
										const int genwidth = (1 << (NUMGEN - 1));

										int f2n = ((flag2 ^ shiftflagmode) & 1);

										int upflagr = upflagit(g, firstpar, genwidth);
										int upflag2r = upflagit(flag2 >> 1, firstpar, genwidth);
										int upshiftr = upflagit(shiftflagmode >> 1, firstpar, genwidth >> (NUMGEN - NUMSHIFTGEN));
										int marker = -(q + 1000);
										int impossibleval = generation * markerposes.size() + marker;

										if (impossible[shiftflagmode & 1][firstpar][f2n][upflagr][upflag2r + 1][upshiftr][marker & 3] == impossibleval)
										{
											goto continueloop;
										}
									}
								}


								val = dous[j]->doanalyze<noneturner>(tb, none, chromstarts[i], chromstarts[i + 1] - 1, classicstop(q, g),
									flag2, true, 0, -15.0 + factor) - factor;

								if (_finite(val) && val > -40.0)
								{
									// shift mode not included, this is the "real" f2n, indicating what value
									// in the marker pair is used, not the strand phase (strand phase is flag2 xored
									// with the other stuff)
									int f2n = ((flag2 /*^ shiftflagmode*/) & 1);


									val = exp(val);
									int marker = -q - 1000;

									double pival = val;
									for (int i = 0; i < 2; i++)
									{
										int index = i * (TYPEBITS / 2);
										int r = dous[j]->n;
										int parindex = ((flag2 >> (index + 1))/* ^ (g >> index)*/) & 1;
										// FOR g shift
										// parindex = !parindex;
										int flag2parindex = (flag2 >> (index + 1)) & 1;

										/*int otherindex = (!i) * (TYPEBITS / 2);
										int oparindex = ((flag2 >> (otherindex + 1)) /*^ (g >> index)*) & 1;*/

										int updateval = f2n ^ i;

										double factor = 1;
										individ* parnow = dous[j]->pars[i];
										
										if (parnow)
										{
											MarkerVal mv = (&(parnow->markerdata[marker].first))[flag2parindex];
											if (mv == UnknownMarkerVal || mv == (&(dous[j]->markerdata[marker].first))[updateval])
											{
												factor = (&(parnow->markersure[marker].first))[flag2parindex];
											}
											else
												factor = 1 - (&(parnow->markersure[marker].first))[flag2parindex];
										}
										factor += 1e-5;
										factor = 1;

										pival *= factor;
									}

									for (int i = 0; i < 2; i++)
									{
										int index = i * (TYPEBITS / 2);
										int r = dous[j]->n;
										int parindex = ((flag2 >> (index + 1))/* ^ (g >> index)*/) & 1;
										// FOR g shift
										// parindex = !parindex;
										int flag2parindex = (flag2 >> (index + 1)) & 1;

										/*int otherindex = (!i) * (TYPEBITS / 2);
										int oparindex = ((flag2 >> (otherindex + 1)) /*^ (g >> index)*) & 1;*/

										int updateval = f2n ^ i;

										for (int evil = 0; evil < 2; evil++)
										{
											individ* parnow = dous[j]->pars[i];
											double factor = 1;

											if (parnow && evil)
											{
												MarkerVal mv = (&(parnow->markerdata[marker].first))[flag2parindex];
												if (mv == UnknownMarkerVal || mv == (&(dous[j]->markerdata[marker].first))[updateval])
												{
													factor = (&(parnow->markersure[marker].first))[flag2parindex];
												}
												else
													factor = 1 - (&(parnow->markersure[marker].first))[flag2parindex];

												factor += 1e-3;
											}
											pival = val;
											//										    pival *= factor;

											dous[j]->parinfprobs[marker][i][parindex][updateval][evil] += pival;

										}

										/*if ((!dous[j]->pars[+i] || (&(dous[j]->pars[+i]->markerdata[marker].first))[+parindex] == UnknownMarkerVal) &&
										(!dous[j]->pars[!i] || (&(dous[j]->pars[!i]->markerdata[marker].first))[oparindex] == UnknownMarkerVal) &&
										(dous[j]->markerdata[marker].first != dous[j]->markerdata[marker].second))
										{
										evil = true;
										}*/										

									}
									int mapval = 0;
									/*if (HAPLOTYPING)
									{
										for (int lev1 = 0; lev1 < 2; lev1++)
										{
											individ* lev1i = dous[j]->pars[lev1];
											if (!lev1i) continue;
											int flag2base = 1 << (1 + lev1 * ((1 << (NUMFLAG2GEN - 1)) - 1));
											int genbase = lev1 * ((1 << (TYPEBITS / 2)));
											if (NUMGEN > 2)
											{
												int f2what = 
												for (int lev2 = 0; lev2 < 2; lev2++)
												{
													individ* lev2i = lev1i->pars[lev2];
													if (!lev2i) continue;

													int f2ninner = (bool) (flag2 & (flag2base << (lev2 + 1)));
													mapval += (&(lev2i->markerdata[-q-1000].first))[f2ninner] == (2 * MarkerValue);
												}
											}
										}
									}*/
									int g3 = (g & 1) + ((bool)(g & 8)) * 2;
									probs[g3] += val;
									if (!full && HAPLOTYPING) dous[j]->updatehaplo(tb, -q - 1000, g, flag2, val);
								}
continueloop:;
							}
						}
					}

					if (q <= -1000 && false)
					{
						double colsums[NUMTYPES] = {0};
						double rowsums[NUMTYPES] = {0};
						double acc3 = 0;
						double acc4 = 0;
						for (int g = 0; g < NUMTYPES; g++)
						{
							double acc1 = 0;
							double acc2 = 0;
							for (int g2 = 0; g2 < NUMTYPES; g2++)
							{
								acc1 += mwvals[g][g2];
								rowsums[g] += mwvals[g][g2] * mwvals[g][g2];

								acc2 += mwvals[g2][g];
								colsums[g] += mwvals[g2][g] * mwvals[g2][g];
							}
							acc3 += acc1;
							acc4 += mwfvals[g];

							colsums[g] -= (acc2 * acc2) / NUMTYPES;
							rowsums[g] -= (acc1 * acc1) / NUMTYPES;
						}
						//												  printf("\t\t%d:%d\t%lf\t\t%lf\n", -q - 1000, dous[j]->n, acc3, acc4);

						double relinfo1 = 0;
						double relinfo2 = 0;
						double infosum[2] = {0};
						double recprob[2];
						double dist = markerposes[ - q - 1000 + 1] - markerposes[-q - 1000];

						for (int k = 0; k < 2; k++)
						{
							recprob[k] = 0.5 * (1.0 - exp(actrec[k][-q - 1000 + 1] * (dist)));
							recprob[k] = max(1e-8, recprob[k]);
							//					if (iter == tofind) recprob[k] = max(recprob[k], 1e-5);
						}

						for (int g = 0; g < NUMTYPES; g++)
						{
							relinfo1 += colsums[g];
							relinfo2 += rowsums[g];
							for (int g2 = 0; g2 < NUMTYPES; g2++)
							{
								double infofactor = /*rowsums[g] * colsums[g2]*/ 1;

								double val = mwvals[g][g2] * infofactor;

								int mask = g ^ g2;
								for (int i = 0; i < TYPEBITS; i++)
								{
									double corr = 1;
									if (NUMGEN > 2)
									{
										if (i % (TYPEBITS >> 1))
										{
											// TODO: handle null
											corr /= dous[j]->pars[i / (TYPEBITS) >> 1]->children;
										}
									}	
									bool switched = ((bool) (mask & (1 << i)));
									double expected = !switched ? 1.0 - recprob[TYPESEXES[i]] : recprob[TYPESEXES[i]];

									mwval[TYPESEXES[i] * 2 + switched] += (val / expected) * corr;
									infosum[TYPESEXES[i]] += corr;
								}
							}
						}

						double summw = (mwval[0] + mwval[1] + mwval[2] + mwval[3]);
						//						if (summw > 1e-4)
						//												printf("%.3lf\t%.3lf\n", infosum[0], infosum[1]);
						for (int z = 0; z < 2; z++)
						{
							infosum[z] /= NUMTYPES * NUMTYPES;
						}
						if (acc3 == 0)
						{
							acc3 = 1;
							infosum[0] = 0;
							infosum[1] = 0;
							//						    printf("Zero sum: %d\n", dous[j]->n);
						}


						double delta = 0;
#pragma omp critical(markerweights)					       
						{
							//						  summw /= relinfo1 * relinfo2;
							//summw /= infosum;
							for (int t = 0; t < 4; t++)
							{
								int tbase = (t >> 1) << 1;

								double dval = (mwval[t] / acc3 - infosum[t / 2]) /*/ summw*/
									* fabs(((t & 1) ? 0.0 : 1.0) - recprob[t / 2])
									/** (-1 + (t & 1) * 2)*/;

								markerweight[-q - 1000][t] += dval;
								if (t == 0) delta = dval;
							}
							//														printf("%.3lf\t%.3lf\t\t%.4lf\t%.4lf\n", markerweight[-q - 1000][0], markerweight[-q - 1000][1], mwval[0], mwval[1]);
							markerweight[-q - 1000][4] = min(markerweight[-q - 1000][4], acc4);
						}

						if (delta < -0.1 && q > qend + 2) fprintf(out, "RECOMB:\t%d\t%d\t%lf\n", dous[j]->n, q, delta);
					}

					//					fprintf(stderr, "marker: %d\n", -q - 1000);
					if (!full)
					{
						int marker = -q - 1000;
						double pinfsum[2] = {0, 0};
						for (int a = 0; a < 2; a++)
						{
							for (int c = 0; c < 2; c++)							
							{
								for (int b = 0; b < 2; b++)
								{
									for (int d = 0; d < 1; d++)
									{
										pinfsum[a] += dous[j]->parinfprobs[marker][a][c][b][d];
									}
								}
							}

							if (pinfsum[a] < 1e-10) pinfsum[a] = 1e-10;
						}


						if (DOINFPROBS) {
#pragma omp critical(infprobs)
						{

							// Add some diffusion to avoid numerical runaway scenarios for essentially symmetrical cases
							/*if (dous[j]->pars[0] && dous[j]->pars[1])
							for (int a = 0; a < 2; a++)
							{
							for (int b = 0; b < 2; b++)
							{
							if ((&(dous[j]->pars[0]->markerdata[marker].first))[a] !=
							(&(dous[j]->pars[1]->markerdata[marker].first))[b]) continue;

							double diff = (&(dous[j]->pars[0]->markersure[marker].first))[a] -
							(&(dous[j]->pars[1]->markersure[marker].first))f[b];

							if (fabs(diff) < 0.001)
							{
							diff *= 0.001 - fabs(diff);
							diff /= 0.001;

							(&(dous[j]->pars[0]->markersure[marker].first))[a] -= diff;
							(&(dous[j]->pars[1]->markersure[marker].first))[b] += diff;
							}
							}
							}*/


							for (int a = 0; a < 2; a++)
							{
								if (!dous[j]->pars[a]) continue;

								for (int c = 0; c < 2; c++)
								{
									double maxval = 0;
									double partsum = 0;
									for (int b = 0; b < 2; b++)
									{
										//		  if ((&(dous[j]->markerdata[marker].first))[b] == UnknownMarkerVal) continue;

										// pinfsum calculated for d, or evil, equalling 0, 
										for (int d = 0; d < 2; d++)
										{
											dous[j]->parinfprobs[marker][a][c][b][d] /= pinfsum[a];

											if ((&(dous[j]->markerdata[marker].first))[b] == UnknownMarkerVal)
											{
												dous[j]->pars[a]->unknowninfprobs[marker][c] += dous[j]->parinfprobs[marker][a][c][b][d];
												continue;
											}


											dous[j]->pars[a]->infprobs[marker][c][make_pair((&(dous[j]->markerdata[marker].first))[b], (&(dous[j]->markerdata[marker].first))[!b])]
											+= dous[j]->parinfprobs[marker][a][c][b][d];

											/*if ((dous[j]->pars[a]->n == 1633 || dous[j]->pars[a]->n == 1726))
											{
											fprintf(out, "%d %d contributes %lf to %d for %d, phase %d, evil %d\n", dous[j]->n, marker, dous[j]->parinfprobs[marker][a][c][b][d],
											(&(dous[j]->markerdata[marker].first))[b], dous[j]->pars[a]->n, c, d);
											}*/
											partsum += dous[j]->parinfprobs[marker][a][c][b][d];
										}

										if ((&(dous[j]->markerdata[marker].first))[b] == UnknownMarkerVal) continue;

										//										if ((&(dous[j]->markerdata[marker].first))[b] == 1 * MarkerValue && dous[j]->pars[a]->n == 1726 && marker == 276 && out) fprintf(out, "Oddone: %d %lf\n", dous[j]->n, dous[j]->parinfprobs[marker][a][c][b][0]);

										if (dous[j]->parinfprobs[marker][a][c][b][0] > maxval) maxval = dous[j]->parinfprobs[marker][a][c][b][0];
									}

									partsum -= maxval;

									double surelimit = 0.99;
									bool homo = false;
									if (dous[j]->markerdata[marker].first == dous[j]->markerdata[marker].second)
									{
										surelimit /= 2;
										homo = true;
										partsum = 0;
									}

									//partsum = 0.5 - maxval;

									//								if (maxval > surelimit)
									{
										for (int b = 0; b < 2; b++)
										{
											if ((&(dous[j]->markerdata[marker].first))[b] == UnknownMarkerVal) continue;

											double toadd = dous[j]->parinfprobs[marker][a][c][b][0] - partsum;
											if (toadd < 0) continue;

											dous[j]->pars[a]->sureinfprobs[marker][c][make_pair((&(dous[j]->markerdata[marker].first))[b], (&(dous[j]->markerdata[marker].first))[!b])]
											+= toadd;
										}
									}
								}
							}
						}
					}				
					}
					

					// TODO: NEGSHIFT DOESN'T TAKE RELMAP FLAG2 RESTRICTIONS INTO ACCOUNT
					// Consider doing haplotype reversal from a specific position and all the way down.
					if (HAPLOTYPING && !early && !full && dous[j]->gen >= 1)
					{
						const int NUMTURNS = 1 << (TYPEBITS + 1);
						double rawvals[NUMTURNS][NUMSHIFTS];
						double rawervals[NUMTURNS][NUMSHIFTS];
						double sumnegval[TYPEBITS + 1] = {0};
						for (int g = 0; g < NUMTURNS; g++)
						{
							for (int s = 0; s < NUMSHIFTS; s++)
							{
								rawvals[g][s] = 0;
								rawervals[g][s] = 0;
							}
						}

						for (int g = 0; g < NUMTURNS; g++)
						{		
							if (g & (flag2ignore >> 1)) continue;

							int c = 0;
							for (int p = 0; p < TYPEBITS + 1; p++)
							{
								if (g & (1 << p)) c++;
							}

														if (c > 1) continue;

							aroundturner turn(g);
							for (shiftflagmode = shifts; shiftflagmode < shiftend; shiftflagmode++)
							{
													if (shiftflagmode & shiftignore) continue;
								// If we are above this limit, we are shifting shift moden
								// within the range and cannot use this heuristic of the
								// aggregated probability to know anything... anything at all!
								// (In other cases, our simulated crossover event amounts to 
								// far below a exp(20) change in probability.)
								int g2 = g;
								if (!g) g2 = (1 << 15) - 1;

								int oldshift = shiftflagmode;								
								rawervals[g][oldshift] = exp(dous[j]->doanalyze<aroundturner>(tb, turn, chromstarts[i],
									chromstarts[i + 1] - 1, classicstop(q, -1), -1, true, 0, -5000 + factor) - factor);
								shiftflagmode = oldshift;

								/*								if (c > 1) continue;*/
								if (rawervals[g][oldshift] < 0) continue;

								rawvals[g][oldshift] = rawervals[g][oldshift];

								for (int t = 0; t < TYPEBITS + 1; t++)
								{
																		if (g2 & (1 << t))
									{
										sumnegval[t] += rawvals[g][shiftflagmode];
									}
								}

							}

						}

						for (int t = 0; t < TYPEBITS + 1; t++)
						{
							if (!sumnegval[t]) sumnegval[t] = 1;
						}


#pragma omp critical(negshifts)
						{
							for (int g = 0; g < NUMTURNS; g++)
							{
								for (int s = shifts; s < shiftend; s++)
								{
														if (s & shiftignore) continue;
									int marker = -q - 1000;
									double val = rawvals[g][s];
									//fprintf(out, "rawvals: %d %d %d %d %lf\n", dous[j]->n, marker, g, s, rawervals[g][s] / rawervals[0][s] -1);
									// Consider switching to all-log
									if (!_finite(val) || val < 1e-10) val = 1e-10;

									//									if (_finite(val) && val > 1e-10)
									{

										{

										}
										int g2 = g;										
										if (!g) g2 = (1 << 15) - 1;							      

										// This is hardcoded for the generation count of 3.
										if (NUMGEN == 3)
										{
											dous[j]->negshift[marker] += val * (1.0 - ((g >> 6) & 1) * 2) * ((g2 >> 6) & 1) / sumnegval[6];

											if (dous[j]->pars[0])
												dous[j]->pars[0]->negshift[marker] += val * (1.0 - ((g >> 0) & 1) * 2) * ((g2 >> 0) & 1) / sumnegval[0];

											if (dous[j]->pars[1])
												dous[j]->pars[1]->negshift[marker] += val * (1.0 - ((g >> 3) & 1) * 2) * ((g2 >> 3) & 1) / sumnegval[3];
											if (dous[j]->gen >= 2)
											{
												if (dous[j]->pars[0] && dous[j]->pars[0]->pars[0])
													dous[j]->pars[0]->pars[0]->negshift[marker] += val * (1.0 - ((g >> 1) & 1) * 2) * ((g2 >> 1) & 1) / sumnegval[1] / dous[j]->pars[0]->children;

												if (dous[j]->pars[0] && dous[j]->pars[0]->pars[1])
													dous[j]->pars[0]->pars[1]->negshift[marker] += val * (1.0 - ((g >> 2) & 1) * 2) * ((g2 >> 2) & 1) / sumnegval[2] / dous[j]->pars[0]->children;

												if (dous[j]->pars[1] && dous[j]->pars[1]->pars[0])
													dous[j]->pars[1]->pars[0]->negshift[marker] += val * (1.0 - ((g >> 4) & 1) * 2) * ((g2 >> 4) & 1) / sumnegval[4] / dous[j]->pars[1]->children;

												if (dous[j]->pars[1] && dous[j]->pars[1]->pars[1])
													dous[j]->pars[1]->pars[1]->negshift[marker] += val * (1.0 - ((g >> 5) & 1) * 2) * ((g2 >> 5) & 1) / sumnegval[5] / dous[j]->pars[1]->children;
											}
										}
										else
										{
											dous[j]->negshift[marker] += val * (1.0 - ((g >> 2) & 1) * 2) /* * ((g2 >> 2) & 1) */ / (sumnegval[2]);

											if (dous[j]->pars[0])
												dous[j]->pars[0]->negshift[marker] += log(val) * (1.0 - ((g >> 0) & 1) * 2) * ((g2 >> 0) & 1) /*/ (sumnegval[0])*/;

											if (dous[j]->pars[1])
												dous[j]->pars[1]->negshift[marker] += log(val) * (1.0 - ((g >> 1) & 1) * 2) * ((g2 >> 1) & 1) /*/ (sumnegval[1])*/;
											if (dous[j]->pars[0] && dous[j]->pars[1])
											{
												nsm[make_pair(dous[j]->pars[0], dous[j]->pars[1])][marker][g] +=
													log(val) - log(rawvals[0][s]);
											}

											if (g == 1)
												nsm[make_pair(dous[j]->pars[0], dous[j]->pars[0])][marker][2] +=
												log(val) - log(rawvals[0][s]);

											if (g == 4)
												nsm[make_pair(dous[j]->pars[0], dous[j])][marker][2] +=
												log(val) - log(rawvals[0][s]);
										}
									}
								}
							}
						}
					}					




					double probsum = 0;
					for (int i = 0; i < 3; i++)
					{
						probsum += probs[i];
					}
					probsum = 1 / probsum;
					for (int i = 0; i < 4; i++)
					{
						char string[255];
						int val;
						sprintf(string, "%.5lf%c%n", probs[i] * probsum, i == 3 ? '\n' : '\t', &val);
						for (int k = 0; k < val; k++)
						{
							outqueue[j].push_back(string[k]);
						}
					}



					//						if (oqp[j] > 50000) printf("%d\t%d\n", j, oqp[j]);
					/*						strcpy(&outqueue[j][oqp[j]], lineout);
					oqp[j] += strlen(lineout);*/
					if (!full)
				{
					int marker = -q - 1000;

					// Contribute haplotype data, but truncate it, i.e. a 50/50 contribution for either interpretation is not added.
					// Instead, we have a cap later on at the maximum change at any iteration.
					// critical section outside the loop to maintain symmetry in floating point ops
#pragma omp critical(update)
					{
#pragma ivdep
						for (int k = 0; k < (int) reltree.size(); k++)
						{
							int i = reltree[k]->n;
/*							for (int side = 0; side < 2; side++)
							{
								for (map<MarkerVal, float>::iterator i = infprobs[i][side].begin(); i != infprobs[i][side].end(); i++)
								{
									reltree[k]->infprobs[marker]
								}*/
							if (haplos[i][0] || haplos[i][1])
							{
								float base;

								bool zerobase = false;
								/*for (int z = 0; z < 2; z++)
								{
								if ((&(reltree[k]->markerdata[marker].first))[z] == UnknownMarkerVal)
								{
								haplos[i][z] -= haplos[i][!z];
								if (haplos[i][z] < 0) haplos[i][z] = 0;
								//zerobase = true;
								}
								}*/

								//									if (/*reltree[k] != dous[j] || true*/ !zerobase && fabs(reltree[k]->haploweight[marker] - 0.5) < 0.49999)
								//									{
								//										base = min(haplos[i][0] / reltree[k]->haploweight[marker], haplos[i][1] / (1.0f - reltree[k]->haploweight[marker]));
								//									}
								//									else
								//										base = 0;
								//
								//									if (i == 1633 && marker < 10)
								//									{
								//										fprintf(out, "HAPLOS: %02d %lf %lf %lf %lf %lf %lf\n", marker, (double) reltree[k]->haploweight[marker], (double) haplos[i][0], (double) haplos[i][1], (double) (haplos[i][0] - base * reltree[k]->haploweight[marker]), (double) (haplos[i][1] + haplos[i][0] - base), (double) base);
								//									}
								//#pragma omp critical(update)
								//									{
								//										getind(i)->haplobase[marker] += haplos[i][0] - base * reltree[k]->haploweight[marker];
								//										getind(i)->haplocount[marker] += haplos[i][1] + haplos[i][0] - base;
								//									}


								if (fabs(reltree[k]->haploweight[marker] - 0.5) < 0.49999)
								{
									double b1 =  (haplos[i][0] + maxdiff * maxdiff * 0.5) /*/ reltree[k]->haploweight[marker] /** (1 - reltree[k]->haploweight[marker])*/;
									double b2 = (haplos[i][1] + maxdiff * maxdiff * 0.5) /*/ (1 - reltree[k]->haploweight[marker]) /** reltree[k]->haploweight[marker]*/;

									double intended = (b1 - b2) / min(reltree[k]->haploweight[marker], 1 - reltree[k]->haploweight[marker]);
									//intended -= reltree[k]->haploweight[marker];

									bool neg = intended < 0;

									//intended /= sqrt(fabs(intended) + 1.0);
									// if (neg) intended = -intended;

									{
										reltree[k]->haplobase[marker] += log(b1/b2);
										reltree[k]->haplocount[marker] += 1;
									}
								}
								haplos[i][0] = 0;
								haplos[i][1] = 0;
							}
						}
					}
				}
				}

			}
		}


				for (unsigned int j = 0; j < outqueue.size(); j++)
		{
			outqueue[j].push_back(0);
		if (out && printalot) fprintf(out, "%s:%d\n", dous[j]->name.c_str(), i + 1);
		if (out && printalot) fprintf(out, "%s\n", &outqueue[j].front());
		}
			     
				if (out) fflush(out);

#ifdef F2MPI
		if (false) reduce(world, markerweight, markerweight2, vectorplus<MWTYPE>(), 0);
#endif
		if (false
#ifdef F2MPI
			&& world.rank() == 0
#endif
			&& !full)
		{
			markerweight = markerweight2;
			for (int q = chromstarts[i]; q < chromstarts[i + 1] - 2; q++)
			{
				double dist = markerposes[q + 1] - markerposes[q];
				if (dist > 0)
				{
					for (int t = 0; t < 2; t++)
					{
						/*double prob = markerweight[q][t * 2 + 1] / (markerweight[q][t * 2 + 1] + markerweight[q][t * 2]);
						prob -= 0.5;
						prob *= 2;

						if (prob > 0)
						{
						const double newpart = 0.99;
						prob = log(prob);
						prob /= dist;
						actrec[t][q + 1] = prob * newpart + actrec[t][q + 1] * (1 - newpart);
						}
						else
						{
						fprintf(out, "Strange prob in marker %d : %lf     %lf:%lf\n", q, prob, markerweight[q][t*2], markerweight[q][t*2+1]);
						}*/
						/*				  double prob2 = (markerweight[q][t * 2 + 1] - markerweight[q][t * 2]) / dous.size() / 2;
						prob2 += 1;
						if (prob2 < 0.5) prob2 = 0.5;
						if (prob2 > 3) prob2 = 3;
						actrec[t][q + 1] *= prob2;*/

						double prob2 = - (markerweight[q][t * 2 + 1] - markerweight[q][t * 2]) / dous.size() / dist;
						if (prob2 > 0.3) prob2 = 0.3;
						if (prob2 < -0.3) prob2 = -0.3;

						//				  actrec[t][q + 1] *= 2;
						actrec[t][q + 1] += prob2;
						if (actrec[t][q + 1] > -1e-5) actrec[t][q + 1] = -1e-4;
						if (actrec[t][q + 1] < -20) actrec[t][q + 1] = -10;

						if (q != chromstarts[i] && markerweight[q - 1][4] < dous.size() * 1e-3 && false)
						{
							double sum = actrec[t][q + 1] + actrec[t][q];
							sum /= 2;

							actrec[t][q] = sum;
							actrec[t][q + 1] = sum;
						}
					}
				}
				fprintf(out, "actrec: %0.3lf\t%0.3lf\n", actrec[0][q + 1] * 50, actrec[1][q + 1] * 50);
			}
		}


		if (!full && HAPLOTYPING)
		{
			vector<set<negshiftcand> > negshiftcands;
			negshiftcands.resize(chromstarts.size());

			for (unsigned int i = 0; i < 1000000; i++)
			{
				individ* ind = getind(i);
				if (!ind || !ind->haplocount.size()) continue;

				{
#ifdef F2MPI
					vector<float> templist;
					templist.resize(ind->haplocount.size());

					reduce(world, ind->haplocount, templist, vectorplus<vector<float> >(), 0);
					ind->haplocount = templist;
					world.barrier();

					reduce(world, ind->haplobase, templist, vectorplus<vector<float > >(), 0);
					ind->haplobase = templist;
					world.barrier();

					reduce(world, ind->negshift, templist, vectorplus<vector<float > >(), 0);
					ind->negshift = templist;
					world.barrier();
#endif


					//		fprintf(out, "%d:%d\n", world.rank(), i);
					//		fflush(out);
				}

#ifdef F2MPI
				if (world.rank()) continue;
#endif	      	      		 

				if (/*ind->pars[0] || ind->pars[1] || */!ind->haplocount.size()) continue;		  

				// Perform the inversions indicated by the negshift data, at most a single one per individual
				// and chromosome, maybe 0.
				for (int c = 0; c < (int) chromstarts.size() - 1; c++)
				{
					int minstart = chromstarts[c + 1];
					double minval = -1e-5;
					bool prevlow = false;

					for (int p = chromstarts[c]; p < (int) chromstarts[c + 1]; p++)
					{
						if (ind->negshift[p] < minval)
						{
							minstart = p;
							minval = ind->negshift[p];
						}
						if (ind->negshift[p] < -1e-5)
						{
							if (!prevlow)
							{
								//						fprintf(stdout, "prevlow: %d %d %lf\n", ind->n, p, ind->negshift[p]);
								negshiftcand ourtuple(ind, 0, p);
								negshiftcands[c].insert(ourtuple);
							}
							prevlow = true;
						}
						else
							prevlow = false;
					}
					if (ind->lastinved[c] == minstart)
					{
						ind->lastinved[c] = -1;
						//continue;
					}

					ind->lastinved[c] = -1;

					if (minstart + 1 < chromstarts[c + 1])
					{
						negshiftcand ourtuple(ind, minval, minstart);
						inferiorrelated pred(ourtuple);
						for (set<negshiftcand>::iterator iter = negshiftcands[c].begin(); iter != negshiftcands[c].end(); )
						{
							if (pred(*iter))
							{
								negshiftcands[c].erase(iter++);
							}
							else
							{
								iter++;
							}

						}
						if (!pred.anymatch)
						{
							negshiftcands[c].insert(ourtuple);
						}
					}

				}
			}

#ifdef F2MPI
			if (!world.rank())
#endif
			{

				for (unsigned int i = 0; i < 1000000; i++)
				{
					individ* ind = getind(i, false);
					if (!ind || !ind->haplocount.size()) continue;

					int cno = 0;
					for (unsigned int j = 0; j < ind->haplocount.size(); j++)
					{
						while (cno + 1 < chromstarts.size() && j >= chromstarts[cno + 1]) cno++;

						MarkerVal bestvals[2] = {UnknownMarkerVal, UnknownMarkerVal};
						double bestsure[2] = {ind->markersure[j].first, ind->markersure[j].second};
						bool foundbest = true;
						bool surefound = false;

						map<MarkerVal, double> surenesses;

						if (DOINFPROBS)
						{
							map<MarkerVal, double> sums[2];
							for (int a = 0; a < 2; a++)
							{
								for (map<pair<MarkerVal, MarkerVal>, double>::iterator i = ind->infprobs[j][a].begin(); i != ind->infprobs[j][a].end(); i++)
								{
									sums[a][i->first.second] += i->second;
								}
							}

							sums[0][UnknownMarkerVal] = 1;
							sums[1][UnknownMarkerVal] = 1;

							for (int a = 0; a < 2; a++)
							{
								// Clear old infprobs
								for (int b = 0; b < 2; b++)
								{
									for (int c = 0; c < 2; c++)
									{
										for (int d = 0; d < 2; d++)
										{
											ind->parinfprobs[j][a][b][c][d] = 0;
										}
									}
								}

								double sum = 0;
								MarkerVal bestmarker;
								double bestval = 0;
								double bestval2 = 0;
								map<MarkerVal, double> infprobs;


								for (map<pair<MarkerVal, MarkerVal>, double>::iterator i = ind->infprobs[j][a].begin(); i != ind->infprobs[j][a].end(); i++)
								{
									double factor = sums[a][i->first.second] / (sums[0][i->first.second] + sums[1][i->first.second] + 1e-10);
									fprintf(out, "Factor: %lf %lf %d %d %d %d %lf %lf\n", i->second, factor, ind->n, j, i->first.first, i->first.second, sums[0][i->first.second], sums[1][i->first.second]);
									factor += 1e-10;

									factor = 1;

									sum += i->second / factor;
									infprobs[i->first.first] += i->second / factor;  
								}
								if (sum <= 1e-12)
								{
									sum = 500;
									if ((&(ind->markerdata[j].first))[a] != UnknownMarkerVal)
									{
										bestval = 1 - (&(ind->markersure[j].first))[a];
										bestval2 = 1 - (&(ind->markersure[j].first))[a];
										//									if (bestval2 < 0.99) bestval2 += 0.01;
										bestmarker = (&(ind->markerdata[j].first))[a];
									}
								}

								for (map<MarkerVal, double>::iterator i = infprobs.begin(); i != infprobs.end(); i++)
								{
									//double factor = sums[a][i->first.second] / (sums[0][i->first.second] + sums[1][i->first.second]);
									double factor = 1;

									double sureness = i->second / factor / sum;
									printf("SURENESS: %lf %d\n", sureness, (int)i->first.value());
									// Add extra uncertainty
									/*								double extra = 1.0 - min(0.5 / (sum), 0.5);

									if (sureness > extra)
									{
									sureness -= (sureness - extra) * 0.5;
									}*/
									double origsureness = sureness;

									if ((&(ind->markerdata[j].first))[a] == UnknownMarkerVal)
									{
										/*if ((&(ind->markersure[j].first))[a] > 1e-3)
										{
										sureness -= (1 - (&(ind->markersure[j].first))[a]);
										sureness /= (&(ind->markersure[j].first))[a];
										}*/
									}
									else
									if (i->first == (&(ind->markerdata[j].first))[a])
									{
										/*						  double bigdenom = 1.0 / (((sum - i->second) / ((&(ind->markersure[j].first))[a] + 1e-5)) + i->second);
										origsureness = i->second * bigdenom;
										sureness = origsureness;*/
									}
									else
									{
										/*						  double bigdenom = 1.0 / (i->second / ((&(ind->markersure[j].first))[a] + 1e-5) + sum - i->second);
										sureness = i->second / ((&(ind->markersure[j].first))[a] + 1e-5) * bigdenom;
										//sureness = i->second / ((&(ind->markersure[j].first))[a]) / sum;
										origsureness = sureness;*/

										//sureness = origsureness;
										if (sureness > 0.9999)
										{
											//									fprintf(out, "Was 1: %lf %d %d\n", sureness, ind->n, j);
											sureness = 0.9999;
										}
										//							origsureness = sureness;
										//if (origsureness < 0.5) origsureness = 0.5;
									}

									if (/*ind->unknowninfprobs[j][a] > 0.001 &&*/ origsureness > 0.9999) origsureness = 0.9999;

									/*origsureness -= 0.5;
									origsureness *= 2;*/

									//origsureness = 1 - ((1 - origsureness) / origsureness);


									//									if ((ind->n == 1633 || ind->n == 1726) && j < 600) fprintf(out, "Sureness: %lf %d %d %d %lf\n", sureness, i->first, j, a, sum);
									surenesses[i->first] += 1 - origsureness;

									if (sureness > bestval)
									{
										bestval = sureness;
										bestval2 = origsureness;
										bestmarker = i->first;
									}
								}


								double origsum = sum;
								double sureness = bestval;

								if (sureness > ((iter % 30 == 19 && false) ? 0.99 : 0.49))
								{
									bestvals[a] = bestmarker;
									bestsure[a] = 1 - bestval2;
								}
								else
								{
									printf("Foundbest now false, with sureness %lf, marker %d, pair-half %d for ind %d\n", sureness, j, a, ind->n);
									foundbest = false;
								}								
							}
						}

						if (DOINFPROBS)
						{
							//if (!ind->sex) foundbest = false;

							if (!foundbest && ind->markerdata[j].first != UnknownMarkerVal && ind->markerdata[j].second != UnknownMarkerVal && (bestvals[0] != UnknownMarkerVal || bestvals[1] != UnknownMarkerVal))
							{
								foundbest = true;
								bestvals[0] = UnknownMarkerVal;
								bestvals[1] = UnknownMarkerVal;

								bestsure[0] = 0;
								bestsure[1] = 0;
							}

							if (ind->lastinved[cno] == -1 || true)
							{
								if (foundbest)
								{
									/*bestsure[0] *= 0.99;
									bestsure[1] *= 0.99;*/
									if (bestvals[0] == bestvals[1] && bestvals[0] != UnknownMarkerVal)
									{
										if (bestsure[0] + bestsure[1] > 0.5 && false)
										{
											double bsum = bestsure[0] + bestsure[1];

											int index = bestsure[0] > bestsure[1];
											bestsure[index] /= 2;

											//if (bestsure[index] < 0.01) bestsure[index] = 0.01;

											MarkerVal bestwhat = UnknownMarkerVal;
											double bestsury = 4;

											for (map<MarkerVal, double>::iterator i = surenesses.begin(); i != surenesses.end(); i++)
											{
												if (i->first == bestvals[0]) continue;

												if (i->second < bestsury)
												{
													bestsury = i->second;
													bestwhat = i->first;
												}
											}

											if (bestwhat == UnknownMarkerVal) bestsury = 0;
											/*else
											{
											bestsury /= 2;
											bestsury = 1 / (1 + bestsury);
											bestsury *= 2;
											bestsury = (1 - bestsury) / bestsury;
											}*/

											if (bestsury < 0 || bestsury > 1) fprintf(out, "Bestsury problem %d %d %lf\n", ind->n, j, bestsury);

											bestsure[!index] = (2 - bestsury) / 2;
											bestvals[!index] = bestwhat;
										}

										double sum = (bestsure[0] + bestsure[1]) / 2;
										double diff = fabs(bestsure[0] / sum - 1) + fabs(bestsure[1] / sum - 1);

										/*if (diff < 0.05)
										{
										bestsure[0] = sum * 0.9;
										bestsure[1] = sum * 1.1;
										}*/

										/*bestsure[0] += 0.03;
										bestsure[1] -= 0.03;
										bestsure[0] = min(0.99999, bestsure[0]);
										bestsure[1] = max(0.00001, bestsure[1]);*/
										/*int index = bestsure[0] > bestsure[1];
										float sumsure = (1 - bestsure[0]) + (1 - bestsure[1]);

										if (sumsure < 1.4)
										{
										bestsure[index] /= 2;
										bestsure[!index] = 0;
										bestvals[!index] = UnknownMarkerVal;
										}*/

										/*bestsure[index] = sumsure;
										if (bestsure[index] > 0.999) bestsure[index] = 0.999;

										sumsure -= bestsure[index];
										bestsure[!index] = sumsure;

										bestsure[0] = 1 - bestsure[0];
										bestsure[1] = 1 - bestsure[1];*/
									}
									if (fabs(0.5 - ind->haploweight[j]) == 0.5)
									{							 
										int min = 0;
										if (bestsure[1] < bestsure[0]) min = 1;
										if (bestvals[0] == bestvals[1] && bestvals[0] != UnknownMarkerVal)
										{
											bestsure[min] = 0;
										}
										else
										{
											bestsure[!min] /*-= bestsure[min]*/ = 0;
											bestsure[min] = 0;
										}
									}
									/*				  if (bestvals[0] == bestvals[1] && bestvals[0] != UnknownMarkerVal)
									{
									double delta = -0.1;

									if (ind->haploweight[j] > 0.5) delta = -delta;

									ind->haploweight[j] += delta;

									if (ind->haploweight[j] < 0.001) ind->haploweight[j] = 0.001;
									if (ind->haploweight[j] > 0.999) ind->haploweight[j] = 0.999;
									}*/

									if ((bestvals[0] == ind->markerdata[j].second || bestvals[1] == ind->markerdata[j].first) && bestvals[0] != bestvals[1])
									{
										/*swap(bestvals[0], bestvals[1]);
										swap(bestsure[0], bestsure[1]);*/
										//ind->haploweight[j] = 1 - ind->haploweight[j];
									}
									bool nochange = true;
									for (int i = 0; i < 2; i++)
									{
										if ((&(ind->markerdata[j].first))[i] != bestvals[i]) nochange = false;
									}

									for (int i = 0; i < 2; i++)
									{
										double old = (&(ind->markersure[j].first))[i];

										if (bestvals[i] == UnknownMarkerVal)
										{
											bestsure[i] = 0;
											continue;
										}



										if ((&(ind->markerdata[j].first))[i] != bestvals[i])
										{
											if (bestvals[i] == bestvals[!i] && (&(ind->markerdata[j].first))[i] == UnknownMarkerVal) old = 0.5;
											else
												old = 1 - old;
										}

										if (old == 0)
										{
											bestsure[i] = 0;
										}

										/*old *= 2;
										bestsure[i] *= old;

										if (bestsure[i] < old * 0.35) bestsure[i] = old * 0.35;*/

										//if (bestsure[i] <= old * 1.021) bestsure[i] *= 0.98;
										//								double factor = fabs(ind->haploweight[j] - 0.5) + 0.01;
										//								bestsure[i] *= factor;
										bestsure[i] += old/* * (0.51 - factor)*/;
										bestsure[i] /= 2 /* 0.51 */;


										/*bestsure[i] *= bestsure[i];
										bestsure[i] *= 2;*/
									}

									/*							if (nochange && bestvals[0] != bestvals[1])
									{
									bestsure[0] *= 0.99;
									bestsure[1] *= 0.99;
									}*/

									if (nochange && bestvals[0] == bestvals[1] && bestvals[0] != UnknownMarkerVal)
									{
										double ratio = bestsure[0] / bestsure[1];
										if (ratio < 1) ratio = 1 / ratio;

										if (ratio < 1.0001)
										{
											double val = bestsure[0] * 0.01;
											bestsure[0] -= val;
											bestsure[1] += val;
											if (bestsure[0] > 0.01)
											{
												fprintf(out, "Bringing aside %d %d %lf\n", ind->n, j, bestsure[0]);
											}
										}
									}

									/*							if (bestvals[0] != ind->markerdata[j].first || bestvals[1] != ind->markerdata[j].second || bestsure[0] + bestsure[1] + ind->markersure[j].first + ind->markersure[j].second > 0.01 || bestvals[0] == UnknownMarkerVal || ind->n == 1633)
									{
									fprintf(out, "Individual %d stochastic fix at marker %d, %d:%d, was %d:%d %c %lf %lf\n", ind->n, j, bestvals[0], bestvals[1], ind->markerdata[j].first, ind->markerdata[j].second, surefound ? 'S' : 'I', bestsure[0], bestsure[1]);
									}*/

									ind->markerdata[j] = make_pair(bestvals[0], bestvals[1]);
									ind->markersure[j] = make_pair(bestsure[0], bestsure[1]);
								}
								/*else if (ind->markerdata[j].first == ind->markerdata[j].second && ind->markerdata[j].first == UnknownMarkerVal)
								{				    
								for (int a = 0; a < 2 && !foundbest; a++)
								{
								map<MarkerVal, double> sums;
								double sum = 0;

								for (map<pair<MarkerVal>, double>::iterator i = ind->infprobs[j][a].begin(); i != ind->infprobs[j][a].end(); i++)
								{
								sums[i->first] += i->second;
								sum += i->second;
								}
								if (sum < 0.1) sum = 0.1;

								for (map<pair<MarkerVal>, double>::iterator i = sums.begin(); i != sums.end() && !foundbest; i++)
								{
								//							  fprintf(out, "Individual %d stochastic tryout at marker %d, :%d:, was %d:%d %lf M\n", ind->n, j, i->first, i->first, ind->markerdata[j].second, i->second / sum);
								if (i->second >= sum * 0.5)
								{
								ind->markerdata[j] = make_pair(i->first, UnknownMarkerVal);
								ind->markersure[j] = make_pair(0.999, 0);
								fprintf(out, "Individual %d stochastic fix at marker %d, :%d:, was %d:%d %lf M\n", ind->n, j, i->first, ind->markerdata[j].first, ind->markerdata[j].second, i->second / sum);
								foundbest = true;
								}
								}
								}					  
								}*/
							}
							for (int a = 0; a < 2; a++)
							{
								ind->sureinfprobs[j][a].clear();
								ind->infprobs[j][a].clear();
								ind->unknowninfprobs[j][a] = 0;
							}
						}

						{
							vector<bool> allhalf;
							vector<bool> anyinfo;
							vector<bool> cleared;
							vector<int> nudgeme;


							anyinfo.resize(chromstarts.size());
							allhalf.resize(chromstarts.size());
							cleared.resize(chromstarts.size());
							nudgeme.resize(chromstarts.size());
							for (int k = 0; k < chromstarts.size(); k++)
							{
								anyinfo[k] = false;
								allhalf[k] = true;
								cleared[k] = false;
								nudgeme[k] = -1;
							}

							cno = 0;
							for (unsigned int j = 0; j < ind->haplocount.size(); j++)
							{

								anyinfo[cno] = true;

								if (!(ind->haploweight[j] && ind->haploweight[j] != 1) && false)
								{
									double b1 = ind->haplobase[j];
									double b2 = ind->haplocount[j] - ind->haplobase[j];
									/*if (!allhalf[cno])
									{
									ind->haploweight[j] = 0.5 - (0.5 - ind->haploweight[j]) * 0.01;
									}*/

									if (!early && fabs(ind->negshift[j]) < 1e-5)
									{
										fprintf(out, "Clearing: %d %lf %lf %lf\n", ind->n, b1, b2, ind->negshift[j]);
										//lockhaplos(ind, cno);
										cleared[cno] = true;
										ind->haploweight[j] = 0.5 /*- (0.5 - ind->haploweight[j]) * 0.99*/;
										ind->haplocount[j] = 0;
									}
									else
										allhalf[cno] = false;
								}


								if (ind->haplocount[j] && ind->haploweight[j] && ind->haploweight[j] != 1)
								{							
									/*double b1 = ind->haplobase[j];
									double b2 = ind->haplocount[j] - ind->haplobase[j];							*/
									double b1 = 1;
									double b2 = 1;

									/*b1 /= ind->haploweight[j];
									b2 /= (1.0 - ind->haploweight[j]);*/

									/*if (ind->markerdata[j].first == UnknownMarkerVal) b1 /= 1.5;
									if (ind->markerdata[j].second == UnknownMarkerVal) b1 *= 1.5;*/

									/*b1 += 1e-10;
									b2 += 1e-10;*/
									double val = exp(ind->haplobase[j] / ind->haplocount[j]);
									val *= (1 - ind->haploweight[j]) / ind->haploweight[j];


									double intended = exp(log(val) * 0.1 + log(ind->haploweight[j] / (1 - ind->haploweight[j])));
									intended = intended / (intended + 1.0);

									if (fabs(b1 + b2) < 1e-4)
									{
										intended = 0.500;
										//if (ind->haploweight[j] == 0.5) ind->negshift[j] -= 1e-4;
									}

									if (!early && allhalf[cno] && fabs(intended - 0.5) > 0.1 &&
										ind->markerdata[j].first != UnknownMarkerVal && ind->markerdata[j].second != UnknownMarkerVal &&
										cleared[cno])
									{
										allhalf[cno] = false;
										fprintf(out, "Locking: %d %d %lf\n", ind->n, j, ind->negshift[j]);
										ind->haploweight[j] = (intended < 0.5) ? 0 : 1;
									}
									else
										//							if ((j == 5 || j == 6 || j == 4) && (ind->n == 2334)) fprintf(out, "B1B2: %d %d %lf %lf\n", ind->n, j, b1, b2);

										//printf("---\n");
										/*for (int q = 0; q < 20; q++)
										{
										double b12 = b1 * intended;
										double b22 = b2 * (1.0 - intended);
										intended = b12 / (b12 + b22);
										printf("%lf\n", intended);
										}*/

										//				double intended = ind->haplobase[j] / ind->haplocount[j];
										//double intended2 = intended * (1.0 - ind->haploweight[j]) +
										//(1.0 - intended) * ind->haploweight[j];

										//if (early)
									{
										// Cap the change if the net difference is small/miniscule
										/*double nnn = 1 + 0.5 * (b1 + b2);						*/
										double nnn = 1.6;
										//								if (nnn > 1.3) nnn = 1.3;
										/*if (ind->markerdata[j].first == UnknownMarkerVal ||
										ind->markerdata[j].second == UnknownMarkerVal) nnn = (nnn - 1.0) / 5 + 1.0;*/
										if (nnn < 1.0) nnn = 1.0;

										double limn = (nnn - 1.0) * ind->haploweight[j] * (-1 + ind->haploweight[j]);

										double limd1 = -1 - (nnn - 1.0) * ind->haploweight[j];
										double limd2 = (nnn - 1.0) * ind->haploweight[j] - nnn;

										double lim = min(limn/limd1, limn/limd2);

										double diff = intended - ind->haploweight[j];
										/*			if (ind->haploweight[j] > 0.5 && i == 2827)
										{
										fprintf(stderr, "2827: %lf %lf %lf %lf %lf %lf %lf\n", limn, limd1, limd2, lim, intended, b1, b2);
										}*/									

										/*if (fabs(diff) > lim)
										{
										intended = ind->haploweight[j] + diff / fabs(diff) * lim;
										}*/

										if (diff > limn/limd1)
										{
											intended = ind->haploweight[j] + limn/limd1;
										}

										if (diff < -limn/limd2)
										{
											intended = ind->haploweight[j] - limn/limd2;
										}

										//								if ((ind->haploweight[j] - 0.5) * (intended - 0.5) < 0) intended = 0.5;

										intended = min((float) intended, 1.0f - maxdiff);
										if ((ind->lastinved[cno] == -1 || true) /*&& !ind->pars[0] && !ind->pars[1]*/)
										{
											ind->haploweight[j] = max((float) intended, maxdiff);

											if ((nudgeme[cno] == -1 || fabs(ind->haploweight[nudgeme[cno]] - 0.5) < fabs(ind->haploweight[j] - 0.5)) && ind->haploweight[j] > maxdiff && ind->haploweight[j] < 1 - maxdiff &&
												fabs(b1 * 2 - b1 + b2) < 0.001 )
											{
												nudgeme[cno] = j;
											}
										}							
									}

									/*							if (ind->haploweight[j] != 0.5)
									{
									allhalf[cno] = false;
									}*/
								}
								else
								{
									if (ind->haploweight[j] && ind->haploweight[j] != 0.5)
									{
										/*								fprintf(stderr, "H c w: %d %d\n", i, j);
										fflush(stderr);*/
									}
								}


							}

							//										if (ind->n < 2327)
							for (int k = 0; k < chromstarts.size() - 1; k++)
							{
								/*if (nudgeme[k] != -1 && ind->haploweight[nudgeme[k]] != 0.5)
								{
								fprintf(out, "Nudging %04d:%02d\n", ind->n, nudgeme[k]);
								ind->haploweight[nudgeme[k]] += ind->haploweight[nudgeme[k]] < 0.5 ? (-0.10) : 0.10;
								ind->haploweight[nudgeme[k]] = min(ind->haploweight[nudgeme[k]], (float) (1.0f - maxdiff));
								ind->haploweight[nudgeme[k]] = max((float) maxdiff, ind->haploweight[nudgeme[k]]);
								}*/
								/*					  						if (allhalf[k] && anyinfo[k])
								{
								fprintf(out, "Locking haplos %d %d\n", ind->n, k);
								lockhaplos(ind, k);
								}*/
							}
						}
					}
					vector<pair<double, boost::tuple<individ*, individ*, int, int> > > allnegshifts;
					map<individ*, double> bestshift;
					for (map<pair<individ*, individ*>, map<int, boost::array<double, 8> > >::iterator i = nsm.begin(); i != nsm.end(); i++)
					{
						for (map<int, boost::array<double, 8> >::iterator j = i->second.begin(); j != i->second.end(); j++)
						{
							// 1-3 allows shifts, but not genotype switches
							// 2 only allows shifts for paren 2, e.g. assumption that paren 1 is part of several half sibships
							for (int k = 1; k <= 4; k++)
							{
								/*		      		      */
								if (k == 2/* || k == 4*/)
									allnegshifts.push_back(make_pair(j->second[k] - 1e-5, boost::make_tuple(i->first.first, i->first.second, k, j->first)));
								/*		      if (k == 4) printf("ANS: %lf %d %d %d %d\n",
								j->second[k], i->first.first->n, i->first.second->n, k, j->first);*/

								/*		      if (k == 4)
								allnegshifts.push_back(make_pair(j->second[k] - 1e-5, make_tuple(i->first.first, i->first.first, 1, j->first)));*/
							}
						}
					}

					sort(allnegshifts.begin(), allnegshifts.end());

					for (int k = allnegshifts.size() - 1; k >= 0; k--)
					{
						if (bestshift[allnegshifts[k].second.get<0>()] < allnegshifts[k].first &&
							bestshift[allnegshifts[k].second.get<1>()] < allnegshifts[k].first)
						{
							bestshift[allnegshifts[k].second.get<0>()] = allnegshifts[k].first;
							bestshift[allnegshifts[k].second.get<1>()] = allnegshifts[k].first;
							int c = 0;
							while (c < chromstarts.size() && (chromstarts[c] <= allnegshifts[k].second.get<3>()))
							{
								c++;
							}

							int phase = allnegshifts[k].second.get<2>();
							individ* inds[2] = {allnegshifts[k].second.get<0>(), allnegshifts[k].second.get<1>()};
							if (rand() > RAND_MAX / 10) continue;

							printf("Inv: %d %d %d %d %lf\n", inds[0]->n, inds[1]->n, phase, allnegshifts[k].second.get<3>(), allnegshifts[k].first);
							for (int m = allnegshifts[k].second.get<3>() + 1; m < chromstarts[c]; m++)
							{
								for (int z = 0; z < 2; z++)
								{			  
									if (phase & 4)
									{
										// These days, we are just emulating the shifts
										if (z == 0 && false)
										{
											swap(inds[0]->haploweight[m], inds[1]->haploweight[m]);
											swap(inds[0]->markerdata[m], inds[1]->markerdata[m]);
											swap(inds[0]->markersure[m], inds[1]->markersure[m]);
										}

										for (int k = 0; k < inds[z]->kids.size(); k++)
										{
											// ONLY inverting those that share both parents
											// This used to be that we would invert every child, but avoid inverting those that share both parents
											// twice, i.e. not inverting them at all
											if (z && inds[z]->kids[k]->pars[0] == inds[0]) {
												//printf("I2: %d %d\n", inds[z]->kids[k]->n, m);
												inds[z]->kids[k]->haploweight[m] = 1.0 - inds[z]->kids[k]->haploweight[m];
											}
										}
									}

									if (phase & (1 << z))
									{
										//			      printf("Inv %d at %d, was %lf\n", inds[z]->n, m, inds[z]->haploweight[m]);
										inds[z]->haploweight[m] = 1.0 - inds[z]->haploweight[m];  
									}

								}
							}
						}
					}




					for (int c = 0; c < (int) chromstarts.size() - 1; c++)
					{
						for_each(negshiftcands[c].begin(), negshiftcands[c].end(), negshifter(c));
					}
				}
			}
		}
	}

	/*	FILE* binut = fopen("form_data/set009/ad_dataset_009_chr_01.bin", "wb");
	for (int pos = (int) markerposes[chromstarts[0]]; pos != (int) markerposes[chromstarts[1] - 1] + 1; pos++)
	{
	for (int half = 0; half < 2; half++)
	{
	for (int ind = 0; ind < dous.size(); ind++)
	{
	fwrite(&realgeno[ind][pos][half], sizeof(float), 1, binut);
	}
	}
	}*/

	generation++;
}

// Handle the fuss of trailing \r, \n characters when combining scanf and gets and stuff.
void clean(char* tlf)
{
	int i = strlen(tlf);
	while (i && tlf[i - 1] < 32)
	{
		tlf[i - 1] = 0;
		i--;
	}
}

void readhaplodata(int innum)
{
	char tlf[255];
	sprintf(tlf, "qtlmas14haplo_%d", innum);
	FILE* in = fopen(tlf, "r");
	bool inind = false;
	bool inact = false;
	individ* ind = 0;
	int mcc = 0;
	int mnum;

	while (fgets(tlf, 255, in))
	{
		float w;
		int a, b;
		float n;

		if (sscanf(tlf, "%f %d %d %f", &w, &a, &b, &n) == 4)
		{
			if (ind)
			{
				ind->haploweight[mnum++] = w;
			}
		}
		else
		{
			ind = 0;
			float a1, a2;

			if (sscanf(tlf, "%d", &a) == 1)
			{
				ind = getind(a);
				printf("Reading ind %d, was at marker %d\n", a, mnum);
				mnum = innum;
			}
			else
				if (sscanf(tlf, "actrec: %f %f", &a1, &a2) == 2)
				{
					if (!inact)
					{
						inact = true;
						mnum = innum;
					}

					mnum++;
					actrec[0][mnum] = a1 / 10000.0;
					actrec[1][mnum] = a2 / 10000.0;
					printf("Reading actrec %d : %f, %f\n", mnum, a1, a2);
				}
				else
					inact = false;
		}
	}
}

// Map individual names to inds
map<string, individ*> indmap;
individ* getind(string name)
{
	while (indmap.find(name) == indmap.end())
	{
		int origindex = indmap.size();	  
	  individ* ind = indmap[name] = getind(origindex, true);
	  if (ind) ind->name= name;
	}

	return indmap[name];
}


void readalphaped(FILE* in)
{
	getind("0");
	char me[255], father[255], mother[255];
	while (fscanf(in, "%s %s %s", me, father, mother) == 3)
	{	
		char line[255];
		fgets(line, 255, in);
		int gen = 0;
		if (sscanf(line, "%d", &gen) != 1) gen = 0;

		// We are all empty until we get some data, not implemented in getind yet due
		// to possible regressions.
		individ* ime = getind(me);	       
		individ* ifounderf = getind(father);
		individ* ifounderm = getind(mother);
		
		if (ime) ime->empty = true;
		if (ifounderf) ifounderf->empty = true;
		if (ifounderm) ifounderm->empty = true;

		if (gen >= 2 && !ifounderf->gen && !ifounderm->gen)
		{
			individ* realpars[2] = {getind(me + (string) "_aux_realf"), getind(me + (string) "_aux_realm")};
			for (int k = 0; k < 2; k++)
			{
				ime->pars[k] = realpars[k];
				realpars[k]->gen = 1;
				realpars[k]->pars[0] = ifounderf;
				realpars[k]->pars[1] = ifounderm;
				realpars[k]->empty = true;
			}
			ime->gen = gen;
		}
		else
		{
			ime->gen = gen;
			ime->pars[0] = ifounderf;
			ime->pars[1] = ifounderm;
		}

		if (gen >= 2)
		{
		  dous.push_back(ime);
		}
	}      
}

void readalphadata(FILE* in)
{
	char me[255];
	while (fscanf(in, "%s", me) == 1)
	{
		individ* ime = getind(me);
		ime->empty = false;
		int data;
		for (int x = 0; x < markerposes.size(); x++)
		{
			char datastr[255];
			fscanf(in, "%s", datastr);
			int data2;
			int numread = sscanf(datastr, "%d/%d", &data, &data2);
			ime->haploweight[x] = 0.5;
			if (numread == 1)
			{
				pair<MarkerVal, MarkerVal> marker;
				switch (data) 
				{
				case 0:
					marker = make_pair(1 * MarkerValue, 1 * MarkerValue);
					break;
				case 1:
					marker = make_pair(1 * MarkerValue, 2 * MarkerValue);
					break;
				case 2:
					marker = make_pair(2 * MarkerValue, 2 * MarkerValue);
					break;
				default:
					marker = make_pair(UnknownMarkerVal, UnknownMarkerVal);
					break;
				}
				ime->markerdata[x] = marker;
				ime->markersure[x] = make_pair(0, 0);
			}
			else
			{
			  printf("%s %d\n", datastr, numread);
				if (data == data2 && !data)
				{
					ime->markerdata[x] = make_pair(UnknownMarkerVal, UnknownMarkerVal);
					ime->markersure[x] = make_pair(0, 0);
				}
				else
				{
					boost::math::binomial_distribution<double> dist1(data);
					boost::math::binomial_distribution<double> dist2(data2);
					double sure1 = 0;
					double sure2 = 0;
					double probsum = 0;

					for (int rl1 = 0; rl1 <= data; rl1++)
					{
						for (int rl2 = 0; rl2 <= data2; rl2++)
						{
							int l1 = rl1;
							int l2 = rl2;
							double overallprob = (data ? boost::math::pdf(dist1, l1) : 1) * (data2 ? boost::math::pdf(dist2, l2) : 1);													

							double sureb1;
							double sureb2;
							while (true)
							{
								sureb1 = 0.5;
								sureb2 = 0.5;
								if (l1 + l2)
									sureb1 = l1 / (double) (l1 + l2);
								
								if (data + data2 - l1 - l2) 
									sureb2 = (data2 - l2) / (double) (data + data2 - l1 - l2);

								if (sureb1 + 1e-9 > 1 - sureb2) break;
								l1 = data - l1;
								l2 = data2 - l2;
							}

							overallprob *= pow(sureb1, l1) * pow(1 - sureb1, l2) * pow(sureb2, (data2 - l2)) * pow(1 - sureb2, data - l1);

							sure1 += sureb1 * overallprob;
							sure2 += sureb2 * overallprob;
							probsum += overallprob;
						}
					}

					sure1 /= probsum;
					sure2 /= probsum;

					pair<MarkerVal, MarkerVal> marker = make_pair(2 * MarkerValue, 1 * MarkerValue);
					pair<double, double> markersure = make_pair(sure1, sure2);

					// Oh, error rate over 0.5, let's invert to make it more passable
					for (int k = 0; k < 2; k++)
					{
						if ((&(markersure.first))[k] > 0.5)
						{
							(&(markersure.first))[k] = 1 - (&(markersure.first))[k];
							(&(marker.first))[k] = (k + 1) * MarkerValue;
						}
					}

					ime->markerdata[x] = marker;
					ime->markersure[x] = markersure;
					printf("%d/%d turned into %d %d with %lf;%lf\n", data, data2, marker.first, marker.second, markersure.first, markersure.second);
				}
			}
		}
	}
}

void readalphamap(FILE* in)
{
	double prevval = 1e30;
	double value;
	int nowmarker = 0;
	while (fscanf(in, "%lf", &value) == 1)
	{
		if (value < prevval)
		{
			chromstarts.push_back(nowmarker);
		}
		nowmarker++;
		markerposes.push_back(value);
		prevval = value;
	}
	chromstarts.push_back(markerposes.size());
}

void readmerlinmap(FILE* mapfile)
{
	// Assuming sex-averaged map
	// Assuming chromosomes and markers to be presented in order
	int chromo;
	int oldchromo = -1;
	double cmpos; // All positions will be shifted making first marker occurs at 0
	double cmbase = -1;
	int bppos;
	char markername[256];
	while (fscanf(mapfile, "%d %255s %lf %d", &chromo, markername, &cmpos, &bppos) == 4)
	{
		chromo--;
		if (chromo != oldchromo)
		{
			chromstarts.push_back(markerposes.size());
			oldchromo = chromo;
			cmbase = cmpos;
		}

		markerposes.push_back(cmpos - cmbase);
		for (int t = 0; t < 2; t++)
	    {
	      actrec[t].push_back(baserec[t]);
	    }

	}

	chromstarts.push_back(markerposes.size());
}

void readmerlinped(FILE* pedfile)
{
	indmap["0"] = 0;
	char famname[256];
	char indname[256];
	char pname[256];
	char mname[256];
	int sex;
	double pheno;
	while (fscanf(pedfile, "%255s %255s %255s %255s %d %lf", famname, indname, pname, mname, &sex, &pheno) == 6)
	{
	  printf("Reading main line for individual %s\n", indname);
		individ* ind = getind(indname);
		ind->pars[0] = getind(pname);
		ind->pars[1] = getind(mname);
		ind->gen = 0;
		ind->sex = --sex;

		if (ind->pars[0] || ind->pars[1])
		{
			dous.push_back(ind);
			ind->gen++;
		}

			ind->markerdata.resize(markerposes.size());
      ind->haplobase.resize(markerposes.size());
      ind->haplocount.resize(markerposes.size());
      ind->haploweight.resize(markerposes.size());
      ind->negshift.resize(markerposes.size());
      
      ind->infprobs.resize(markerposes.size());
      ind->sureinfprobs.resize(markerposes.size());
      ind->unknowninfprobs.resize(markerposes.size());
      ind->parinfprobs.resize(markerposes.size());
      ind->markersure.resize(markerposes.size());
			//		ind->semishift.resize(5000);
      ind->lastinved.resize(chromstarts.size());
      ind->lockstart.resize(chromstarts.size());
      
      for (int i = 0; i < chromstarts.size(); i++)
		{
			ind->lastinved[i] = -1;
			ind->lockstart[i] = 0;
		}
      
	  int a, b;
      for (int k = 0; k < markerposes.size(); k++)
		{
		  ind->haploweight[k] = 0.5;
		  fscanf(pedfile, "%d %d", &a, &b);
	  
		  // Assume 1e-7 allele error rate
			  ind->markersure[k] = make_pair(a ? 1e-7 : 0.0, b ? 1e-7 : 0.0);
	      
		  ind->markerdata[k] = make_pair(a * MarkerValue, b * MarkerValue);

		}
	}
}

int family = 1;

void domerlinind(FILE* pedfile, individ* ind)
{
  int pn[2] = {0};

  for (int j = 0; j < 2; j++)
    {
      if (ind->pars[j]) pn[j] = ind->pars[j]->n;
    }

  fprintf(pedfile, "%d\t%d\t%d\t%d\t%d", family, ind->n, pn[0], pn[1], ind->sex + 1);
  for (int k = 0; k < chromstarts[1]; k++)
    {
      fprintf(pedfile, "\t%d\t%d", ind->markerdata[k].first.value(), ind->markerdata[k].second.value());
    }

  fprintf(pedfile, "\n");
}



void readhaplodata(FILE* in, int swap)
{
  char tlf[255];
	bool inind = false;
	bool inact = false;
	individ* ind = 0;
	int mnum = 0;
	int mcc = 0;
	int ncc = 0;
	
	int mixes[] = {
	  0, 0,
	  1215, 1739,
	  1241, 2254,
	  1264, 1886,
	  1314, 1905,
	  1332, 1485,
	  1351, 1430,
	  1391, 1758,
	  1457, 1576,
	  1425, 1860,
	  1482, 1823,
	  1527, 1564,
	  1538, 1855,
	  1552, 2259,
	  1590, 2057,
	  1694, 2101,
	  1772, 1931,
	  1773, 2096,
	  1805, 2020,
	  1849, 1900,
	  1876, 1935,
	  1877, 2076,
	  1936, 2029,
	  1953, 2288,
	  2012, 2233,
	  2034, 2293,
	  2118, 2237,		       		       
	  0, 0};

	printf("Reading haplos\n");
	while (fgets(tlf, 255, in))
	{
		float w;
		int a, b;
		float n;

		if (sscanf(tlf, "%f %d %d %f", &w, &a, &b, &n) == 4)
		{
			if (ind)
			{
			  int i = mnum;
			  if ((!((a == ind->markerdata[i].first.value() && b == ind->markerdata[i].second.value()) ||
				 (b == ind->markerdata[i].first.value() && a == ind->markerdata[i].second.value()))) && (a || b))
			  {
			    printf("Mismatch: %d %d\t%d %d %d %d\n", ind->n, mnum, a, b, ind->markerdata[i].first.value(), ind->markerdata[i].second.value());
			    mcc++;
			  }
			  if (!(a || b)) ncc++;

				ind->haploweight[mnum++] = w;
			}
		}
		else
		  {
			float a1, a2;

			if (sscanf(tlf, "%d", &a) == 1)
			{
			  int i = 0;
			  int n = a;
			  /*      if (n <= 220 && swap)
	{
	  if (n > 20)
	{
	  n -= 21;
	  n /= 10;
	  n++;
	}
      else
	{
	  n--;
	  n *= 10;
	  n += 21;
	}
	}*/
      a = n;
			  while (mixes[i])
			    {
			      for (int j = 0; j < 2; j++)
				{
				  if (a == mixes[i + j])
				    {
				      a = mixes[i + (!j)];
				      break;
				    }
				}
			      i += 2;
			    }

			  if (mnum)
			    {
			      printf("%d Ind %d, %d mismatches (%d)\n", swap, ind->n, mcc, ncc);
			    }

				ind = getind(a);
				printf("Reading ind %d, was at marker %d\n", a, mnum);
				if (a > 220) ind = 0;

				mnum = 0;
				mcc = 0;
				ncc = 0;
			}
			else
				if (sscanf(tlf, "actrec: %f %f", &a1, &a2) == 2)
				{
					if (!inact)
					{
						inact = true;
						mnum = 0;
					}

					mnum++;
					actrec[0][mnum] = a1 / 10000.0;
					actrec[1][mnum] = a2 / 10000.0;
					printf("Reading actrec %d : %f, %f\n", mnum, a1, a2);
				}
				else
					inact = false;
		}
	}
}

int main(int argc, char* argv[])
{
#ifdef _MSC_VER
	// Performant?
	// This turns off some aspects of IEEE precision, but we do not really need it, as we avoid denormalized values
	// anyway.
	_controlfp(_DN_FLUSH, _MCW_DN);
	_set_printf_count_output(1);
#endif
#ifdef F2MPI
	mpi::environment env(argc, argv);
	mpi::communicator world;
#endif


	//	scanf("%lf", &discstep);

	printf("Number of sexes in map: ");	
	//	scanf("%d", &sexc);
	discstep = 1;
	sexc = 2;

	// Not really related to the number of generations, but doing it like this makes the
	// estimates similar for the two cases. Whether the old 4-state, 2-gen model was
	// actually correct is another issue entirely.
	//
	// / 50 everywhere is probably more appropriate
	if (NUMGEN == 3)
	{
		baserec[0] = -discstep / 50.0;
	}
	else
	{
		baserec[0] = -discstep / 50.0;
	}
	baserec[1] = baserec[0];
	for (int gen = 0; gen < 2; gen++)
	{
		genrec[gen] = baserec[0];
	}
	//selfgen = 1;
	genrec[2] = baserec[0] /** selfgen*/;

	char tlf[255];
	//	fgets(tlf, 255, stdin);

	/*	printf("Marker info file: ");
	fgets(tlf, 255, stdin);
	clean(tlf);
	FILE* in = fopen(tlf, "r");
	readmarkerinfo(in);
	fclose(in);

	printf("Pedigree file: ");
	fgets(tlf, 255, stdin);
	clean(tlf);
	in = fopen(tlf, "r");
	readped(in);
	fclose(in);

	printf("Marker data file: ");
	fgets(tlf, 255, stdin);
	clean(tlf);
	in = fopen(tlf, "r");
	readmarkerdata(in);
	fclose(in);	*/

	if (argc < 4)
	  {
	    printf("Three args expected: map, ped files in Merlin format, followed by output filename.\n");
	    return -1;
	  }

	FILE* mapfile = fopen(argv[1], "rt");
	readmerlinmap(mapfile);
	FILE* pedfile = fopen(argv[2], "rt");
	readmerlinped(pedfile);
	//	return 0;
	CORRECTIONINFERENCE = false;
	postmarkerdata();
	CORRECTIONINFERENCE = false;
	int chromnum;

	/*	sscanf(argv[1], "%d", &chromnum);
	chromstarts[0] = chromstarts[chromnum - 1];
	chromstarts[1] = chromstarts[chromnum];

	chromstarts.resize(2);
	dous.resize(2400);*/
	//sprintf(tlf, "/glob/nettel/qtlmas14_2gen.%d", world.rank());
	//	sprintf(tlf, "qtlmas14aug03");
	/*	long long seed;
	scanf("%lld", &seed);
	rng.seed((boost::mt19937::result_type) seed);

	in = fopen("sexcorrhapend", "r");
	readhaploweights(in);
	fclose(in);

	postmarkerdata();*/

	/*	readhaplodata(0);
	readhaplodata(1976);
	readhaplodata(4034);
	readhaplodata(6082);
	readhaplodata(8073);*/


	//	sscanf(argv[5], "%d", &chromstarts[1]);

	FILE* out = fopen(argv[3], "w");
	int COUNT = 1;
	if (argc == 5) sscanf(argv[4], "%d", &COUNT);

	if (HAPLOTYPING || true)
		for (int i = 0; i < COUNT; i++)
		{
			//		  	  	{
			early = (i < 1);
			doit<false>((i == COUNT - 1) ? out : stdout, /*i == COUNT - 1*/ true
#ifdef F2MPI
				, world
#endif
				);
			//		}

			//	doit<true>(out);
			//		if (i == COUNT - 1)
			{
				//		    fprintf(out, "\n");
				for (unsigned int j = 0; j < chromstarts[1]; j++)
				{
					for (int t = 0; t < 5; t++)
					{
						//			    fprintf(out, "%0.5lf\t", markerweight[j][t]);
						markerweight[j][t] = 0;
					}
					//			fprintf(out, "\n");		
				}
				//		    fprintf(out, "\n");
			}
			fflush(stdout);
			fflush(out);

			for (unsigned int i2 = 0; i2 < 1000000; i2++)
			{
				individ* ind = getind(i2);
				if (!ind) continue;

				if (ind->haplocount.size())
				{
#ifdef F2MPI
					if (!world.rank())
#endif
					  /*if (i == COUNT - 1)*/						fprintf(stdout, "%d\n", i2);
					// Printing of haplotype data for each iteration
					for (unsigned int c = 0; c < chromstarts.size() - 1; c++)

					{
						for (unsigned int j = chromstarts[c]; j < chromstarts[c + 1]; j++)
						{

#ifdef F2MPI
							if (!world.rank())
#endif
							  /*if (i == COUNT - 1)*/ fprintf(stdout, "%f\t%d\t%d\t\t%f\t%lf %lf\n", ind->haploweight[j], ind->markerdata[j].first.value(), ind->markerdata[j].second.value(), ind->negshift[j],
									    ind->markersure[j].first, ind->markersure[j].second);
							ind->negshift[j] = 0;
						}
						//if (!world.rank()) fprintf(out, "\n");
					}
				}
			}
			fflush(stdout);
			fflush(out);
		}

		return 0;
}

