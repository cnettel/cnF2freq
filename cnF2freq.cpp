// cnF2freq, (c) Carl Nettelblad, Department of Information Technology, Uppsala University
// 2008-2019
//
// PlantImpute 1.5, with support for forward-backward and old "true" tree-style cnF2freq
// algorithm. This reduces mamoery requirements and speeds up imputation. A lot. Still
// som e kinks to work out. Please contact the author regarding proper references.
//
// carl.nettelblad@it.uu.se
//
// This code is allowed to be freely used for any commercial or research purpose. If the code is
// used or integrated into another project largely unchanged, attribution to the original author
// should be included. No warranties are given.


#define NDEBUG
#define BOOST_NO_EXCEPTIONS
#undef __cpp_exceptions
#define __cpp_exceptions 0
#include <exception>
#include <cstdlib>
namespace boost
{
	void throw_exception( std::exception const & e )
	{std::abort();}
}

#include <vector>
#include <string.h>
#include <stdio.h>
#include <omp.h>
#include <limits>

const int ANALYZE_FLAG_FORWARD = 0;
const int ANALYZE_FLAG_BACKWARD = 16;
const int ANALYZE_FLAG_STORE = 32;

float templgeno[8] = { -1, -0.5,
	0,  0.5,
	0,  0.5,
	1, -0.5 };

const int NO_EQUIVALENCE = -1;
const int ZERO_PROPAGATE = 1;
const long long WEIGHT_DISCRETIZER = 1000000;

#include <array>
#include <ranges>

#include <boost/math/distributions/binomial.hpp>
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
#include <memory>
#include <string>

#include <boost/iostreams/device/mapped_file.hpp>
#include <boost/container/flat_map.hpp>
#include <boost/container/flat_set.hpp>
#include <boost/container/small_vector.hpp>
#include <boost/container/static_vector.hpp>
#include <boost/program_options.hpp>
#include <boost/interprocess/file_mapping.hpp>
#include <boost/interprocess/mapped_region.hpp>
#include <boost/spirit/home/x3.hpp>
#include <boost/spirit/home/x3/support/ast/variant.hpp>
#include <boost/fusion/container/vector.hpp>
#include <boost/fusion/sequence.hpp>
#include <boost/fusion/include/vector.hpp>
#include <boost/fusion/container/vector/vector_fwd.hpp>
#include <boost/fusion/include/vector_fwd.hpp>

#include <boost/fusion/include/adapt_struct.hpp>
#include <boost/fusion/include/io.hpp>
#include <boost/fusion/adapted/std_pair.hpp>
#include <boost/fusion/include/std_pair.hpp>
#include <boost/fusion/adapted/std_tuple.hpp>
#include <boost/fusion/include/std_tuple.hpp>
#include <boost/spirit/include/support_istream_iterator.hpp>
#include <boost/fusion/include/at_c.hpp>
//#include <boost/numeric/odeint.hpp>
#include <boost/math/quadrature/gauss_kronrod.hpp>
#include <iostream>
#include <fstream>

#if LIBSTATGEN
#include <VcfRecordGenotype.h>
#include <VcfFileReader.h>
#include <VcfFileWriter.h>
#endif

#include <vector>

using namespace boost::iostreams;
using namespace boost::spirit;
namespace po = boost::program_options;
//namespace ode = boost::numeric::odeint;

#include <errno.h>
#include <assert.h>
#include <stdlib.h>
#include <unordered_set>
#include <set>
#include <algorithm>
#include <math.h>
#include <type_traits>
#include <map>
#include <float.h>
#include <numeric>
#include <type_traits>
#include <spawn.h>
#include <sys/types.h>
#include <sys/wait.h>


using namespace std; // use functions that are part of the standard library

//using namespace boost::mpi;
using namespace boost::random;
using boost::container::static_vector;
using boost::container::small_vector;
using boost::container::flat_map;
using boost::container::flat_set;

#define none cnF2freqNONE


#ifndef _MSC_VER
#define _isnan isnan
#define _finite isfinite
#endif

#include "settings.h"

#if XSTDBITSET
#include <xstd/bit_set.hpp>
#endif

#define __assume(cond) do { if (!(cond)) __builtin_unreachable(); } while (0)
#define restrict __restrict__
EXTERNFORGCC boost::random::mt19937 rng;

#pragma omp threadprivate(rng)

#ifdef DOEXTERNFORGCC
boost::random::mt19937 rng;
#endif

struct spectype
{

};
spectype globspec;

int myrand(int max)
{
  boost::uniform_int<> orig(0, max - 1);

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

struct MarkerValueType {};

constexpr MarkerValueType MarkerValue;

struct MarkerVal
{
private:
	int val;

public:
	constexpr explicit MarkerVal(int val) : val(val) {}
	constexpr MarkerVal() : val(0) {}
	constexpr int value() const
	{
		return val;
	}

	constexpr bool operator == (const MarkerVal& rhs) const
	{
		return val == rhs.val;
	}

	constexpr bool operator != (const MarkerVal& rhs) const
	{
		return val != rhs.val;
	}

	constexpr bool operator < (const MarkerVal& rhs) const
	{
		return val < rhs.val;
	}
};

constexpr const MarkerVal operator* (const int& val, const MarkerValueType& rhs)
{
	return MarkerVal(val);
}

typedef pair<MarkerVal, MarkerVal> MarkerValPair;


const MarkerVal UnknownMarkerVal = (MarkerVal)0;
const MarkerVal sexmarkerval = 9 * MarkerValue;

const float maxdiff = 0.000005;


bool early = false;

vector<double> markerposes;;
vector<double> actrec[2];
map<string, int> markernames;

//int selfgen = 0;
double genrec[3];
vector<unsigned int> chromstarts;

vector<int> markertranslation;
typedef vector<array<double, 5 > > MWTYPE;

template<class T> class vectorplus;

template<> class vectorplus<MWTYPE >
{
public:
	MWTYPE operator () (const MWTYPE& l, const MWTYPE& r)
	{
		MWTYPE result;
		result.resize(l.size());

		for (size_t i = 0; i < l.size(); i++)
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

		for (size_t i = 0; i < l.size(); i++)
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
template<int zeropropagate> bool markermiss(MarkerVal& a, const MarkerVal b)
{
	// This is the real logic; we do not bind anything at all when zeropropagate is true
	if (zeropropagate == ZERO_PROPAGATE) return false;

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
int upflagit(int flag, int parnum, unsigned int genwidth)
{

	if (flag < 0) return flag;
	flag >>= parnum * (genwidth - 1);
	flag &= ((1 << (genwidth - 1)) - 1);

	return flag;
}

struct individ;

std::string tmppath{ "." };

int generation = 1;
int shiftflagmode;
int lockpos[NUMSHIFTS];
int quickmark[NUMSHIFTS];
int quickgen[NUMSHIFTS];

template<class T2> class PerStateArray
{
public:
	typedef array<T2, NUMTYPES> T;
};

template<class T2> class StateToStateMatrix
{
public:
	typedef typename PerStateArray<
		typename PerStateArray<T2>::T >
		::T T;
};

#if !DOFB
// The quick prefixes are caches that retain the last invocation.
// Only used fully when we stop at marker positions exactly, i.e. not a general grid search.
double quickfactor[NUMSHIFTS];
PerStateArray<int>::T quickendmarker[NUMSHIFTS];
PerStateArray<double>::T quickendfactor[NUMSHIFTS];
StateToStateMatrix<double>::T quickendprobs[NUMSHIFTS];
PerStateArray<double>::T quickmem[NUMSHIFTS];
#endif

small_vector<std::pair<MarkerVal, double>, 5> testlista;

template<class K, class T, int N = 2> using small_map = flat_map<K, T, less<K>, small_vector<std::pair<K, T>, N>>;

// A hashed store of inheritance pathway branches that are known to be impossible.
// Since we can track the two branches that make up the state in the F_2 individual independently,
// this optimization can reduce part of the cost by sqrt(number of states).
typedef array<array<array<array<array<array<array<int, 4>, HALFNUMSHIFTS>, HALFNUMPATHS + 1>, HALFNUMTYPES>, 2>, 2>, 2> IAT;

EXTERNFORGCC IAT impossible;

// A memory structure storing haplo information for later update.
// By keeping essentially thread-independent copies, no critical sections have to
// be acquired during the updates.
EXTERNFORGCC array<array<double, 2>, INDCOUNT> haplos;
EXTERNFORGCC array<array<small_map<MarkerVal, double>, 2>, INDCOUNT> infprobs;

#if !DOFB
// done, factors and cacheprobs all keep track of the same data
// done indicates that a specific index (in the binary tree of blocks of multi-step transitions) is done
// with a "generation id" that's semi-unique, meaning no active clearing of the data structure is performed
EXTERNFORGCC vector<int> done[NUMSHIFTS];
// factors contain the mantissas of the extended floating-point representation
EXTERNFORGCC vector<PerStateArray<double>::T > factors[NUMSHIFTS];
// cacheprobs contain actual transitions from every possible state to every possible other state
EXTERNFORGCC vector<StateToStateMatrix<double>::T > cacheprobs[NUMSHIFTS];
#else
EXTERNFORGCC vector<array<PerStateArray<double>::T, 3> > fwbw[NUMSHIFTS];
EXTERNFORGCC vector<array<double, 3> > fwbwfactors[NUMSHIFTS];
int fwbwdone[NUMSHIFTS];
#endif

EXTERNFORGCC vector<individ*> reltree;
EXTERNFORGCC vector<individ*> reltreeordered;
EXTERNFORGCC flat_map<individ*, int> relmap;
EXTERNFORGCC flat_map<individ*, int> relmapshift;

//#pragma omp threadprivate(realdone, realfactors, realcacheprobs)
#pragma omp threadprivate(generation, shiftflagmode, impossible, haplos, lockpos, reltree, reltreeordered, relmap, relmapshift, infprobs)
#if !DOFB
#pragma omp threadprivate(quickmark, quickgen, quickmem, quickfactor, quickendfactor, quickendprobs, done, factors, cacheprobs)
#else
#pragma omp threadprivate(fwbw,fwbwfactors,fwbwdone)
#endif

#ifdef DOEXTERNFORGCC
IAT impossible;
array<array<double, 2>, INDCOUNT> haplos;
vector<PerStateArray<double>::T > factors[NUMSHIFTS];
vector<individ*> reltree;
vector<individ*> reltreeordered;
flat_map<individ*, int> relmap; //containing flag2 indices
flat_map<individ*, int> relmapshift; //containing flag2 indices
std::array<std::array<small_map<MarkerVal, double>, 2>, INDCOUNT> infprobs;
#if !DOFB
vector<int> done[NUMSHIFTS];
vector<StateToStateMatrix<double>::T > cacheprobs[NUMSHIFTS];
#else
vector<std::array<PerStateArray<double>::T, 3> > fwbw[NUMSHIFTS];
vector<std::array<double, 3> > fwbwfactors[NUMSHIFTS];
#endif
#endif




// We put all thread-local structures in our own separate struct. This is because many compilers implementing OpenMP
// use a relatively expensive call for determining thread-local global data, which is not cached over function calls.
// The simple pointer-arithmetic for a lookup within a struct is cheaper.
struct threadblock
{
	int* const generation;
	int* const shiftflagmode;
#if !DOFB
	int* const quickmark;
	int* const quickgen;
#endif
	int* const lockpos;
#if !DOFB
	double* const quickfactor;
	PerStateArray<double>::T* const quickmem;
	PerStateArray<int>::T* const quickendmarker;
#endif
	IAT* const impossible;
	std::array<std::array<double, 2>, INDCOUNT>* const haplos;
	std::array<std::array<small_map<MarkerVal, double>, 2>, INDCOUNT>* infprobs;
#if !DOFB
	vector<int>* const done;
	vector<PerStateArray<double>::T >* const factors;
	vector<StateToStateMatrix<double>::T >* const cacheprobs;
	PerStateArray<double>::T* const quickendfactor;
	StateToStateMatrix<double>::T* const quickendprobs;
#else
	vector<std::array<PerStateArray<double>::T, 3> >* fwbw;
	vector<std::array<double, 3> >* fwbwfactors;
	int* fwbwdone;
#endif	

	threadblock() : generation(&::generation), shiftflagmode(&::shiftflagmode), impossible(&::impossible),
		haplos(&::haplos), lockpos(::lockpos), infprobs(&::infprobs),
#if !DOFB
		done(::done), factors(::factors), cacheprobs(::cacheprobs),
		quickmark(::quickmark), quickgen(::quickgen), quickmem(::quickmem),
		quickfactor(::quickfactor), quickendfactor(::quickendfactor), quickendprobs(::quickendprobs),
		quickendmarker(::quickendmarker),
#else
		fwbw(::fwbw), fwbwfactors(::fwbwfactors), fwbwdone(::fwbwdone)
#endif
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
			// Mis-guided (?) attempt to correct for relskews
			if (RELSKEWSTATES && false)
			{
				this->turn |= (flagmodeshift & 1) << BITS_W_SELF;
			}
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
	const double updateval;
	int* const gstr;

	trackpossibleparams() : updateval(0.0), gstr(0)
	{
	}

	trackpossibleparams(double updateval, int* gstr, MarkerVal markerval = UnknownMarkerVal) : updateval(updateval), gstr(gstr)
	{
	}
} tpdefault;



struct classicstop
{
	int lockpos;
	int genotype;

	classicstop(int lockpos, int genotype) : lockpos(lockpos), genotype(genotype)
	{}

	explicit operator int() const
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
	typedef std::array<std::array<float, 2>, 2> miniactrecT;
	miniactrecT actrec;

	stopmodpair(int lockpos, miniactrecT actrec) : tssmcommon(lockpos), actrec(actrec)
	{
		genotype[0] = -1;
		genotype[1] = -1;
	}
};

template<class G> const bool canquickend(int startmark, const G& stopdata)
{
	return false;
}

/*template<> const bool canquickend<smnonecommon>(int startmark, const smnonecommon &stopdata)
{
return false;
}*/

template<> const bool canquickend<classicstop>(int startmark, const classicstop& stopdata)
{
	return stopdata.lockpos <= -1000 && stopdata.genotype >= 0;
}

template<> const bool canquickend<twicestop>(int startmark, const twicestop& stopdata)
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

	if (index == 0 || index == 1) return stopdata.actrec[k][index];

	return actrec[k][j];
}

const int HAPLOS = 1;
const int GENOS = 2;
const int HOMOZYGOUS = 4;
const int GENOSPROBE = 8;

struct clause
{
	long long weight;
	static_vector<int, 8> cinds;
	//vector<individ*> individuals;

	string toString() const {
		string s = boost::lexical_cast<std::string>(weight);
		//std::stringstream ints;
		//std::copy(cInds.begin(), cInds.end(), std::ostream_iterator<int>(ints, " "));
		for (size_t i = 0; i < cinds.size(); i++) {
			if (cinds[i]) {
				s = s + " " + boost::lexical_cast<std::string>(cinds[i]);
			}
		}
		//s = s + " " + ints.str();
		return s;// note each line ends in a space
	}

	stringstream& clausetostringstream() const {
		static thread_local stringstream s;
		s.str("");
		s.clear();
		for (size_t i = 0; i < cinds.size(); i++) {
			if (cinds[i]) {
			  s << " " << std::to_string(cinds[i]);
			}
		}
		return s;
	}

	string weighttostring() const {
	  return std::to_string(weight);
	}

};

template<int skipper, int update, int zeropropagate>
struct vectortrackpossible
{
	const threadblock& tb;
	MarkerVal inmarkerval[TURNBITS][NUMPATHS / skipper];
	double secondval[TURNBITS][NUMPATHS / skipper];
	int geno[TURNBITS][NUMPATHS / skipper]; // Used to be flag
	int flag2[TURNBITS][NUMPATHS / skipper]; // Used to be flag99 and flag2 everywhere else
	int localshift[TURNBITS][NUMPATHS / skipper];
	double result[NUMPATHS / skipper];
	const trackpossibleparams& extparams;
	const unsigned int marker;
	int meid;
	int menow;

	vectortrackpossible(unsigned int marker, const trackpossibleparams& extparams, int flag2);
};

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
	int descendants;
	// Am I explicitly created to not contain genotypes?
	bool empty;
	// Do I lack (relevant) parents?
	bool founder;
	// Parents.
	individ* pars[2];
	//
	vector<individ*> kids;
	// Sex.
	bool sex;
	// Line or strain of origin, should only exist in founders.
	int strain;
	// Marker data as a list of pairs. No specific ordering assumed.
	vector<MarkerValPair> markerdata;
	vector<double> variances;
	vector<pair<double, double>> markersure;

	vector<MarkerValPair> priormarkerdata;
	vector<pair<double, double>> priormarkersure;
	vector<array<double, 3>> priorgenotypes;

	// Temporary storage of all possible marker values, used in fixparents.
	vector<small_map<MarkerVal, pair<int, double> > > markervals;
	// The haplotype weight, or skewness. Introducing an actual ordering of the value in markerdata.
	vector<double> haploweight;
	// Relative skewness, i.e. shifts between adjacent markers.
	vector<double> relhaplo;
	// The cost-benefit value of inverting the haplotype assignment from an arbitrary marker point on.
	vector<double> negshift;
	vector<int> lastinved;
	vector<unsigned int> lockstart;

	vector<array<small_map<MarkerVal, double>, 2> > infprobs;
	vector<array<double, 2>> homozyg;

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
		// False is the safe option, other code will be take shortcuts if founder is true.
		founder = false;
		descendants = 0;
	}

	bool arerelated(individ* b, vector<individ*> stack = vector<individ*>(), int gens = 0)
	{
		if (gens > 2) return false;
		if (!this) abort();

		if (b == this) return true;
		if (find(stack.begin(), stack.end(), this) != stack.end())
		{
			return false;
		}

		stack.push_back(this);
		if (stack.size() == 1)
		{
			if (b && b->arerelated(this, stack, gens + 1)) return true;
		}

		if (pars[0] && pars[0]->arerelated(b, stack, gens + 1)) return true;
		if (pars[1] && pars[1]->arerelated(b, stack, gens + 1)) return true;

		for (individ* i : kids)
		{
			for (individ* j : b->kids)
			{
				if (i == j) return true;
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
	template<int update, int zeropropagate, class extra=const spectype, int templgenwidth = -2, int templflag2 = -2> struct recursetrackpossible
	{
		individ* mother;
		int upflagr;
		int upflag2r;
		int upshiftr;
		double secondval;
		const trackpossibleparams& extparams;
		const unsigned int genwidth;
		const MarkerVal markerval;
		const int marker;
		const threadblock& tb;
		int firstpar;
		int* impossibleref;
		int impossibleval;
		bool prelok;
		extra& data;

		recursetrackpossible(individ* mother, const threadblock& tb, MarkerVal markerval, double secondval, int marker, int upflag, int upflag2, int upshift, unsigned int genwidth, int f2n, int firstpar, int numrealpar, const trackpossibleparams& extparams, extra& data) noexcept :
			genwidth(genwidth), extparams(extparams), markerval(markerval), marker(marker), tb(tb), firstpar(firstpar), mother(mother), secondval(secondval), data(data)

		{
			if constexpr (templgenwidth >= 0) __assume(genwidth == templgenwidth);
			if constexpr (templflag2 > -2)
			{
				if constexpr (templflag2 == -1) __assume(upflag2 == -1);
				else
				__assume(upflag2 > -1);
			}
			upflagr = upflagit(upflag, firstpar, genwidth);
			upflag2r = upflagit(upflag2, firstpar, genwidth >> (NUMGEN - NUMFLAG2GEN));
			upshiftr = upflagit(upshift, firstpar, genwidth >> (NUMGEN - NUMSHIFTGEN));
			prelok = true;

			
			if constexpr (DOIMPOSSIBLE)
			{
			if (zeropropagate != ZERO_PROPAGATE && genwidth == (1 << (NUMGEN - 1)))
			{
				if (DOIMPOSSIBLE)
				{
					impossibleref = &(*tb.impossible)[*(tb.shiftflagmode) & 1][firstpar][f2n][upflagr][upflag2r + 1][upshiftr][marker & 3];
					impossibleval = (*tb.generation) * markerposes.size() + marker;
				}
				else
				{
					impossibleref = 0;
				}

				if (DOIMPOSSIBLE && impossibleref && *impossibleref == impossibleval)
				{
					prelok = false;
				}
			}
		}
		}

		double compute(const spectype& data) noexcept
		{
			if (templgenwidth >= 0) __assume(genwidth == templgenwidth);
			return mother->pars[firstpar]->trackpossible<update, zeropropagate, const spectype, (templgenwidth >> 1), templflag2>(tb, markerval, secondval, marker,
					upflagr,
					upflag2r,
					upshiftr, extparams, genwidth >> 1);
		}

		template<int skipper, int update2>
		double compute(vectortrackpossible<skipper, update2, zeropropagate>& data) noexcept
			{
			int meid2 = data.meid + 1 + (genwidth >> 1) * firstpar;

			data.inmarkerval[meid2][data.menow] = markerval;
			data.secondval[meid2][data.menow] = secondval;
			data.geno[meid2][data.menow] = upflagr;
			data.flag2[meid2][data.menow] = upflag2r;
			data.localshift[meid2][data.menow] = upshiftr;

			return 1.0;
		}

		operator double() noexcept
		{
			if (templgenwidth >= 0) __assume(genwidth == templgenwidth);
			if constexpr (DOIMPOSSIBLE) if (!prelok)
			{
				return 0;
			}

			if (!mother || !mother->pars[firstpar])
			{
				return 1 + secondval;
			}

			double baseval = compute(data);


			if (DOIMPOSSIBLE && impossibleref && zeropropagate != ZERO_PROPAGATE && !update && genwidth == (1 << (NUMGEN - 1)) && !baseval)
			{
				*impossibleref = impossibleval;
			}

			return baseval;
		}
	};

	// zeropropagate also implies that there is a gstr value, we propagate zeros to find any possible source strain
	// The main logic of tracking a specific inheritance pathway, computing haplotype weights, and overall feasibility.
	// update: Should haplotype weights be updated?
	// zeropropagate: Should zero marker values be kept, or changed to the actual values they are matched against. Is it ok
		// treat equivalent values as one.
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
	template<int update, int zeropropagate, class extra=const spectype, int templgenwidth=-2, int templflag2=-2> double trackpossible(const threadblock& tb, MarkerVal inmarkerval, double secondval, const unsigned int marker,
		const unsigned int flag, const int flag99, int localshift = 0, const trackpossibleparams& extparams = tpdefault,
		const unsigned int genwidth = 1 << (NUMGEN - 1), extra& data = globspec) noexcept /*const*/
	{
		if constexpr (templgenwidth == -2)
		{
			// TODO NICER FORWARDING
			if (genwidth == 1 << (NUMGEN - 1))
			{
				return trackpossible<update, zeropropagate, extra, (1 << (NUMGEN - 1))>(tb, inmarkerval, secondval, marker, flag, flag99, localshift, extparams, genwidth, globspec);
			}
			else
			{
				return trackpossible<update, zeropropagate, extra, -1>(tb, inmarkerval, secondval, marker, flag, flag99, localshift, extparams, genwidth, globspec);
			}
		}
		if constexpr (templflag2 == -2)
		{
			// TODO NICER FORWARDING
			if (flag99 >= 0)
			{
				return trackpossible<update, zeropropagate, extra, templgenwidth, 0>(tb, inmarkerval, secondval, marker, flag, flag99, localshift, extparams, genwidth, globspec);
			}
			else
			{
				return trackpossible<update, zeropropagate, extra, templgenwidth, -1>(tb, inmarkerval, secondval, marker, flag, flag99, localshift, extparams, genwidth, globspec);
			}
		}
		// This used to be a nice silent null check. Compilers don't like that, so we abort in order to try to detect those cases.
		//if (this == NULL) abort();
		if (templgenwidth >= 0) __assume(genwidth == templgenwidth);
		__assume(localshift >= 0);

		if constexpr (templflag2 > -2)
		{
			if constexpr (templflag2 == -1) __assume(flag99 == -1);
			else
			__assume(flag99 > -1);
		}

		// TYPEBITS are the ordinary TYPEBITS. Anything set beyond those indicates selfing. g is multiplied by 2 to become flag, hence TYPEBITS + 1
		const bool rootgen = (genwidth == (1 << (NUMGEN - 1)));
		bool selfingNOW = false;
		bool relskewingNOW = false;

		const bool attopnow = !(update & HOMOZYGOUS) && ((genwidth == HAPLOTYPING) || founder);

		const int selfval = (flag >> (TYPEBITS + 1))& SELFMASK;

		if (rootgen)
		{
			selfingNOW = SELFING && selfval;
			relskewingNOW = RELSKEWSTATES;
		}

		MarkerVal selfmarker[2];
		double selfsure[2];

		int upflag2 = -1;
		const int upflag = flag >> 1;
		const int upshift = localshift >> 1;
		int f2s = 0;
		const MarkerVal* restrict const themarker = selfingNOW ? selfmarker : &markerdata[marker].first;
		const double* const themarkersure = selfingNOW ? selfsure : &markersure[marker].first;
		int f2end = 2;

		if (flag99 != -1 && genwidth >> (NUMGEN - NUMFLAG2GEN) > 0)
		{
			upflag2 = flag99 >> 1;
			f2s = (flag99 & 1);
			f2end = (flag99 & 1) + 1;
		}

		if (relskewingNOW)
		{
			const int relskewval = (flag >> (BITS_W_SELF + 1));
			//std::cout << relskewval << std::endl;
			f2s = max(f2s, relskewval);
			f2end = min(f2end, relskewval + 1);
		}

		int firstpar = flag & 1;
		double ok = 0;

		//#pragma ivdep
		// This ivdep is quite complicated, since we actually change markerval, but not in
		// conditions where this alters the results of the other iteration.

		// flag2 determines which value in the tuple to check. flag and localshift determine the haplotype weight value
		// assigned to that test, and which parent to match against that value.
#pragma ivdep		
		for (int flag2 = f2s; flag2 < f2end && (HAPLOTYPING || !ok); flag2++)
		{
			int f2n = (flag2 & 1);
			if constexpr (SELFING) if (selfingNOW)
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
				baseval = themarkersure[f2n];
				if (themarkersure[f2n] && secondval) mainsecondval = (1.0 - themarkersure[f2n]) * secondval;
			}
			else
			{
				double effectivesecondval = (inmarkerval == UnknownMarkerVal && markerval != UnknownMarkerVal) ? 1 : secondval;
				baseval = 1.0 - themarkersure[f2n];

				double effectivemarkersure = (themarker[f2n] == UnknownMarkerVal ? 1 : themarkersure[f2n]);
				mainsecondval = effectivemarkersure * effectivesecondval;
			}

			// Include it all in one big thing, because it doesn't matter, or we want to enforce the genotype
			if (attopnow || (update & (GENOS || GENOSPROBE)))
			{
				baseval += mainsecondval;
				mainsecondval = 0;
			}
			else
			{
				if (mainsecondval) mainsecondval /= baseval;
			}

			//if (!baseval) continue;
			bool doupdatehaplo = true;

			// Normalize, in some sense.
			f2n ^= ((firstpar ^ localshift) & 1);

			if (zeropropagate || !genwidth)
			{
				baseval *= 0.5;
				doupdatehaplo = false;
			}
			//			else if (/*!empty &&*/ (allthesame && (CORRECTIONINFERENCE) || (themarker[0] == UnknownMarkerVal && themarker[1] == UnknownMarkerVal && themarkersure[0] + themarkersure[1] == 0)))
			else if (/*!empty &&*/ ((!relskewingNOW) && allthesame && ((CORRECTIONINFERENCE) || (themarkersure[0] == themarkersure[1]))) || selfingNOW)
			{
				baseval *= ((f2n) ? 1.0 : 0.0);
				doupdatehaplo = false;
			}
			else
			{
				if (HAPLOTYPING)
				{
					// TODO interaction relskew/selfing
					baseval *= fabs((f2n ? 1.0 : 0.0) - haploweight[marker]);
				}
				else
				{
					// No haplotype weights, all interpretations allowed.
					baseval *= 0.5;
				}
			}

			/*if (!baseval)
			{
				continue;
			}*/

			// If we are at maximum depth, by depth limit or by lack of ancestor
			if (baseval && (attopnow || !pars[firstpar]))
			{
				// TODO: If pars[firstpar] exists and is empty, things are messy
				// The empty one could, in turn, have parents with genotypes.
				if (zeropropagate && extparams.gstr)
				{
					*(extparams.gstr) += (themarker[realf2n] == (2 * MarkerValue));
				}
			}

			// There should be some other flag for the actual search depth
			if (!baseval || attopnow)
			{
			}
			else
			{
				// Track to the parent generation, creating evaluation objects first
				// These do a lookup in a special hash for combinations known to be 0, to avoid unnecessary calls
				// Both are checked before either branch of the pedigree is actually traced.
				auto subtrack1 =
					recursetrackpossible<update & ~HOMOZYGOUS, zeropropagate, extra, templgenwidth, templflag2>(this, tb, markerval, mainsecondval, marker,
						upflag,
						upflag2,
						upshift,
						genwidth,
						f2n,
						firstpar,
						0,
						extparams,
						data);

				if (subtrack1.prelok && (!zeropropagate || rootgen) && !(update & GENOS))
				{
					double secsecondval = 0;
					MarkerVal secmark = themarker[!realf2n];

					if (!(update & HOMOZYGOUS))
					{
						if (themarkersure[!realf2n])
						{
							baseval *= (1 - themarkersure[!realf2n]);
							secsecondval = themarkersure[!realf2n] / (1 - themarkersure[!realf2n]);
						}
					}
					else
					{
						if (markerval != secmark)
						{
							if (secmark != UnknownMarkerVal)
							{
								baseval *= themarkersure[!realf2n];
							}
							secmark = markerval;
							secsecondval = 0;
						}
						else
						{
							baseval *= (1 - themarkersure[!realf2n]);
							secsecondval = 0;
						}
					}

					baseval *= recursetrackpossible<update & ~HOMOZYGOUS, zeropropagate, extra, templgenwidth, templflag2>(this, tb, secmark, secsecondval,
						marker,
						upflag,
						upflag2,
						upshift,
						genwidth,
						f2n,
						!firstpar,
						1,
						extparams,
						data);

					if (selfingNOW && extparams.gstr && !(selfval & (1 << !firstpar))) *extparams.gstr = 0;
				}
				if (baseval)
				{
				if (!(selfingNOW && extparams.gstr && !(selfval & (1 << firstpar))))
					baseval *= subtrack1;
			}
			}

			if (baseval)
			{
			ok += baseval;

			if ((update & HAPLOS) /*&& !allthesame*/ && doupdatehaplo)
			{
				(*tb.haplos)[n][f2n] += extparams.updateval;
			}
			if (update & GENOS)
			{
				(*tb.infprobs)[n][realf2n][markerval] += extparams.updateval;
			}
		}
		}
		if (selfingNOW && extparams.gstr) *extparams.gstr *= 2;
		return ok;
	}

	template<int skipper, int update, int zeropropagate> void trackpossiblevector(vectortrackpossible<skipper, update, zeropropagate>& data)
	{
		int meid = data.meid;
		for (int geno = 0, i = 0; geno < NUMTYPES * 2; geno += skipper, i++)
		{
			data.menow = i;
			int flag2 = data.flag2[meid][i];
			__assume(flag2 >= 0 && flag2 < NUMPATHS);
			//if (meid == 0)
				data.result[i] *= trackpossible<update, zeropropagate>(data.tb, data.inmarkerval[meid][i], data.secondval[meid][i], data.marker, data.geno[meid][i],
					flag2, data.localshift[meid][i], data.extparams, GENWIDTHS[meid], data);
			/*else
				data.result[i] *= trackpossible<update & ~HOMOZYGOUS, zeropropagate>(data.tb, data.inmarkerval[meid][i], data.secondval[meid][i], data.marker, data.geno[meid][i],
					data.flag2[meid][i], data.localshift[meid][i], data.extparams, GENWIDTHS[meid], data);*/
		}
	}


	// calltrackpossible is a slight wrapper that hides at least some of the internal parameters needed for the recursion from outside callers
	template<int update, bool zeropropagate> double calltrackpossible(const threadblock& tb, const MarkerVal* const markervals, const unsigned int marker,
		const int genotype, const unsigned int offset, const int flag2, const double updateval = 0.0)
	{
		return trackpossible<update, zeropropagate, const spectype, 1 << (NUMGEN - 1)>(tb, UnknownMarkerVal, 0,
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
		individ* parp[3] = { pars[0], pars[1], this };

#pragma omp critical (parmarkerval)
		for (int i = 0; i < 3; i++)
		{
			if (parp[i])
			{
				while (parp[i]->markervals.size() <= marker)
				{
					parp[i]->markervals.resize(markerdata.size());
				}
			}
		}

		double okvals[2] = { 0 };
		// We only accept an interpretation when it is by exclusion the only possible one. As soon as one interpretation has gained acceptance,
		// no need to retest it.
		for (shiftflagmode = 0; shiftflagmode < 1; shiftflagmode += 1)
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
			fprintf(stderr, "Clearing %d:%d (was %d,%d)\n", this->n, marker, markerdata[marker].first.value(), markerdata[marker].second.value());
			markerdata[marker] = make_pair(UnknownMarkerVal, UnknownMarkerVal);
			markersure[marker] = make_pair(0.0, 0.0);
		}

		if ((((bool)okvals[0]) ^ ((bool)okvals[1])) || latephase)
		{
			for (int flag2 = 0; flag2 < 2; flag2++)
			{
				if (!okvals[flag2]) continue;

				for (int k = 0; k < 2; k++)
				{
					if (pars[k])
					{
						int u = ((k ^ flag2) & 1);
						if (latephase || (&themarker.first)[u] != UnknownMarkerVal)
#pragma omp critical (parmarkerval)
						{
							double old1 = 1;
							int old2 = 0;
							auto i = pars[k]->markervals[marker].find((&themarker.first)[u]);
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

	void fixkid(unsigned int marker)
	{
		MarkerValPair& themarker = markerdata[marker];
		if (themarker != make_pair(UnknownMarkerVal, UnknownMarkerVal)) return;

		for (int p = 0; p < 2; p++)
		  {
		    if (pars[p] && pars[p]->markerdata.size() > marker)
		      {
			MarkerValPair& parmarker = pars[p]->markerdata[marker];
			if (parmarker.first == UnknownMarkerVal) continue;
			if (parmarker.first != parmarker.second) continue;

			(&themarker.first)[p] = parmarker.first;
			(&markersure[marker].first)[p] = 0.5;
		      }
		  }

	}

	void addvariance(unsigned int marker, int flag2ignore = 0)
	{
		MarkerValPair& themarker = markerdata[marker];
		auto& thesure = markersure[marker];

		// We only accept an interpretation when it is by exclusion the only possible one. As soon as one interpretation has gained acceptance,
		// no need to retest it.
		double sum = 0;
		double sqsum = 0;
		individ* relp[7] = { pars[0], nullptr, nullptr, pars[1], nullptr, nullptr, this };
		for (int i = 0; i < 2; i++)
		{
			if (relp[i * 3])
			{
				relp[i * 3 + 1] = relp[i * 3]->pars[0];
				relp[i * 3 + 2] = relp[i * 3]->pars[1];
			}
		}
		bool debugtime = (n == 68) && (marker == 0 || marker == 15);

		int count = 0;
			for (shiftflagmode = 0; shiftflagmode < 2; shiftflagmode += 1)
			{
				threadblock tb;
				for (int majori = 0; majori < 2; majori++)
				{
					//int flag2 = -1;
					for (int majorflag2 = 0; majorflag2 < 2; majorflag2++)
					{
						double fullok = 0;
							  double ok = 0;
						for (int subi = 0, i = majori; subi < NUMTYPES*2; subi+=2, i+=2)
						{
							  
							for (int flag2 = majorflag2; flag2 < NUMPATHS; flag2+=2)
					{
						if (flag2 & flag2ignore) {} else {

								  for (int allele = 0; allele < 2; allele++)
								    {
								      double term = trackpossible<false, NO_EQUIVALENCE>(tb, (&themarker.first)[allele], (&thesure.first)[allele],
							marker, i, flag2, *(tb.shiftflagmode), trackpossibleparams());
								      ok += term * (allele ? 1 : -1);
								      fullok += term;
								    }
								count++;
						//				double ok = calltrackpossible<false, false>(tb, 0, marker, i, 0, flag2);
								}
							}
							

						}
						ok = fabs(ok);
							sum += fullok;
						sqsum += ok * ok;
							if (debugtime && fullok) fprintf(stderr, "%02d at %02d: %d %d %03d %03d %lf %lf\n", n, marker, 0, shiftflagmode, majori, majorflag2, ok, fullok);
				}
			}
		}


		if (!sum) return;
		double contrib = (sqsum /*- (sum * sum / count)*/) /* / (sum * sum) * count*/;
		/*		for (individ* rel : relp)
		{
		if (!rel) continue;*/
		variances[marker] = contrib;
		if (debugtime) fprintf(stderr, "Variance for %d at %d is %lf (%lf:%lf)\n", n, marker, contrib, sqsum, sum);
		//		}
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

	// Adjust the probability, i.e. filter all probability values based on the haplotype weights and overall admissibility
	// for the different states.
	void adjustprobs(const threadblock& tb, PerStateArray<double>::T& probs, const unsigned int marker, double& factor,
		const bool oldruleout, int flag99)
	{
		double sum = 0;
		//PerStateArray<double>::T probs2;

		const MarkerValPair& themarker = markerdata[marker];
		const bool ruleout = true; // TODO

		for (int q = 0; q <= (int)!ruleout; q++)
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
				double val = probs[i];

				// We will multiply this already small number with an even smaller number... let's assume it's zero and be done with it.
				if (val < 1e-300)
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
					val *= (double)realok;
				}
				else
				{
					val *= (bool)realok;
				}
				sum += val;
				probs[i] = val;
			}

			if (!ruleout && sum == 0)
			{
				assert(false && "Bring back probs2");
				for (int i = 0; i < NUMTYPES; i++)
				{
					//					probs[i] = probs2[i];
				}

				// The code sees: This doesn't make sense, ignore this marker!
				// NOTE: If the parent-grandparent genotypes are completely incompatible, the problem
				// is NOT eliminated.

				// Current caching scheme does not execute full elimination. Efficiency gains possible by not allowing
				// for genotyping error in every single marker. A small epsilon should otherwise be introduced whenever
				// any probs[i] is zero.
				markerdata[marker] = make_pair(UnknownMarkerVal, UnknownMarkerVal);
				sum = 1;
			}
			else
				break;
		}

		if (sum <= 0)
		{
		  //		  fprintf(stderr, "Error in %d, marker %d, impossible path\n", this->n, marker);
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
#if !DOFB
	double fillortake(const threadblock& tb, const int index, const unsigned int startmark, const unsigned int endmark, PerStateArray<double>::T& probs)
	{
		if ((tb.done[*(tb.shiftflagmode)])[index] != (*tb.generation))
		{

			for (unsigned int i = 0; i < NUMTYPES; i++)
			{
				PerStateArray<double>::T probs2 = { {0} };
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
#pragma GCC ivdep
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
					sum = 1 / sum;
				}
				double* probdest = &(tb.cacheprobs[*tb.shiftflagmode])[index][i][0];

#pragma ivdep
#pragma GCC ivdep
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

		PerStateArray<double>::T probs2 = { {0} };

		for (int i = 0; i < NUMTYPES; i++)
		{
			double step = (tb.factors[*tb.shiftflagmode])[index][i] - factor;
			if (probs[i] == 0.0 || step <= -700.0f) continue;
			double basef = exp((double)step) * probs[i];
			if (basef == 0.0) continue;
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
		const unsigned int endmark, const G& stopdata, const int flag2, bool ruleout, PerStateArray<double>::T& probs,
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

			// A single step, either due to the fixated genotypee being within this range, or simply because we've gone all the
			// way down the tree.
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
					tb.lockpos[*tb.shiftflagmode] = (int)stopdata;
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
					factor += realanalyze<0, T>(tb, turner, startmark, startmark + stepsize, stopdata, flag2, ruleout, &probs);
				}

				if (!_finite(factor) || factor <= minfactor)
				{
					if (!_finite(factor)) fprintf(stderr, "Non-finiteA %d %d\n", n, startmark);
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
				if (!_finite(factor)) fprintf(stderr, "Non-finiteB %d %d\n", n, startmark);
				if (!frommem && !stopdata.okstep(startmark, endmark))
				{
					tb.quickgen[*tb.shiftflagmode] = *tb.generation;
					tb.lockpos[*tb.shiftflagmode] = (int)stopdata;
					tb.quickfactor[*tb.shiftflagmode] = MINFACTOR;
				}

				return MINFACTOR;
			}
		}

		return factor;
	}
#else
// Analyze for a specific range, including a possible fixed specific state at some position (determined by stopdata)
	template<bool inclusive, class T, class G> double quickanalyze(const threadblock& tb, const T& turner, unsigned int startmark,
		const unsigned int endmark, const G& stopdata, const int flag2, bool ruleout, PerStateArray<double>::T& probs,
		float minfactor = MINFACTOR)
	{
		// Probe the distance to test
		int newstart = startmark;
		int origstart = startmark;
		bool allowfull = inclusive;

		while (stopdata.okstep(startmark, newstart + 1))
		{
			unsigned int stepsize;
			for (stepsize = 1; stepsize < (endmark - newstart + allowfull) &&
				stopdata.okstep(newstart, newstart + stepsize); stepsize *= 2)
			{
			}

			if (stepsize == 1) break;

			stepsize /= 2;
			newstart += stepsize;
		}

		startmark = newstart;
		int genotype = stopdata.getgenotype(startmark);
		int pad = genotype == -1 ? 2 : 0;

		double factor = tb.fwbwfactors[*tb.shiftflagmode][startmark][pad];
		probs = tb.fwbw[*tb.shiftflagmode][startmark][pad];
		

		if (factor < minfactor) return factor;

		while (startmark < endmark + inclusive)
		{
			int stepsize = 1;

			bool willquickend = /*(turner.canquickend() && canquickend(startmark, stopdata))*/ stopdata.okstep(startmark + 1, endmark);

			if (willquickend)
			{
				// If we are doing a quick end
				if (pad == 0)
				{
				factor += realanalyze<4, T>(tb, turner, startmark, startmark + stepsize, stopdata, flag2, ruleout, &probs);
				}
				else
				{
					factor += realanalyze<4 | 2, T>(tb, turner, startmark, startmark + stepsize, stopdata, flag2, ruleout, &probs);
				}

				// We might have turned, this value might not exist
				initfwbw(tb, origstart, endmark, 2);

				double sum = 0;

				if (factor < minfactor) return factor;

				for (int k = 0; k < NUMTYPES; k++)
				{
					probs[k] *= tb.fwbw[*tb.shiftflagmode][startmark][1][k];
					sum += probs[k];
				}

				factor += tb.fwbwfactors[*tb.shiftflagmode][startmark][1];

				if (sum <= 0)
				{
				  fprintf(stderr, "No non-zero prob for ind %d at marker %d, factor %lf\n", this->n, startmark, factor);
					factor = MINFACTOR;
				}
				else
				{
					// Normalize, update auxiliary exponent
					sum = 1 / sum;
					for (int i = 0; i < NUMTYPES; i++)
					{
						probs[i] *= sum;
					}
					factor -= log(sum);
				}

				return factor;
			}
			else
			{
				std::cerr << "Incomplete at " << startmark << std::endl;
				factor += realanalyze<0, T>(tb, turner, startmark, startmark + stepsize, stopdata, flag2, ruleout, &probs);
				abort();
			}

			startmark += stepsize;
			allowfull |= true;
		}

		return factor;
	}
#endif

#if !DOFB
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
#else
	void initfwbw(const threadblock& tb, const int startmark, const int endmark, int domask = 3)
	{
		if (tb.fwbwdone[*(tb.shiftflagmode)] != (*(tb.generation) << 2) + domask)
		{
			int donemask = 0;
			if (tb.fwbwdone[*(tb.shiftflagmode)] >> 2 == *(tb.generation))
			{
				donemask = tb.fwbwdone[*(tb.shiftflagmode)] & 3;
			}
			domask &= ~donemask;
			PerStateArray<double>::T probs;

			if (domask & 1)	
			{
			// Initialize forward-backward matrices in one big go.
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
				probs[i] = EVENGEN * (SELFING ?
					selfingfactors[(i >> TYPEBITS)& SELFMASK] : 1.0);
			}

			realanalyze<ANALYZE_FLAG_STORE | ANALYZE_FLAG_FORWARD | 1, noneturner>(tb, noneturner(), startmark, endmark, NONESTOP, -1, false, &probs);
			}

			if (domask & 2)
			{
			for (int i = 0; i < NUMTYPES; i++)
			{
				probs[i] = 1.0;
			}
			realanalyze<ANALYZE_FLAG_STORE | ANALYZE_FLAG_BACKWARD | 1, noneturner>(tb, noneturner(), startmark, endmark, NONESTOP, -1, false, &probs);
			}

			tb.fwbwdone[*(tb.shiftflagmode)] = (*(tb.generation) << 2) | donemask | domask;
		}
	}

	template<class T, class G> double doanalyze(const threadblock& tb, const T& turner, const int startmark, const int endmark, const G& stopdata,
		const int flag2, bool ruleout = false, PerStateArray<double>::T* realprobs = 0, float minfactor = MINFACTOR)
	{
		if (realprobs != 0) fprintf(stderr, "THIS WAS NOT EXPECTED. WHEN IS THIS USED?\n");
		PerStateArray<double>::T probs;

		initfwbw(tb, startmark, endmark);

		return quickanalyze<true, T>(tb, turner, startmark, endmark, stopdata, flag2, ruleout, probs, minfactor);
	}
#endif

	// This is the actual analyzing code. It works with no caches, and can function independently, but is generally only used to patch in those
	// elements not found in caches by quickanalyze and fillortake.
	//
	// Both transition and emission (through adjustprobs) probabilities are handled here.
	//
	// first bit in updateend signals whether the interval is end-inclusive at endmark
	// the second bit in updateend will be 0 if the interval is end-inclusive at startmark, and 1 IF NOT
	// the third bit will cause the code to quit early, after processing the genotype and turner condition
	// the fourth bit will cause the code to go backwards, if set
	// the fifth bit will cause whe code to store probs in fwbw structures
	// Be careful how things switch meaning when going backwards! Could be refactored!
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
		int d = 1;
		int firstmark = startmark;
		int lastmark = endmark;
#if DOFB		
		auto savefwbw = [&tb, &probs, &factor] (int m, int pad = 0)
		{
				copy(probs.begin(), probs.end(),
					tb.fwbw[*tb.shiftflagmode][m][pad + (bool)(updateend & ANALYZE_FLAG_BACKWARD)].begin());

				tb.fwbwfactors[*tb.shiftflagmode][m][pad + (bool)(updateend & ANALYZE_FLAG_BACKWARD)] = factor;
		};
#endif

		if (updateend & ANALYZE_FLAG_BACKWARD)
		{
			d = -1;
			swap(firstmark, lastmark);
#if DOFB
			if ((updateend & ANALYZE_FLAG_BACKWARD) && (updateend & ANALYZE_FLAG_STORE))
			{
				savefwbw(firstmark);
			}
#endif
		}

		for (int j = firstmark + d; j != lastmark + d; j += d)
		{
			double startpos = markerposes[j - d];
			double endpos = markerposes[j];
			int genotype = -1;

			if (updateend & ANALYZE_FLAG_BACKWARD) swap(startpos, endpos);

			bool tofind = stopdata.fixtofind(genotype, startpos, endpos, j);

			int f2use = -1;

			if (tofind)
			{
				f2use = flag2;
			}

#if DOFB
			if (!(updateend & ANALYZE_FLAG_BACKWARD) && (updateend & ANALYZE_FLAG_STORE))
			{
				savefwbw(j - d);				
			}
#endif

			auto filtergenos = [&]
			{
				if (genotype >= 0)
				{
					for (int i = 0; i < NUMTYPES; i++)
					{
						probs[i] *= (i == genotype);
					}
				}
			};

			// Speed up if our turner allows it, this makes the following adjustprobs cheaper
			if (tofind && startpos == endpos && turner.canquickend())
			{
				filtergenos();
			}

			if (genotype != -2)
			{
				// If we are at the very first position, and the specific flag was set, include the emission probabilities for
				// the previous marker. Used to maximize the caching.
				if (!((updateend & 2) && (j == firstmark + d))) adjustprobs(tb, probs, j - d, factor, ruleout, f2use);
			}
			else
			{
				// We could do some stuff here to avoid excessive adjustprobs calls
				// a -2 genotype does not only mean that all genotypes are allowed, but indeed that the marker data at this marker
				// is ignored!
			}

			#if DOFB
			if (!(updateend & ANALYZE_FLAG_BACKWARD) && (updateend & ANALYZE_FLAG_STORE))
			{
				savefwbw(j - d, 2);
			}
			#endif

			// For a specific intra-marker region, we have two cases: the case of a fixated position between the two markers,
			// and the simple case of no fixated position.
			for (int iter = 0; iter <= (int)tofind; iter++)
			{
				// If iter is 1, we have currently handled the transition all the way to the fixated position. Now filter to keep
				// only a single state value positive.
				if (iter)
				{
					turner(probs);
					filtergenos();

					// Were we asked to stop at this very state, in the middle of things?
					if (updateend & 4) return factor;
				}


				double dist = (endpos - startpos);

				// Compute transitions, assuming there is any distance to transfer over.
				if (dist > 0)
				{
					PerStateArray<double>::T probs2 = { {0} };
					double recprob[2 + SELFING][2];
					// Compute recombination probabilities for this specific distance, for the two sexes.
					// (as the sex-dependent marker distance might not be a simple transformation, the actrec
					// data comes into play).
					const int selfgen = gen - 2;
#pragma ivdep
					for (int gen = 0; gen < 2 + SELFING; gen++)
					{
						for (int k = 0; k < 2; k++)
						{
							recprob[gen][k] = 0.5 * (1.0 - exp((gen == 2 ? selfgen : 1) * getactrec(stopdata, startpos, endpos, k, j - d + 1, gen) * (dist)));
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

					array<double, NONSELFNUMTYPES> recombprec;

#pragma ivdep
					for (int index = 0; index < NONSELFNUMTYPES; index++)
					{
						recombprec[index] = 1;
					}

					double selfprec[4 * SELFING + 1][4 * SELFING + 1] = { 0 };
					if constexpr (SELFING)
					{
						selfprec[0][1] = selfprec[0][2] = recprob[2][0];
						selfprec[0][0] = 1 - 2 * recprob[2][0];
						selfprec[1][0] = selfgen ? selfprec[0][1] * 2.0 / ((1 << selfgen) - 1) : 1;
						selfprec[1][2] = selfprec[1][0] * selfprec[0][1]; // TODO, UNDERESTIMATED!
						selfprec[2][0] = selfprec[1][0];
						selfprec[2][1] = selfprec[1][2];
						selfprec[2][2] = selfprec[1][1] = 1 - selfprec[1][0] - selfprec[1][2];
					}

					// Compute probabilities for arbitrary xor values of current and future state
					for (int t = 0; t < TYPEBITS; t++)
					{
						int sex = TYPESEXES[t];
						int gen = TYPEGENS[t];
#pragma ivdep
#pragma GCC ivdep
						for (int index = 0; index < NONSELFNUMTYPES; index++)
						{
							int val = !((index >> t) & 1);
							recombprec[index] *= other[gen][val][sex];
						}
					}

					float relscore[2] = { 1, 1 };
					if (RELSKEWSTATES && iter == tofind)
					{
						relscore[0] = relhaplo[min(j, j - d)];
						relscore[1] = 1 - relhaplo[min(j, j - d)];
					}

					// Use those xor values
					// For the 4-state model, this is an inefficient way to go about it, but it is quite a bit more efficient for
					// the 64-state model (or beyond).
					for (int from = 0; from < VALIDSELFNUMTYPES; from++) // SELFING condition removes the double-bit set case, which is not allowed
					{
						double fromval = probs[from];
						if (fromval <= 0) continue;
#pragma ivdep
#pragma GCC ivdep
						for (int to = 0; to < VALIDSELFNUMTYPES; to++)
						{
							int xored = from ^ to;
							probs2[to] += fromval * recombprec[SELFING ? xored & (NONSELFNUMTYPES - 1) : xored] * (SELFING ? selfprec[(from >> TYPEBITS)& SELFMASK][(to >> TYPEBITS)& SELFMASK] : 1)*
								(RELSKEWSTATES ? relscore[(xored >> BITS_W_SELF) & 1] : 1);
						}
					}

					copy(probs2.begin(), probs2.end(), probs.begin());
				}
				//				else
				{
					//				if (iter == tofind)
					{
						/*for (int c = 0; c < VALIDSELFNUMTYPES; c++)
						{
							if (probs[c] < 1e-200) probs[c] = 1e-200;
						}*/
					}
				}

				// startpos and endpos are always defined in the forward sense, i.e. startpos being upstream of endpos
				if (updateend & ANALYZE_FLAG_BACKWARD)
				{
					endpos = startpos;
					startpos = markerposes[j];
				}
				else
				{
					startpos = endpos;
					endpos = markerposes[j];
				}
			}

#if DOFB
			if ((updateend & ANALYZE_FLAG_BACKWARD) && (updateend & ANALYZE_FLAG_STORE))
			{
				savefwbw(j);				
			}
#endif
		}

		if (updateend & 1)
		{
#if DOFB
			if (!(updateend & ANALYZE_FLAG_BACKWARD) && (updateend & ANALYZE_FLAG_STORE))
			{
				savefwbw(lastmark);				
			}
#endif
			adjustprobs(tb, probs, lastmark, factor, ruleout, -1); // TODO
#if DOFB
			if (!(updateend & ANALYZE_FLAG_BACKWARD) && (updateend & ANALYZE_FLAG_STORE))
			{
				savefwbw(lastmark, 2);
			}
#endif			
		}

		return factor;
	}
};

template<int skipper, int update, int zeropropagate>
vectortrackpossible<skipper, update, zeropropagate>::vectortrackpossible(unsigned int marker, const trackpossibleparams& extparams, int flag2in) : marker(marker), extparams(extparams), tb(tb)
{
	for (double& val : result)
	{
		val = 1.0;
	}
	for (int i = 0; i < NUMPATHS / skipper; i++)
	{
		inmarkerval[0][i] = UnknownMarkerVal;
		secondval[0][i] = 0.;
		geno[0][i] = i * skipper;
		flag2[0][i] = flag2in;
		localshift[0][i] = *tb.shiftflagmode;
	}

	for (int j = 0; j < TURNBITS; j++)
	{
		if (reltreeordered[j] != nullptr)
		{
			meid = j;
			reltreeordered[j]->trackpossiblevector<skipper, update, zeropropagate>(*this);
		}
	}
}	

// Oh, how we waste memory, in a pseudo-O(1) manner
individ* individer[INDCOUNT];
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
		ind->homozyg.resize(markerposes.size());
		ind->markersure.resize(markerposes.size());

		ind->variances.resize(markerposes.size());

		if (RELSKEWS)
		{
			ind->relhaplo.resize(markerposes.size());
		}

		for (size_t x = 0; x < markerposes.size(); x++)
		{
			ind->markerdata[x] = make_pair(UnknownMarkerVal, UnknownMarkerVal);
			ind->haplobase[x] = 0;
			ind->haplocount[x] = 0;
			ind->haploweight[x] = 0.5;
			ind->negshift[x] = 0;
			ind->markersure[x] = make_pair(0, 0);
			if (RELSKEWS)
			{
				ind->relhaplo[x] = 0.5;
			}
		}

		//		ind->semishift.resize(5000);
		ind->lastinved.resize(chromstarts.size());
		ind->lockstart.resize(chromstarts.size());

		for (size_t i = 0; i < chromstarts.size(); i++)
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
			chromstarts.push_back(n - 1);
			n2 = 0;
			oldc = c;
		}
		markerposes.push_back(bppos / 1000000.0);
		for (int t = 0; t < 2; t++)
		{
			actrec[t].push_back(baserec[t]);
		}
		n2++;
	}
	chromstarts.push_back(n);
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

		if (ind->gen > 0) dous.push_back(ind);

		ind->sex = (sex[0] == 'F');
		ind->strain = 1;
		ind->markerdata.resize(markerposes.size());

		ind->haplobase.resize(markerposes.size());
		ind->haplocount.resize(markerposes.size());
		ind->haploweight.resize(markerposes.size());
		ind->negshift.resize(markerposes.size());

		ind->infprobs.resize(markerposes.size());
		ind->markersure.resize(markerposes.size());
		//		ind->semishift.resize(5000);
		ind->lastinved.resize(chromstarts.size());
		ind->lockstart.resize(chromstarts.size());

		for (size_t i = 0; i < chromstarts.size(); i++)
		{
			ind->lastinved[i] = -1;
			ind->lockstart[i] = 0;
		}

		for (size_t i = 0; i < markerposes.size(); i++)
		{
			ind->haploweight[i] = 0.5;
		}
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

		for (unsigned int i = 0; i < /*markerposes.size() - 2*/ 10031; i++)
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

		for (size_t i = 0; i < markerposes.size(); i++)
		{
			ind->haploweight[i] = 0.5;
		}

		for (size_t i = 0; i < markertranslation.size(); i++)
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
		for (unsigned int i = chromstarts[0]; i < chromstarts[1]; i++)
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
		for (int i = 0; i < INDCOUNT; i++)
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
					for (unsigned int marker = 0; marker < chromstarts[1]; marker++)
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

							gval |= mask * ((bool)(ind->pars[p]->genotypegrid[i] & (1 << parhalf)));

							for (; markerposes[marker] < i + 1 && marker < chromstarts[1]; marker++)
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
				for (unsigned int marker = 0; marker < chromstarts[1]; marker++)
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

	if (ind->lockstart[i] >= chromstarts[i + 1])
	{
		ind->lockstart[i] = 0;
	}
	int bestpos = -1;
	double bestvar = 0;

	for (j = max(chromstarts[i], ind->lockstart[i]); j != chromstarts[i + 1]; j++)
	{
		if (ind->variances[j] > bestvar)
		{
			bestpos = j;
			bestvar = ind->variances[j];
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
		ind->haploweight[j] = ind->haploweight[j] <= 0.5 ? 0 : 1;
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

pair<int, int> fixtrees(individ* ind)
{
	int flag2ignore = 0;
	int shiftignore = 0;

	reltree.clear();
	relmap.clear();
	relmapshift.clear();
	reltreeordered.clear();

	reltree.push_back(ind);
	reltreeordered.resize(TURNBITS); // Right constant? Lots of equalities with some settings.
	reltreeordered[0] = ind;
	relmap[ind] = 1;
	relmapshift[ind] = 1;

	if (HAPLOTYPING)
	{
		flag2ignore = 1;
		shiftignore = 0;
		bool anylev1 = false;
		for (int lev1 = 0; lev1 < 2; lev1++)
		{
			individ* lev1i = ind->pars[lev1];
			if (!lev1i) continue;
			int flag2index = 1 + lev1 * ((1 << (NUMFLAG2GEN - 1)) - 1);
			int shiftval = (NUMGEN == 3) ? (2 << lev1) : 0;

			if (!lev1i->empty)
			{
				flag2ignore |= 1 << flag2index;
				relmap[lev1i] |= 1 << flag2index;
				relmapshift[lev1i] |= shiftval;
				reltreeordered[flag2index] = lev1i;
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
						flag2ignore |= 1 << (flag2index + lev2 + 1);
						relmap[lev2i] |= 1 << (flag2index + lev2 + 1);
						relmapshift[lev2i] |= 0;
						reltreeordered[flag2index + lev2 + 1] = lev2i;
						anypars = true;
					}
					reltree.push_back(lev2i);
				}
			}

			if (anypars)
			{
				shiftignore |= shiftval;
			}
			else
			{
				//lev1i->founder = true;
			}
			// Any information of relevance in parents
			if (anypars || !lev1i->empty)
			{
				anylev1 = true;
			}
		}
		if (anylev1)
		{
			shiftignore |= 1;
		}
		else
		{
			ind->founder = true;
		}
		flag2ignore ^= (NUMPATHS - 1);
		shiftignore ^= (NUMSHIFTS - 1);
	}

	//reltreeordered = reltree;
	sort(reltree.begin(), reltree.end());
	reltree.resize(unique(reltree.begin(), reltree.end()) - reltree.begin());

	return make_pair(shiftignore, flag2ignore);
}


// Some operations performed when marker data has been read, independent of format.
void postmarkerdata(int indcount = INDCOUNT)
{
	int any, anyrem;
	bool latephase = false;


	markerweight.resize(markerposes.size());
	// If inference is active, add "new" marker data until all data has been found.
	if (CORRECTIONINFERENCE) do
	{
#pragma omp parallel for schedule(dynamic,32)
		// all haploweights MUST be non-zero at this point, as we do not explore all shiftflagmode values
		for (int i = 1; i < indcount; i++)
		{
			individ* ind = getind(i);
			if (ind) ind->children = 0;
		}
#pragma omp parallel for schedule(dynamic,32)		
		for (int i = 1; i < indcount; i++)
		{
			individ* ind = getind(i);
			if (!ind) continue;

			if (ind->markerdata.size())
			{
				generation++;

				for (size_t g = 0; g < ind->markerdata.size(); g++)
				{
				  ind->fixkid(g);
				}
			}
		}
		map<individ*, int> upsent;
		bool anychange;
		do
		{
			anychange = false;
			for (int i = 1; i < indcount; i++)
			{
				individ* ind = getind(i);
				if (!ind) continue;
				auto& upsentval = upsent.insert({ind, 0}).first->second;
				int nowdesc = ind->descendants;
				if (!nowdesc) nowdesc = 1;
				nowdesc -= upsentval;
				if (nowdesc > 0)
				{
					for (individ* par : ind->pars)
					{
						if (par)
						{
							par->descendants += nowdesc;
						}
					}
					upsentval += nowdesc;
					anychange = true;
				}
			}
		} while (anychange);
		for (int i = 1; i < indcount; i++)
		{
			individ* ind = getind(i);
			if (!ind) continue;
			if (ind->descendants == 0) ind->descendants = 1;
		}
#pragma omp parallel for schedule(dynamic,32)
		for (int i = 1; i < indcount; i++)
		{
			individ* ind = getind(i);
			if (!ind) continue;

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

				for (size_t g = 0; g < ind->markerdata.size(); g++)
				{
					ind->fixparents(g, latephase);
				}
			}
		}

		any = 0;
		anyrem = 0;
		for (int i = 1; i < indcount; i++)
		{
			individ* ind = getind(i);
			if (!ind) continue;

			for (size_t g = 0; g < ind->markervals.size(); g++)
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
				if (ind->markerdata[g].first != UnknownMarkerVal) ind->markervals[g].insert(make_pair(ind->markerdata[g].first, make_pair(ind->children, ind->markersure[g].first)));
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
					ind->markersure[g] = make_pair(dosureval(ind->markervals[g].begin()->second.first, ind->markervals[g].begin()->second),
								       dosureval(ind->markervals[g].begin()->second.first, ind->markervals[g].begin()->second));
				} // DANGEROUS ASSUMPTIONS
				else if (!latephase && ind->markervals[g].size() == 1 && known == 0)
				{
					any++;
					ind->markerdata[g] = make_pair(ind->markervals[g].begin()->first, UnknownMarkerVal);
					ind->markersure[g] = make_pair(dosureval(ind->markervals[g].begin()->second.first, ind->markervals[g].begin()->second), 0.0);
				}


				if (any != oldany) printf("Correction at %d, marker %d (%d;%d) (%lf;%lf)\n", i, g,
					ind->markerdata[g].first.value(), ind->markerdata[g].second.value(),
					ind->markersure[g].first, ind->markersure[g].second);

			}

			for (size_t g = 0; g < ind->markervals.size(); g++)
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
	} while (any > anyrem);

#pragma omp parallel for schedule(dynamic,32)
	for (int i = 1; i < indcount; i++)
	{
		individ* ind = getind(i);

		// Lock the first position
		if (HAPLOTYPING && ind && ind->haploweight.size())
			// These days we do the locking in all generations
		{
			auto [shiftignore, flag2ignore] = fixtrees(ind);
			for (int j = 0; j < markerposes.size(); j++)
			{
				// TODO: Variance could take relskew into account...
				ind->addvariance(j, flag2ignore);
			}
		}
	}


#pragma omp parallel for schedule(dynamic,32)
	for (int i = 1; i < indcount; i++)
	{
		individ* ind = getind(i);

		// Lock the first position
		if (HAPLOTYPING && ind && ind->haploweight.size())
			// These days we do the locking in all generations
		{
			ind->markervals.clear();
			// Lock the haplotype (in an arbitrary manner) for the first marker in each linkage groupfl
			// Propagation would be more efficient if we locked the mid-position (which should really be determined in cM)
			// Doing so would be much more opaque, though...
			for (size_t i = 0; i < chromstarts.size() - 1; i++)
			{
				//if (!ind->pars[0] && !ind->pars[1])
				lockhaplos(ind, i);
			}
		}
	}
}

typedef boost::tuple<individ*, double, unsigned int> negshiftcand;
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
		for (int p = minstart + 1; p < (int)chromstarts[c + 1]; p++)
		{
			if (p == minstart + 1) fprintf(stdout, "Inv: %d %d\n", ind->n, p);
			ind->haploweight[p] = 1.0f - ind->haploweight[p];
			if (ind->n == 24 && p > 8755 && p < 8765)
			{
				printf("CHECK: %lf\n", ind->haploweight[p]);
			}
			ind->haplobase[p] = ind->haplocount[p] - ind->haplobase[p];
		}
	}
};

bool ignoreflag2(int flag2, int g, int shiftflagmode, int q, int flag2ignore, const flat_map<individ*, int>& relmap, const flat_map<individ*, int>& relmapshift)
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

	int marker = -q - 1000;
	for (auto i = as_const(relmap).begin(), j = as_const(relmapshift).begin(); i != relmap.end(); i++, j++)
	{
		int currfilter = (i->second & flag2filter);
		int filtered = ((flag2 ^ (g * 2)) & currfilter);
		// Require ALL bits in the flag to be set, if at least one is set
		if (filtered && filtered != currfilter) return true;
		//if (marker >= 0 && i->first->markerdata[marker].first == UnknownMarkerVal && i->first->markerdata[marker].second == UnknownMarkerVal && (!(flag2 & i->second)))
		if (marker >= 0 && i->first->markerdata[marker].first == i->first->markerdata[marker].second && i->first->markersure[marker].first == i->first->markersure[marker].second && !(((bool)filtered) ^ ((bool)(shiftflagmode & j->second))) &&
			((!RELSKEWSTATES || currfilter != 1) && (!SELFING/* || selfgen == 0*/)))
		{
			//			return false;
			return true;
		}
	}
	return false;
}


template<int N> struct valuereporter
{
	array<double, N> probs;

	valuereporter()
	{
		for (int k = 0; k < N; k++)
		{
			probs[k] = 0;
		}
	}

	void report(vector<char>& outqueue)
	{
		double probsum = 0;
		for (int i = 0; i < N; i++)
		{
			probsum += probs[i];
		}
		probsum = 1 / probsum;
		for (int i = 0; i < N; i++)
		{
			char string[255];
			int val;
			sprintf(string, "%.5lf%c%n", probs[i] /** probsum*/, i == N - 1 ? '\n' : '\t', &val);
			for (int k = 0; k < val; k++)
			{
				outqueue.push_back(string[k]);
			}
		}
	}
};

struct genotypereporter : valuereporter<3>
{
	void addval(int q, int mapval, int g, int flag2, double val)
	{
		probs[mapval] += val;
	}
};

struct statereporter : valuereporter<NUMTYPES>
{
	void addval(int q, int mapval, int g, int flag2, double val)
	{
		probs[g] += val;
	}
};

template<typename T> T square(T v)
{
	return v * v;
}

void resizecaches()
{
	// Some heaps are not properly synchronized. Putting a critical section here makes the operations not safe,
	// but *safer*.
#pragma omp critical(uglynewhack)
	for (unsigned int t = 0; t < NUMSHIFTS; t++)
	{
#if !DOFB 
		factors[t].resize(markerposes.size());
		cacheprobs[t].resize(markerposes.size());
		done[t].resize(markerposes.size());
#else
		fwbwfactors[t].resize(markerposes.size());
		fwbw[t].resize(markerposes.size());
		fwbwdone[t] = 0;
#endif
}
}

// Global scale factor, 1.0 meaning "use unscaled gradient".
double scalefactor = 0.013;
double entropyfactor = 1;


void moveinfprobs(int i, int k, int marker, double norm, double descfactor)
{
	// TODO: Sanitize for zero, inf, nan.
	// Compensate for duplicates
	// TODO: Does this compensation make sense when infprobs is cleared below?
	norm *= 2;
	for (individ* ind : reltreeordered)
	{
		if (ind == reltree[k]) norm /= 2;
	}
	norm *= descfactor;

	for (int side = 0; side < 2; side++)
	{
		for (auto infval : infprobs[i][side])
		{
			reltree[k]->infprobs[marker][side][infval.first] += infval.second * norm /* * sum*/;
		}
		infprobs[i][side].clear();
	}
}

void movehaplos(int i, int k, int marker, double descfactor)
{
	if (haplos[i][0] || haplos[i][1])
	{
		if (fabs(reltree[k]->haploweight[marker] - 0.5) < 0.5 - 1e-12)
		{
			double b1 = (haplos[i][0] + exp(-400) * maxdiff * maxdiff * 0.5) /*/ reltree[k]->haploweight[marker] /** (1 - reltree[k]->haploweight[marker])*/;// * (1 + 1e-10 - rhfactor);
			double b2 = (haplos[i][1] + exp(-400) * maxdiff * maxdiff * 0.5) /*/ (1 - reltree[k]->haploweight[marker]) /** reltree[k]->haploweight[marker]*/;// * (rhfactor + 1e-10);
			//if (i == 89 || i == 90) fprintf(stderr, "IND %d K %d MARKER %d %lf %lf\n", i, k, marker, b1, b2);
			{
				reltree[k]->haplobase[marker] += /*log(b1 / b2)*/b1 / (b1 + b2) * descfactor;
				reltree[k]->haplocount[marker] += descfactor;
			}
		}
		haplos[i][0] = 0;
		haplos[i][1] = 0;
	}
}

void calcdistancecolrowsums(double mwvals[1][1], double rowsums[NUMTYPES], double colsums[NUMTYPES], double& acc3, double& acc4, double mwfvals[1])
{
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
}

void updatenegshifts(const bool validg[NUMTURNS], int shifts, int shiftend, int shiftignore, const double rawvals[NUMTURNS][NUMSHIFTS], int j, int marker)
{
	for (int g = 0; g < NUMTURNS; g++)
	{
		if (!validg[g]) continue;

		double val = 0;

		for (int s = shifts; s < shiftend; s++)
		{
			if (s & shiftignore) continue;
			if (rawvals[g][s] < 0) continue;
			val += rawvals[g][s];
		}

		// TODO: Capping val here, but not in sumnegval, leads to some
		// strange results, indeed. Maybe?
		if (!_finite(val) || val < 1e-174) val = 1e-174;
		{
			int g2 = g;
			if (!g) g2 = (1 << 15) - 1;

			// This is hardcoded for the generation count of 3.
			if (NUMGEN == 3)
			{
				dous[j]->negshift[marker] += log(val) * (1.0 - ((g >> 6) & 1) * 2)* ((g2 >> 6) & 1) /*/ sumnegval[6]*/;

				if (dous[j]->pars[0])
					dous[j]->pars[0]->negshift[marker] += log(val) * (1.0 - ((g >> 0) & 1) * 2)* ((g2 >> 0) & 1) /* / sumnegval[0] */;

				if (dous[j]->pars[1])
					dous[j]->pars[1]->negshift[marker] += log(val) * (1.0 - ((g >> 3) & 1) * 2)* ((g2 >> 3) & 1) /* / sumnegval[3] */;

				if (dous[j]->gen >= 2)
				{
					if (dous[j]->pars[0] && dous[j]->pars[0]->pars[0])
						dous[j]->pars[0]->pars[0]->negshift[marker] += log(val) * (1.0 - ((g >> 1) & 1) * 2)* ((g2 >> 1) & 1) /* / sumnegval[1]*/ / dous[j]->pars[0]->children;

					if (dous[j]->pars[0] && dous[j]->pars[0]->pars[1])
						dous[j]->pars[0]->pars[1]->negshift[marker] += log(val) * (1.0 - ((g >> 2) & 1) * 2)* ((g2 >> 2) & 1) /*/ sumnegval[2]*/ / dous[j]->pars[0]->children;

					if (dous[j]->pars[1] && dous[j]->pars[1]->pars[0])
						dous[j]->pars[1]->pars[0]->negshift[marker] += log(val) * (1.0 - ((g >> 4) & 1) * 2)* ((g2 >> 4) & 1) /*/ sumnegval[4]*/ / dous[j]->pars[1]->children;

					if (dous[j]->pars[1] && dous[j]->pars[1]->pars[1])
						dous[j]->pars[1]->pars[1]->negshift[marker] += log(val) * (1.0 - ((g >> 5) & 1) * 2)* ((g2 >> 5) & 1) /*/ sumnegval[5]*/ / dous[j]->pars[1]->children;
				}
			}
#if false
			else
			{
				dous[j]->negshift[marker] += val * (1.0 - ((g >> 2) & 1) * 2) /* * ((g2 >> 2) & 1) */ / (sumnegval[2]);

				if (dous[j]->pars[0])
					dous[j]->pars[0]->negshift[marker] += log(val) * (1.0 - ((g >> 0) & 1) * 2)* ((g2 >> 0) & 1) /*/ (sumnegval[0])*/;

				if (dous[j]->pars[1])
					dous[j]->pars[1]->negshift[marker] += log(val) * (1.0 - ((g >> 1) & 1) * 2)* ((g2 >> 1) & 1) /*/ (sumnegval[1])*/;
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
#endif
		}
	}
}

#if 0
void oldinfprobslogic(individ* ind, unsigned int j, int iter, int cno, FILE* out)
{
	MarkerVal bestvals[2] = { UnknownMarkerVal, UnknownMarkerVal };
	double bestsure[2] = { ind->markersure[j].first, ind->markersure[j].second };
	bool foundbest = true;
	bool surefound = false;

	flat_map<MarkerVal, double> surenesses;

	if (DOINFPROBS)
	{
		flat_map<MarkerVal, double> sums[2];
		for (int a = 0; a < 2; a++)
		{
			for (auto i = ind->infprobs[j][a].begin(); i != ind->infprobs[j][a].end(); i++)
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
			flat_map<MarkerVal, double> infprobs;


			for (flat_map<pair<MarkerVal, MarkerVal>, double>::iterator i = ind->infprobs[j][a].begin(); i != ind->infprobs[j][a].end(); i++)
			{
				sum += i->second;
				infprobs[i->first.first] += i->second;
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

			for (flat_map<MarkerVal, double>::iterator i = infprobs.begin(); i != infprobs.end(); i++)
			{
				double sureness = i->second / sum;
				double origsureness = sureness;

				if (!((&(ind->markerdata[j].first))[a] == UnknownMarkerVal))
				{
					if (i->first != (&(ind->markerdata[j].first))[a])
					{
						if (sureness > 0.9999)
						{
							sureness = 0.9999;
						}
					}
				}

				if (origsureness > 0.9999) origsureness = 0.9999;
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
		if (ind->sex) foundbest = false;

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

						for (auto i = surenesses.begin(); i != surenesses.end(); i++)
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
		}
		for (int a = 0; a < 2; a++)
		{
			ind->sureinfprobs[j][a].clear();
			ind->infprobs[j][a].clear();
			ind->unknowninfprobs[j][a] = 0;
		}
	}
}
#endif

int oldhitnnn = 0;
int oldhitnnn2 = 0;
double caplogitchange(double intended, double orig, double epsilon, std::atomic_int& hitnnn, bool breakathalf)
{
	double nnn = 3;
	if (nnn < 1.0) nnn = 1.0;

	double limn = (nnn - 1.0) * orig * (-1 + orig);

	double limd1 = -1 - (nnn - 1.0) * orig;
	double limd2 = (nnn - 1.0) * orig - nnn;

	intended = min((double)intended, 1.0 - epsilon);
	intended = max((double)intended, epsilon);
	double diff = intended - orig;


	if (diff > limn / limd1)
	{
		//		if (!hitnnn)  fprintf(stderr, "CAP: Exceeding limit %lf > %lf, intended %lf, orig %lf\n", diff, limn / limd1, intended, orig);
		intended = orig + limn / limd1;
		if (intended < 0.5) hitnnn++;
	}

	if (diff < -limn / limd2)
	{
		//if (!hitnnn) fprintf(stderr, "CAP: Underflowing limit %lf < %lf, intended %lf, orig %lf\n", diff, -limn / limd2, intended, orig);
		intended = orig - limn / limd2;
		if (intended > 0.5) hitnnn++;
	}

	if (breakathalf && (intended - 0.5) * (orig - 0.5) < 0) intended = 0.5 * (0.5 + orig);

	return intended;
}

template<class T> double cappedgd(T& gradient, double orig, double epsilon, std::atomic_int& hitnnn, bool breakathalf = false)
{    
  double randomdrift = 0 /*boost::random::normal_distribution(0., 1e-0)(rng)*/;
#if false
  auto prelcompute = [&] (double startval, double accuracy, double starttime, double endtime) -> double
  {
	  array<double, 1> state{startval - 0.5};
  	  array<double, 1> endstate{startval - 0.5};
  	int hitcount = 0;
  
	auto adaptive = ode::make_adaptive_time_range(ode::make_controlled<ode::runge_kutta_cash_karp54< array<double, 1> >>(accuracy, accuracy/*epsilon * 0.1, 0.1*/),  	
				[&] (array<double, 1>& in,
					array<double, 1>& out, double time)
				{
					double val = std::clamp(in[0] + 0.5, epsilon, 1 - epsilon);
					
				//if (val < epsilon || val > 1 - epsilon) out[0] = 0;
				//else
				{
				gradient(val, out[0], time);
				out[0] += randomdrift;
				out[0] += log((orig / (1-orig)) / (val / (1 - val)));
				if (!isfinite(out[0])) printf("Invalid gradient %lf for %lf, starting from %lf\n", out[0], in[0] + 0.5, orig);
				}
				}, state,
				starttime, endtime, endtime * 0.01);
	for (auto [nowstate, time] : boost::make_iterator_range(adaptive.first, adaptive.second))
		{
			if (nowstate[0] < epsilon - 0.5 || nowstate[0] > 0.5 - epsilon) hitcount++;
			endstate = nowstate;
			if (hitcount > 96) break;		
		}
		return endstate[0] + 0.5;
  };
  
  if (false)
  {
	double res = orig;
	double newaccuracy = 1e-2;
	double starttime = 0;
	int step = 0;
	do
	{
		double prel = prelcompute(res, newaccuracy, starttime, scalefactor);
		prel = std::clamp(prel, epsilon, 1 - epsilon);
		newaccuracy = fabs(res - prel) * 1e-2;
		if (newaccuracy <= 1e-6)
		{
			newaccuracy = 1e-6;
			break;
		}
		if (step++ >= 4) break;
		double endtime = starttime + (scalefactor - starttime) * 0.9;
		res = prelcompute(res, newaccuracy, starttime, endtime);
		res = std::clamp(res, epsilon, 1 - epsilon);

		starttime = endtime;
	} while (true);

	res = prelcompute(res, newaccuracy, starttime, scalefactor);
				

	return caplogitchange(res, orig, epsilon, hitnnn, breakathalf);
  }
  else
#endif
  {
	  std::atomic_int dumpval;	  
	  auto actualgradient = [orig, epsilon, randomdrift, &gradient] (double val) -> double
		{
			val = std::clamp(val, epsilon, 1 - epsilon);
			double toret;
			gradient(val, toret, -1);
			return 1. / (toret + randomdrift/*+ log((orig / (1-orig)) / (val / (1 - val)))*/);
		};
	  // We make the binary search limits slightly too large for our final caplogitchange to adjust hitnnn				
	  double lolim = caplogitchange(epsilon, orig, epsilon, dumpval, breakathalf);
	  double lo = lolim - epsilon * 0.125;
	  double hilim = caplogitchange(1 - epsilon, orig, epsilon, dumpval, breakathalf);
	  double hi = hilim + epsilon * 0.125;

	  orig = caplogitchange(orig, orig, epsilon, dumpval, breakathalf);

	  double gradval = actualgradient(orig);
	  if (!isfinite(gradval) || !scalefactor)
	    {
	      lo = orig;
	      hi = orig;
	    }
	  bool lowside = gradval < 0;
	  (lowside ? hi : lo) = orig;
	  for (int i = 0; i < 51 && scalefactor; i++)
	  {
		  // We're outside true bounds
		  if (lo > hilim || hi < lolim) break;

		  double mid = (lo + hi) / 2;
		  double prel = 0;
		  double gradval = actualgradient(mid);
		  if (gradval < 0 ^ lowside || !isfinite(gradval))
		  {
			  prel = (scalefactor + 0.1) * 1.1;
		  }
		  else
		  {
			double start = orig;
			double end = mid;		  
			if (start > end) swap(start, end);
			if (end - start < 1e-10) break;
			//double prel = boost::math::quadrature::gauss_kronrod<double, 15>::integrate(actualgradient, start, end, 10, 1e-3);		
			prel = boost::math::quadrature::gauss<double, 15>::integrate(actualgradient, start, end);		
			if (end != mid) prel = -prel;
			if (!isfinite(prel)) prel = (scalefactor + 0.1) * 1.1;
		  }
		if (fabs(prel - scalefactor) < scalefactor * 1e-3)
		{
			break;
		}

		if ((prel < scalefactor) ^ lowside)
		{
			lo = mid;
		}
		else
		{
			hi = mid;
		}
	  }

	  if (!scalefactor)
	  {
		  lo = orig;
		  hi = orig;
	  }

	  return caplogitchange((lo + hi) / 2, orig, epsilon, hitnnn, breakathalf);
  }
}

void processinfprobs(individ* ind, const unsigned int j, const int side, int iter, std::atomic_int& hitnnn)
{
	double bestprob = 0;
	MarkerVal bestmarker = UnknownMarkerVal;
	double sum = 0;
	bool doprint = false;


	for (auto probpair : ind->infprobs[j][side])
	{
		sum += probpair.second;
	}

	MarkerVal priorval = UnknownMarkerVal;
	if (ind->priormarkerdata.size() > j)
	{
		priorval = (&ind->priormarkerdata[j].first)[side];
	}

	for (int i = 0; i < 2; i++)
	{
		if (doprint) fprintf(stdout, "PROBHZYG  : %d %d %d   %lf\n", ind->n, j, i, ind->homozyg[j][i]);
	}

	double hzygcorrsum = 0;
	for (auto probpair : ind->infprobs[j][side])
	{
		if (probpair.first.value() >= 1 && probpair.first.value() <= 2)
		{
			double hzygcorrterm = 0;
			MarkerVal copymv = probpair.first;
			double otherside = fabs((!markermiss<false>(copymv, (&ind->markerdata[j].first)[!side])
				? 1.0 : 0.0) - (&ind->markersure[j].first)[!side]);
			double hzygval = ind->homozyg[j][probpair.first.value() - 1];
			hzygcorrterm -= hzygval / probpair.second * sum;
			hzygcorrterm += hzygval / otherside;

			hzygcorrsum += hzygcorrterm * (probpair.first.value() == 1 ? 1 : -1);
		}
	}

	double ef = exp(0 * -0.01 * iter) * entropyfactor;

	for (auto probpair : ind->infprobs[j][side])
	{
		if (doprint) fprintf(stdout, "PROBPAIR A: %d %d %d %d %lf\n", ind->n, j, side, probpair.first.value(), probpair.second);
		double curprob = 0.5;
		auto curmarker = (&ind->markerdata[j].first)[side];

		if (curmarker != UnknownMarkerVal)
		{
			curprob = fabs((curmarker == probpair.first ? 1 : 0) - (&ind->markersure[j].first)[side]);
		}

		double hzygcorred = probpair.second;
		/*				if (probpair.first.value() >= 1 && probpair.first.value() <= 2)
		{
			hzygcorred += hzygcorrsum * (probpair.first.value() == 1 ? 1 : -1);
			}
		*/		if (doprint) fprintf(stdout, "PROBPAIR a: %d %d %d %d %lf\n", ind->n, j, side, probpair.first.value(), hzygcorred);

		double hw = ind->haploweight[j];
		double etf = 1 * ef/*+ (side ? 0 : -1) * 4 * (hw - hw * hw)*/;
		double epsilon = maxdiff / (ind->children + 1);
		double priord = 0;
		double priorprob = 0.5;
			if (priorval != UnknownMarkerVal)
			{
				priorprob = 1.0 - (&ind->priormarkersure[j].first)[side];

				MarkerVal nowval = probpair.first;
				if (nowval != priorval)
				{
					priorprob = 1.0 - priorprob;
				}
				if (priorprob == 0)
				{
					priord -= 10000;
				}
				else if (priorprob == 1)
				{
					priord += 10000;
				}
				else
				{
					priorprob = std::clamp(priorprob, (double) 1e-14, (double) (1 - 1e-14));

			priord += log(priorprob) - log(1 - priorprob);
			}
			}

		auto gradient = [&](const double& x, double& out, const double)
		{
			//double d = (hzygcorred - sum * x) / (x - x *x);
			// expr6 = ((h)*(1-x)*log((1-x))+g*x*log(x))/(h*(1-x)+g*x)
			// str(simplify(simplify(factor(simplify(-val.subs(g,g/y).subs(h,(h-g)/(1-y)).subs(y, curprob).subs(g, hzygcorred).subs(h, sum))))))
			double d = -(-square(curprob*hzygcorred)*log(x) + square(curprob*hzygcorred)*log(1 - x) + square(curprob)*hzygcorred*sum*log(x) - square(curprob)*hzygcorred*sum*log(1 - x) - square(curprob)*hzygcorred*sum - square(curprob*sum)*x + square(curprob*sum) + curprob*square(hzygcorred)*log(x) - curprob*square(hzygcorred)*log(1 - x) + curprob*square(hzygcorred) + 2*curprob*hzygcorred*sum*x - curprob*hzygcorred*sum*log(x) + curprob*hzygcorred*sum*log(1 - x) - curprob*hzygcorred*sum - square(hzygcorred)*x)/square(curprob*hzygcorred + curprob*sum*x - curprob*sum - hzygcorred*x);
			/*double d = (hzygcorred - sum * curprob + curprob * curprob * log(curprob/priorprob)
					   -curprob * curprob * log((curprob - 1)/(priorprob - 1)) - curprob * log(curprob/priorprob)
					   + curprob * log((curprob - 1)/(priorprob - 1)));*/

			double et = log(1 / x - 1);
			//double et = log((1-curprob) * priorprob / ((1 - priorprob) * curprob));
			d += etf * et; // Entropy term
			d += etf * priord;

			out = d;
			//if (fabs(out) > 1e8) printf("Large grad %lf, marker %d, ind %d, in %lf, hzygcorred %lf, sum %lf\n", out, j, ind->n, in, hzygcorred, sum);
		};

		double intended = cappedgd(gradient, curprob, epsilon, hitnnn);
		ind->infprobs[j][side][probpair.first] = intended;
	}

	for (auto probpair : ind->infprobs[j][side])
	{
		if (probpair.second > bestprob - (side ? 1e-30 : 0))
		{
			bestmarker = probpair.first;
			bestprob = probpair.second;
		}
		if (doprint) fprintf(stdout, "PROBPAIR B: %d %d %d %d %lf\n", ind->n, j, side, probpair.first.value(), probpair.second);
	}

	// We might have stats for empty inds, but those stats are incomplete
	// Putting them in will make the ind non-empty, breaking assumptions
	if (!ind->empty && (bestmarker != UnknownMarkerVal || bestprob > 0))
	{
		//if (iter <= 3)
		if (ind->priormarkerdata.size() > j)
		{
			(&ind->markerdata[j].first)[side] = bestmarker;
			double intended = 1.0 - bestprob;
			(&ind->markersure[j].first)[side] = intended;
		}
	}
	ind->infprobs[j][side].clear();
	if (side == 1)
	{
		for (int i = 0; i < 2; i++)
		{
			ind->homozyg[j][i] = 0;
		}
	}
}

struct relskewhmm
{
	static constexpr bool realhmm = true;
	typedef array<double, 2> halfstate;
	typedef array<halfstate, 2> state;
	vector<state> relskewfwbw;
	vector<double> ratio;

	const int firstmarker;
	const individ* ind;

	relskewhmm() : firstmarker(-1), ind(nullptr) {}

	relskewhmm(int firstmarker, int endmarker, const individ* ind) : firstmarker(firstmarker), ind(ind)
	{
		relskewfwbw.resize(endmarker - firstmarker);
		ratio.resize(endmarker - firstmarker);
		halfstate s = { 0.5 * realhmm, 0.5 * realhmm };

		const auto doemissions = [&s, &ind](int m)
		{
			double w;
			if (false && ind->haplocount[m]) w = ind->haplobase[m] / ind->haplocount[m];
			else
			{
				//		      fprintf(stderr, "FILLING IN %d %d\n", ind->n, m);
				w = ind->haploweight[m];
			}
			for (int k = 0; k < 2; k++)
			{
				auto val = fabs(!k - w);
				if (realhmm) s[k] *= val;
				else
					s[k] += val;
			}
		};
		const auto dotransitions = [&s, &ind](int m)
		{
			double n = ind->relhaplo[m];
			double nb = 1 - n;

			if (!realhmm)
			{
				double base = min(n, nb);
				nb -= base;
				n -= base;
			}
			halfstate nexts;
			for (int k = 0; k < 2; k++)
			{
				nexts[k] = s[k] * n + s[!k] * nb;
			}

			if (!realhmm)
			{
				double sum = 0;
				for (int k = 0; k < 2; k++)
				{
					sum += nexts[k];
				}

				sum *= (1.0 - n - nb) / max((double)(n + nb), (double)maxdiff);
				sum = 1.0 / max((double)sum, (double)maxdiff);

				for (int k = 0; k < 2; k++)
				{
					nexts[k] *= sum;
				}
			}
			s = nexts;
		};
		const auto renormalizes = [&s]
		{
			if (!realhmm) return; // HACK: Disable normalization
			double sum = s[0] + s[1];
			if (sum < 1e-10)
			{
				s[0] *= 1e20;
				s[1] *= 1e20;
			}
			/*if (sum == 0)
			{
				fprintf(stderr, "Renormalization problem.\n");
			}*/
		};

		// FW
		for (int m = firstmarker; m < endmarker; m++)
		{
			doemissions(m);
			relskewfwbw[m - firstmarker][0] = s;

			renormalizes();
			dotransitions(m);
		}

		// BW, but not really, emissions included everywhere
		s = { 0.5 * realhmm, 0.5 * realhmm };
		int m = endmarker - firstmarker - 1;
		ratio[m] = relskewfwbw[m][0][1] / (relskewfwbw[m][0][0] + relskewfwbw[m][0][1]);
		for (int m = endmarker - 2; m >= firstmarker; m--)
		{
			doemissions(m + 1);
			relskewfwbw[m - firstmarker + 1][1] = s;

			dotransitions(m);
			renormalizes();
			double ratiofactors[2] = { 0 };
			for (int k = 0; k < 2; k++)
			{
				//for (int i = 0; i < 2; i++)
				{
					ratiofactors[k] += s[k] * relskewfwbw[m - firstmarker][0][k];
				}
			}

			ratio[m - firstmarker] = ratiofactors[1] / (ratiofactors[0] + ratiofactors[1]);
		}
	}
	double getratio(int m)
	{
		return ratio[m - firstmarker];
	}

	double getweight(int m, int dir)
	{
		int realm = m - firstmarker;
		halfstate s = relskewfwbw[realm][dir];
		/*const halfstate& bws = relskewfwbw[realm][1];

		for (int k = 0; k < 2; k++)
		{
			s[k] *= bws[k];
		}*/

		double sum = s[0] + s[1];
		if (sum == 0)
		{
			fprintf(stderr, "Renormalization problem for ind %d at %d, dir %d.\n", ind->n, m, dir);
		}
		return s[1] / sum;
	}
};

array<double, TURNBITS> calcskewterms(int marker, relskewhmm* relskews)
{
	array<double, TURNBITS> skewterms;

	for (auto& i : skewterms)
	{
		i = 0;
	}
	for (size_t i = 0; i < skewterms.size() && i < 1; i++)
	{
		individ* ind = reltreeordered[i];
		if (!ind) continue;

		int truei = i;
		// A different ordering regime for flag2 values vs. turn values
		if (i == 0)
		{
			truei = TURNBITS - 1;
		}
		else
		{
			truei -= 1;
		}

		double vals[2] = { relskews[i].getweight(marker, 0), relskews[i].getweight(marker + 1, 1) };

		for (int ix = 0; ix < 2; ix++)
		{
			double hw = ind->haploweight[marker + !ix];
			double hwo = ind->haploweight[marker + ix];
			double val = vals[ix];
			double rh = ind->relhaplo[marker];
			val = hwo; // TODO skip relskewhmm for skewterm logic if we keep this line
			
			double now =
				hw * val * (log(rh) + log(hw) + log(hwo)) +
				(1 - hw) * (1 - val) * (log(rh) + log(1 - hw) + log(1 - hwo)) +
				hw * (1 - val) * (log(1 - rh) + log(hw) + log(1 - hwo)) +
				(1 - hw) * val * (log(1 - rh) + log(1 - hw) + log(hwo));
			// This version assumed the inv was happening on the relskewhmm computed side
			/*double then =
				hw * (1 - val) * (log(rh) + log(hw)) +
				(1 - hw) * val * (log(rh) + log(1 - hw)) +
				hw * val * (log(1 - rh) + log(hw)) +
				(1 - hw) * (1 - val) * (log(1 - rh) + log(1 - hw));*/
			double then =
				(1 - hw) * val * (log(rh) + log(1 - hw) + log(hwo)) +
				hw * (1 - val) * (log(rh) + log(hw) + log(1 - hwo)) +
				(1 - hw) * (1 - val) * (log(1 - rh) + log(1 - hw) + log(1 - hwo)) +
				hw * val * (log(1 - rh) + log(hw) + log(hwo));

			skewterms[truei] -= then - now;
			if (ind->haplocount[marker + ix])
			{
			double gonext = ind->haplobase[marker + ix] / ind->haplocount[marker + ix];
			skewterms[truei] += (gonext - hw) * (hw - 0.5) < 0 ? 25000 : 0;
		}		
		}		
		if (ind->n == 66 && marker >= 4865 && marker <= 4875) printf("CALCSKEWTERMS %d: %d %lf\n", ind->n, marker, skewterms[truei]);
	}

	return skewterms;
}

void updatehaploweights(individ* ind, FILE* out, int iter, std::atomic_int& hitnnn)
{
	vector<bool> allhalf;
	vector<bool> anyinfo;
	vector<bool> cleared;
	vector<int> nudgeme;

	anyinfo.resize(chromstarts.size());
	allhalf.resize(chromstarts.size());
	cleared.resize(chromstarts.size());
	nudgeme.resize(chromstarts.size());
	for (size_t k = 0; k < chromstarts.size(); k++)
	{
		anyinfo[k] = false;
		allhalf[k] = true;
		cleared[k] = false;
		nudgeme[k] = -1;
	}

	unsigned int cno = (unsigned int)-1;
	std::unique_ptr<relskewhmm> relskews;

	for (size_t j = 0; j < ind->haplocount.size(); j++)
	{
		while (cno + 1 < chromstarts.size() && j >= chromstarts[cno + 1])
		{
			cno++;
			for (int k = j; k < ind->haplocount.size() && k < chromstarts[cno + 1] && !anyinfo[cno]; k++)
			{
				if (ind->haplocount[k]) anyinfo[cno] = true;
		}
			if (anyinfo[cno]) relskews = std::make_unique<relskewhmm>(chromstarts[cno], chromstarts[cno + 1], ind);
		}		

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


		if ((anyinfo[cno]) && ind->haploweight[j] && ind->haploweight[j] != 1 /*&& (!ind->founder || ind->children)*/)
		{
			double val;
			if (ind->haplocount[j])
			{
				//val = exp(ind->haplobase[j] / ind->haplocount[j]);
				val = ind->haplobase[j] / ind->haplocount[j];
				val *= (1 - ind->haploweight[j]) / ind->haploweight[j];
			}
			else
			{
				val = 1;
			}

				double relskewterm = 0;
				if (RELSKEWS)
				{
					relskewterm = relskews->getratio(j);
				}
				if (RELSKEWS && false)
				{
					for (int d = -1; d < 1; d++)
					{
						double otherval;
						if (d == -1)
						{
							if (j == chromstarts[cno]) continue;
							otherval = relskews->getweight(j - 1, 0);
						}
						else
						{
							if (j + 1 >= chromstarts[cno + 1]) continue;
							otherval = relskews->getweight(j + 1, 1);
						}
						// arctanh arises from log(1-x) - log(x)
						relskewterm += 2 * atanh((2 * ind->relhaplo[j + d] - 1) * (2 * otherval - 1));
						//if (!isfinite(relskewterm)) printf("Invalid relskewterm %lf for ind %d, marker %d, d %d, relahplo %lf, otherval %lf\n", relskewterm, ind->n, j, d, ind->relhaplo[j + d], otherval);
					}
					// Each direction is counted twice, for two different markers
					if (j > chromstarts[cno] && j + 1 < chromstarts[cno + 1]) relskewterm *= 0.5;
					/*static double minrelskewterm = 0;
					static double maxrelskewterm = 0;
					if (minrelskewterm > relskewterm)
					{
						minrelskewterm = relskewterm;
						printf("New ever lo relskewterm %lf for ind %d, marker %d\n", relskewterm, ind->n, j);
					}
					if (maxrelskewterm < relskewterm)
					{
						maxrelskewterm = relskewterm;
						printf("New ever hi relskewterm %lf for ind %d, marker %d\n", relskewterm, ind->n, j);
					}*/
					// relskewterm *= 0.5;
				}

			double scorea = 1.0 - ind->markersure[j].first;
			double scoreb = 1.0 - ind->markersure[j].second;
			if (ind->markerdata[j].first != ind->markerdata[j].second) scoreb = 1 - scoreb;

			/*double similarity = (scorea * scoreb + (1 - scorea) * (1 - scoreb)) / sqrt((scorea * scorea + (1-scorea) * (1-scorea)) *
			  (scoreb * scoreb + (1-scoreb) * (1-scoreb)));*/
			double sims[2];
			for (int k = 0; k < 2; k++)
			{
				sims[k] = fabs(k - scorea) / fabs(k - scoreb);
				if (sims[k] > 1) sims[k] = 1 / sims[k];
			}
			double similarity = min(sims[0], sims[1]);
			similarity = 0;
			similarity = (scorea * scoreb + (1 - scorea) * (1 - scoreb));

			if (!ind->haplocount[j] || similarity == 1.0)
			{
				ind->haplocount[j] = max(1.0, (double)ind->haplocount[j]);
				ind->haplobase[j] = ind->haploweight[j] * ind->haplocount[j];
			}
			else
			{
				if (similarity >= 1 - maxdiff) similarity = 1 - maxdiff;

				double count = ind->haplocount[j];
				ind->haplobase[j] -= count * ind->haploweight[j];
				count = count - similarity * count;
				ind->haplobase[j] += count * ind->haploweight[j];
				//ind->haplocount[j] = count; // OR use following line:
				ind->haplobase[j] *= ind->haplocount[j] / count;
				if (ind->haplobase[j] < 0) ind->haplobase[j] = 0;
				if (ind->haplobase[j] >= ind->haplocount[j]) ind->haplobase[j] = ind->haplocount[j];
			}

			double ef = exp(0 * -0.01 * iter) * entropyfactor;
			if (j == chromstarts[1] - 1 || j == chromstarts[1] - 2)
			{
				printf("HAPLO DATA ind %3d: %lf %lf %lf %lf\n", ind->n, ind->haploweight[j], ind->haplobase[j], ind->haplocount[j], relskewterm);
			}
			auto gradient = [&](const double& in, double& out, const double)
			{
				double x = in;
				//out = 0;
				out = -(-square(ind->haploweight[j]*ind->haplobase[j])*log(x) + square(ind->haploweight[j]*ind->haplobase[j])*log(1 - x) + square(ind->haploweight[j])*ind->haplobase[j]*ind->haplocount[j]*log(x) - square(ind->haploweight[j])*ind->haplobase[j]*ind->haplocount[j]*log(1 - x) - square(ind->haploweight[j])*ind->haplobase[j]*ind->haplocount[j] - square(ind->haploweight[j]*ind->haplocount[j])*x + square(ind->haploweight[j]*ind->haplocount[j]) + ind->haploweight[j]*square(ind->haplobase[j])*log(x) - ind->haploweight[j]*square(ind->haplobase[j])*log(1 - x) + ind->haploweight[j]*square(ind->haplobase[j]) + 2*ind->haploweight[j]*ind->haplobase[j]*ind->haplocount[j]*x - ind->haploweight[j]*ind->haplobase[j]*ind->haplocount[j]*log(x) + ind->haploweight[j]*ind->haplobase[j]*ind->haplocount[j]*log(1 - x) - ind->haploweight[j]*ind->haplobase[j]*ind->haplocount[j] - square(ind->haplobase[j])*x)/square(ind->haploweight[j]*ind->haplobase[j] + ind->haploweight[j]*ind->haplocount[j]*x - ind->haploweight[j]*ind->haplocount[j] - ind->haplobase[j]*x);
				out +=
					//(ind->haplobase[j] - in * (ind->haplocount[j])) / (in - in * in) +
					 ((1 - similarity) * 1 * (ef * log(1 / in - 1)) + // Entropy term
								     (relskewterm - in) / (in - in * in) * ind->/*haplocount[j]*/descendants);
				/*out =
					((ind->haplobase[j] - in * (ind->haplocount[j]) + (in * in - in)) / (in - in * in) +
					 (1 - similarity) * 1 * ef * (log((1-in)/in*relskewterm/(1-relskewterm)))); // Entropy term*/
				//				if (fabs(out) > 1e8) printf("Large grad %lf, marker %d, ind %d, haplobase %lf, in %lf, relskewterm %lf\n ", out, j, ind->n, ind->haplobase[j], in, relskewterm);
			};


			if (/*ind->children &&*/ (ind->lastinved[cno] == -1 || true) /*&& !ind->pars[0] && !ind->pars[1]*/)
			{
				// Cap the change if the net difference is small/miniscule
				double intended = cappedgd(gradient, ind->haploweight[j], maxdiff / (ind->children + 1), hitnnn, ind->lastinved[cno] != -1);
				//				fprintf(stderr, "HAPLOS: %d %d %lf %lf %lf %lf\n", ind->n, j, intended, ind->haplobase[j], ind->haplocount[j], ind->haploweight[j]);

				//								if ((ind->haploweight[j] - 0.5) * (intended - 0.5) < 0) intended = 0.5;


				/*										if (!(intended < 0.5) && ind->haploweight[j] < 0.5)
				{
				cout << "CROSSOVER " << ind->name << " " << ind->n << " " << j << " " << intended << " " << ind->haploweight[j] << " " << limn << " " << limd1 << std::endl;
				}*/
				//if (similarity < 1 - maxdiff && fabs(intended - 0.5) < 1e-5) intended += 1e-5; 
				/*if (iter <= 3)*/ ind->haploweight[j] = intended;
				if (j > 8755 && j < 8765 && ind->n == 24)
				{
					printf("FLECK: %lf\n", ind->haploweight[j]);
				}


				// Nudging flag currently not respected
				if ((nudgeme[cno] == -1 || fabs(ind->haploweight[nudgeme[cno]] - 0.5) < fabs(ind->haploweight[j] - 0.5)) && ind->haploweight[j] > maxdiff&& ind->haploweight[j] < 1 - maxdiff)
				{
					nudgeme[cno] = j;
				}
			}

			/*							if (ind->haploweight[j] != 0.5)
			{
			allhalf[cno] = false;
			}*/
		}
	}
}

#if XSTDBITSET
using covertype = xstd::bit_set<TOULBARINDCOUNT + 1>;
#else
using covertype = set<int>;
#endif

struct canddata
{
	long long score;
	covertype cover;
	vector<negshiftcand> cands;
};

auto operator < (const canddata& a, const canddata& b) noexcept {
	return std::tie(a.score, a.cover) < std::tie(b.score, b.cover);
}

void fillcandsexists(individ* ind, array<int, 7>& cands, array<bool, 7>& exists)
{
  flat_set<int, less<int>, static_vector<int, 7>> family;
	int temp = ind->n;
	for (int i = 0; i < 7; i++)
	  {
	    exists[i] = 0;
	    cands[i] = 0;
	  }
	cands[6] = temp;
	exists[6] = true;
	family.insert(temp);



	// If incest, we preted the person did not sire anyone the second time they show up in the focus tree

	//test << "Mark: " << mark << " Individ: " << ind->n;

	if (ind->pars[0]) {
		temp = ind->pars[0]->n;
		if (family.insert(temp).second) {//if family member is unique
			cands[0] = temp;
			exists[0] = true;
		}
		//test << " Parent1: " << temp;

		if (ind->pars[0]->pars[0]) {
			temp = ind->pars[0]->pars[0]->n;
			if (family.insert(temp).second) {
				cands[1] = temp;
				exists[1] = true;
			}
			//test << " Parent1's parent1: " << ind->pars[0]->pars[0]->n;
		}
		if (ind->pars[0]->pars[1]) {
			temp = ind->pars[0]->pars[1]->n;
			if (family.insert(temp).second) {
				cands[2] = temp;
				exists[2] = true;
			}
			//test << " Parent1's parent2: " << ind->pars[0]->pars[1]->n;
		}
	}

	if (ind->pars[1]) {
		temp = ind->pars[1]->n;
		if (family.insert(temp).second) {
			cands[3] = temp;
			exists[3] = true;
		}
		//test << " Parent2: " << ind->pars[1]->n;
		if (ind->pars[1]->pars[0]) {
			temp = ind->pars[1]->pars[0]->n;
			if (family.insert(temp).second) {
				cands[4] = temp;
				exists[4] = true;
			}
			//test << " Parent2: " << ind->pars[1]->pars[0]->n;
		}
		if (ind->pars[1]->pars[1]) {
			temp = ind->pars[1]->pars[1]->n;
			if (family.insert(temp).second) {
				cands[5] = temp;
				exists[5] = true;
			}
			//test << " Parent2: " << ind->pars[1]->pars[1]->n;
		}
	}
}

long long computesumweight(const int m, const vector<int>& tf, const vector<vector<clause>>& toulinput, covertype& cover)
{
	long long sumweight = 0;
	int numviol = 0;
	for (const clause& c : toulinput[m])
	{
		bool viol = true;
		bool anyswitch = false;
		for (int val : c.cinds)
		{
			int ind = val < 0 ? -val : val;
			if (val < 0)
			{
				anyswitch = true;
			}
			if (tf[ind - 1] == (val > 0))
			{
				viol = false;
			}
		}
		if (viol)
		{
			numviol++;
			sumweight += c.weight;
			if (anyswitch) for (int val : c.cinds)
			{
				int ind = val < 0 ? -val : val;
				cover.insert(ind);
			}
		}
	}

	if (numviol != dous.size())
	{
		fprintf(stderr, "Wrong number of violated clauses %d/%d at %d\n", numviol, toulinput[m].size(), m);
	}
	return sumweight;
}

template<class C1, class C2>
bool smartincludes(const C1& set1, const C2& set2)
{
	if (set1.size() < set2.size()) return false;
	return std::includes(set1.begin(), set1.end(), set2.begin(), set2.end());
}

#if XSTDBITSET
template<>
bool smartincludes<covertype, covertype>(const covertype& set1, const covertype& set2)
{
	return set2.is_subset_of(set1);
}
#endif



vector<canddata> computecandcliques(const int m, const vector<int>& tf, const vector<vector<clause>>& toulinput, long long bias)
{
	vector<canddata> result;

	int numviol = 0;
	for (const clause& c : toulinput[m])
	{
		bool viol = true;
		bool anyswitch = false;
		for (int val : c.cinds)
		{
			int ind = val < 0 ? -val : val;
			if (val < 0)
			{
				anyswitch = true;
			}
			if (tf[ind - 1] == (val > 0))
			{
				viol = false;
			}
		}
		if (viol)
		{
			numviol++;
						
			if (anyswitch)
			{
				int useindex = -1;
				for (int i = 0; i < result.size(); i++)
				{
					for (int val : c.cinds)
					{
						int ind = val < 0 ? -val : val;
						if (result[i].cover.find(ind) != result[i].cover.end())
						{
							if (useindex == -1)
							{
								useindex = i;
							}
							else
							{
#if XSTDBITSET								
								result[useindex].cover |= result[i].cover;
#else							
								result[useindex].cover.merge(result[i].cover);
#endif								
								result[useindex].cands.insert(result[useindex].cands.end(), result[i].cands.begin(), result[i].cands.end());
								result[useindex].score += result[i].score;								
								//								printf("Merging %d and %d\n", useindex, i);
								result.erase(result.begin() + i);
								i--;
							}
							break;
						}
					}
				}
				if (useindex == -1)
				{
					result.emplace_back();
					useindex = result.size() - 1;
					result[useindex].score = 0;
				}
				result[useindex].score -= bias - c.weight;
				for (int val : c.cinds)
				{
					int ind = val < 0 ? -val : val;
					//					printf("Doing insert of %d into %d\n", ind, useindex);
					//					fflush(stdout);
#if XSTDBITSET								
					bool isnew = !result[useindex].cover.contains(ind);
					result[useindex].cover.insert(ind);
#else							
					bool isnew = result[useindex].cover.insert(ind).second;
#endif														
					if (isnew && tf[ind - 1])
					{
						result[useindex].cands.emplace_back(getind(ind), c.weight, m);
					}
				}				
			}
		}
	}

	if (numviol != dous.size())
	{
		fprintf(stderr, "Wrong number of violated clauses %d/%d at %d\n", numviol, toulinput[m].size(), m);
	}

	return result;
}

void createtoulbarfile(const string toulin, long long maxweight, vector<clause>& clauses)
{
	std::fstream infile(toulin, ios::out | ios::in | ios::trunc);
	if (!infile) {
		perror("Toulbars input file failed to open to be written to because: ");
	}

	infile << "c In Weigthed Partial Max-SAT, the parameters line is 'p wcnf nbvar nbclauses top'\n";
	infile << "c Where p is the weight (true maxweight" << maxweight << ")\n";
	infile << "c nbvar is the number of a variables appearing in the file\n";
	infile << "c nbclauses is the exact number of clauses contained in the file\n";
	infile << "c see http://maxsat.ia.udl.cat/requirements/\n";

	int nbclauses = (int)clauses.size();
	int nbc = nbclauses;
	//cout<<"nbvar: " <<nbvar<< "\n"; // problem solving
	//cout<<"nbclauses: " <<nbc<< "\n"; // problem solving
	infile << "p wcnf " << TOULBARINDCOUNT << " " << nbc << "\n"; //" " <<std::numeric_limits<int>::max()<<"\n";

	for (clause& c : clauses) {
		if (c.weight <= 0)
		{
			fprintf(stderr, "Non-positive weight, weight %lld, maxweight %lld\n", c.weight, maxweight);
			abort();
		}
		infile << c.weighttostring() << c.clausetostringstream().rdbuf() << " 0\n";
		//infile<< toulInput[m][g].toString() << "\n";
		//cout<<"TEST " <<toulInput[m][g].toString()<< "\n"; // problem solving
	}
}

typedef map<pair<individ*, individ*>, map<int, array<double, 8> > > nsmtype;

void parentswapnegshifts(nsmtype& nsm)
{
	vector<pair<double, boost::tuple<individ*, individ*, int, int> > > allnegshifts;
	flat_map<individ*, double> bestshift;
	for (auto i = nsm.begin(); i != nsm.end(); i++)
	{
		for (auto j = i->second.begin(); j != i->second.end(); j++)
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
		if (bestshift[allnegshifts[k].second.get<0>()] < allnegshifts[k].first&&
			bestshift[allnegshifts[k].second.get<1>()] < allnegshifts[k].first)
		{
			bestshift[allnegshifts[k].second.get<0>()] = allnegshifts[k].first;
			bestshift[allnegshifts[k].second.get<1>()] = allnegshifts[k].first;
			unsigned int c = 0;
			while (c < chromstarts.size() && (chromstarts[c] <= allnegshifts[k].second.get<3>()))
			{
				c++;
			}

			int phase = allnegshifts[k].second.get<2>();
			individ* inds[2] = { allnegshifts[k].second.get<0>(), allnegshifts[k].second.get<1>() };
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

						for (size_t k = 0; k < inds[z]->kids.size(); k++)
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
}


void cheat_system(const char* line)
{
  char* const argv[4] = {"sh", "-c", (char*) line, 0};
  pid_t respid;
  int wstatus;
  posix_spawn(&respid, "/bin/sh", nullptr, nullptr, argv, environ);

  while (waitpid(respid, &wstatus, 0) == -1 && errno == EINTR) {};
}

void mergebestcands(std::set<canddata>& bestcands, int ceiling, int clearto)
{
	bool toolarge;
	bool delprev = false;
	do
	{
		toolarge = false;
		vector<remove_reference<decltype(bestcands)>::type::iterator> toremove;
		for (auto i = bestcands.rbegin(); i != bestcands.rend() && !toolarge; i++)
		{
			if (delprev)
			{
				// base is the forward iterator one step later, i.e. our previous entry
				toremove.push_back(i.base());
			}

			delprev = false;				
			for (auto j = bestcands.begin(); *j < *i && !toolarge; j++)
			{
				int covered = 0;
				bool fullcover = false;
#if XSTDBITSET
				if (i->cover.intersects(j->cover))
				{
					covered = 1;
					if (i->cover == j->cover) fullcover = true;
				}
#else
				auto jj = j->cover.begin();
				for (auto ind : i->cover)
				{
					while (jj != j->cover.end() && *jj < ind)
					{
						jj++;
					}
					if (jj == j->cover.end()) break;

					if (*jj == ind)
					{
						covered++;
						if (i->cover.size() != j->cover.size()) break;
					}
				}
				fullcover = covered == j->cover.size() && covered == i->cover.size();
#endif				

				if (!covered)
				{
					canddata newcand{ i->score + j->score, i->cover, i->cands };
#if XSTDBITSET
					newcand.cover |= j->cover;
#else		
					newcand.cover.insert(j->cover.begin(), j->cover.end());
#endif							
					newcand.cands.insert(newcand.cands.end(), j->cands.begin(), j->cands.end());
					bestcands.insert(std::move(newcand));
					// The greedy part, replaced by maximum size limit
					// // break;
				}
				if (fullcover)
				{
					delprev = true;
					break;
				}

				if (bestcands.size() > ceiling) break;
			}
			if (bestcands.size() > ceiling)
			{
				toolarge = true;
				break;
			}
		}
		for (auto i : toremove)
		{
			bestcands.erase(i);
		}
		while (bestcands.size() > ceiling)
		{
			bestcands.erase(--bestcands.end());
		}
	} while (toolarge);
	while (bestcands.size() > clearto)
	{
			bestcands.erase(--bestcands.end());
	}
}

// The actual walking over all chromosomes for all individuals in "dous"
// If "full" is set to false, we assume that haplotype inference should be done, over marker positions.
// A full scan is thus not the iteration that takes the most time, but the scan that goes over the full genome grid, not only
// marker positions.
template<bool full, typename reporterclass> void doit(FILE* out, bool printalot
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
	vector<vector<array<float, 2> > > realgeno;

	realgeno.resize(dous.size());

	for (size_t j = 0; j < dous.size(); j++)
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

	nsmtype nsm;

	for (int i = 0; i < INDCOUNT; i++)
	{
		individ* ind = getind(i);
		if (!ind) continue;
		ind->children = 0;
		ind->kids.clear();
		if (!ind || !ind->haplocount.size()) continue;

		for (size_t j = 0; j < ind->haplocount.size(); j++)
		{
			ind->haplocount[j] = 0.0f;
			ind->haplobase[j] = 0.0f;
	}

#ifdef F2MPI
		broadcast(world, ind->haploweight, 0);

		//		fprintf(out, "C%d:%d\n", world.rank(), i);

		world.barrier();
#endif

}

	for (size_t j = 0; j < dous.size(); j++)
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

	vector<set<negshiftcand>> negshiftcands;
	vector<covertype> negshiftcovers;
	negshiftcands.resize(chromstarts.size());
	negshiftcovers.resize(chromstarts.size());
	vector<omp_lock_t> markerlocks;
	markerlocks.resize(markerposes.size());
	for (auto& lock : markerlocks)
	  {
	    omp_init_lock(&lock);
	  }
	// Create a vector where each element corresponds to a marker and
	//contains a referense to a vector containing all the clauses for said marker
	//EBBA also: Here starts the parallels, investigate names
#if DOTOULBAR
	static vector<vector<clause>> toulInput;
	toulInput.resize(markerposes.size());
	for (auto& subv : toulInput)
	  {
	    subv.clear();
	  }
#endif

	for (size_t i = 0; i < chromstarts.size() - 1; i++)
	{
		//printf("Chromosome %d\n", i + 1);

		// The output for all individuals in a specific iteration is stored, as we have parallelized the logic and 
		// want the output to be done in order.
		vector<vector<char> > outqueue;
		outqueue.resize(dous.size());

		long long maxweight = 0;

#pragma omp parallel for schedule(dynamic,1)
		for (size_t j = 0; j < dous.size(); j++)
		{
#ifdef F2MPI
			if (j % world.size() != world.rank()) continue;
#endif
			//			realgeno[j].resize(markerposes[chromstarts[1] - 1] + 1);

			generation++;
			threadblock tborig;
			threadblock tb = tborig;

			resizecaches();

			if (dous[j]->markerdata.size())
			{
				dous[j]->lastinved[i] = -1;
				int qstart = -1000 - chromstarts[i];
				int qend = -1000 - chromstarts[i + 1];
				//TODO: WRONG RANGE
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

				// Special optimization hardcoded for this population structure, eagerly skipping flags that do not
				// correspond to any inheritance, i.e. if not the full pedigree of 6 individuals back is present.
				auto [shiftignore, flag2ignore] = fixtrees(dous[j]);

				bool skipsome = false;
				for (size_t u = 0; u < reltree.size() && !skipsome; u++)
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
					qstart = (int)markerposes[chromstarts[i]];
					qend = (int)markerposes[chromstarts[i + 1] - 1] + 1;
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
				/*if (!anygood)
				{
					shiftignore = 6;
					flag2ignore = 0;
					}*/

				double factor = -1e15;
				double factors[NUMSHIFTS];
				for (shiftflagmode = shifts; shiftflagmode < shiftend; shiftflagmode++)
				{
					if (shiftflagmode & shiftignore) //continue;
						factors[shiftflagmode] = -1e30;
					else
						factors[shiftflagmode] = dous[j]->doanalyze<noneturner>(tb, none, chromstarts[i], chromstarts[i + 1] - 1, NONESTOP, -1, false, 0, -40000 + factor);
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
				//factor -= log(1 << countf2i);

				// This output can get ugly due to race conditions. One shouldn't rely on it.
				printf("%d,%03d,%03d: %lf\t", dous[j]->n, flag2ignore, shiftignore, factor);
				factor += log(realfactor);
				printf("%lf %d\n", factor, shiftend);
				fflush(stdout);
				if (std::isnan(factor) || factor < MINFACTOR) continue;

				// Walk over all chromosome positions, whether it be markers (negative q values <= -1000) or grid positions
				for (int q = qstart; q != qend; q += qd)
				{
					reporterclass reporter;
					//double mwvals[NUMTYPES][NUMTYPES] = {0};
					//double mwfvals[NUMTYPES] = {0};
					double mwvals[1][1];
					double mwfvals[1];
					double mwval[4] = { 0 };
					constexpr double unusualstate = -200;

					for (int g = 0; g < NUMTYPES; g++)
					{
						for (shiftflagmode = shifts; shiftflagmode < shiftend; shiftflagmode++)
						{
							if (shiftflagmode & shiftignore) continue;
							if (factor - factors[shiftflagmode] > 40) continue;
							if (q <= -1000 && DOREMAPDISTANCES)
							{
								double val;
								for (int g2 = 0; g2 < NUMTYPES; g2++)
								{
									val = dous[j]->doanalyze<noneturner>(tb, none, chromstarts[i], chromstarts[i + 1] - 1, twicestop(q, g, g2),
										-1, true, 0, -5000.0 + factor) - factor;

									val = exp(val);
									mwvals[g][g2] += val;
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
								/*vectortrackpossible<2, 0, 0> vtp(-q-1000, tpdefault, flag2);
								printf("%lf\n", vtp.result[0]);*/
								if (ignoreflag2(flag2, g, shiftflagmode, q, flag2ignore, relmap, relmapshift)) continue;
								//if (flag2 & (flag2ignore)) continue;

								int firstpar = 0;
								double val;

								if (DOIMPOSSIBLE && q <= -1000 && flag2 != -1)
								{
									// Do a lookup in the impossible structure. If we are at a marker position, those branches that are impossible will have an exact
									// score of 0.
									// If we are NOT at a marker, this information can still be used in trackpossible, but not in here, as there is some recombination
									// added between the marker and fixated position we're checking.
									for (int a = 0; a < 2; a++)
									{
										if (a) firstpar = !firstpar;
										const unsigned int genwidth = (1 << (NUMGEN - 1));

										int f2n = ((flag2 ^ shiftflagmode) & 1);

										int upflagr = upflagit(g, firstpar, genwidth);
										int upflag2r = upflagit(flag2 >> 1, firstpar, genwidth);
										int upshiftr = upflagit(shiftflagmode >> 1, firstpar, genwidth >> (NUMGEN - NUMSHIFTGEN));
										int marker = -(q + 1000);
										int impossibleval = generation * markerposes.size() + marker;

										if (impossible[shiftflagmode & 1][firstpar][f2n][upflagr][upflag2r + 1][upshiftr][marker & 3] == impossibleval)
										{
											std::cerr << "IMPOSSIBLE FACILITY USED" << std::endl;
											goto continueloop;
										}
									}
								}


								// This is the big main call								
								val = dous[j]->doanalyze<noneturner>(tb, none, chromstarts[i], chromstarts[i + 1] - 1, classicstop(q, g),
									flag2, true, 0, unusualstate + factor) - factor;

								if (_finite(val) && val > unusualstate)
								{
									// shift mode not included, this is the "real" f2n, indicating what value
									// in the marker pair is used, not the strand phase (strand phase is flag2 xored
									// with the other stuff)

									val = exp(val);
									int marker = -q - 1000;

									int mapval = 0;
									double outmapval = dous[j]->trackpossible<false, true>(tb, UnknownMarkerVal, 0, marker, g * 2, flag2, *(tb.shiftflagmode), trackpossibleparams(0, &mapval));
									double sidevals[2][2] = { 0 };
									double sidevalsums[2] = { 0 };
									double homozyg[2] = { 0 };

									if (DOINFPROBS)
									{
										for (int side = 0; side < 2; side++)
										{
											#pragma ivdep
											for (int i = 1; i <= 2; i++)
											{
												MarkerVal markervalcopy = i * MarkerValue;
												double sideval = dous[j]->trackpossible<GENOSPROBE, false>(tb, markervalcopy, 0, marker, g * 2 + side, flag2 ^ side, *(tb.shiftflagmode), trackpossibleparams(0, nullptr));
												sidevals[side][markervalcopy.value() - 1] += sideval;
												sidevalsums[side] += sideval;
											}
										}

										#pragma ivdep
										for (auto markerval : { 1 * MarkerValue, 2 * MarkerValue })
										{
											auto markervalcopy = markerval;
											double sideval = dous[j]->trackpossible<HOMOZYGOUS, false>(tb, markervalcopy, 0, marker, g * 2, flag2, *(tb.shiftflagmode), trackpossibleparams(0, nullptr));
											homozyg[markerval.value() - 1] += sideval;
										}
									}
									//									double pairvalsum = accumulate(begin(pairvals), end(pairvals), 0.0);

									/*									for (int i = 0; &pairvals[i] != end(pairvals); i++)
																		{
																		  reporter.addval(q, i, g, flag2, val * pairvals[i] / pairvalsum);
																		}*/
																		/*/*if (!outmapval && val > 1e-3)
																		{
																			//std::cerr << "ERROR TERROR " << -q - 1000 << " " << g * 2 << " " << flag2 << " " << *(tb.shiftflagmode) << "\t" << val << "\n";
																		}
																		if (mapval < 0 || mapval > 2)
																		  {
																		  std::cerr << "Incorrect mapval " << -q - 1000 << " " << dous[j]->n << std::endl;
																		  }*/
																		  //reporter.addval(q, mapval, g, flag2, val);
									if (!full && HAPLOTYPING)
									{
										dous[j]->updatehaplo(tb, marker, g, flag2, val);

										if (DOINFPROBS)
										{
											for (int side = 0; side < 2; side++)
											{
												for (auto markerval : { 1 * MarkerValue, 2 * MarkerValue })
												{
													//std::cout << "EXTREME VETTING IND " << dous[j]->n << " MARKER " << marker << ":" << markerval.value() << ", fl2" << flag2 << ", sfm " << *(tb.shiftflagmode) << ", VAL: " << val << " SIDEVAL " << sidevals[side][markerval.value() - 1] << ", SIDEVALSUM " << sidevals[side][markerval.value() - 1] << std::endl;
													double updateval = val * sidevals[side][markerval.value() - 1] / sidevalsums[side];
													dous[j]->trackpossible<GENOS, false>(tb, markerval, 0, marker, g * 2 + side, flag2 ^ side, *(tb.shiftflagmode), trackpossibleparams(updateval, nullptr));
												}
											}

											double hzs = 0;
											for (auto markerval : { 1 * MarkerValue, 2 * MarkerValue })
											{
												dous[j]->homozyg[marker][markerval.value() - 1] += val * homozyg[markerval.value() - 1] / sidevalsums[0];
												hzs += homozyg[markerval.value() - 1];
											}
											if (hzs - 1e-7 > sidevalsums[0] * (1 + 1.e-7)) fprintf(stdout, "HZS MISMATCH %d %d %d %d %lf %lf %lf\n", marker, g * 2, flag2, *(tb.shiftflagmode), hzs, sidevalsums[0], sidevalsums[1]);
										}
									}
								}
							continueloop:;
							}
						}
					}

					// Coordinate estimated distances from all individuals.
					if (q <= -1000 && DOREMAPDISTANCES)
					{
						double colsums[NUMTYPES] = { 0 };
						double rowsums[NUMTYPES] = { 0 };
						double acc3 = 0;
						double acc4 = 0;
						calcdistancecolrowsums(mwvals, rowsums, colsums, acc3, acc4, mwfvals);

						double relinfo1 = 0;
						double relinfo2 = 0;
						double infosum[2] = { 0 };
						double recprob[2];
						double dist = markerposes[-q - 1000 + 1] - markerposes[-q - 1000];

						for (int k = 0; k < 2; k++)
						{
							recprob[k] = 0.5 * (1.0 - exp(actrec[k][-q - 1000 + 1] * (dist)));
							recprob[k] = max(1e-8, recprob[k]);
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
									bool switched = ((bool)(mask & (1 << i)));
									double expected = !switched ? 1.0 - recprob[TYPESEXES[i]] : recprob[TYPESEXES[i]];

									mwval[TYPESEXES[i] * 2 + switched] += (val / expected) * corr;
									infosum[TYPESEXES[i]] += corr;
								}
							}
						}

						double summw = (mwval[0] + mwval[1] + mwval[2] + mwval[3]);
						for (int z = 0; z < 2; z++)
						{
							infosum[z] /= NUMTYPES * NUMTYPES;
						}
						if (acc3 == 0)
						{
							acc3 = 1;
							infosum[0] = 0;
							infosum[1] = 0;
						}


						double delta = 0;
#pragma omp critical(markerweights)					       
						{
							for (int t = 0; t < 4; t++)
							{
								double dval = (mwval[t] / acc3 - infosum[t / 2]) /*/ summw*/
									* fabs(((t & 1) ? 0.0 : 1.0) - recprob[t / 2])
									/** (-1 + (t & 1) * 2)*/;

								markerweight[-q - 1000][t] += dval;
								if (t == 0) delta = dval;
							}
							markerweight[-q - 1000][4] = min(markerweight[-q - 1000][4], acc4);
						}
						if (delta < -0.1 && q > qend + 2) fprintf(out, "RECOMB:\t%d\t%d\t%lf\n", dous[j]->n, q, delta);
					}

					// TODO: NEGSHIFT DOESN'T TAKE RELMAP FLAG2 RESTRICTIONS INTO ACCOUNT
					// Consider doing haplotype reversal from a specific position and all the way down.
					if (HAPLOTYPING && !early && !full /*&& dous[j]->gen >= 0*/)
					{
						int marker = -q - 1000;

						double rawvals[NUMTURNS][NUMSHIFTS];
						double rawervals[NUMTURNS][NUMSHIFTS];
						double sumnegval[TURNBITS] = { 0 };
						bool validg[NUMTURNS] = { 0 };

						for (int g = 0; g < NUMTURNS; g++)
						{
							for (unsigned int s = 0; s < NUMSHIFTS; s++)
							{
								rawvals[g][s] = -1;
								rawervals[g][s] = std::numeric_limits<double>::quiet_NaN();
							}
						}

						for (int g = 0; g < NUMTURNS; g++)
						{
							if (g & (flag2ignore >> 1)) continue;

#if !DOTOULBAR
							int c = 0;
							for (int p = 0; p < TYPEBITS + 1; p++)
							{
								if (g & (1 << p)) c++;
							}

							if (c > 1) continue;
#endif

							double skewterm = 0;
							for (int i = 0; i < TURNBITS; i++)
							{
								if (g & (1 << i))
								{
									//									skewterm += skewterms[i];
								}
							}
							aroundturner turn(g);
							for (shiftflagmode = shifts; shiftflagmode < shiftend; shiftflagmode++)
							{
								if (shiftflagmode & shiftignore) continue;
								// If we are above this limit, we are shifting shift mode
								// within the range and cannot use this heuristic of the
								// aggregated probability to know anything... anything at all!
								// (In other cases, our simulated crossover event amounts to 
								// far below a exp(20) change in probability.)
								int g2 = g;
								if (!g) g2 = (1 << 15) - 1;

								int oldshift = shiftflagmode;
								rawervals[g][oldshift] = dous[j]->doanalyze<aroundturner>(tb, turn, chromstarts[i],
									chromstarts[i + 1] - 1, classicstop(q, -1), -1, true, 0, -50000 + factor) - factor;

								shiftflagmode = oldshift;

								/*								if (c > 1) continue;*/
								validg[g] = true;

								double expval = exp(rawervals[g][oldshift]);

								if (!isfinite(expval))
								{
									if (expval > 1e300)
										rawvals[g][oldshift] = 1e300;
									else
										rawvals[g][oldshift] = 0;
									continue;
								}

								rawvals[g][oldshift] = expval;

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

#if DOTOULBAR
						// Changed for EBBA's degree project
						// Reminder: NUMTURNS is a bitflag with relevant individuals that should be changed
						// remaining individuals should not be changed (= negative numbers)
						// All individuals ancestors (assuming max 3 gens) should be included in clause
						// The weight is the sum of rawvals (?)
						// The current individual is dous[j]
						// structure of container: vector<vector<clause*>> toulInput;

						{
							//Markers stored in q (see line 2844 in original cnF2freq)
							//such that
							int mark = -q - 1000;

							//Do we have parents and grand parents?
							//Store their identifying numbers in an array.
							//Hard coded for max 3 gens.
							//std::fstream test("test2.txt", ios::out | ios::in | ios::trunc);//TEST//TEST
							array<int, 7> cands;
							array<bool, 7> exists;
							fillcandsexists(dous[j], cands, exists);

							//test << "Number of individuals:  " << numbind << " End of input into cands \n ";//remember, incesters only counted once

																											//Use structure clause to store weight and values.

																											//loop over NUMTURNS
																											// get relevant weights and individuals, insert to container

																											//for (int g = 0; g < NUMTURNS; g++) {
							double normfactor = 0;
							const int shiftignore2 = shiftignore;
							auto computew = [&rawervals, &normfactor, &j, shiftignore = shiftignore2, shiftend, shifts] (int g)
							{
								double w = MINFACTOR;
								int count = g > 0 ? std::popcount((unsigned int) g) : 1;

								for (int s = shifts; s < shiftend; s++) {
									if (s & shiftignore) continue;
									w = max(w, rawervals[g][s]);
								}
								double normsum = 0;
								for (int s = shifts; s < shiftend; s++) {
									if (s & shiftignore) continue;
									normsum += exp(rawervals[g][s] - w);
								}
								w += log(normsum);
								w -= normfactor;

								return w * dous[j]->descendants/*/ count*/;
							};
							normfactor = computew(0);
							static_vector<clause, NUMTURNS> subInput;
							long long submax = 0;

							for (int g = 0; g < NUMTURNS; g++) {
								if (g & (flag2ignore >> 1)) continue;

								static_vector<int, NUMTURNS> claus;
								int shiftcount = 0;
								for (int b = 0; b < TURNBITS; b++) {
									if (exists[b]) {
										//if (find(claus.begin(), claus.end(), cind | -cind) == claus.end()){// avoid inbreeding results.
										if (g & (1 << b)) {
											claus.push_back(-cands[b]);
											shiftcount++;
										}
										else {
											claus.push_back(cands[b]);
										}
										//}
									}
								}
								double w = computew(g);
								//Now simply construct a clause type and send it to the right marker
								clause c;
								if (isfinite(w))
								{
								}
								else
								{
									if (w < 0)
									{
										w = -1000000;
									}
									else
									{
										w = 25000;
									}
								}
								w = std::clamp(w, -1000000., 25000.);
								c.weight = w * WEIGHT_DISCRETIZER;
								c.weight -= g;
								c.cinds = claus;
								//test << "Mark: " << mark << "ClausToString: " << c.toString() << " Current maxweight: " << maxweight << endl;//TEST
								if (c.weight > submax) {
									submax = c.weight;

								}
								subInput.push_back(c);
							}
							{
							  omp_set_lock(&markerlocks[mark]);
							  toulInput[mark].insert(toulInput[mark].end(), subInput.begin(), subInput.end());
							  omp_unset_lock(&markerlocks[mark]);
#pragma critical negshifts
							  if (submax > maxweight) maxweight = submax;
							}
							}
#else
						updatenegshifts(validg, shifts, shiftend, shiftignore, rawvals, j, marker);
#endif
						}

					//reporter.report(outqueue[j]);

					if (!full)
					{
						int marker = -q - 1000;

						double sum = 0;
#pragma ivdep						
						for (auto p : infprobs[dous[j]->n][0])
						{
							sum += p.second;
						}
						sum = 1 / sum;

						for (int i = 0; i < 2; i++)
						{
							dous[j]->homozyg[marker][i] *= sum;
						}

						{
						  omp_set_lock(&markerlocks[marker]);
#pragma ivdep
#pragma GCC ivdep
							for (size_t k = 0; k < reltree.size(); k++)
							{
								int i = reltree[k]->n;
								moveinfprobs(i, k, marker, sum, dous[j]->descendants);
								movehaplos(i, k, marker, dous[j]->descendants);
							}
						  omp_unset_lock(&markerlocks[marker]);
						}
					}
						}

					}
			fprintf(stderr, "LAST: %d %s\n", dous[j]->n, dous[j]->name.c_str());
			fflush(stderr);

			if (RELSKEWS && !RELSKEWSTATES)
			{
				vector<relskewhmm> relskews;
				for (individ* indr : reltreeordered)
				{
					if (indr != nullptr)
					{
						relskews.emplace_back(chromstarts[i], chromstarts[i + 1], dous[j]);
					}
					else
					{
						relskews.emplace_back();
					}
				}


				array<double, TURNBITS> skewterms;

#if DOTOULBAR
				long long submax = 0;
				for (int marker = chromstarts[i]; marker < chromstarts[i + 1] - 1; marker++)
				{
					skewterms = calcskewterms(marker, &relskews[0]);
					double w = skewterms[TURNBITS - 1] * 0.5;	
					if (!isfinite(w) || fabs(w) > 25000)
					{
						if (w < -25000)
						{
							w = -25000;
						}
						else
						{
							w = 25000;
						}
					}

					omp_set_lock(&markerlocks[marker]);
					for (clause& c : toulInput[marker])
					{
						if (*(--c.cinds.end()) == -dous[j]->n)
						{
							c.weight -= w * dous[j]->descendants * WEIGHT_DISCRETIZER;
							if (c.weight > submax) {
								submax = c.weight;
							}
						}
					}
					omp_unset_lock(&markerlocks[marker]);
				}
#pragma critical negshifts
					if (submax > maxweight) maxweight = submax;
#endif
				}
			}

		//Print all information to seperate files EBBA
		// stored by marker in toulIn  vector<vector<clause>>
		//Then run toulbar and save best solution in relevant negshift vectors
		//Remember: typedef boost::tuple<individ*, double, int> negshiftcand;

#if DOTOULBAR
		std::set<canddata> bestcands;
		const int maxcandcount = 1000;
		//for (int m=0; m < (int) toulInput.size(); m++ ){//TODO change so that it is valid for more than one chromosome
		int donext = 1;
		bool solexists = false;

#pragma omp parallel for schedule(dynamic,1) firstprivate(donext)
		for (unsigned int m = chromstarts[i]; m < chromstarts[i + 1] - 1; m++) {
			if ((m % 1) == (iter % 1)) donext++;
			if (!donext) continue;
			donext--;
			std::string tid = boost::lexical_cast<std::string>(omp_get_thread_num());
			std::string toulin(tmppath + "/" + std::string("toul_in") + tid + ".wcnf");
			std::string toulout(/*tmppath + "/" + std::string("toul_out") + tid + ".txt"*/ "/dev/null");
			std::string sol(tmppath + "/" + std::string("sol") + tid + ".sol");

			long long fakegain = 0;
			long long fakegainterm = 0;
			long long prevlast = -1;
#if XSTDBITSET			
			covertype uofakecover;
#else			
			unordered_set<int> uofakecover;
#endif			
			for (clause& c : toulInput[m])
			{
				int mainind = abs(*(--c.cinds.end()));
				if (mainind != prevlast)
				{
					fakegain += fakegainterm;
					fakegainterm = 0;
					prevlast = mainind;
				}

				if (c.weight > fakegainterm)
				{
					fakegainterm = c.weight;
					for (int i : c.cinds)
					{
						uofakecover.insert(abs(i));
					}
				}
				c.weight = maxweight - c.weight + 1;
			}
			fakegain += fakegainterm;
			if (!fakegain)
			{
				donext++;
				continue;
			}

#if XSTDBITSET
			covertype& fakecover = uofakecover;
#else			
			vector<int> fakecover;
			fakecover.reserve(uofakecover.size());
			std::copy(uofakecover.begin(), uofakecover.end(), std::back_inserter(fakecover));
			std::sort(fakecover.begin(), fakecover.end());
#endif			

			fakegain = -fakegain;
			bool skippable;

#pragma omp critical(negshifts)
			{
				skippable = !fakegain || (bestcands.size() >= maxcandcount && (--bestcands.end())->score < fakegain);
				for (const canddata& elem : bestcands)
				{
					if (skippable) break;
					/*					int anym = -1;
					for (negshiftcand cand : elem.cands)
					  {
					    int m = cand.get<2>();
					    if (anym == -1)
					      {
						anym = m;
					      }
					    if (m != anym)
					      {
						goto nextcand;
					      }
					      }*/
					if (elem.score > fakegain)
					{
						// Sorted order, no use in checking further, we're done here
						break;
					}

					if (smartincludes(fakecover, elem.cover))
					{
							skippable = true;
						}
				nextcand:;
				}
			}

			if (skippable)
			{
				donext++;
				continue;
			}

			createtoulbarfile(toulin, maxweight, toulInput[m]);

			// Run toulbar2 with partitioning rules that match the fact that we use small pedigrees
			// Feed previous sol file as certificate/starting point
			string str = "toulbar2 " + toulin /*+ (solexists ? " " + sol + " " : "") */+ " -p=8 -O=-1 -m=1 -w=" + sol + " -s > " + toulout; //works as in it runs, not as in it actually does what we want
																			 //string str = "toulbar2 brock200_4.clq.wcnf -m=1 -w -s";//TEST


																			 // Convert string to const char * as system requires
			const char* command = str.c_str();
			cheat_system(command);
			solexists = true;

			//Read outfile and store best result in negshift
			std::fstream touloutput(sol, ios::in);
			//read from file to string of 0s and 1s
			int rawinput;
			vector<int> tf;
			while (touloutput >> rawinput) {
				tf.push_back(rawinput);
				if (rawinput)
				{
					fprintf(stdout, "Ind %d marker %d shift indicated\n", tf.size(), m);
				}
			}

			//Identify all violated clauses, elimination step means optimum cost data from toulbar not usable.
			/*set<int> cover;
			long long sumweight = computesumweight(m, tf, toulInput, cover);

			canddata data;
			data.score = sumweight - (maxweight + 1) * (long long)dous.size();
			data.cover = std::move(cover);
			for (size_t g = 0; g < tf.size(); g++) {
				if (tf[g])
				{
					data.cands.emplace(getind(g + 1), sumweight, m);
				}
			}
			if (data.cover.size() && data.score >= 0)
			{
				fprintf(stderr, "ERROR: POSITIVE SCORE %lld %lld %lld %d %d\n", data.score, sumweight, maxweight, data.cover.size(), data.cands.size());
			}*/

			
			for (canddata data : computecandcliques(m, tf, toulInput, maxweight + 1))
			{
				if (data.cover.size() && data.score < 0)
				{
					fprintf(stdout, "Candidate at marker %d with score %lld (%lld)\n", m, data.score, maxweight);
	#pragma omp critical(negshifts)
				{				
					vector<decltype(bestcands)::iterator> toremove;
					bool addme = true;
					for (const canddata& elem : bestcands)
					{
						bool firstmatch = false;
						if (smartincludes(elem.cover, data.cover))
						{
							if (data.score <= elem.score)
							{
								// Stupid to search for ourselves...
								toremove.push_back(bestcands.find(elem));
							}
							firstmatch = true;
						}
						if (smartincludes(data.cover, elem.cover))
						{
							if (elem.score <= data.score)
							{
								addme = false;
							}
							if (firstmatch) break;
						}
					}

					if (addme)
					{
						for (auto i : toremove)
						{
							bestcands.erase(i);
						}
						bestcands.insert(std::move(data));
					}
					if (bestcands.size() > maxcandcount)
					{
						// Do some mergings, then make room for new ones
						mergebestcands(bestcands, maxcandcount * 2, maxcandcount * 0.5);
					}
				}
			  }
			}
		}

		mergebestcands(bestcands, maxcandcount, maxcandcount);

		if (bestcands.size())
		{
			negshiftcands[i].clear();
			negshiftcands[i].insert(bestcands.begin()->cands.begin(), bestcands.begin()->cands.end());
			negshiftcovers[i] = bestcands.begin()->cover;
		}
		else
		{
			negshiftcands[i].clear();
			negshiftcovers[i].clear();
		}
#endif

		for (size_t j = 0; j < outqueue.size(); j++)
		{
			outqueue[j].push_back(0);
			if (out && printalot) fprintf(out, "%s:%d\n", dous[j]->name.c_str(), i + 1);
			if (out && printalot) fprintf(out, "%s\n", &outqueue[j].front());
			}

		if (out) fflush(out);

#ifdef F2MPI
		// TODO: Shouldn't this be turned on if F2MPI is enabled?
		if (false) reduce(world, markerweight, markerweight2, vectorplus<MWTYPE>(), 0);
#endif
		if (DOREMAPDISTANCES
#ifdef F2MPI
			&& world.rank() == 0
#endif
			&& !full)
		{
			markerweight = markerweight2;
			for (unsigned int q = chromstarts[i]; q < chromstarts[i + 1] - 2; q++)
			{
				double dist = markerposes[q + 1] - markerposes[q];
				if (dist > 0)
				{
					for (int t = 0; t < 2; t++)
					{
						double prob2 = -(markerweight[q][t * 2 + 1] - markerweight[q][t * 2]) / dous.size() / dist;
						if (prob2 > 0.3) prob2 = 0.3;
						if (prob2 < -0.3) prob2 = -0.3;

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
			for (unsigned int i = 0; i < INDCOUNT; i++)
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


					fprintf(out, "FIRST PASS: %d\n", i);
					fflush(out);
				}

#ifdef F2MPI
				if (world.rank()) continue;
#endif	      	      		 

				if (/*ind->pars[0] || ind->pars[1] || */!ind->haplocount.size()) continue;

#if !DOTOULBAR				
				// Perform the inversions indicated by the negshift data, at most a single one per individual
				// and chromosome, maybe 0.
				for (size_t c = 0; c < chromstarts.size() - 1; c++)
				{
					int minstart = chromstarts[c + 1];
					double minval = -1e-10;

					for (int p = chromstarts[c]; p < (int)chromstarts[c + 1]; p++)
					{
						if (ind->negshift[p] < minval)
						{
							minstart = p;
							minval = ind->negshift[p];
						}
						// TODO: Was a seed element ever needed here?
						/*						if (ind->negshift[p] < -1e-10)
						{
							if (!prevlow)
							{
								negshiftcand ourtuple(ind, 0, p);
								negshiftcands[c].insert(ourtuple);
							}
							prevlow = true;
							}*/
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
						if (!pred.anymatch && rand() / (RAND_MAX / 5))
						{
							negshiftcands[c].insert(ourtuple);
						}
						}

	}
#endif
			}

#ifdef F2MPI
			if (!world.rank())
#endif
			{
				std::atomic_int hitnnn(0);
				double oldscalefactor = scalefactor;
				bool any = false;
				for (size_t c = 0; c < chromstarts.size() - 1; c++)
				{		
				  if (negshiftcands[c].size()) any = true;
					for_each(negshiftcands[c].begin(), negshiftcands[c].end(), negshifter(c));
				}
								if (any)
								scalefactor = 0;

#pragma omp parallel for schedule(dynamic,1)
				for (unsigned int i = 0; i < INDCOUNT; i++)
				{

					individ* ind = getind(i, false);
					if (!ind || !ind->haplocount.size()) continue;

					fprintf(out, "SKEWNESS PASS: %d\n", i);
					fflush(out);

					unsigned int cno = 0;
					if (DOINFPROBS)
						for (size_t j = 0; j < ind->haplocount.size(); j++)
						{
							while (cno + 1 < chromstarts.size() && j >= chromstarts[cno + 1]) cno++;

							for (int side = 0; side < 2; side++)
							{
							  processinfprobs(ind, j, side, iter, hitnnn);
							}
							// oldinfprobslogic(ind, j, iter, cno, out);
						}

					updatehaploweights(ind, out, iter, hitnnn);
				}
#if !DOTOULBAR				
				parentswapnegshifts(nsm);
#endif				

				bool badhit = hitnnn > max(oldhitnnn, oldhitnnn2);
				if (badhit)
				{
					scalefactor /= 1.1;
				}
				bool goodhit = hitnnn < max(min(oldhitnnn, oldhitnnn2), (int) dous.size() / TURNBITS) * 0.99;
				if (goodhit)
				{
					scalefactor *= 1.21;
				}
				scalefactor *= 0.997;			
				//if (scalefactor < 0.01) scalefactor = 0.01;
								if (any) scalefactor = oldscalefactor;
								else
				{
					oldhitnnn2 = oldhitnnn;
					oldhitnnn = hitnnn;
					//if (goodhit) entropyfactor *= 0.99;
				}
				fprintf(stdout, "Scale factor now %lf, entropy %lf, hitnnn %d\n", scalefactor, entropyfactor, oldhitnnn);
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
	bool inact = false;
	individ* ind = 0;
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
individ* getind(string name, bool create = true)
{
	decltype(indmap.begin()) i;
	while ((i = indmap.find(name)) == indmap.end() && create)
	{
		int origindex = indmap.size();
		individ* ind = indmap[name] = getind(origindex, true);
		if (ind) ind->name = name;
	}

	return i == indmap.end() ? nullptr : i->second;
}

const individ* zeroguy = getind("0");

void readalphaped(FILE* in)
{
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
			individ* realpars[2] = { getind(me + (string)"_aux_realf"), getind(me + (string)"_aux_realm") };
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
	individ* haplo = getind("haplo");
	for (size_t x = 0; x < markerposes.size(); x++)
	{
		haplo->markerdata[x] = make_pair(sexmarkerval, sexmarkerval);
		haplo->markersure[x] = make_pair(0, 0);
	}

	char me[255];
	while (fscanf(in, "%s", me) == 1)
	{
		individ* ime = getind(me);
		bool doublehaplo = (ime->pars[1] == getind("haplo"));
		ime->empty = false;
		int data;
		for (size_t x = 0; x < markerposes.size(); x++)
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
				if (marker.first != UnknownMarkerVal)
				{
					ime->markersure[x] = make_pair(0.02, 0.02);
				}
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
									sureb1 = l1 / (double)(l1 + l2);

								if (data + data2 - l1 - l2)
									sureb2 = (data2 - l2) / (double)(data + data2 - l1 - l2);

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
					printf("%d/%d turned into %d %d with %lf;%lf\n", data, data2, marker.first.value(), marker.second.value(), markersure.first, markersure.second);
				}
				if (doublehaplo)
				{
					ime->markerdata[x] = make_pair(ime->markerdata[x].first, sexmarkerval);
				}
			}
		}
		ime->priormarkerdata = ime->markerdata;
		ime->priormarkersure = ime->markersure;
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
		ind->markersure.resize(markerposes.size());
		//		ind->semishift.resize(5000);
		ind->lastinved.resize(chromstarts.size());
		ind->lockstart.resize(chromstarts.size());

		for (size_t i = 0; i < chromstarts.size(); i++)
		{
			ind->lastinved[i] = -1;
			ind->lockstart[i] = 0;
		}

		int a, b;
		for (size_t k = 0; k < markerposes.size(); k++)
		{
			ind->haploweight[k] = 0.5;
			fscanf(pedfile, "%d %d", &a, &b);

			// Assume 1 percent allele error rate
			ind->markersure[k] = make_pair(a ? 1e-7 : 0.0, b ? 1e-7 : 0.0);

			ind->markerdata[k] = make_pair(a * MarkerValue, b * MarkerValue);

		}
	}
}

int family = 1;

void domerlinind(FILE* pedfile, individ* ind)
{
	int pn[2] = { 0 };

	for (int j = 0; j < 2; j++)
	{
		if (ind->pars[j]) pn[j] = ind->pars[j]->n;
	}

	fprintf(pedfile, "%d\t%d\t%d\t%d\t%d", family, ind->n, pn[0], pn[1], ind->sex + 1);
	for (unsigned int k = 0; k < chromstarts[1]; k++)
	{
		fprintf(pedfile, "\t%d\t%d", ind->markerdata[k].first.value(), ind->markerdata[k].second.value());
	}

	fprintf(pedfile, "\n");
}

template<class RuleType, class AttrType> void parseToEndWithError(mapped_file_source& file, const RuleType& rule, AttrType& target)
{
	bool res = phrase_parse(file.begin(), file.end(), rule, x3::space - x3::eol, target);

	if (!res)
	{
		std::string val;
		//file >> val;
		throw logic_error("Parsing failed. " + (std::string) __func__ + " " + val);
	}

	// TODO
	/*if (!file.eof())
	{
		throw logic_error("Not reaching end of file in parser. " + (std::string) __func__);
	}*/
}

template<class RuleType> void parseToEndWithError(mapped_file_source& file, const RuleType& rule)
{
	bool res = phrase_parse(file.begin(), file.end(), rule, x3::space - x3::eol);

	if (!res)
	{
		std::string val;
		/*while (!file.eof())
		  {
			file >> val;
			std::cout << val;
		  }
		std::cout << std::endl;*/
		throw logic_error("Parsing failed. " + (std::string) __func__ + " # " + val);

	}

	// TODO
	/*if (!file.eof())
	{
		throw logic_error("Not reaching end of file in parser. " + (std::string) __func__);
	}*/
}

#ifdef READHAPSSAMPLE
std::string filterExisting(const set<std::string>& names, std::string name)
{
	if (names.find(name) == names.end())
	{
		return "0";
	}

	return name;
}

vector<int> mapIndices;
vector<bool> hapmonomorphs;

auto word_ = x3::lexeme[+(x3::char_ - x3::space)];

typedef std::vector<std::tuple<std::string, std::string, std::string>> sampletype;

const auto marker_ = (x3::int_ > word_);
const auto hapsLine = (marker_ > x3::omit[x3::float_] > word_ > word_ > (+x3::int_));
const auto hapsLineIgnoreGenotypes = (marker_ > x3::long_long > word_ > word_ > x3::omit[(+x3::int_)]);

const auto sampleHeader = x3::omit[
	x3::repeat(7)[word_] > x3::eol >
		x3::int_ > x3::int_ > x3::int_ > x3::repeat(4)[word_] > x3::eol];
const auto sampleLine = (x3::omit[x3::int_] >> word_ >> x3::omit[x3::int_] >> word_ >> word_ >> x3::omit[x3::int_] >> x3::omit[x3::int_]);

x3::rule<class samplesRule, sampletype, true> const samplesRule = "samplesRule";
auto const samplesRule_def = sampleHeader >> (sampleLine% x3::eol);

BOOST_SPIRIT_DEFINE(samplesRule)

struct samplereader
{
	sampletype samples;

	void read(mapped_file_source& sampleFile)
	{
		using namespace x3;

		try
		{
			parseToEndWithError(sampleFile, samplesRule, samples);
			std::cout << samples.size() << " samples read." << std::endl;
		}
		catch (expectation_failure<mapped_file_source::iterator> const& x)
		{
			std::cerr << "expected: " << x.which();
			std::cerr << "got: \"" << x.where() << '"' << std::endl;

			throw x;
		}
	}
};

using SnpDataType = std::vector<std::tuple<int, std::string, std::string, std::string, std::vector<int>>>;

template<class T1>
void priorGenotypesFromHaps(const SnpDataType& snpData, const vector<individ*>&
	sampleInds, T1 indexconv, double unit)
{
	for (size_t j = 0; j < sampleInds.size(); j++)
	{
		sampleInds[j]->priorgenotypes.resize(snpData.size());
	}
	for (size_t i = 0; i < snpData.size(); i++)
	{		
		const vector<int>& markers = get<4>(snpData[i]);
		for (size_t j = 0; j < sampleInds.size(); j++)
		{
		  int geno = indexconv(markers[j * 2], i).value() + indexconv(markers[j * 2 + 1], i).value() - 2;
		  if (geno >= 0 && geno <= 2)
			sampleInds[j]->priorgenotypes[i][geno] += unit;
		}
	}
}

template<class T1, class T2>
void readFirstHaps(const SnpDataType& snpData, const vector<individ*>&
	sampleInds, T1 dohaploweight, T2 indexconv)
{
	for (size_t i = 0; i < snpData.size(); i++)
	{
		const vector<int>& markers = get<4>(snpData[i]);
		for (size_t j = 0; j < sampleInds.size(); j++)
		{
			float sureVal = 0;
			/*if (sampleInds[j]->gen == 2)*/ sureVal = 0;
			auto straightmarker = make_pair(indexconv(markers[j * 2], i), indexconv(markers[j * 2 + 1], i));
			for (int p = 1; p <= 2; p++)
			{
				auto marker = (p == 1) ? straightmarker :
					make_pair(indexconv(markers[j * 2 + 1], i), indexconv(markers[j * 2], i));
				if (marker == sampleInds[j]->markerdata[i])
				{
					sampleInds[j]->priormarkerdata[i] = straightmarker;
					break;
				}
			}
			sampleInds[j]->markerdata[i] = straightmarker;
			if (dohaploweight(sampleInds[j])) sampleInds[j]->haploweight[i] = 1e-3;
			sampleInds[j]->markersure[i] = { sureVal, sureVal };
			if (RELSKEWS)
			{
				if (i != markerposes.size() - 1)
				{
					//sampleInds[j]->relhaplo[i] = 0.51;//(markers[j * 2] == markers[j * 2 + 1]) ? 1.0 : 0.51;
					// Assume that our local bp/cM model is 1e-6 and the population-level rho used by ShapeIT is 0.0004
					// And that a good proxy for haplotype accuracy as a function of length is the global rho...
					// Most important point is to try to get negshift inversion boundaries located on actual marker map gaps.
					sampleInds[j]->relhaplo[i] = 0.5 + 0.5 * exp(-(markerposes[i + 1] - markerposes[i]) /* * 1e6 * 0.0004 */);
				}
			}
		}
	}
}

template<class T1, class T2>
void readOtherHaps(const SnpDataType& snpData,
		   const vector<individ*>& sampleInds, double unit, double genounit, T1 dohaploweight, T2 indexconv)
{
	vector<int> phases;
	vector<int> origPhases;
	phases.resize(sampleInds.size());
	origPhases.resize(sampleInds.size());
	auto findMatch = [&indexconv, &sampleInds, &snpData](int i, int j) -> pair<int, int> {
		const vector<int>& markers = get<4>(snpData[i]);
		int matchNum;
		int numMatches = 0;

		for (int p = 1; p <= 2; p++)
		{
			auto marker = (p == 1) ? make_pair(indexconv(markers[j * 2], i), indexconv(markers[j * 2 + 1], i)) :
				make_pair(indexconv(markers[j * 2 + 1], i), indexconv(markers[j * 2], i));
			if (marker == sampleInds[j]->markerdata[i])
			{
				matchNum = p;
				numMatches++;
			}
		}
		return make_pair(matchNum, numMatches);

	};

	for (size_t i = 0; i < snpData.size(); i++)
	{
		for (size_t j = 0; j < sampleInds.size(); j++)
		{
			if (!origPhases[j])
			{
				auto [matchNum, numMatches] = findMatch(i, j);
				if (numMatches == 1)
				{
					origPhases[j] = matchNum;
					phases[j] = origPhases[j];
				}
			}
		}
	}

	for (size_t j = 0; j < sampleInds.size(); j++)
	{
		if (!origPhases[j])
		{
			origPhases[j] = 1;
			phases[j] = 1;
		}
	}

	for (size_t i = 0; i < snpData.size(); i++)
	{
		const vector<int>& markers = get<4>(snpData[i]);


		for (size_t j = 0; j < sampleInds.size(); j++)
		{
			int oldPhase = phases[j];
			auto [matchNum, numMatches] = findMatch(i, j);

			// Not conclusive, assume same phase
			if (numMatches == 2 || numMatches == 0)
			{
				matchNum = oldPhase;
			}

			phases[j] = matchNum;
			if (dohaploweight(sampleInds[j]) && origPhases[j] != phases[j])
			{
				sampleInds[j]->haploweight[i] += unit;
			}
			if (RELSKEWS && i)
			{
				// The definition of relhaplo is that marker i
				// defines the skewness shift to the NEXT marker.
			  sampleInds[j]->relhaplo[i - 1] += unit * ((oldPhase == 0 || phases[j] == oldPhase) * 1.0 + (numMatches == 0) * (0 * -0.5));
			}
			if (!numMatches)
			{
				MarkerVal ms[2] = { indexconv(markers[j * 2], i), indexconv(markers[j * 2 + 1], i) };
				if (phases[j] == 2) swap(ms[0], ms[1]);

				bool nomatch[2] = { true, true };
				for (int p = 0; p < 2; p++)
				{
					for (int z = p; z < p + 1; z++)
					{
						if (ms[z] == (&sampleInds[j]->markerdata[i].first)[p])
						{
							nomatch[p] = false;
						}
					}
				}
				if (!nomatch[0] && !nomatch[1]) nomatch[0] = nomatch[1] = true;
				sampleInds[j]->markersure[i] = make_pair(min(sampleInds[j]->markersure[i].first + genounit * nomatch[0], (1 - unit)), min(sampleInds[j]->markersure[i].second + genounit * nomatch[1], (1 - unit)));
			}
		}
	}
}

template<class T1>
double initPadding(const vector<individ*>& sampleInds, int count, T1 dohaploweight)
{
	const double padding = 1e-2;
	double unit = 1.0 / (count + padding);
	for (size_t j = 0; j < sampleInds.size(); j++)
	{
		for (size_t i = 0; i < markerposes.size(); i++)
		{
			if (RELSKEWS)
			{
				sampleInds[j]->relhaplo[i] = unit;
			}
			if (dohaploweight(sampleInds[j])) sampleInds[j]->haploweight[i] = unit * padding * 0.5;
			sampleInds[j]->markersure[i] = make_pair(padding * unit, padding * unit);
		}
	}

	return unit;
}


MarkerVal mvFromIndex(int index, size_t)
{
	return (index + 1) * MarkerValue;
}

void readhapsfull(const sampletype& samples, mapped_file_source& bimFile, vector<mapped_file_source*>& hapsFile)
{
	using namespace x3;

	SnpDataType snpData;
	map<std::pair<int, std::string>, pair<int, int> > geneMap;

	auto alleles_ = word_ > word_;
	auto bimLine = (marker_ > omit[float_] > int_ > omit[alleles_]);

	try
	{
		parseToEndWithError(bimFile, (bimLine[([&](auto& context)
			{
				using namespace boost::fusion;
				auto attr = _attr(context);
				int index = geneMap.size();
				geneMap[make_pair(at_c<0>(attr), at_c<1>(attr))] = make_pair(at_c<2>(attr), index);
			})])
			% eol);
		std::cout << geneMap.size() << " entries read in map." << std::endl;

		parseToEndWithError(*hapsFile[0], hapsLine % eol, snpData);
		std::cout << snpData.size() << " SNPs read." << std::endl;
	}
	catch (expectation_failure<mapped_file_source::iterator> const& x)
	{
		std::cerr << "expected: " << x.which();
		std::cerr << "got: \"" << x.where() << '"' << std::endl;

		throw x;
	}

	int lastchrom = -1;
	double basepos = 0;

	// Walk over SNPs to create map info
	for (auto snp : snpData)
	{
		int chrom = get<0>(snp);
		int bppos;
		int index;

		std::tie(bppos, index) = geneMap[make_pair(chrom, get<1>(snp))];

		double pos = bppos * 1e-6;

		if (chrom != lastchrom)
		{
			chromstarts.push_back(markerposes.size());
			basepos = pos;
		}

		markerposes.push_back(pos - basepos);
		mapIndices.push_back(index);
		lastchrom = chrom;
		hapmonomorphs.push_back(get<2>(snp) == get<3>(snp));
	}

	chromstarts.push_back(markerposes.size());

	vector<individ*> sampleInds;
	set<std::string> sampleNames;

	// Collect all names for which we really have samples
	for (std::tuple<std::string, std::string, std::string> sample : samples)
	{
		sampleNames.insert(get<0>(sample));
	}

	for (std::tuple<std::string, std::string, std::string> sample : samples)
	{
		individ* me = getind(get<0>(sample));
		// We don't care about sex? Is this OK?
		me->sex = 0;
		me->pars[0] = getind(filterExisting(sampleNames, get<1>(sample)));
		me->pars[1] = getind(filterExisting(sampleNames, get<2>(sample)));

		printf("%x %x\n", me->pars[0], me->pars[1]);

		// Hack the generation to make non-founders full citizens
		me->gen = 2 * (me->pars[0] || me->pars[1]);
		dous.push_back(me);

		sampleInds.push_back(me);
	}

	auto dohaploweight = [](individ* ind) { return (ind->gen < 2); };

	readFirstHaps(snpData, sampleInds, dohaploweight, mvFromIndex);

	double unit = initPadding(sampleInds, hapsFile.size(), dohaploweight);

	for (size_t k = 1; k < hapsFile.size(); k++)
	{
		snpData.clear();
		parseToEndWithError(*hapsFile[k], hapsLine % eol, snpData);
		std::cout << snpData.size() << " SNPs read." << std::endl;

		readOtherHaps(snpData, sampleInds, unit, unit, dohaploweight, mvFromIndex);
	}

	for (size_t j = 0; j < sampleInds.size(); j++)
	{
		sampleInds[j]->priormarkerdata = sampleInds[j]->markerdata;
		sampleInds[j]->priormarkersure = sampleInds[j]->markersure;
	}
}

void readhapsonly(vector<mapped_file_source*>& hapsFile)
{
	using namespace x3;

	SnpDataType snpData;
	bool DH = true;
	auto dohaploweight = [&DH](individ* ind) { /*return (ind->gen < 2);*/ return DH; };

	auto mapToSnpGeno = [&snpData](int index, size_t snp)
	{
		switch (index)
		{
		case 0:
			return (std::get<2>(snpData[snp])[0] - 48) * MarkerValue;
		case 1:
			return (std::get<3>(snpData[snp])[0] - 48) * MarkerValue;
		default:
			fprintf(stderr, "Encountered index %d at snp %d\n", index, snp);
			abort();
		}
	};	
	
	/*for (auto file : hapsFile)
	{
		snpData.clear();
		parseToEndWithError(*file, hapsLine % eol, snpData);
		priorGenotypesFromHaps(snpData, dous, mapToSnpGeno, unit);
		}

	for (size_t j = 0; j < dous.size(); j++)
	{
		for (size_t i = 0; i < markerposes.size(); i++)
		{
			// TODO: Only one marker missing.
			MarkerValPair& marker = dous[j]->markerdata[i];
			switch (std::max_element(dous[j]->priorgenotypes[i].begin(), dous[j]->priorgenotypes[i].end()) - dous[j]->priorgenotypes[i].begin())
			{
			case 0:
				marker = { 1 * MarkerValue, 1 * MarkerValue };
				break;
			case 1:
			  // Preserve opposite heterzygote if present
			  if (marker != pair{2 * MarkerValue, 1 * MarkerValue})
				marker = { 1 * MarkerValue, 2 * MarkerValue };
				break;
			case 2:
				marker = { 2 * MarkerValue, 2 * MarkerValue };
				break;
			}
		}
		}*/

	parseToEndWithError(*hapsFile[0], hapsLine % eol, snpData);
	readFirstHaps(snpData, dous, dohaploweight, mapToSnpGeno);
	
	double unit = initPadding(dous, hapsFile.size(), dohaploweight);


	for (size_t k = 1; k < hapsFile.size(); k++)
	{
		snpData.clear();
		parseToEndWithError(*hapsFile[k], hapsLine % eol, snpData);
		std::cout << snpData.size() << " SNPs read." << std::endl;

		readOtherHaps(snpData, dous, unit, unit /* * 0.5 */, dohaploweight, mapToSnpGeno);
	}

	for (size_t j = 0; j < dous.size(); j++)
	{
		for (size_t i = 0; i < markerposes.size(); i++)
		{
			if (dous[j]->priormarkerdata.size() > i&& dous[j]->priormarkerdata[i] == pair(UnknownMarkerVal, UnknownMarkerVal))
			{
				dous[j]->priormarkerdata[i] = dous[j]->markerdata[i];
				dous[j]->priormarkersure[i] = dous[j]->markersure[i];
			}
		}
	}

	



	//return; // NOTE


}

void createhapfile(const sampletype& samples, mapped_file_source& oldhapfile, ostream& newhapfile)
{
	vector<individ*> sampleInds;
	for (auto sample : samples)
	{
		sampleInds.push_back(getind(get<0>(sample), false));
	}

	std::vector<std::tuple<int, std::string, long long, std::string, std::string>> snpData;

	using namespace x3;

	try
	{
		parseToEndWithError(oldhapfile, hapsLineIgnoreGenotypes % eol, snpData);
		std::cout << snpData.size() << " SNPs read." << std::endl;
	}
	catch (expectation_failure<mapped_file_source::iterator> const& x)
	{
		std::cerr << "expected: " << x.which();
		std::cerr << "got: \"" << x.where() << '"' << std::endl;

		throw x;
	}

	auto translateMarker = [](MarkerVal x) -> std::string
	{
		if (x == UnknownMarkerVal) return "?";
		if (x == 1 * MarkerValue) return "0";
		if (x == 2 * MarkerValue) return "1";
		assert(false);
	};

	int i = 0;
	for (auto marker : snpData)
	{
		newhapfile << get<0>(marker) << " " << get<1>(marker) << " " << get<2>(marker) << " " << get<3>(marker) << " " << get<4>(marker);
		for (auto ind : sampleInds)
		{
			pair<MarkerVal, MarkerVal> data = ind->markerdata[i];
			if (ind->haploweight[i] > 0.5 && (((ind->pars[0] && !ind->pars[0]->empty) || (ind->pars[1] && !ind->pars[1]->empty)) || ind->children))
			{
				swap(data.first, data.second);
			}

			newhapfile << " " << translateMarker(data.first) << " " << translateMarker(data.second);
		}
		newhapfile << "\n";
		i++;
	}
}

void readfambed(std::string famFileName, std::string bedFileName, bool readall = true, bool dooverwrite = false)
{
	using namespace boost::interprocess;
	using namespace x3;

	auto word = lexeme[+(char_ - space)];
	std::ifstream file(famFileName);
	file >> std::noskipws;
	auto parseriter = boost::spirit::istream_iterator(file);
	boost::spirit::istream_iterator end;
	map<std::string, int> indNums;

	phrase_parse(parseriter, end, (omit[int_] > word > omit[word] > omit[word] > omit[int_] > omit[int_])
		[([&](auto& ctx)
			{
				int index = indNums.size();
				indNums[_attr(ctx)] = index;
			})] % eol, space - eol);
	cout << indNums.size() << " individuals found." << std::endl;

	int blocksize = (indNums.size() + 3) / 4;
	cout << blocksize << " bytes per SNP block." << std::endl;

	file_mapping bedFile(bedFileName.c_str(), read_only);

	mapped_region snpRegion(bedFile, read_only, 3); // Skip header
	unsigned char* snpdata = (unsigned char*)snpRegion.get_address();
	size_t size = snpRegion.get_size();
	cout << size << " bytes, " << size / blocksize << " SNPs." << std::endl;

	vector<int> indArray;
	for (individ* ind : dous)
	{
		indArray.push_back(indNums[ind->name]);
		cout << ind->name << " " << indArray[indArray.size() - 1] << std::endl;
	}

	for (size_t i = 0; i < markerposes.size(); i++)
	{
		unsigned char* thisSnp = &snpdata[mapIndices[i] * blocksize];
		for (size_t j = 0; j < dous.size(); j++)
		{
			/*if (!(
			  dous[j]->pars[0] && !dous[j]->pars[0]->empty &&
			  dous[j]->pars[1] && !dous[j]->pars[1]->empty)) continue;*/

			int index = indArray[j];
			int thisval = (thisSnp[index / 4] >> (2 * (index % 4))) & 3;
			pair<MarkerVal, MarkerVal> marker;

			bool isachange = false;
			/*			if (rand() / (RAND_MAX / 10) == 0)
			  {
				cout << "Masking out " << i << "at " << j << "\n";
				thisval = 1;
				}*/
			switch (thisval)
			{
			case 0:
				marker = make_pair(1 * MarkerValue, 1 * MarkerValue);
				isachange = marker != dous[j]->priormarkerdata[i];
				break;
			case 1:
				marker = make_pair(UnknownMarkerVal, UnknownMarkerVal);
				cout << "MISSING BED " << dous[j]->name << " " << i << std::endl;
				break;
			case 2:
				marker = make_pair(1 * MarkerValue, 2 * MarkerValue);
				isachange = dous[j]->priormarkerdata[i].first == dous[j]->priormarkerdata[i].second;
				break;
			case 3:
				// ShapeIT will turn A A to 0 A, making all genotypes homozygotes for the second allele, rather than the first
				int val = 2 - hapmonomorphs[i];
				marker = make_pair(val * MarkerValue, val * MarkerValue);
				isachange = marker != dous[j]->priormarkerdata[i];
				break;
			}
			if (isachange)
			{
				cout << "!!! " << dous[j]->name << " " << i << " " << marker.first.value() << marker.second.value() << " " << std::endl;
			}

			/*if (readall || marker.first == UnknownMarkerVal || isachange)
			  {
				//cout << "/// " << dous[j]->name << " " << i << std::endl;
				dous[j]->priormarkerdata[i] = marker;
				if (marker.first == UnknownMarkerVal)
				  {
				/*if (RELSKEWS)
				  {
					dous[j]->relhaplo[i] = 1;
					}*/
					/*dous[j]->priormarkersure[i] = make_pair(0.f, 0.f);
					  }
				  }*/
			if (isachange || marker.first == UnknownMarkerVal)
			{
				dous[j]->priormarkersure[i] = make_pair(
					0.5 * (0.5 + dous[j]->priormarkersure[i].first),
					0.5 * (0.5 + dous[j]->priormarkersure[i].second));
				cout << "Increasing prior uncertainty individual " << dous[j]->n << ", marker " << i << std::endl;
				if (dooverwrite)
				{
					dous[j]->markerdata[i] = marker;
					// TODO MARKERSURE
				}
			}
			/*			if (dous[j]->priormarkerdata[i].first != UnknownMarkerVal)
			  {
				dous[j]->priormarkersure[i] = make_pair(
									max(1e-3, dous[j]->priormarkersure[i].first),
									max(1e-3, dous[j]->priormarkersure[i].second));
									}*/
		}
	}
}
#endif

auto mapline = x3::omit[x3::int_] > ::word_ > x3::double_ > x3::omit[x3::int_];

void readgigidata(mapped_file_source&& map, mapped_file_source&& ped)
{
	using namespace x3;
	chromstarts.push_back(0);
	parseToEndWithError(map, mapline
		[([&](auto& ctx)
			{
				using namespace boost::fusion;

				auto attr = _attr(ctx);
				auto [name, pos] = make_pair(at_c<0>(attr), at_c<1>(attr));
				markernames[name] = markerposes.size();
				markerposes.push_back(pos);
			})] % eol);
	chromstarts.push_back(markerposes.size());

	individ* ind;
	int nowpar = 0;
	int nowmarker = 0;

	auto indid = [&ind, &nowpar, &nowmarker](auto& ctx)
	{
		ind = getind(_attr(ctx));
		nowpar = 0;
		nowmarker = 0;
		dous.push_back(ind);
	};

	auto parid = [&ind, &nowpar](auto& ctx) -> int
	{
		ind->pars[nowpar++] = getind(_attr(ctx));
		return nowpar;
	};

	auto setSex = [&ind](auto& ctx)
	{
		ind->sex = _attr(ctx) - 1;
	};

	auto setMarker = [&ind, &nowmarker](auto& ctx) -> int
	{
		using namespace boost::fusion;

		auto attr = _attr(ctx);
		std::pair<int, int> data = make_pair(at_c<0>(attr), at_c<1>(attr));
		ind->markerdata[nowmarker] = make_pair(data.first * MarkerValue, data.second * MarkerValue);
		ind->markersure[nowmarker] = make_pair(0, 0);
		nowmarker++;

		return nowmarker;
	};

	auto transferPrior = [&ind](auto& ctx) -> int
	{
		ind->priormarkerdata = ind->markerdata;
		ind->priormarkersure = ind->markersure;

		return 0;
	};

	parseToEndWithError(ped,
		(omit[word_ > (word_[indid]) > repeat(2)[word_[parid]] > int_[setSex] > omit[word_] >
			repeat(markerposes.size())[(int_ > int_)[setMarker]]])[transferPrior] % eol);
}


void addprotmarkers(set<double>& protmarkers, mapped_file_source&& source)
{
	using namespace x3;
	phrase_parse(source.begin(), source.end(), 
		     omit[x3::lit("map") > x3::lit("marker") > x3::lit("positions")] >
			  
			  +(x3::double_[([&](auto& ctx)
			{
				using namespace boost::fusion;

				double pos = _attr(ctx);

				protmarkers.insert(pos);
			})]),
		     x3::space - x3::eol);
}

void addprotinds(set<individ*>& protinds, mapped_file_source&& source)
{
	using namespace x3;
	parseToEndWithError(source, (omit[word_] > word_
		[([&](auto& ctx)
			{
				std::string name = _attr(ctx);
				individ* ind = getind(name, false);
				if (ind == 0) // TODO: Correct null ind?
				{
					fprintf(stderr, "Incorrect individual name. %s\n", name.c_str());
				}
				protinds.insert(ind);
			})]) % eol);
}

void clearunprotected(set<individ*>& protinds, set<double>& protmarkers)
{
	for (individ* ind : dous)
	{
		if (protinds.find(ind) != protinds.end()) continue;

		int lastmarker = 0;
		for (double i : protmarkers)
		{
			for (; markerposes[lastmarker] < i; lastmarker++)
			{
				ind->markerdata[lastmarker] = { UnknownMarkerVal, UnknownMarkerVal };
				ind->priormarkerdata[lastmarker] = { UnknownMarkerVal, UnknownMarkerVal };
				ind->priormarkersure[lastmarker] = make_pair(0, 0);
				ind->markersure[lastmarker] = make_pair(0, 0);
			}
			lastmarker++;
		}
	}
}

void compareimputedoutput(istream& filteredOutput)
{
	int j = 0;
	while (!filteredOutput.eof())
	{
		std::cout << j++ << std::endl;
		for (individ* ind : dous)
		{
			std::string name;
			do
				filteredOutput >> name;
			while (name == "--");

			for (unsigned int i = 0; i < chromstarts[1]; i++)
			{
				double val[3];
				int maxval = 0;
				for (int k = 0; k < 3; k++)
				{
					std::string temp;
					filteredOutput >> temp;
					if (sscanf(temp.c_str(), "%lf", &val[k]) != 1)
					{
						val[k] = -1;
					}
					if (val[k] > val[maxval]) maxval = k;
				}

				int data = (ind->markerdata[i].first == 2 * MarkerValue) + (ind->markerdata[i].second == 2 * MarkerValue);
				int origmaxval = maxval;
				int oridata = data;

				// If reference allele is not aligned in shapit output:
				/*				  maxval = abs(1-maxval);
								  data = abs(1-data);*/
								  // That will only look at hetero vs. homo, not what hetero
				if (maxval != data &&
					ind->pars[0] && !ind->pars[0]->empty &&
					ind->pars[1] && !ind->pars[1]->empty &&
					i != chromstarts[1] - 1 && val[origmaxval] >= 0 && ind->markerdata[i].first != UnknownMarkerVal)
				{
					std::cout << ind->name << " " << j << ":" << i << " " << oridata << "\t";
					for (double oneVal : val)
					{
						std::cout << oneVal << "\t";
					}
					std::cout << std::endl;
				}
			}
		}
	}
}

void readhaplodata(FILE* in, int swap)
{
	char tlf[255];
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
	  0, 0 };

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

void deserialize(istream& stream)
{
	/*
	This is what we're unserializing
	fprintf(stdout, "%f\t%d\t%d\t\t%f\t%lf %lf %lf\n", ind->haploweight[j], ind->markerdata[j].first.value(), ind->markerdata[j].second.value(), ind->negshift[j],
												ind->markersure[j].first, ind->markersure[j].second, RELSKEWS ? ind->relhaplo[j] : 0);*/
	auto haploline = x3::double_ >> x3::int_ >> x3::int_ >> x3::double_ >> x3::double_ >> x3::double_;

	while (!stream.eof())
	{
		string line;
		std::getline(stream, line);

		pair<int, string> data;

		if (x3::phrase_parse(line.begin(), line.end(), x3::int_ >> word_ >> x3::eoi, x3::space, data))
		{
			int n;
			string name;
			std::tie(n, name) = data;

			individ* ind = getind(n, false);
			individ* indcheck = getind(name, false);

			if (ind && ind == indcheck)
			{
				int oldphase = 0;
				int switches = 0;
				for (size_t i = 0; i < markerposes.size(); i++)
				{
					std::getline(stream, line);
					std::tuple<double, int, int, double, double, double> output;
					if (!x3::phrase_parse(line.begin(), line.end(), haploline, x3::space, output))
					{
						std::cerr << "Reading haplotype for marker " << i << " for individual " << ind->name << " failed: " << line << std::endl;
					}
					else
					{
						ind->haploweight[i] = std::get<0>(output);

						pair<MarkerVal, MarkerVal> pmv = make_pair(std::get<1>(output) * MarkerValue, std::get<2>(output) * MarkerValue);
						pair<MarkerVal, MarkerVal> rmv = make_pair(std::get<2>(output) * MarkerValue, std::get<1>(output) * MarkerValue);
						bool inv = false;
						bool match = true;
						if (pmv != ind->markerdata[i])
						{
						  if (rmv != ind->markerdata[i])
						    {
							std::cerr << "Genotype mismatch for marker " << i << " for individual " << ind->name << " (" << ind->markerdata[i].first.value() << "," << ind->markerdata[i].second.value() << ") to " <<
								" (" << pmv.first.value() << "," << pmv.second.value() << ")" << std::endl;
							match = false;
						    }
						  else inv = true;
						}
						ind->markerdata[i] = pmv;
						ind->markersure[i] = make_pair(std::get<4>(output), std::get<5>(output));
						if (ind->haploweight[i] == 0.5) continue;
						if (pmv == rmv) continue;
						if (!match) continue;

						int newphase = 1 + ((ind->haploweight[i] > 0.5) ^ inv);
						if (oldphase && oldphase != newphase) switches++;

						oldphase = newphase;
						}
					}

				if (ind->children || (ind->pars[0] && !ind->pars[0]->empty) || (ind->pars[1] && !ind->pars[1]->empty)) std::cout << "Switches " << ind->n << " " << ind->name << "\t" << switches << std::endl;
				}
			else
			{
				std::cerr << "Supposed individual header not a header: " << line << std::endl;
			}
			}
		}
	}

std::string getname(individ* ind)
{
	if (ind)
	{
		return ind->name;
	}
	else
	{
		return "0";
	}
}

void outputped(const std::string& filename)
{
	FILE* f = fopen(filename.c_str(), "w");
	for (individ* ind : dous)
	{

		fprintf(f, "%d %s %s %s %d %d", 1, getname(ind).c_str(), getname(ind->pars[0]).c_str(), getname(ind->pars[1]).c_str(), ind->sex + 1, -9);
		for (int j = 0; auto [a, b] : ind->markerdata)
		{
			if (ind->haploweight[j++] > 0.5) swap(a, b);
			fprintf(f, " %d %d", a.value(), b.value());
		}
		fprintf(f, "\n");
	}
	fclose(f);
}

#if LIBSTATGEN
void outputvcf(const std::string& templatefilename, const std::string& outputfilename)
{
	VcfFileReader reader;
	VcfHeader header_read;

	reader.open(templatefilename.c_str(), header_read);	
	
	VcfFileWriter writer;
	writer.open(outputfilename.c_str(), header_read, InputFile::BGZF);

	VcfRecord record_template;
	while (reader.readRecord(record_template))
	{
		const char* markername = record_template.getIDStr();		
		int refnum = boost::lexical_cast<int>(record_template.getRefStr());
		auto pos = markernames.at(markername);

		for (int i = 0; i < header_read.getNumSamples(); i++)
		{
			const char* sampleName = header_read.getSampleName(i);
			individ* ind = getind(&sampleName[0], false);
			if (!ind) ind = getind(&sampleName[2], false); // Removing extra 1_ prefix...
			if (!ind)
			{
				fprintf(stderr, "ERROR getting %s\n", sampleName);
				abort();
			}

			auto [a, b] = ind->markerdata[pos];
			auto markerVal2Str = [refnum] (MarkerVal m) -> string
			{
				if (m == UnknownMarkerVal) return ".";
				return std::to_string(m.value() != refnum);
			};

			if (ind->haploweight[pos] > 0.5)
			{
				swap(a, b);
			}

			// TODO: Change to std::format when supported
			record_template.getGenotypeInfo().setString("GT", i, markerVal2Str(a) + "|" + markerVal2Str(b));
		}
		writer.writeRecord(record_template);
	}
	reader.close();
	writer.close();
}
#endif

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
	discstep = 1;
	sexc = 2;
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

#ifdef READHAPSSAMPLE
	po::options_description desc;
	po::variables_map inOptions;
	string impoutput, famfilename, bedfilename, deserializefilename, outputfilename, outputhapfilename, genfilename, pedfilename, mapfilename, samplefilename, protmarkersfn, protindsfn, gigimapfilename, gigipedfilename, outputpedfilename, 
		outputvcffilename, templatevcffilename;
	bool clear;
	int COUNT;
	int limit;

	desc.add_options()("samplefile", po::value<string>(&samplefilename), "ShapeIT-style .sample file")
		("bimfile", po::value<string>(), "BIM file")
		("hapfiles", po::value<vector<string> >()->multitoken(), "One or more HAP files, maximum realization followed by others.")
		("deserialize", po::value<string>(&deserializefilename), "Load existing Chaplink output as starting point, with reporting on number of inversions.")
		("impoutput", po::value<string>(&impoutput), "Imputed genotype output from previous run.")
		("famfile", po::value<string>(&famfilename), "Original PLINK fam file. Use with bedfile.")
		("bedfile", po::value<string>(&bedfilename), "Original PLINK bed file. Use with famfile.")
		("count", po::value<int>(&COUNT)->default_value(3), "Number of iterations")
		("limit", po::value<int>(&limit)->default_value(INDCOUNT), "Maximum number of individuals")
		("output", po::value<string>(&outputfilename), "Output file name")
		("tmppath", po::value<string>(&tmppath), "Directory for toulbar temp files")
		("capmarker", po::value<int>()->notifier([&](int cap)
			{
				markerposes.resize(cap);
				chromstarts[1] = min(cap, (int)chromstarts[1]);
			}), "Limit to marker count.")
			("mapfile", po::value<string>(&mapfilename), "map file in original PlantImpute format, similar to AlphaImpute.")
				("pedfile", po::value<string>(&pedfilename), "ped file in original PlantImpute format, similar to AlphaImpute.")
				("genfile", po::value<string>(&genfilename), "Genotype file in original PlantImpute format, similar to AlphaImpute.")
				("gigimapfile", po::value<string>(&gigimapfilename), "map file in format compatible with Gigi.")
				("gigipedfile", po::value<string>(&gigipedfilename), "ped file in format compatible with Gigi (including genotypes).")
				("outputpedfile", po::value<string>(&outputpedfilename), "Output ped file based on read data (including genotypes).")
#if LIBSTATGGEN				
				("templatevcffile", po::value<string>(&templatevcffilename), "Template vcf file to use.")
				("outputvcffile", po::value<string>(&outputvcffilename), "Output vcf file based on read data (including genotypes, template required).")
#endif				
				("protmarkers", po::value<string>(&protmarkersfn), "File of mapping distances for protected markers. Used with --clear.")
				("protinds", po::value<string>(&protindsfn), "File of individuals to not clear. Used with --clear.")
				("clear", po::bool_switch(&clear), "Clear all non-protected markers in all non-protected individuals.")
				("createhapfile", po::value<string>(&outputhapfilename), "Output a hapfile based on input haplotypes.");

			auto parser = po::command_line_parser(argc, argv);
			parser.options(desc);
			po::store(parser.run(), inOptions);
			po::notify(inOptions);

			if (mapfilename != "")
			{
				FILE* mapfile = fopen(mapfilename.c_str(), "rt");
				cerr << "Reading map file " << mapfilename << "\n";
				readalphamap(mapfile);
			}

			if (pedfilename != "")
			{
				FILE* pedfile = fopen(pedfilename.c_str(), "rt");
				cerr << "Reading pedigree file " << pedfilename << "\n";
				readalphaped(pedfile);
			}

			if (genfilename != "")
			{
				FILE* genofile = fopen(genfilename.c_str(), "rt");
				cerr << "Reading genotype file " << genfilename << "\n";
				readalphadata(genofile);
			}

			if (gigimapfilename != "" && gigipedfilename != "")
			{
				readgigidata(mapped_file_source(gigimapfilename), mapped_file_source(gigipedfilename));
			}


			auto clearer = [&] ()
			{
				set<double> protmarkers;
				set<individ*> protinds;
				if (!protmarkersfn.empty())
				{
					addprotmarkers(protmarkers, mapped_file_source(protmarkersfn));
				}
				if (!protindsfn.empty())
				{
					addprotinds(protinds, mapped_file_source(protindsfn));
				}
				clearunprotected(protinds, protmarkers);
			};

			// TODO: Make sets of required params.
			if (clear && deserializefilename == "") clearer();
			

			samplereader samples;
			vector<mapped_file_source*> hapFiles;
			vector<string> hapsfileOption;
			if (inOptions.count("hapfiles")) hapsfileOption = inOptions["hapfiles"].as<vector<string>>();

			for (string filename : hapsfileOption)
			{
				hapFiles.push_back(new mapped_file_source(filename));
			}

			if (samplefilename != "")
			{
				mapped_file_source sampleFile(samplefilename);
				samples.read(sampleFile);

				mapped_file_source bimFile(inOptions["bimfile"].as<string>());


				readhapsfull(samples.samples, bimFile, hapFiles);
				std::cout << "readhapssample finished." << std::endl;
			}
			else
			{
				if (hapFiles.size()) readhapsonly(hapFiles);
}

			bool docompare = (impoutput != "");
			if (inOptions.count("famfile") + inOptions.count("bedfile") == 2)
			{
				std::cout << "readfambed started." << std::endl;
				readfambed(famfilename, bedfilename, docompare);
			}

			// Put generation 2 first, since those are more complex to analyze, avoiding a few threads
			// getting stuck towards the end.
			//	stable_sort(dous.begin(), dous.end(), [] (individ* a, individ* b) { return a->gen > b->gen; } );

			if (docompare)
			{
				std::ifstream filteredOutput(impoutput);
				compareimputedoutput(filteredOutput);

				return 0;
			}
#endif

			//	return 0;
			CORRECTIONINFERENCE = true;
			postmarkerdata(limit);
			CORRECTIONINFERENCE = false;

			if (deserializefilename != "")
			{
				std::ifstream deserializationFile(deserializefilename);

				std::cout << "deserialize started." << std::endl;
				deserialize(deserializationFile);
				std::cout << "deserialize finished." << std::endl;

				if (clear) clearer();
			}

			if (outputpedfilename != "")
			{
				outputped(outputpedfilename);
			}

#if LIBSTATGEN
			if (outputvcffilename != "")
			{
				outputvcf(templatevcffilename, outputvcffilename);
			}
#endif			


			if (samplefilename != "" && outputhapfilename != "")
			{
				std::ofstream outputhapfile(outputhapfilename);
				createhapfile(samples.samples, *hapFiles[0], outputhapfile);

				return 0;
}

			FILE* out = stdout;
			if (outputfilename != "")
			{
				out = fopen(outputfilename.c_str(), "w");
			}
			if (limit < dous.size()) dous.resize(limit);
			//chromstarts[1] = 100;

			if (HAPLOTYPING || true)
				for (int i = 0; i < COUNT; i++)
				{
					//		  	  	{
					early = (i < 1);
					if (!early) doit<false, genotypereporter>((i == COUNT - 1) ? out : stdout, /*i == COUNT - 1*/ true
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

					for (unsigned int i2 = 0; i2 < INDCOUNT; i2++)
					{
					  	if (i2 > limit) continue;
						individ* ind = getind(i2);
						if (!ind) continue;

						if (ind->haplocount.size())
						{
#ifdef F2MPI
							if (!world.rank())
#endif
								/*if (i == COUNT - 1)*/						fprintf(out, "%d %s\n", i2, ind->name.c_str());
							// Printing of haplotype data for each iteration
							for (unsigned int c = 0; c < chromstarts.size() - 1; c++)

							{
								for (unsigned int j = chromstarts[c]; j < chromstarts[c + 1]; j++)
								{

#ifdef F2MPI
									if (!world.rank())
#endif
										/*if (i == COUNT - 1)*/ 
										if (ind->priormarkerdata.size() > j) fprintf(out, "%f\t%d\t%d\t\t%f\t%lf %lf %lf\t%d\t%d\t%lf\t%lf\n", ind->haploweight[j], ind->markerdata[j].first.value(), ind->markerdata[j].second.value(), ind->negshift[j],
											ind->markersure[j].first, ind->markersure[j].second, RELSKEWS ? ind->relhaplo[j] : 0,
											ind->priormarkerdata[j].first.value(), ind->priormarkerdata[j].second.value(),
											ind->priormarkersure[j].first, ind->priormarkersure[j].second);
											else
											fprintf(out, "%f\t%d\t%d\t\t%f\t%lf %lf %lf\n", ind->haploweight[j], ind->markerdata[j].first.value(), ind->markerdata[j].second.value(), ind->negshift[j],
											ind->markersure[j].first, ind->markersure[j].second, RELSKEWS ? ind->relhaplo[j] : 0);
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

