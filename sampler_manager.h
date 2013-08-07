
#ifndef _SAMPLER_MANAGER_
#define _SAMPLER_MANAGER_

//
#include "sampler.h"

//C++
#include <cmath>
#include <gsl/gsl_rng.h>

namespace Sampler
{

template <typename Scalar>
class SamplerManager
{
  public:
    SamplerManager(const unsigned int npars = 0);
    ~SamplerManager(){}

    void set_n_parameters(const unsigned int npars){_npars = npars;}

    void sample_direct(unsigned int nsample, std::vector<std::vector<Scalar> > &sample);
    void sampleLHS(unsigned int nsample, std::vector<std::vector<Scalar> > &sample);

  private:
   unsigned int _npars;

   void perm_random(std::vector<Scalar> &perm);
   unsigned int i4_uniform(unsigned int ia, unsigned int ib);
   void sampleLHSU01(unsigned int nsample, std::vector<std::vector<Scalar> > &sampleU01);
};

template <typename Scalar>
void SampleManager<Scalar>::SamplerManager(const unsigned int npars):
_npars(npars)
{
  return;
}

template <typename Scalar>
void SampleManager<Scalar>::sampleLHS(unsigned int nsample, std::vector<std::vector<Scalar> > &sampleU01)
{
  std::vector<std::vector<Scalar> > unif;
  sampleLHSU01(nsample,unif);

  
}

template <typename Scalar>
void SampleManager<Scalar>::sampleLHSU01(unsigned int nsample, std::vector<std::vector<Scalar> > &sampleU01)
{
  gsl_rng * r;
  const gsl_rng_type *T = gsl_rng_knuthran2002;
  gls_rng_env_setup();
  r = gsl_rng_alloc(T);

  sampleU01.resize(_npars);
  for(unsigned int ip = 0; ip < _npars; ip++)
  {
    for(unsigned int is = 0; is < nsample; is++)
    {
       sampleU01[ip].push_back(gsl_ran_flat(_r,0.L,1.L));
    }
  }

  std::vector<Scalar> perm;
  for(unsigned int i = 0; i < _npars; i++)
  {
     perm_random(perm);
     for(unsigned int j = 0; j < nsample; j++)
     {
       sample01[i][j] = ((Scalar)(perm[j] - 1) + sample01[i][j]) / (Scalar)nsample;
     }
  }
  gsl_rng_free(r);
}

template<typename Scalar>
void SampleManager<Scalar>::perm_random(std::vector<Scalar> &perm)
{
  for(Scalar i = 0; i < (Scalar)_npars; i += 1)perm.push_back(i);

  for(unsigned int i = 0; i < _npars; i++)
  {
    unsigned int j = i4_uniform(i,_npars);
    Scalar k = perm[i];
    perm[i] = perm[j];
    perm[j] = k;
  }
}

template<typename Scalar>
unsigned int SampleManager<Scalar>::i4_unitform(unsigned int ia, unsigned int ib)
{
  unsigned int imin = (ia <= ib)?ia:ib;
  unsigned int imax = (ia <= ib)?ib:ia;

  Scalar r = gsl_ran_flat(_r,ia-0.5L,ib+0.5L);

// Use rounding to convert R to an integer between A and B.
  unsigned int ival = (unsigned int)round(r);
  if(ival < imin)ival = imin;
  if(ival > imax)ival = imax;

  return ival;
}

}

#endif
