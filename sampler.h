#ifndef _SAMPLER_FROM_GSL_
#define _SAMPLER_FROM_GSL_

//C++
#include <cmath>
#include <vector>
#include <iostream>

//GSL
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_sort.h>

namespace Sampler
{

template<typename Scalar = double>
class Distribution
{
  public:
    Distribution(const gsl_rng_type * rT = gsl_rng_knuthran2002);
    virtual ~Distribution();

    virtual void set_parameters(const std::vector<Scalar> &parameters)              {return;}
    virtual void set_parameters(const std::vector<std::vector<Scalar> > &parameters){return;}
    virtual void set_parameters(const unsigned int npars)                           {return;}

    virtual void sample(unsigned int nsample, std::vector<Scalar> &sample)              const {return;}
    virtual void sample(unsigned int nsample, std::vector<std::vector<Scalar> >&sample) const {return;}

    virtual void sample_indirect(unsigned int nsample, std::vector<Scalar> &U01, std::vector<Scalar> &sample) const {return;}
    virtual void sample_indirect(unsigned int nsample, std::vector<std::vector<Scalar> > &U01, 
                                                       std::vector<std::vector<Scalar> > &sample) const {return;}

    virtual unsigned int npars() const  = 0;

    Distribution<Scalar> * mom() const {return _mother;}

    void set_mother(Distribution<Scalar> *mommy, unsigned int ipar) {mommy->set_daughter(this,ipar);}
    void set_mother(Distribution<Scalar> *mommy) {_mother = mommy;}
    void set_daughter(Distribution<Scalar> * daught, unsigned int ipar);

    void set_daughters(std::vector<Distribution<Scalar> *> daught);

    void sample_tree(unsigned int nsample,std::vector<std::vector<Scalar> > &sample);

    std::vector<Distribution<Scalar>*> daughters()     const {return _daughters;}
    Distribution<Scalar> * daughter(unsigned int idau) const;

    const gsl_rng * random_gsl() const {return _r;}

    virtual void print(std::ostream &out = std::cout) const {out << "nothing to do here" << std::endl;}

  protected:
   const gsl_rng_type * _T;
   gsl_rng * _r;
   Distribution<Scalar> * _mother;
   std::vector<Distribution<Scalar> * > _daughters;

   void sample_branch(unsigned int nsample, std::vector<Scalar> &sample, Distribution<Scalar> *mother,
                      std::vector<std::vector<Scalar> > &sampleFill, unsigned int idaughter);

};

template<typename Scalar>
Distribution<Scalar>::Distribution(const gsl_rng_type * rT):
_T(rT),
_mother(NULL)
{
  gsl_rng_env_setup();
  _r = gsl_rng_alloc(_T);
}

template<typename Scalar>
Distribution<Scalar>::~Distribution()
{
  gsl_rng_free(_r);
}

template<typename Scalar>
Distribution<Scalar> * Distribution<Scalar>::daughter(unsigned int idau) const
{
  if(_daughters.empty())return NULL;
  if(idau >= _daughters.size())return NULL;
  return _daughters[idau];
}

template<typename Scalar>
void Distribution<Scalar>::set_daughter(Distribution<Scalar> * daught, unsigned int ipar)
{
  if(_daughters.empty())_daughters.resize(npars(),NULL);
  _daughters[ipar] = daught;
  _daughters[ipar]->set_mother(this);
}

template<typename Scalar>
void Distribution<Scalar>::set_daughters(std::vector<Distribution<Scalar> *> daught)
{
  _daughters = daught;
}

template<typename Scalar>
void Distribution<Scalar>::sample_tree(unsigned int nsample, std::vector<std::vector<Scalar> > &sample)
{
  sample.clear();
  Distribution<Scalar> *cur = this;
  while(cur->mom() != NULL)cur = cur->mom();
  std::vector<std::vector<Scalar> > sampletmp;
  cur->sample(nsample,sampletmp);//init
  for(unsigned int ibranch = 0; ibranch < cur->npars(); ibranch++)
  {
     std::vector<std::vector<Scalar> > sampleDaughter;
     this->sample_branch(nsample,sampletmp[ibranch],cur,sampleDaughter,ibranch);
     for(unsigned int id = 0; id < sampleDaughter.size(); id++)
     {
       sample.push_back(sampleDaughter[id]);
     }
  }
}

template<typename Scalar>
void Distribution<Scalar>::sample_branch(unsigned int nsample, std::vector<Scalar> &sample, 
                                        Distribution<Scalar> *mother,std::vector<std::vector<Scalar> > &sampleFill, unsigned int idaughter)
{
  Distribution<Scalar> *daughter = mother->daughter(idaughter);
  if(daughter == NULL)
  {
    sampleFill.push_back(sample);
    return;
  }
  std::vector<std::vector<Scalar> > sampletmp;
  daughter->sample(nsample,sampletmp);//init
  for(unsigned int ibranch = 0; ibranch < daughter->npars(); ibranch++)
  {
     std::vector<std::vector<Scalar> > sampleDaughter;
     this->sample_branch(nsample,sampletmp[ibranch],daughter,sampleDaughter,ibranch);
     for(unsigned int ibranch = 0; ibranch < sampleDaughter.size(); ibranch++)
     {
       for(unsigned int is = 0; is < nsample; is++)
       {
         sampleDaughter[ibranch][is] *= sample[is];
       }
       sampleFill.push_back(sampleDaughter[ibranch]);
     }
  }

  return;
}

}
#endif
