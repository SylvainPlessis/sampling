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


//!gaussian distribution
template <typename T = double>
class Norm:public Distribution<T>
{
  public:
    Norm():_npars(2){}
    Norm(const std::vector<std::vector<T> >&pars):_npars(2) {this->set_parameters(pars);}
    Norm(const std::vector<T> & v, const std::vector<T> & s):_mean(v),_sig(s),_npars(2){}
    Norm(const T & v, const T & s):_npars(2){set_parameters(v, s);}
    Norm(const Norm<T> &rhs):_npars(2){if(this == &rhs)return;_mean = rhs.mean();_sig=rhs.sigma();this->_r=rhs.random_gsl();}
    ~Norm(){}

   void set_parameters(const std::vector<std::vector<T> > & parameters){_mean = parameters[0]; _sig = parameters[1];}
   void set_parameters(const T & v, const T & s){_mean.push_back(v);_sig.push_back(s);}

   void sample(unsigned int nsample, std::vector<std::vector<T> >&sample) const;
   void sample_indirect(unsigned int nsample, const std::vector<std::vector<T> >&U01, std::vector<std::vector<T> >&sample) const;

   void sample(unsigned int nsample, std::vector<T> &sample) const {this->sample_intern(nsample,_mean[0],_sig[0],sample);}
   void sample_indirect(unsigned int nsample, const std::vector<T> &U01, std::vector<T> &sample) const 
                                                {this->sample_indirect_intern(nsample,_mean[0],_sig[0],U01,sample);}

   unsigned int npars()          const {return _mean.size();}
   const std::vector<T> mean()   const {return _mean;}
   const std::vector<T> sigma()  const {return _sig;}

   void print(std::ostream &out = std::cout) const;
   friend std::ostream & operator << (std::ostream &out, const Norm<T> &normal)
    {
        normal.print(out);
        return out;
    }

  private:
   std::vector<T> _mean;
   std::vector<T> _sig;
   const unsigned int _npars;

   void sample_intern(unsigned int nsample, const T &v, const T &s, std::vector<T> &sample) const;
   void sample_indirect_intern(unsigned int nsample, const T& v, const T &s, const std::vector<T> &U01, std::vector<T> &sample) const;
};

template <typename T>
void Norm<T>::print(std::ostream &out) const
{
  out << "Norm{";
  for(unsigned int i = 0; i < this->npars() - 1; i++)
  {
        out << _mean[i] << ",";
  }
  out << _mean.back() << " ; ";
  for(unsigned int i = 0; i < this->npars() - 1; i++)
  {
        out << _sig[i] << ",";
  }
  out << _sig.back() << "}";
}

template <typename T>
void Norm<T>::sample_indirect(unsigned int nsample, const std::vector<std::vector<T> >&U01, std::vector<std::vector<T> >&sample) const
{
  sample.clear();
  sample.resize(this->npars());
  for(unsigned int ipar = 0; ipar < this->npars(); ipar++)
  {
    this->sample_indirect(nsample,_mean[ipar],_sig[ipar],U01[ipar],sample[ipar]);
  }
  return;
}


template <typename T>
void Norm<T>::sample_indirect_intern(unsigned int nsample, const T &v, const T &s, const std::vector<T> &U01, std::vector<T> &sample) const
{
  sample.clear();
  for(unsigned int i = 0; i < nsample; i++)
  {
    sample.push_back(gsl_cdf_gaussian_Pinv(U01[i],s) + v);
  }
  return;
}

template <typename T>
void Norm<T>::sample_intern(unsigned int nsample, const T &v, const T &s, std::vector<T> &sample) const
{
  sample.clear();
  for(unsigned int i = 0; i < nsample; i++)
  {
    sample.push_back(gsl_ran_gaussian(this->_r,s) + v);
  }
  return;
}

template <typename T>
void Norm<T>::sample(unsigned int nsample, std::vector<std::vector<T> >&sample) const
{
  sample.clear();
  sample.resize(this->npars());
  for(unsigned int ipar = 0; ipar < this->npars(); ipar++)
  {
    this->sample_intern(nsample,_mean[ipar],_sig[ipar],sample[ipar]);
  }
  return;
}

//!unif distribution
template <typename T = double>
class Unif:public Distribution<T>
{
  public:
    Unif():_npars(2){}
    Unif(const std::vector<std::vector<T> > &pars):_npars(2) {this->set_parameters(pars);}
    Unif(const std::vector<T> & min, const std::vector<T> & max):_min(min),_max(max),_npars(2){}
    Unif(const T & min, const T & max):_npars(2){set_parameters(min,max);}
    Unif(const Unif<T> &rhs):_npars(2){if(this == &rhs)return;_min = rhs.min();_max=rhs.max();this->_r=rhs.random_gsl();}
    ~Unif(){}

   void set_parameters(const std::vector<std::vector<T> >&parameters){_min = parameters[0]; _max = parameters[1];}
   void set_parameters(const T & min, const T & max){_min.push_back(min);_max.push_back(max);}

   void sample(unsigned int nsample, std::vector<std::vector<T> > &sample) const;
   void sample_indirect(unsigned int nsample, const std::vector<std::vector<T> >&U01, std::vector<std::vector<T> >&sample) const;

   void sample(unsigned int nsample, std::vector<T> &sample) const {this->sample_intern(nsample,_min[0],_max[0],sample);}
   void sample_indirect(unsigned int nsample, const std::vector<T> &U01, std::vector<T> &sample) const
                                                 {this->sample_indirect_intern(nsample,_min[0],_max[0],U01,sample);}
   
   unsigned int npars()       const {return _min.size();}
   const std::vector<T> min() const {return _min;}
   const std::vector<T> max() const {return _max;}

   void print(std::ostream &out = std::cout) const;
   friend std::ostream & operator << (std::ostream &out, const Unif<T> &uniform)
    {
        uniform.print(out);
        return out;
    }


  private:
   std::vector<T> _min;
   std::vector<T> _max;
   const unsigned int _npars;

   void sample_intern(unsigned int nsample, const T &mi, const T& max, std::vector<T> &sample) const;
   void sample_indirect_intern(unsigned int nsample, const T &mi, const T &ma, const std::vector<T> &U01, std::vector<T> &sample) const;

};

template <typename T>
void Unif<T>::print(std::ostream &out) const
{
  out << "Unif{";
  for(unsigned int i = 0; i < this->npars() - 1; i++)
  {
        out << "[" << _min[i] << "," << _max[i] << "],";
  }
  out << "[" << _min.back() << "," << _max.back() << "]}";
}

template <typename T>
void Unif<T>::sample_intern(unsigned int nsample, const T& mi, const T& ma, std::vector<T> &sample) const
{
  sample.clear();
  for(unsigned int i = 0; i < nsample; i++)
  {
    sample.push_back(gsl_ran_flat(this->_r,mi,ma));
  }
  return;
}

template <typename T>
void Unif<T>::sample(unsigned int nsample, std::vector<std::vector<T> >&sample) const
{
  sample.clear();
  sample.resize(this->npars());
  for(unsigned int ipar = 0; ipar < this->npars(); ipar++)
  {
    this->sample_intern(nsample,_min[ipar],_max[ipar],sample[ipar]);
  }
  return;
}


template <typename T>
void Unif<T>::sample_indirect_intern(unsigned int nsample, const T &mi, const T& ma, const std::vector<T> &U01, std::vector<T> &sample) const
{
  sample.clear();
  for(unsigned int i = 0; i < nsample; i++)
  {
    sample.push_back(gsl_cdf_flat_Pinv(U01[i],mi,ma));
  }
  return;
}

template <typename T>
void Unif<T>::sample_indirect(unsigned int nsample, const std::vector<std::vector<T> >&U01, std::vector<std::vector<T> > &sample) const
{
  sample.clear();
  sample.resize(this->npars());
  for(unsigned int ipar =0; ipar < this->npars(); ipar++)
  {
    this->sample_indirect_intern(nsample,_min[ipar],_max[ipar],U01[ipar],sample[ipar]);
  }
  return;
}


//!logn distribution
template <typename T = double>
class LogN:public Distribution<T>
{
  public:
    LogN():_npars(2){}
    LogN(const std::vector<std::vector<T> >&pars):_npars(2) {this->set_parameters(pars);}
    LogN(const std::vector<T> & v, const std::vector<T> & s):_mean(v),_F(s),_npars(2){}
    LogN(const T & v, const T & s):_npars(2){set_parameters(v,s);}
    LogN(const LogN<T> &rhs):_npars(2){if(this == &rhs)return;_mean = rhs.mean();_F=rhs.F();this->_r=rhs.random_gsl();}
    ~LogN(){}

   void set_parameters(const std::vector<std::vector<T> >&parameters){_mean = parameters[0]; _F = parameters[1];}
   void set_parameters(const T & v, const T & s){_mean.push_back(v);_F.push_back(s);}

   void sample(unsigned int nsample, std::vector<std::vector<T> >&sample) const;
   void sample_indirect(unsigned int nsample, const std::vector<std::vector<T> >&U01, std::vector<std::vector<T> >&sample) const;

   void sample(unsigned int nsample, std::vector<T> &sample) const {this->sample_intern(nsample,_mean[0],_F[0],sample);}
   void sample_indirect(unsigned int nsample, const std::vector<T> &U01, std::vector<T> &sample) const
                                                 {this->sample_indirect_intern(nsample,_mean[0],_F[0],U01,sample);}
   
   unsigned int npars()        const {return _mean.size();}
   const std::vector<T> mean() const {return _mean;}
   const std::vector<T> F()    const {return _F;}

   void print(std::ostream &out = std::cout) const;
   friend std::ostream & operator << (std::ostream &out, const LogN<T> &lognormal)
    {
        lognormal.print(out);
        return out;
    }

  private:
   std::vector<T> _mean;
   std::vector<T> _F;
   const unsigned int _npars;

   void sample_intern(unsigned int nsample, const T& v, const T& F, std::vector<T> &sample) const;
   void sample_indirect_intern(unsigned int nsample, const T& v, const T& F, const std::vector<T> &U01, std::vector<T> &sample) const;

};

template <typename T>
void LogN<T>::print(std::ostream &out) const
{
  out << "LogN{";
  for(unsigned int i = 0; i < this->npars() - 1; i++)
  {
        out << _mean[i] << ",";
  }
  out << _mean.back() << " ; ";
  for(unsigned int i = 0; i < this->npars() - 1; i++)
  {
        out << _F[i] << ",";
  }
  out << _F.back() << "}";
}


template <typename T>
void LogN<T>::sample_intern(unsigned int nsample, const T& v, const T& F, std::vector<T> &sample) const
{
  sample.clear();
  double sigma = std::log(F);
  double zeta  = std::log(v);
  for(unsigned int i = 0; i < nsample; i++)
  {
    sample.push_back(gsl_ran_lognormal(this->_r,zeta,sigma));
  }
  return;
}

template <typename T>
void LogN<T>::sample(unsigned int nsample, std::vector<std::vector<T> >&sample) const
{
  sample.clear();
  sample.resize(this->npars());
  for(unsigned int ipar = 0; ipar < this->npars(); ipar++)
  {
    this->sample_intern(nsample,_mean[ipar],_F[ipar],sample[ipar]);
  }
  return;
}


template <typename T>
void LogN<T>::sample_indirect_intern(unsigned int nsample, const T& v, const T &F, const std::vector<T> &U01, std::vector<T> &sample) const
{
  sample.clear();
  double sigma = std::log(F);
  double zeta  = std::log(v);
  for(unsigned int i = 0; i < nsample; i++)
  {
    sample.push_back(gsl_cdf_lognormal_Pinv(U01[i],zeta,sigma));
  }
  return;
}

template <typename T>
void LogN<T>::sample_indirect(unsigned int nsample, const std::vector<std::vector<T> >&U01, std::vector<std::vector<T> >&sample) const
{
  sample.clear();
  sample.resize(this->npars());
  for(unsigned int ipar = 0; ipar < this->npars(); ipar++)
  {
    this->sample_indirect_intern(nsample,_mean[ipar],_F[ipar],U01[ipar],sample[ipar]);
  }
  return;
}

//!logu distribution
template <typename T = double>
class LogU:public Distribution<T>
{
  public:
    LogU():_npars(2){}
    LogU(const std::vector<std::vector<T> >&pars):_npars(2) {this->set_parameters(pars);}
    LogU(const std::vector<T> & mi, const std::vector<T> & ma):_min(mi),_max(ma),_npars(2){}
    LogU(const T & mi, const T & ma):_npars(2){set_parameters(mi,ma);}
    LogU(const LogU<T> &rhs):_npars(2){if(this == &rhs)return;_min = rhs.min();_max=rhs.max();this->_r=rhs.random_gsl();}
    ~LogU(){}

   void set_parameters(const std::vector<std::vector<T> >&parameters){_min = parameters[0]; _max = parameters[1];}
   void set_parameters(const T & mi, const T & ma){_min.push_back(mi);_max.push_back(ma);}


   void sample(unsigned int nsample, std::vector<T> &sample) const {sample_intern(nsample,_min[0],_max[0],sample);}
   void sample_indirect(unsigned int nsample, const std::vector<T> &U01, std::vector<T> &sample) const
                                                 {sample_indirect_inter(nsample,_min[0],_max[0],U01,sample);}

   void sample(unsigned int nsample, std::vector<std::vector<T> >&sample) const;
   void sample_indirect(unsigned int nsample, const std::vector<std::vector<T> >&U01, std::vector<std::vector<T> >&sample) const;
   
   unsigned int npars()       const {return _min.size();}
   const std::vector<T> min() const {return _min;}
   const std::vector<T> max() const {return _max;}

   void print(std::ostream &out = std::cout) const;
   friend std::ostream & operator << (std::ostream &out, const LogU<T> &loguniform)
    {
        loguniform.print(out);
        return out;
    }

  private:
   std::vector<T> _min;
   std::vector<T> _max;
   const unsigned int _npars;

   void sample_intern(unsigned int nsample, const T& min, const T& max, std::vector<T> &sample) const;
   void sample_indirect_intern(unsigned int nsample, const T &min, const T &max, const std::vector<T> &U01, std::vector<T> &sample) const;
};

template <typename T>
void LogU<T>::print(std::ostream &out) const
{
  out << "LogU{";
  for(unsigned int i = 0; i < this->npars() - 1; i++)
  {
        out << "[" << _min[i] << "," << _max[i] << "],";
  }
  out << "[" << _min.back() << "," << _max.back() << "]}";
}

template <typename T>
void LogU<T>::sample_intern(unsigned int nsample, const T& min, const T& max, std::vector<T> &sample) const
{
  sample.clear();
  double lnmin = std::log(min);
  double lnmax = std::log(max);
  for(unsigned int i = 0; i < nsample; i++)
  {
    sample.push_back(std::exp(gsl_ran_flat(this->_r,lnmin,lnmax)));
  }
  return;
}

template <typename T>
void LogU<T>::sample(unsigned int nsample, std::vector<std::vector<T> >&sample) const
{
  sample.clear();
  sample.resize(this->npars());
  for(unsigned int ipar = 0; ipar < this->npars(); ipar++)
  {
    this->sample_intern(nsample,_min[ipar],_max[ipar],sample[ipar]);
  }
  return;
}


template <typename T>
void LogU<T>::sample_indirect_intern(unsigned int nsample, const T& min, const T & max, const std::vector<T> &U01, std::vector<T> &sample) const
{
  sample.clear();
  double lnmin = std::log(min);
  double lnmax = std::log(max);
  for(unsigned int i = 0; i < nsample; i++)
  {
    sample.push_back(std::exp(gsl_cdf_flat_Pinv(U01[i],lnmin,lnmax)));
  }
  return;
}


template <typename T>
void LogU<T>::sample_indirect(unsigned int nsample, const std::vector<std::vector<T> >&U01, std::vector<std::vector<T> >&sample) const
{
  sample.clear();
  sample.resize(this->npars());
  for(unsigned int ipar = 0; ipar < this->npars(); ipar++)
  {
    this->sample_indirect_intern(nsample,_min[ipar],_max[ipar],U01[ipar],sample[ipar]);
  }
  return;
}


//!nort distribution
template <typename T = double>
class NorT:public Distribution<T>
{
  public:
    NorT():_npars(4){}
    NorT(const std::vector<std::vector<T> >&pars):_npars(4) {this->set_parameters(pars);}
    NorT(const std::vector<T> & v, const std::vector<T> & s, const std::vector<T> & min, const std::vector<T> & max):
                _mean(v),_sig(s),_min(min),_max(max),_npars(4){}
    NorT(const T & v, const T & s, const T & min, const T & max):_npars(4){set_parameters(v,s,min,max);}
    NorT(const NorT<T> &rhs):_npars(4){if(this == &rhs)return;_mean = rhs.mean();_sig=rhs.sigma();_min = rhs.min();_max=rhs.max();this->_r=rhs.random_gsl();}
    ~NorT(){}

   void set_parameters(const std::vector<std::vector<T> >&parameters){_mean = parameters[0]; _sig = parameters[1];_min = parameters[2]; _max = parameters[3];}
   void set_parameters(const T & v, const T & s, const T & min, const T & max)
                {_mean.push_back(v); _sig.push_back(s);_min.push_back(min); _max.push_back(max);}

   void sample(unsigned int nsample, std::vector<T> &sample) const {this->sample_intern(nsample,_mean[0],_sig[0],_min[0],_max[0],sample);}
   void sample_indirect(unsigned int nsample, const std::vector<T> &U01, std::vector<T> &sample) const
                                  {this->sample_indirect_intern(nsample,_mean[0],_sig[0],_min[0],_max[0],U01,sample);}

   void sample(unsigned int nsample, std::vector<std::vector<T> >&sample) const;
   void sample_indirect(unsigned int nsample, const std::vector<std::vector<T> > &U01, std::vector<std::vector<T> > &sample) const;

   unsigned int npars()         const {return _mean.size();}
   const std::vector<T> mean()  const {return _mean;}
   const std::vector<T> sigma() const {return _sig;}
   const std::vector<T> min()   const {return _min;}
   const std::vector<T> max()   const {return _max;}

   void print(std::ostream &out = std::cout) const;
   friend std::ostream & operator << (std::ostream &out, const NorT<T> &normaltrunc)
    {
        normaltrunc.print(out);
        return out;
    }

  private:
   std::vector<T> _mean;
   std::vector<T> _sig;
   std::vector<T> _min;
   std::vector<T> _max;
   const unsigned int _npars;

   void sample_intern(unsigned int nsample, const T &v, const T &s, const T& min, const T& max, std::vector<T> &sample) const;
   void sample_indirect_intern(unsigned int nsample, const T &v, const T& s, const T& min, const T& max, 
                                const std::vector<T> &U01, std::vector<T> &sample) const;
};

template <typename T>
void NorT<T>::print(std::ostream &out) const
{
  out << "NorT{";
  for(unsigned int i = 0; i < this->npars() - 1; i++)
  {
        out << _mean[i] << ",";
  }
  out << _mean.back() << " ; ";
  for(unsigned int i = 0; i < this->npars() - 1; i++)
  {
        out << _sig[i] << ",";
  }
  out << _sig.back() << " ; ";
  for(unsigned int i = 0; i < this->npars() - 1; i++)
  {
        out << _min[i] << ",";
  }
  out << _min.back() << " ; ";
  for(unsigned int i = 0; i < this->npars() - 1; i++)
  {
        out << _max[i] << ",";
  }
  out << _max.back() << "}";
}

template <typename T>
void NorT<T>::sample_intern(unsigned int nsample, const T& v, const T& s, const T& min, const T& max, std::vector<T> &sample) const
{
  sample.clear();
  double minCDF = gsl_cdf_gaussian_P(min - v,s);
  double maxCDF = gsl_cdf_gaussian_P(max - v,s);
  for(unsigned int i = 0; i < nsample; i++)
  {
     sample.push_back(gsl_cdf_gaussian_Pinv(gsl_ran_flat(this->_r,minCDF,maxCDF),s) + v);
  }
  return;
}

template <typename T>
void NorT<T>::sample(unsigned int nsample, std::vector<std::vector<T> >&sample) const
{
  sample.clear();
  sample.resize(this->npars());
  for(unsigned int ipar = 0; ipar < this->npars(); ipar++)
  {
    this->sample_intern(nsample,_mean[ipar],_sig[ipar],_min[ipar],_max[ipar],sample[ipar]);
  }
  return;
}

template <typename T>
void NorT<T>::sample_indirect_intern(unsigned int nsample, const T& v, const T& s, const T& min, const T& max, const std::vector<T> &U01, std::vector<T> &sample) const
{
  sample.clear();
  double minCDF = gsl_cdf_gaussian_P(min - v,s);
  double maxCDF = gsl_cdf_gaussian_P(max - v,s);
  for(unsigned int i = 0; i < nsample; i++)
  {
    sample.push_back(gsl_cdf_gaussian_Pinv(gsl_cdf_flat_Pinv(U01[i],minCDF,maxCDF),s) + v);
  }
  return;
}

template <typename T>
void NorT<T>::sample_indirect(unsigned int nsample, const std::vector<std::vector<T> >&U01, std::vector<std::vector<T> >&sample) const
{
  sample.clear();
  sample.resize(this->npars());
  for(unsigned int ipar = 0; ipar < this->npars(); ipar++)
  {
    this->sample_indirect_intern(nsample,_mean[ipar],_sig[ipar],_min[ipar],_max[ipar],U01[ipar],sample[ipar]);
  }
  return;
}


//!diri distribution
template <typename T = double>
class Diri:public Distribution<T>
{
  public:
    Diri():_mean(0.),_sig(0.),_npars(2){}
    Diri(const std::vector<std::vector<T> > &pars):_npars(2) {this->set_parameters(pars);}
    Diri(const Diri<T> &rhs):_npars(4)
    {
        if(this == &rhs)return;
        _mean = rhs.mean();
        _sig=rhs.sigma();
        this->compute_gamma();
        this->_r=rhs.random_gsl();
    }
    ~Diri(){}

   void set_parameters(const std::vector<std::vector<T> > &parameters);

   void sample(unsigned int nsample,std::vector<std::vector<T> > &sample) const;
   
   unsigned int npars()         const {return _mean.size();}
   const std::vector<T> mean()  const {return _mean;}
   const std::vector<T> sigma() const {return _sig;}

   void print(std::ostream &out = std::cout) const;
   friend std::ostream & operator << (std::ostream &out, const Diri<T> &dirichlet)
    {
        dirichlet.print(out);
        return out;
    }

  private:
   std::vector<T> _mean;
   std::vector<T> _sig;
   T _gamma;
   const unsigned int _npars;

   void compute_gamma();
};

template <typename T>
void Diri<T>::print(std::ostream &out) const
{
  out << "Diri{";
  for(unsigned int i = 0; i < this->npars() - 1; i++)
  {
        out << _mean[i] << ",";
  }
  out << _mean.back() << " ; ";
  for(unsigned int i = 0; i < this->npars() - 1; i++)
  {
        out << _sig[i] << ",";
  }
  out << _sig.back() << "}";
}


template <typename T>
void Diri<T>::set_parameters(const std::vector<std::vector<T> > &parameters)
{
  _mean = parameters[0];
  _sig  = parameters[1];
  this->compute_gamma();
}

template <typename T>
void Diri<T>::compute_gamma()
{
  double a(0.L),b(0.L),mingamma(1e303);
  for(unsigned int i = 0; i < _mean.size(); i++)
  {
     double locmax = (_mean[i] > 0.5L)?_mean[i]:1.L - _mean[i];
     if(locmax < mingamma)mingamma = locmax;

     a += _mean[i] * (1.L - _mean[i]);
     b += _sig[i] * std::sqrt(_mean[i] * (1. - _mean[i]) );
  }
  _gamma = std::pow(a/b,2) - 1.L;
  if(_gamma < 1.L/mingamma)_gamma = 1.L/mingamma;
  return;
}

template <typename T>
void Diri<T>::sample(unsigned int nsample, std::vector<std::vector<T> > &sample) const
{
  sample.clear();
  sample.resize(this->npars());
  double  diri[_mean.size()];
  double alpha[_mean.size()];

  for(unsigned int i = 0; i < _mean.size(); i++)alpha[i] = _gamma * _mean[i];

  for(unsigned int i = 0; i < nsample; i++)
  {
     gsl_ran_dirichlet(this->_r,_mean.size(),alpha,diri);
     for(unsigned int j = 0; j < this->npars(); j++)
     {
       sample[j].push_back(diri[j]);
     }
  }

  return;
}

//!diun distribution
template <typename T = double>
class DiUn:public Distribution<T>
{
  public:
    DiUn():_npars(0){}
    DiUn(const unsigned int npar):_npars(0) {_ndiri = npar;}
    DiUn(const DiUn<T> &rhs):_npars(0){if(this == &rhs)return;_ndiri = rhs.ndiri();}
    ~DiUn(){}

   void set_parameters(const unsigned int npar) {_ndiri = npar;}

   void sample(unsigned int nsample,std::vector<std::vector<T> > &sample) const;
   
   unsigned int npars()         const {return _ndiri;}

   void print(std::ostream &out = std::cout) const {out << "DiUn{" << this->npars() << "}";}
   void print_long(std::ostream &out = std::cout) const;
   friend std::ostream & operator << (std::ostream &out, const DiUn<T> &dirichlet)
    {
        dirichlet.print(out);
        return out;
    }

  private:
   unsigned int _ndiri;
   const unsigned int _npars;
};

template <typename T>
void DiUn<T>::print_long(std::ostream &out) const
{
  out << "DiUn{";
  for(unsigned int i = 0; i < this->npars() - 1; i++)
  {
    out << "1,";
  }
  out << "1}";
}

template <typename T>
void DiUn<T>::sample(unsigned int nsample, std::vector<std::vector<T> > &sample) const
{
  sample.clear();
  sample.resize(this->npars());
  double alpha[_ndiri];
  double diri[_ndiri];
  for(unsigned int i = 0; i < _ndiri; i++)alpha[i] = 1.L;

  for(unsigned int i = 0; i < nsample; i++)
  {
     gsl_ran_dirichlet(this->_r,_ndiri,alpha,diri);
     std::vector<T> diriv;
     for(unsigned int j = 0; j < this->npars(); j++)
     {
       sample[j].push_back(diri[j]);
     }
  }

  return;
}

//!dirg distribution
template <typename T = double>
class DirG:public Distribution<T>
{
  public:
    DirG():_mean(0.),_sig(0.),_npars(2){}
    DirG(const std::vector<std::vector<T> > &pars):_npars(2) {this->set_parameters(pars);}
    DirG(const DirG<T> &rhs):_npars(4){if(this == &rhs)return;_mean = rhs.mean();_sig=rhs.sigma();this->_r=rhs.random_gsl();}
    ~DirG(){}

   void set_parameters(const std::vector<std::vector<T> > &parameters);

   void sample(unsigned int nsample,std::vector<std::vector<T> > &sample) const;
   void sample_indirect(unsigned int nsample, std::vector<std::vector<T> > &U01, std::vector<std::vector<T> > &sample) const;
   
   unsigned int npars()         const {return _mean.size();}
   const std::vector<T> mean()  const {return _mean;}
   const std::vector<T> sigma() const {return _sig;}

   void print(std::ostream &out = std::cout) const;
   friend std::ostream & operator << (std::ostream &out, const DirG<T> &dirichlet)
    {
        dirichlet.print(out);
        return out;
    }

  private:
   std::vector<T> _mean;
   std::vector<T> _sig;
   const unsigned int _npars;
};

template <typename T>
void DirG<T>::print(std::ostream &out) const
{
  out << "DirG{";
  for(unsigned int i = 0; i < this->npars() - 1; i++)
  {
        out << _mean[i] << ",";
  }
  out << _mean.back() << " ; ";
  for(unsigned int i = 0; i < this->npars() - 1; i++)
  {
        out << _sig[i] << ",";
  }
  out << _sig.back() << "}";
}

template <typename T>
void DirG<T>::set_parameters(const std::vector<std::vector<T> > &parameters)
{
  _mean = parameters[0];
  _sig  = parameters[1];
}

template <typename T>
void DirG<T>::sample(unsigned int nsample, std::vector<std::vector<T> > &sample) const
{
  sample.clear();
  sample.resize(this->npars());
//dirg
  double dirgA[_mean.size()];
  double dirgB[_mean.size()];
  for(unsigned int i = 0; i < _mean.size(); i++)
  {
     dirgA[i] = std::pow(_mean[i]/_sig[i],2);
     dirgB[i] = std::pow(_sig[i],2)/_mean[i];
  }

  for(unsigned int i = 0; i < nsample; i++)
  {
    double x[_mean.size()],sum(0.L);
    std::vector<T> diriv;
    for(unsigned int ibr = 0; ibr < this->npars(); ibr++)
    {
       x[ibr] = gsl_ran_gamma_knuth(this->_r,dirgA[ibr],dirgB[ibr]);
       sum += x[ibr];
    }
    for(unsigned int ibr = 0; ibr < this->npars(); ibr++)
    {
       sample[ibr].push_back(x[ibr]/sum);
    }
  }

  return;
}

template <typename T>
void DirG<T>::sample_indirect(unsigned int nsample, std::vector<std::vector<T> > &U01, std::vector<std::vector<T> > &sample) const
{
  sample.clear();
  sample.resize(this->npars());
//dirg
  double dirgA[_mean.size()];
  double dirgB[_mean.size()];
  for(unsigned int i = 0; i < _mean.size(); i++)
  {
     dirgA[i] = std::pow(_mean[i]/_sig[i],2);
     dirgB[i] = std::pow(_sig[i],2)/_mean[i];
  }

  for(unsigned int i = 0; i < nsample; i++)
  {
    double x[_mean.size()],sum(0.L);
    std::vector<T> diriv;
    for(unsigned int ibr = 0; ibr < this->npars(); ibr++)
    {
       x[ibr] = gsl_cdf_gamma_Pinv(U01[i][ibr],dirgA[ibr],dirgB[ibr]);
       sum += x[ibr];
    }
    for(unsigned int ibr = 0; ibr < _mean.size(); ibr++)
    {
       sample[ibr].push_back(x[ibr]/sum);
    }
  }

  return;
}

//!diut distribution
template <typename T = double>
class DiUT:public Distribution<T>
{
  public:
    DiUT():_npars(2){}
    DiUT(const std::vector<std::vector<T> > &pars):_npars(2) {this->set_parameters(pars);}
    DiUT(const DiUT<T> &rhs):_npars(4){if(this == &rhs)return;_min = rhs.min();_max=rhs.max();this->_r=rhs.random_gsl();}
    ~DiUT(){}

   void set_parameters(const std::vector<std::vector<T> > &parameters);

   void sample(unsigned int nsample,std::vector<std::vector<T> > &sample) const;
   
   unsigned int npars()       const {return _min.size();}
   const std::vector<T> min() const {return _min;}
   const std::vector<T> max() const {return _max;}

   void print(std::ostream &out = std::cout);
   friend std::ostream & operator << (std::ostream &out, const DiUT<T> &dirichlet)
    {
        dirichlet.print(out);
        return out;
    }

  private:
   std::vector<T> _min;
   std::vector<T> _max;
   const unsigned int _npars;

   bool accept(const double *run) const;
};

template <typename T>
void DiUT<T>::print(std::ostream &out)
{
  out << "DiUT{";
  for(unsigned int i = 0; i < this->npars() - 1; i++)
  {
        out << "[" << _min[i] << "," << _max[i] << "],";
  }
  out << "[" << _min.back() << "," << _max.back() << "]}";
}

template <typename T>
bool DiUT<T>::accept(const double *run) const
{
  for(unsigned int i = 0; i < _min.size(); i++)
  {
     if(_min[i] > run[i] || _max[i] < run[i])return false;
  }

  return true;
}

template <typename T>
void DiUT<T>::set_parameters(const std::vector<std::vector<T> > &parameters)
{
  _min = parameters[0];
  _max = parameters[1];
}

template <typename T>
void DiUT<T>::sample(unsigned int nsample, std::vector<std::vector<T> > &sample) const
{
  sample.clear();
  sample.resize(this->npars());
//diut -> troncated diun
  double alpha[_min.size()];
  double diri[_min.size()];
  for(unsigned int i = 0; i < _min.size(); i++)alpha[i] = 1.L;

  unsigned int ns(0);
  while(ns != nsample)
  {
     gsl_ran_dirichlet(this->_r,_min.size(),alpha,diri);
     if(!accept(diri))continue;
     ns++;
     std::vector<T> diriv;
     for(unsigned int j = 0; j < this->npars(); j++)
     {
       sample[j].push_back(diri[j]);
     }
  }

  return;
}

//!dior distribution
template <typename T = double>
class DiOr:public Distribution<T>
{
  public:
    DiOr():_npars(2){}
    DiOr(const std::vector<T> &pars):_npars(2) {this->set_parameters(pars);}
    DiOr(const DiOr<T> &rhs):_npars(4){if(this == &rhs)return;_par = rhs.index_pars();this->_r=rhs.random_gsl();}
    ~DiOr(){}

   void set_parameters(const std::vector<T>  &parameters);

   void sample(unsigned int nsample,std::vector<std::vector<T> > &sample) const;
   
   unsigned int npars()              const {return _par.size();}
   const std::vector<T> index_pars() const {return _par;}

   void print(std::ostream &out = std::cout);
   friend std::ostream & operator << (std::ostream &out, const DiOr<T> &dirichlet)
    {
        dirichlet.print(out);
        return out;
    }

  private:
   std::vector<T> _par;
   const unsigned int _npars;
   unsigned int offset;
};

template <typename T>
void DiOr<T>::print(std::ostream &out)
{
  out << "DiOr{";
  for(unsigned int i = 0; i < this->npars() - 1; i++)
  {
       out << _par[i] << ",";
  }
  out << _par.back() << "}";
}

template <typename T>
void DiOr<T>::set_parameters(const std::vector<T>  &parameters)
{
  _par = parameters;
  offset = _par.size();
  for(unsigned int i = 0; i < this->npars(); i++)
  {
     if(_par[i] < offset)offset = _par[i];
  }
}

template <typename T>
void DiOr<T>::sample(unsigned int nsample, std::vector<std::vector<T> > &sample) const
{
  sample.clear();
  sample.resize(this->npars());
//dior -> sorted diun
  double alpha[_par.size()];
  double diri[_par.size()];
  for(unsigned int i = 0; i < this->npars(); i++)alpha[i] = 1.L;

  for(unsigned int ns = 0; ns < nsample; ns++)
  {
     gsl_ran_dirichlet(this->_r,this->npars(),alpha,diri);
     gsl_sort(diri,1,_par.size());
     std::vector<T> diriv;
     for(unsigned int j = 0; j < this->npars(); j++)
     {
       sample[j].push_back(diri[_par.size() - 1 - (unsigned int)_par[j] + offset]); //reverse order
     }
  }

  return;
}

}
#endif
