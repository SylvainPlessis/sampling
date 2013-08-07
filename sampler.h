#ifndef _SAMPLER_FROM_GSL_
#define _SAMPLER_FROM_GSL_

//C++
#include <cmath>
#include <vector>

//GSL
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_sort.h>

namespace Sampler
{

template<typename DistType, typename Scalar = double>
class Distribution
{
  public:
    Distribution(const gsl_rng_type * rT = gsl_rng_knuthran2002);
    ~Distribution();

    void set_distribution(const std::vector<Scalar> &parameters)              {_sampler.set_parameters(parameters);}
    void set_distribution(const std::vector<std::vector<Scalar> > &parameters){_sampler.set_parameters(parameters);}
    void set_distribution(const unsigned int npars)                           {_sampler.set_parameters(npars);}

    void sample(unsigned int nsample, std::vector<Scalar> &sample)              const {_sampler.sample(nsample,sample);}
    void sample(unsigned int nsample, std::vector<std::vector<Scalar> >&sample) const {_sampler.sample(nsample,sample);}

    void sample_indirect(unsigned int nsample, std::vector<Scalar> &U01, 
                                               std::vector<Scalar> &sample) const {_sampler.sample_indirect(nsample,U01,sample);}
    void sample_indirect(unsigned int nsample, std::vector<std::vector<Scalar> > &U01, 
                                               std::vector<std::vector<Scalar> > &sample)
                                                                            const {_sampler.sample_indirect(nsample,U01,sample);}

    const DistType &sampler() {return _sampler;}

  private:
   const gsl_rng_type * _T;
   gsl_rng * _r;
   DistType _sampler;

};

template<typename DistType, typename Scalar>
Distribution<DistType,Scalar>::Distribution(const gsl_rng_type * rT):
_T(rT)
{
  gsl_rng_env_setup();
  _r = gsl_rng_alloc(_T);
  _sampler.set_r(_r);
}

template<typename DistType, typename Scalar>
Distribution<DistType,Scalar>::~Distribution()
{
  gsl_rng_free(_r);
}
//!gaussian distribution
template <typename T = double>
class Norm
{
  public:
    Norm():_mean(0.),_sig(0.),_r(NULL),_npars(2){}
    Norm(const std::vector<T> &pars):_npars(2) {this->set_parameters(pars);}
    Norm(const T & v, const T & s, const gsl_rng * ra):_mean(v),_sig(s),_r(ra),_npars(2){}
    Norm(const Norm<T> &rhs):_npars(2){if(this == &rhs)return;_mean = rhs.mean();_sig=rhs.sigma();_r=rhs.random_gsl();}
    ~Norm(){}

   void set_parameters(const std::vector<T> &parameters){_mean = parameters[0]; _sig = parameters[1];}

   void sample(unsigned int nsample,std::vector<T> &sample) const;
   void sample_indirect(unsigned int nsample, std::vector<T> &U01, std::vector<T> &sample) const;
   
   unsigned int npars() const {return 1;}
   const T mean()       const {return _mean;}
   const T sigma()      const {return _sig;}

//GSL stuff
   void  set_r(gsl_rng *rD)       {_r = rD;}
   const gsl_rng * random_gsl()   const {return _r;}

  private:
   T _mean;
   T _sig;
   gsl_rng *_r;
   const unsigned int _npars;
};

template <typename T>
void Norm<T>::sample_indirect(unsigned int nsample, std::vector<T> &U01, std::vector<T> &sample) const
{
  sample.clear();
  for(unsigned int i = 0; i < nsample; i++)
  {
    sample.push_back(gsl_cdf_gaussian_Pinv(U01[i],_sig) + _mean);
  }
  return;
}

template <typename T>
void Norm<T>::sample(unsigned int nsample, std::vector<T> &sample) const
{
  sample.clear();
  for(unsigned int i = 0; i < nsample; i++)
  {
     sample.push_back(gsl_ran_gaussian(_r,_sig) + _mean);
  }
  return;
}

//!unif distribution
template <typename T = double>
class Unif
{
  public:
    Unif():_min(0.),_max(0.),_r(NULL),_npars(2){}
    Unif(const std::vector<T> &pars):_npars(2) {this->set_parameters(pars);}
    Unif(const T & v, const T & s, const gsl_rng *ra):_min(v),_max(s),_r(ra),_npars(2){}
    Unif(const Unif<T> &rhs):_npars(2){if(this == &rhs)return;_min = rhs.min();_max=rhs.max();_r=rhs.random_gsl();}
    ~Unif(){}

   void set_parameters(const std::vector<T> &parameters){_min = parameters[0]; _max = parameters[1];}

   void sample(unsigned int nsample,std::vector<T> &sample) const;
   void sample_indirect(unsigned int nsample, std::vector<T> &U01, std::vector<T> &sample) const;
   
   unsigned int npars() const {return 1;}
   const T min()        const {return _min;}
   const T max()        const {return _max;}

//GSL stuff
   void  set_r(gsl_rng *rD)       {_r = rD;}
   const gsl_rng * random_gsl()   const {return _r;}

  private:
   T _min;
   T _max;
   gsl_rng *_r;
   const unsigned int _npars;
};

template <typename T>
void Unif<T>::sample(unsigned int nsample, std::vector<T> &sample) const
{
  sample.clear();
  for(unsigned int i = 0; i < nsample; i++)
  {
     sample.push_back(gsl_ran_flat(_r,_min,_max));
  }
  return;
}

template <typename T>
void Unif<T>::sample_indirect(unsigned int nsample, std::vector<T> &U01, std::vector<T> &sample) const
{
  sample.clear();
  for(unsigned int i = 0; i < nsample; i++)
  {
    sample.push_back(gsl_cdf_flat_Pinv(U01[i],_min,_max));
  }
  return;
}


//!logn distribution
template <typename T = double>
class LogN
{
  public:
    LogN():_mean(0.),_F(0.),_r(NULL),_npars(2){}
    LogN(const std::vector<T> &pars):_npars(2) {this->set_parameters(pars);}
    LogN(const T & v, const T & s, const gsl_rng *ra):_mean(v),_F(s),_r(ra),_npars(2){}
    LogN(const LogN<T> &rhs):_npars(2){if(this == &rhs)return;_mean = rhs.mean();_F=rhs.F();_r=rhs.random_gsl();}
    ~LogN(){}

   void set_parameters(const std::vector<T> &parameters){_mean = parameters[0]; _F = parameters[1];}

   void sample(unsigned int nsample,std::vector<T> &sample) const;
   void sample_indirect(unsigned int nsample, std::vector<T> &U01, std::vector<T> &sample) const;
   
   unsigned int npars() const {return 1;}
   const T mean()       const {return _mean;}
   const T F()          const {return _F;}

//GSL stuff
   void  set_r(gsl_rng *rD)       {_r = rD;}
   const gsl_rng * random_gsl()   const {return _r;}

  private:
   T _mean;
   T _F;
   gsl_rng *_r;
   const unsigned int _npars;
};

template <typename T>
void LogN<T>::sample(unsigned int nsample, std::vector<T> &sample) const
{
  sample.clear();
  double sigma = std::exp(_F);
  for(unsigned int i = 0; i < nsample; i++)
  {
     sample.push_back(gsl_ran_lognormal(_r,_mean,sigma));
  }
  return;
}

template <typename T>
void LogN<T>::sample_indirect(unsigned int nsample, std::vector<T> &U01, std::vector<T> &sample) const
{
  sample.clear();
  double sigma = std::exp(_F);
  for(unsigned int i = 0; i < nsample; i++)
  {
    sample.push_back(gsl_cdf_lognormal_Pinv(U01[i],_mean,sigma));
  }
  return;
}

//!logu distribution
template <typename T = double>
class LogU
{
  public:
    LogU():_min(0.),_max(0.),_r(NULL),_npars(2){}
    LogU(const std::vector<T> &pars):_npars(2) {this->set_parameters(pars);}
    LogU(const T & v, const T & s, const gsl_rng *ra):_min(v),_max(s),_r(ra),_npars(2){}
    LogU(const LogU<T> &rhs):_npars(2){if(this == &rhs)return;_min = rhs.min();_max=rhs.max();_r=rhs.random_gsl();}
    ~LogU(){}

   void set_parameters(const std::vector<T> &parameters){_min = parameters[0]; _max = parameters[1];}

   void sample(unsigned int nsample,std::vector<T> &sample) const;
   void sample_indirect(unsigned int nsample, std::vector<T> &U01, std::vector<T> &sample) const;
   
   unsigned int npars() const {return 1;}
   const T min()        const {return _min;}
   const T max()        const {return _max;}

//GSL stuff
   void  set_r(gsl_rng *rD)       {_r = rD;}
   const gsl_rng * random_gsl()   const {return _r;}

  private:
   T _min;
   T _max;
   gsl_rng *_r;
   const unsigned int _npars;
};

template <typename T>
void LogU<T>::sample(unsigned int nsample, std::vector<T> &sample) const
{
  sample.clear();
  double lnmin = std::log(_min);
  double lnmax = std::log(_max);
  for(unsigned int i = 0; i < nsample; i++)
  {
     sample.push_back(std::exp(gsl_ran_flat(_r,lnmin,lnmax)));
  }
  return;
}

template <typename T>
void LogU<T>::sample_indirect(unsigned int nsample, std::vector<T> &U01, std::vector<T> &sample) const
{
  sample.clear();
  double lnmin = std::log(_min);
  double lnmax = std::log(_max);
  for(unsigned int i = 0; i < nsample; i++)
  {
    sample.push_back(std::exp(gsl_cdf_flat_Pinv(U01[i],lnmin,lnmax)));
  }
  return;
}


//!nort distribution
template <typename T = double>
class NorT
{
  public:
    NorT():_mean(0.),_sig(0.),_min(0.),_max(0.),_r(NULL),_npars(4){}
    NorT(const std::vector<T> &pars):_npars(4) {this->set_parameters(pars);}
    NorT(const T & v, const T & s, const T & min, const T & max, const gsl_rng *ra):
                _min(v),_max(s),_r(ra),_npars(4){}
    NorT(const NorT<T> &rhs):_npars(4){if(this == &rhs)return;_mean = rhs.mean();_sig=rhs.sigma();_min = rhs.min();_max=rhs.max();_r=rhs.random_gsl();}
    ~NorT(){}

   void set_parameters(const std::vector<T> &parameters){_mean = parameters[0]; _sig = parameters[1];_min = parameters[2]; _max = parameters[3];}

   void sample(unsigned int nsample,std::vector<T> &sample) const;
   void sample_indirect(unsigned int nsample, std::vector<T> &U01, std::vector<T> &sample) const;

   unsigned int npars() const {return 1;}
   const T mean()       const {return _mean;}
   const T sigma()      const {return _sig;}
   const T min()        const {return _min;}
   const T max()        const {return _max;}

//GSL stuff
   void  set_r(gsl_rng *rD)       {_r = rD;}
   const gsl_rng * random_gsl()   const {return _r;}

  private:
   T _mean;
   T _sig;
   T _min;
   T _max;
   gsl_rng *_r;
   const unsigned int _npars;
};

template <typename T>
void NorT<T>::sample(unsigned int nsample, std::vector<T> &sample) const
{
  sample.clear();
  double minCDF = gsl_cdf_gaussian_P(_min - _mean,_sig);
  double maxCDF = gsl_cdf_gaussian_P(_max - _mean,_sig);
  for(unsigned int i = 0; i < nsample; i++)
  {
     sample.push_back(gsl_cdf_gaussian_Pinv(gsl_ran_flat(_r,minCDF,maxCDF),_sig) + _mean);
  }
  return;
}

template <typename T>
void NorT<T>::sample_indirect(unsigned int nsample, std::vector<T> &U01, std::vector<T> &sample) const
{
  sample.clear();
  double minCDF = gsl_cdf_gaussian_P(_min - _mean,_sig);
  double maxCDF = gsl_cdf_gaussian_P(_max - _mean,_sig);
  for(unsigned int i = 0; i < nsample; i++)
  {
    sample.push_back(gsl_cdf_gaussian_Pinv(gsl_cdf_flat_Pinv(U01[i],minCDF,maxCDF),_sig) + _mean);
  }
  return;
}


//!diri distribution
template <typename T = double>
class Diri
{
  public:
    Diri():_mean(0.),_sig(0.),_r(NULL),_npars(2){}
    Diri(const std::vector<std::vector<T> > &pars):_npars(2) {this->set_parameters(pars);}
    Diri(const Diri<T> &rhs):_npars(4){if(this == &rhs)return;_mean = rhs.mean();_sig=rhs.sigma();_r=rhs.random_gsl();}
    ~Diri(){}

   void set_parameters(const std::vector<std::vector<T> > &parameters);

   void sample(unsigned int nsample,std::vector<std::vector<T> > &sample) const;
   
   unsigned int npars()         const {return _mean.size();}
   const std::vector<T> mean()  const {return _mean;}
   const std::vector<T> sigma() const {return _sig;}

//GSL stuff
   void  set_r(gsl_rng *rD)       {_r = rD;}
   const gsl_rng * random_gsl()   const {return _r;}

  private:
   std::vector<T> _mean;
   std::vector<T> _sig;
   gsl_rng *_r;
   const unsigned int _npars;
};
template <typename T>
void Diri<T>::set_parameters(const std::vector<std::vector<T> > &parameters)
{
  _mean = parameters[0];
  _sig  = parameters[1];
}

template <typename T>
void Diri<T>::sample(unsigned int nsample, std::vector<std::vector<T> > &sample) const
{
  sample.clear();
  double  diri[_mean.size()];
  double alpha[_mean.size()];
  double a(0.L),b(0.L),mingamma(1e303);
  for(unsigned int i = 0; i < _mean.size(); i++)
  {
     double locmax = (_mean[i] > 0.5L)?_mean[i]:1.L - _mean[i];
     if(locmax < mingamma)mingamma = locmax;

     a += _mean[i] * (1.L - _mean[i]);
     b += _sig[i] * std::sqrt(_mean[i] * (1. - _mean[i]) );
  }
  double gamma = std::pow(a/b,2) - 1.L;
  if(gamma < 1.L/mingamma)gamma = 1.L/mingamma;
  for(unsigned int i = 0; i < _mean.size(); i++)alpha[i] = gamma * _mean[i];

  for(unsigned int i = 0; i < nsample; i++)
  {
     gsl_ran_dirichlet(_r,_mean.size(),alpha,diri);
     std::vector<T> diriv;
     for(unsigned int j = 0; j < _mean.size(); j++)
     {
       diriv.push_back(diri[j]);
     }
     sample.push_back(diriv);
  }

  return;
}

//!diun distribution
template <typename T = double>
class DiUn
{
  public:
    DiUn():_r(NULL),_npars(0){}
    DiUn(const unsigned int npar):_npars(0) {_ndiri = npar;}
    DiUn(const DiUn<T> &rhs):_npars(0){if(this == &rhs)return;_ndiri = rhs.ndiri();}
    ~DiUn(){}

   void set_parameters(const unsigned int npar) {_ndiri = npar;}

   void sample(unsigned int nsample,std::vector<std::vector<T> > &sample) const;
   
   unsigned int npars()         const {return _ndiri;}

//GSL stuff
   void  set_r(gsl_rng *rD)       {_r = rD;}
   const gsl_rng * random_gsl()   const {return _r;}

  private:
   unsigned int _ndiri;
   gsl_rng *_r;
   const unsigned int _npars;
};

template <typename T>
void DiUn<T>::sample(unsigned int nsample, std::vector<std::vector<T> > &sample) const
{
  sample.clear();
  double alpha[_ndiri];
  double diri[_ndiri];
  for(unsigned int i = 0; i < _ndiri; i++)alpha[i] = 1.L;

  for(unsigned int i = 0; i < nsample; i++)
  {
     gsl_ran_dirichlet(_r,_ndiri,alpha,diri);
     std::vector<T> diriv;
     for(unsigned int j = 0; j < _ndiri; j++)
     {
       diriv.push_back(diri[j]);
     }
     sample.push_back(diriv);
  }

  return;
}

//!dirg distribution
template <typename T = double>
class DirG
{
  public:
    DirG():_mean(0.),_sig(0.),_r(NULL),_npars(2){}
    DirG(const std::vector<std::vector<T> > &pars):_npars(2) {this->set_parameters(pars);}
    DirG(const DirG<T> &rhs):_npars(4){if(this == &rhs)return;_mean = rhs.mean();_sig=rhs.sigma();_r=rhs.random_gsl();}
    ~DirG(){}

   void set_parameters(const std::vector<std::vector<T> > &parameters);

   void sample(unsigned int nsample,std::vector<std::vector<T> > &sample) const;
   void sample_indirect(unsigned int nsample, std::vector<std::vector<T> > &U01, std::vector<std::vector<T> > &sample) const;
   
   unsigned int npars()         const {return _mean.size();}
   const std::vector<T> mean()  const {return _mean;}
   const std::vector<T> sigma() const {return _sig;}

//GSL stuff
   void  set_r(gsl_rng *rD)       {_r = rD;}
   const gsl_rng * random_gsl()   const {return _r;}

  private:
   std::vector<T> _mean;
   std::vector<T> _sig;
   gsl_rng *_r;
   const unsigned int _npars;
};

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
    for(unsigned int ibr = 0; ibr < _mean.size(); ibr++)
    {
       x[ibr] = gsl_ran_gamma_knuth(_r,dirgA[ibr],dirgB[ibr]);
       sum += x[ibr];
    }
    for(unsigned int ibr = 0; ibr < _mean.size(); ibr++)
    {
       diriv.push_back(x[ibr]/sum);
    }
    sample.push_back(diriv);
  }

  return;
}

template <typename T>
void DirG<T>::sample_indirect(unsigned int nsample, std::vector<std::vector<T> > &U01, std::vector<std::vector<T> > &sample) const
{
  sample.clear();
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
    for(unsigned int ibr = 0; ibr < _mean.size(); ibr++)
    {
       x[ibr] = gsl_cdf_gamma_Pinv(U01[i][ibr],dirgA[ibr],dirgB[ibr]);
       sum += x[ibr];
    }
    for(unsigned int ibr = 0; ibr < _mean.size(); ibr++)
    {
       diriv.push_back(x[ibr]/sum);
    }
    sample.push_back(diriv);
  }

  return;
}

//!diut distribution
template <typename T = double>
class DiUT
{
  public:
    DiUT():_r(NULL),_npars(2){}
    DiUT(const std::vector<std::vector<T> > &pars):_npars(2) {this->set_parameters(pars);}
    DiUT(const DiUT<T> &rhs):_npars(4){if(this == &rhs)return;_min = rhs.min();_max=rhs.max();_r=rhs.random_gsl();}
    ~DiUT(){}

   void set_parameters(const std::vector<std::vector<T> > &parameters);

   void sample(unsigned int nsample,std::vector<std::vector<T> > &sample) const;
   
   unsigned int npars()       const {return _min.size();}
   const std::vector<T> min() const {return _min;}
   const std::vector<T> max() const {return _max;}

//GSL stuff
   void  set_r(gsl_rng *rD)       {_r = rD;}
   const gsl_rng * random_gsl()   const {return _r;}

  private:
   std::vector<T> _min;
   std::vector<T> _max;
   gsl_rng *_r;
   const unsigned int _npars;

   bool accept(const double *run) const;
};


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
//diut -> troncated diun
  double alpha[_min.size()];
  double diri[_min.size()];
  for(unsigned int i = 0; i < _min.size(); i++)alpha[i] = 1.L;

  unsigned int ns(0);
  while(ns != nsample)
  {
     gsl_ran_dirichlet(_r,_min.size(),alpha,diri);
     if(!accept(diri))continue;
     ns++;
     std::vector<T> diriv;
     for(unsigned int j = 0; j < _min.size(); j++)
     {
       diriv.push_back(diri[j]);
     }
     sample.push_back(diriv);
  }

  return;
}

//!dior distribution
template <typename T = double>
class DiOr
{
  public:
    DiOr():_r(NULL),_npars(2){}
    DiOr(const std::vector<T> &pars):_npars(2) {this->set_parameters(pars);}
    DiOr(const DiOr<T> &rhs):_npars(4){if(this == &rhs)return;_par = rhs.index_pars();_r=rhs.random_gsl();}
    ~DiOr(){}

   void set_parameters(const std::vector<T>  &parameters);

   void sample(unsigned int nsample,std::vector<std::vector<T> > &sample) const;
   
   unsigned int npars()              const {return _par.size();}
   const std::vector<T> index_pars() const {return _par;}

//GSL stuff
   void  set_r(gsl_rng *rD)       {_r = rD;}
   const gsl_rng * random_gsl()   const {return _r;}

  private:
   std::vector<T> _par;
   gsl_rng *_r;
   const unsigned int _npars;
   unsigned int offset;
};

template <typename T>
void DiOr<T>::set_parameters(const std::vector<T>  &parameters)
{
  _par = parameters;
  offset = _par.size();
  for(unsigned int i = 0; i < _par.size(); i++)
  {
     if(_par[i] < offset)offset = _par[i];
  }
}

template <typename T>
void DiOr<T>::sample(unsigned int nsample, std::vector<std::vector<T> > &sample) const
{
  sample.clear();
//dior -> sorted diun
  double alpha[_par.size()];
  double diri[_par.size()];
  for(unsigned int i = 0; i < _par.size(); i++)alpha[i] = 1.L;

  unsigned int ns(0);
  for(unsigned int ns = 0; ns < nsample; ns++)
  {
     gsl_ran_dirichlet(_r,_par.size(),alpha,diri);
     gsl_sort(diri,1,_par.size());
     std::vector<T> diriv;
     for(unsigned int j = 0; j < _par.size(); j++)
     {
       diriv.push_back(diri[_par.size() - 1 - (unsigned int)_par[j] + offset]); //reverse order
     }
     sample.push_back(diriv);
  }

  return;
}

}
#endif
