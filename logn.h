#ifndef _LOGNORMAL_PDF_
#define _LOGNORMAL_PDF_

//
#include "sampler.h"

namespace Sampler{
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


}

#endif
