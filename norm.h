#ifndef _NORMAL_PDF_
#define _NORMAL_PDF_

//
#include "sampler.h"

namespace Sampler{
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


}

#endif
