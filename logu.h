#ifndef _LOGUNIFORM_PDF_
#define _LOGUNIFORM_PDF_

//
#include "sampler.h"

namespace Sampler{
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

}

#endif
