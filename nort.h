#ifndef _NORMAL_TRUNCATED_PDF_
#define _NORMAL_TRUNCATED_PDF_

//
#include "sampler.h"

namespace Sampler{
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

}

#endif
