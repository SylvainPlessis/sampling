#ifndef _DELTA_PDF_
#define _DELTA_PDF_

//
#include "sampler.h"

namespace Sampler{
//!delta distribution
template <typename T = double>
class Delt:public Distribution<T>
{
  public:
    Delt():_npars(1){}
    Delt(const std::vector<T> & val):_val(val),_npars(1){}
    Delt(const T & val):_npars(1){set_parameters(val);}
    Delt(const Delt<T> &rhs):_npars(1){if(this == &rhs)return;_val = rhs.val();this->_r=rhs.random_gsl();}
    ~Delt(){}

   void set_parameters(const std::vector<T> &parameters){_val = parameters[0];}
   void set_parameters(const T & val){_val.push_back(val);}

   void sample(unsigned int nsample, std::vector<std::vector<T> > &sample) const;
   void sample_indirect(unsigned int nsample, const std::vector<std::vector<T> >&U01, std::vector<std::vector<T> >&sample) const;

   void sample(unsigned int nsample, std::vector<T> &sample) const {this->sample_intern(nsample,_val[0],sample);}
   void sample_indirect(unsigned int nsample, const std::vector<T> &U01, std::vector<T> &sample) const
                                                 {this->sample_indirect_intern(nsample,_val[0],U01,sample);}
   
   unsigned int npars()       const {return _val.size();}
   const std::vector<T> val() const {return _val;}

   void print(std::ostream &out = std::cout) const;
   friend std::ostream & operator << (std::ostream &out, const Delt<T> &uniform)
    {
        uniform.print(out);
        return out;
    }


  private:
   std::vector<T> _val;
   const unsigned int _npars;

   void sample_intern(unsigned int nsample, const T &va, std::vector<T> &sample) const;
   void sample_indirect_intern(unsigned int nsample, const T &va, const std::vector<T> &U01, std::vector<T> &sample) const;

};

template <typename T>
void Delt<T>::print(std::ostream &out) const
{
  out << "Delt{";
  for(unsigned int i = 0; i < this->npars() - 1; i++)
  {
        out << _val[i] << ",";
  }
  out << _val.back() <<  "}";
}

template <typename T>
void Delt<T>::sample_intern(unsigned int nsample, const T& val, std::vector<T> &sample) const
{
  sample.clear();
  sample.resize(nsample,val);
  return;
}

template <typename T>
void Delt<T>::sample(unsigned int nsample, std::vector<std::vector<T> >&sample) const
{
  sample.clear();
  sample.resize(this->npars());
  for(unsigned int ipar = 0; ipar < this->npars(); ipar++)
  {
    this->sample_intern(nsample,_val[ipar],sample[ipar]);
  }
  return;
}


template <typename T>
void Delt<T>::sample_indirect_intern(unsigned int nsample, const T &val, const std::vector<T> &U01, std::vector<T> &sample) const
{
  sample.clear();
  sample.resize(nsample,val);
  return;
}

template <typename T>
void Delt<T>::sample_indirect(unsigned int nsample, const std::vector<std::vector<T> >&U01, std::vector<std::vector<T> > &sample) const
{
  sample.clear();
  sample.resize(this->npars());
  for(unsigned int ipar =0; ipar < this->npars(); ipar++)
  {
    this->sample_indirect_intern(nsample,_val[ipar],U01[ipar],sample[ipar]);
  }
  return;
}



}

#endif
