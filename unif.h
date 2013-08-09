#ifndef _UNIFORM_PDF_
#define _UNIFORM_PDF_

//
#include "sampler.h"

namespace Sampler{
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



}

#endif
