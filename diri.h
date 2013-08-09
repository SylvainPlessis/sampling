#ifndef _DIRICHLET_PDF_
#define _DIRICHLET_PDF_

//
#include "sampler.h"

namespace Sampler{
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


}

#endif
