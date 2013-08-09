#ifndef _DIRICHLET_GENERALISE_PDF_
#define _DIRICHLET_GENERALISE_PDF_

//
#include "sampler.h"

namespace Sampler{
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

}

#endif
