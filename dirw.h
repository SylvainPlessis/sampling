#ifndef _DIRICHLET_WONG_PDF_
#define _DIRICHLET_WONG_PDF_

//
#include "sampler.h"

namespace Sampler{
//!dirw distribution
template <typename T = double>
class DirW:public Distribution<T>
{
  public:
    DirW():_mean(0.),_sig(0.),_npars(2){}
    DirW(const std::vector<std::vector<T> > &pars):_npars(2) {this->set_parameters(pars);}
    DirW(const DirW<T> &rhs):_npars(4){if(this == &rhs)return;_mean = rhs.mean();_sig=rhs.sigma();this->_r=rhs.random_gsl();}
    ~DirW(){}

   void set_parameters(const std::vector<std::vector<T> > &parameters);

   void sample(unsigned int nsample,std::vector<std::vector<T> > &sample) const;
   void sample_indirect(unsigned int nsample, std::vector<std::vector<T> > &U01, std::vector<std::vector<T> > &sample) const;
   
   unsigned int npars()         const {return _mean.size();}
   const std::vector<T> mean()  const {return _mean;}
   const std::vector<T> sigma() const {return _sig;}

   void print(std::ostream &out = std::cout) const;
   friend std::ostream & operator << (std::ostream &out, const DirW<T> &dirichlet)
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
void DirW<T>::print(std::ostream &out) const
{
  out << "DirW{";
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
void DirW<T>::set_parameters(const std::vector<std::vector<T> > &parameters)
{
  _mean = parameters[0];
  _sig  = parameters[1];
}

template <typename T>
void DirW<T>::sample(unsigned int nsample, std::vector<std::vector<T> > &sample) const
{
  sample.clear();
  sample.resize(this->npars());
//dirw
  for(unsigned int i = 0; i < nsample; i++)
  {
    double z[_mean.size() - 1];
    double x[_mean.size()];
    std::vector<T> diriv;
    for(unsigned int ibr = 0; ibr < this->npars() - 1; ibr++)
    {
       z[ibr] = gsl_ran_beta(this->_r,_mean[ibr],_sig[ibr]);
    }
    double sum(0.L);
    for(unsigned int ibr = 0; ibr < this->npars()-1; ibr++)
    {
       x[ibr] = z[ibr];
       double coef(1.L);
       for(unsigned int j = 0; j < ibr; j++)
       {
          coef -= x[j];
       }
       x[ibr] *= coef;
       sum += x[ibr];
    }
    x[this->npars() - 1] = 1.L - sum;
    for(unsigned int ibr = 0; ibr < this->npars(); ibr++)
    {
       sample[ibr].push_back(x[ibr]);
    }
  }

  return;
}

template <typename T>
void DirW<T>::sample_indirect(unsigned int nsample, std::vector<std::vector<T> > &U01, std::vector<std::vector<T> > &sample) const
{
  return;
}

}

#endif
