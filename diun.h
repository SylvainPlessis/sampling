#ifndef _DIRICHLET_UNIFORM_PDF_
#define _DIRICHLET_UNIFORM_PDF_

//
#include "sampler.h"

namespace Sampler{
//!diun distribution
template <typename T = double>
class DiUn:public Distribution<T>
{
  public:
    DiUn():_npars(0){}
    DiUn(const unsigned int npar):_npars(0) {_ndiri = npar;}
    DiUn(const DiUn<T> &rhs):_npars(0){if(this == &rhs)return;_ndiri = rhs.ndiri();}
    ~DiUn(){}

   void set_parameters(const unsigned int npar) {_ndiri = npar;}

   void sample(unsigned int nsample,std::vector<std::vector<T> > &sample) const;
   
   unsigned int npars()         const {return _ndiri;}

   void print(std::ostream &out = std::cout) const {out << "DiUn{" << this->npars() << "}";}
   void print_long(std::ostream &out = std::cout) const;
   friend std::ostream & operator << (std::ostream &out, const DiUn<T> &dirichlet)
    {
        dirichlet.print(out);
        return out;
    }

  private:
   unsigned int _ndiri;
   const unsigned int _npars;
};

template <typename T>
void DiUn<T>::print_long(std::ostream &out) const
{
  out << "DiUn{";
  for(unsigned int i = 0; i < this->npars() - 1; i++)
  {
    out << "1,";
  }
  out << "1}";
}

template <typename T>
void DiUn<T>::sample(unsigned int nsample, std::vector<std::vector<T> > &sample) const
{
  sample.clear();
  sample.resize(this->npars());
  double alpha[_ndiri];
  double diri[_ndiri];
  for(unsigned int i = 0; i < _ndiri; i++)alpha[i] = 1.L;

  for(unsigned int i = 0; i < nsample; i++)
  {
     gsl_ran_dirichlet(this->_r,_ndiri,alpha,diri);
     std::vector<T> diriv;
     for(unsigned int j = 0; j < this->npars(); j++)
     {
       sample[j].push_back(diri[j]);
     }
  }

  return;
}

}

#endif
