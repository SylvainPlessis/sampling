#ifndef _DIRICHLET_ORDONNE_PDF_
#define _DIRICHLET_ORDONNE_PDF_

//
#include "sampler.h"

namespace Sampler{
//!dior distribution
template <typename T = double>
class DiOr:public Distribution<T>
{
  public:
    DiOr():_npars(2){}
    DiOr(const std::vector<T> &pars):_npars(2) {this->set_parameters(pars);}
    DiOr(const DiOr<T> &rhs):_npars(4){if(this == &rhs)return;_par = rhs.index_pars();this->_r=rhs.random_gsl();}
    ~DiOr(){}

   void set_parameters(const std::vector<T>  &parameters);

   void sample(unsigned int nsample,std::vector<std::vector<T> > &sample) const;
   
   unsigned int npars()              const {return _par.size();}
   const std::vector<T> index_pars() const {return _par;}

   void print(std::ostream &out = std::cout);
   friend std::ostream & operator << (std::ostream &out, const DiOr<T> &dirichlet)
    {
        dirichlet.print(out);
        return out;
    }

  private:
   std::vector<T> _par;
   const unsigned int _npars;
   unsigned int offset;
};

template <typename T>
void DiOr<T>::print(std::ostream &out)
{
  out << "DiOr{";
  for(unsigned int i = 0; i < this->npars() - 1; i++)
  {
       out << _par[i] << ",";
  }
  out << _par.back() << "}";
}

template <typename T>
void DiOr<T>::set_parameters(const std::vector<T>  &parameters)
{
  _par = parameters;
  offset = _par.size();
  for(unsigned int i = 0; i < this->npars(); i++)
  {
     if(_par[i] < offset)offset = _par[i];
  }
}

template <typename T>
void DiOr<T>::sample(unsigned int nsample, std::vector<std::vector<T> > &sample) const
{
  sample.clear();
  sample.resize(this->npars());
//dior -> sorted diun
  double alpha[_par.size()];
  double diri[_par.size()];
  for(unsigned int i = 0; i < this->npars(); i++)alpha[i] = 1.L;

  for(unsigned int ns = 0; ns < nsample; ns++)
  {
     gsl_ran_dirichlet(this->_r,this->npars(),alpha,diri);
     gsl_sort(diri,1,_par.size());
     std::vector<T> diriv;
     for(unsigned int j = 0; j < this->npars(); j++)
     {
       sample[j].push_back(diri[_par.size() - 1 - (unsigned int)_par[j] + offset]); //reverse order
     }
  }

  return;
}

}

#endif
