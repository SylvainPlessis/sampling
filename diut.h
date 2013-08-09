#ifndef _DIRICHLET_UNIFORM_TRONCATED_PDF_
#define _DIRICHLET_UNIFORM_TRONCATED_PDF_

//
#include "sampler.h"

namespace Sampler{
//!diut distribution
template <typename T = double>
class DiUT:public Distribution<T>
{
  public:
    DiUT():_npars(2){}
    DiUT(const std::vector<std::vector<T> > &pars):_npars(2) {this->set_parameters(pars);}
    DiUT(const DiUT<T> &rhs):_npars(4){if(this == &rhs)return;_min = rhs.min();_max=rhs.max();this->_r=rhs.random_gsl();}
    ~DiUT(){}

   void set_parameters(const std::vector<std::vector<T> > &parameters);

   void sample(unsigned int nsample,std::vector<std::vector<T> > &sample) const;
   
   unsigned int npars()       const {return _min.size();}
   const std::vector<T> min() const {return _min;}
   const std::vector<T> max() const {return _max;}

   void print(std::ostream &out = std::cout);
   friend std::ostream & operator << (std::ostream &out, const DiUT<T> &dirichlet)
    {
        dirichlet.print(out);
        return out;
    }

  private:
   std::vector<T> _min;
   std::vector<T> _max;
   const unsigned int _npars;

   bool accept(const double *run) const;
};

template <typename T>
void DiUT<T>::print(std::ostream &out)
{
  out << "DiUT{";
  for(unsigned int i = 0; i < this->npars() - 1; i++)
  {
        out << "[" << _min[i] << "," << _max[i] << "],";
  }
  out << "[" << _min.back() << "," << _max.back() << "]}";
}

template <typename T>
bool DiUT<T>::accept(const double *run) const
{
  for(unsigned int i = 0; i < _min.size(); i++)
  {
     if(_min[i] > run[i] || _max[i] < run[i])return false;
  }

  return true;
}

template <typename T>
void DiUT<T>::set_parameters(const std::vector<std::vector<T> > &parameters)
{
  _min = parameters[0];
  _max = parameters[1];
}

template <typename T>
void DiUT<T>::sample(unsigned int nsample, std::vector<std::vector<T> > &sample) const
{
  sample.clear();
  sample.resize(this->npars());
//diut -> troncated diun
  double alpha[_min.size()];
  double diri[_min.size()];
  for(unsigned int i = 0; i < _min.size(); i++)alpha[i] = 1.L;

  unsigned int ns(0);
  while(ns != nsample)
  {
     gsl_ran_dirichlet(this->_r,_min.size(),alpha,diri);
     if(!accept(diri))continue;
     ns++;
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
