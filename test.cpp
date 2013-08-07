#include <fstream>
#include <iostream>
#include <iomanip>

#include "sampler.h"

int main()
{
  Sampler::Distribution<Sampler::DirG<double>,double> dirg;
  Sampler::Distribution<Sampler::Diri<double>,double> diri;
  Sampler::Distribution<Sampler::DiUT<double>,double> diut;
  Sampler::Distribution<Sampler::DiOr<double>,double> dior;
  Sampler::Distribution<Sampler::DiUn<double>,double> diun;
  Sampler::Distribution<Sampler::Norm<double>,double> norm;
  Sampler::Distribution<Sampler::NorT<double>,double> nort;
  Sampler::Distribution<Sampler::Unif<double>,double> unif;
  Sampler::Distribution<Sampler::LogN<double>,double> logn;
  Sampler::Distribution<Sampler::LogU<double>,double> logu;



  std::vector<double> means,sigs,min,max,ord;
  means.push_back(0.3); sigs.push_back(0.30);
  means.push_back(0.5); sigs.push_back(0.04);
  means.push_back(0.2); sigs.push_back(0.01);
  min.push_back(0.00); max.push_back(1.00);
  min.push_back(0.38); max.push_back(0.62);
  min.push_back(0.17); max.push_back(0.23);
  ord.push_back(2);
  ord.push_back(1);
  ord.push_back(3);

  std::vector<std::vector<double> > pars;
  pars.push_back(means);
  pars.push_back(sigs);
  std::vector<double> parsNorm,parsNort,parsLogn,parsLogu,parsUnif;
  parsNorm.push_back(5.); parsNorm.push_back(2.);
  parsNort.push_back(5.); parsNort.push_back(2.);parsNort.push_back(1.5);parsNort.push_back(10.8);
  parsLogn.push_back(5e7);parsLogn.push_back(0.5);
  parsLogu.push_back(5e7);parsLogu.push_back(8e15);
  parsUnif.push_back(5.);parsUnif.push_back(12.);


  dirg.set_distribution(pars);
  diri.set_distribution(pars);

  pars.clear();
  pars.push_back(min);
  pars.push_back(max);

  diut.set_distribution(pars);

  dior.set_distribution(ord);
  diun.set_distribution(means.size());

  norm.set_distribution(parsNorm);
  nort.set_distribution(parsNort);
  unif.set_distribution(parsUnif);
  logn.set_distribution(parsLogn);
  logu.set_distribution(parsLogu);

  std::vector<std::vector<double> > sampleDiri;
  std::vector<std::vector<double> > sampleDirg;
  std::vector<std::vector<double> > sampleDiut;
  std::vector<std::vector<double> > sampleDior;
  std::vector<std::vector<double> > sampleDiun;
  std::vector<double> sampleNorm;
  std::vector<double> sampleNort;
  std::vector<double> sampleLogn;
  std::vector<double> sampleLogu;
  std::vector<double> sampleUnif;

  unsigned int ns(5000);

  diri.sample(ns,sampleDiri);
  dirg.sample(ns,sampleDirg);
  diut.sample(ns,sampleDiut);
  dior.sample(ns,sampleDior);
  diun.sample(ns,sampleDiun);
  norm.sample(ns,sampleNorm);
  nort.sample(ns,sampleNort);
  unif.sample(ns,sampleUnif);
  logn.sample(ns,sampleLogn);
  logu.sample(ns,sampleLogu);

  std::ofstream out("data.dat");
  out << "x y xg yg xt yt xo yo xu yu norm nort unif logn logu" << std::endl;
  for(unsigned int i = 0; i < ns; i++)
  {
    out << -(sampleDiri[i][0]-1./3.)/std::sqrt(2.) + (sampleDiri[i][1]-1./3.)/std::sqrt(2.) << " "
        << -(sampleDiri[i][0]-1./3.)/std::sqrt(6.) - (sampleDiri[i][1]-1./3.)/std::sqrt(6.) + 2.*(sampleDiri[i][2]-1./3.)/std::sqrt(6.) << " ";
    out << -(sampleDirg[i][0]-1./3.)/std::sqrt(2.) + (sampleDirg[i][1]-1./3.)/std::sqrt(2.) << " "
        << -(sampleDirg[i][0]-1./3.)/std::sqrt(6.) - (sampleDirg[i][1]-1./3.)/std::sqrt(6.) + 2.*(sampleDirg[i][2]-1./3.)/std::sqrt(6.) << " ";
    out << -(sampleDiut[i][0]-1./3.)/std::sqrt(2.) + (sampleDiut[i][1]-1./3.)/std::sqrt(2.) << " "
        << -(sampleDiut[i][0]-1./3.)/std::sqrt(6.) - (sampleDiut[i][1]-1./3.)/std::sqrt(6.) + 2.*(sampleDiut[i][2]-1./3.)/std::sqrt(6.) << " ";
    out << -(sampleDior[i][0]-1./3.)/std::sqrt(2.) + (sampleDior[i][1]-1./3.)/std::sqrt(2.) << " "
        << -(sampleDior[i][0]-1./3.)/std::sqrt(6.) - (sampleDior[i][1]-1./3.)/std::sqrt(6.) + 2.*(sampleDior[i][2]-1./3.)/std::sqrt(6.) << " ";
    out << -(sampleDiun[i][0]-1./3.)/std::sqrt(2.) + (sampleDiun[i][1]-1./3.)/std::sqrt(2.) << " "
        << -(sampleDiun[i][0]-1./3.)/std::sqrt(6.) - (sampleDiun[i][1]-1./3.)/std::sqrt(6.) + 2.*(sampleDiun[i][2]-1./3.)/std::sqrt(6.) << " ";
    out << sampleNorm[i] << " "
        << sampleNort[i] << " "
        << sampleUnif[i] << " "
        << sampleLogn[i] << " "
        << sampleLogu[i];
    out << std::endl;
  }
  out.close();

  return 0;
}
