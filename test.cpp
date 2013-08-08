#include <fstream>
#include <iostream>
#include <iomanip>

#include "sampler.h"

int main()
{
  Sampler::DirG<double> dirg;
  Sampler::Diri<double> diri;
  Sampler::DiUT<double> diut;
  Sampler::DiOr<double> dior;
  Sampler::DiUn<double> diun;
  Sampler::Norm<double> norm(5.,2.);
  Sampler::NorT<double> nort(5.,2.,1.5,10.8);
  Sampler::Unif<double> unif(5.,12.);
  Sampler::LogN<double> logn(5e7,0.5);
  Sampler::LogU<double> logu(5e7,8e15);


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

  dirg.set_parameters(pars);
  diri.set_parameters(pars);

  pars.clear();
  pars.push_back(min);
  pars.push_back(max);

  diut.set_parameters(pars);

  dior.set_parameters(ord);
  diun.set_parameters(means.size());

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

  unsigned int ns(1000);

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

  Sampler::DiUn<double> tree1(4),tree2(2),tree3(2),tree4(2),tree5(2);

  tree2.set_mother(&tree1,1);
  tree3.set_mother(&tree1,2);
  tree4.set_mother(&tree2,0);
  tree5.set_mother(&tree2,1);

  std::vector<std::vector<double> > reaction;
  tree1.sample_tree(ns,reaction);


  std::ofstream out("data.dat");
  out << "x y xg yg xt yt xo yo xu yu norm nort unif logn logu tr1 tr2 tr3 tr4 tr5 tr6 tr7 tr8" << std::endl;
  for(unsigned int i = 0; i < ns; i++)
  {
    out << -(sampleDiri[0][i]-1./3.)/std::sqrt(2.) + (sampleDiri[1][i]-1./3.)/std::sqrt(2.) << " "
        << -(sampleDiri[0][i]-1./3.)/std::sqrt(6.) - (sampleDiri[1][i]-1./3.)/std::sqrt(6.) + 2.*(sampleDiri[2][i]-1./3.)/std::sqrt(6.) << " ";
    out << -(sampleDirg[0][i]-1./3.)/std::sqrt(2.) + (sampleDirg[1][i]-1./3.)/std::sqrt(2.) << " "
        << -(sampleDirg[0][i]-1./3.)/std::sqrt(6.) - (sampleDirg[1][i]-1./3.)/std::sqrt(6.) + 2.*(sampleDirg[2][i]-1./3.)/std::sqrt(6.) << " ";
    out << -(sampleDiut[0][i]-1./3.)/std::sqrt(2.) + (sampleDiut[1][i]-1./3.)/std::sqrt(2.) << " "
        << -(sampleDiut[0][i]-1./3.)/std::sqrt(6.) - (sampleDiut[1][i]-1./3.)/std::sqrt(6.) + 2.*(sampleDiut[2][i]-1./3.)/std::sqrt(6.) << " ";
    out << -(sampleDior[0][i]-1./3.)/std::sqrt(2.) + (sampleDior[1][i]-1./3.)/std::sqrt(2.) << " "
        << -(sampleDior[0][i]-1./3.)/std::sqrt(6.) - (sampleDior[1][i]-1./3.)/std::sqrt(6.) + 2.*(sampleDior[2][i]-1./3.)/std::sqrt(6.) << " ";
    out << -(sampleDiun[0][i]-1./3.)/std::sqrt(2.) + (sampleDiun[1][i]-1./3.)/std::sqrt(2.) << " "
        << -(sampleDiun[0][i]-1./3.)/std::sqrt(6.) - (sampleDiun[1][i]-1./3.)/std::sqrt(6.) + 2.*(sampleDiun[2][i]-1./3.)/std::sqrt(6.) << " ";
    out << sampleNorm[i] << " "
        << sampleNort[i] << " "
        << sampleUnif[i] << " "
        << sampleLogn[i] << " "
        << sampleLogu[i] << " ";
    for(unsigned int j = 0; j < reaction.size(); j++)
    {
        out << reaction[j][i] << " ";
    }
    out << std::endl;
  }
  out.close();

  return 0;
}
