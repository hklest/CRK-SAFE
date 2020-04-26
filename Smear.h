#ifndef __SMEAR_H__
#define __SMEAR_H__

class Smear
{
public:
  Smear() {}
  virtual ~Smear() {}

  // Supply the initial cherenkov angle...get a new value in return...
  virtual double smr(double ThetaCherenkov, double lambda)=0;
};
	
#endif /* __SMEAR_H__ */
