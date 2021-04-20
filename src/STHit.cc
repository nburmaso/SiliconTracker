#include "STHit.hh"

STHit::STHit()
  : fX(),
    fY(),
    fZ(),
    fDx(),
    fDy(),
    fDz(),
    fLayer(),
    fAmplitude(),
    fDigiIndices()
{
}

STHit::~STHit()
{
  fDigiIndices.clear();
}

void STHit::SetXYZ(double x, double y, double z)
{
  fX = x;
  fY = y;
  fZ = z;
}

void STHit::SetDxDyDz(double dx, double dy, double dz)
{
  fDx = dx;
  fDy = dy;
  fDz = dz;
}
