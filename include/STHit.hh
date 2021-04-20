#ifndef S1Hit_h
#define S1Hit_h 1

#include <vector>

class STHit
{
 public:
  STHit();
  virtual ~STHit();

  void AddDigiIndex(int index) { fDigiIndices.push_back(index); }
  int GetNDigiIndices() { return fDigiIndices.size(); }
  int GetDigiIndex(int i) { return fDigiIndices[i]; }
  double X() { return fX; }
  double Y() { return fY; }
  double Z() { return fZ; }
  double Dx() { return fDx; }
  double Dy() { return fDy; }
  double Dz() { return fDz; }
  int GetLayer() { return fLayer; }
  int GetAmplitude() { return fAmplitude; }

  void SetXYZ(double x, double y, double z);
  void SetDxDyDz(double dx, double dy, double dz);
  void SetLayer(int layer) { fLayer = layer; }
  void SetAmplitude(int adc) { fAmplitude = adc; }

 private:
  double fX;
  double fY;
  double fZ;
  double fDx;
  double fDy;
  double fDz;
  int fLayer;
  int fAmplitude;
  std::vector<int> fDigiIndices;
};

#endif
