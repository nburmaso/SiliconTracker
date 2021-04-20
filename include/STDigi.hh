#ifndef S1Digi_h
#define S1Digi_h 1

static const int gnRowsX = 200;
static const int gnRowsY = 100;
static const int gnPadsPerLayer = gnRowsX * gnRowsY;

class STDigi
{
 public:
  STDigi(int channel);
  virtual ~STDigi() {}

  void SetTime(int t) { fTime = t; }
  void SetAmplitude(int adc) { fAmplitude = adc; }
  int GetChannel() { return fChannel; }
  int GetTime() { return fTime; }
  int GetAmplitude() { return fAmplitude; }
  int GetLayer() { return fChannel / gnPadsPerLayer; }
  int GetRowX() { return (fChannel - (fChannel / gnPadsPerLayer) * gnPadsPerLayer) / gnRowsY; }
  int GetRowY() { return fChannel % gnRowsY; }

 private:
  int fChannel;   // channel index
  int fTime;      // time in clock counts
  int fAmplitude; // amplitude in ADC counts
};

#endif
