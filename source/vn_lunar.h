#ifndef vn_lunar_h
#define vn_lunar_h
#include <Arduino.h>

class vn_lunar
{
    public:
      uint64_t getJulius(uint32_t dd, uint32_t mm, uint32_t yy);
      uint32_t getNewMoonDay(uint32_t k);
      uint32_t getSunLongitude(uint32_t jdn);
      uint32_t getLunarMonthll(uint32_t yy);
      uint32_t getLeapMonthOffset(int64_t a11);
      uint32_t convertSolar2Lunar(uint32_t dd, uint32_t mm, uint32_t yy);
      uint32_t convert2CanChi(uint32_t dd, uint32_t mm, uint32_t yy, uint32_t lunar_dd, uint32_t lunar_mm, uint32_t lunar_yy);
      uint32_t get_lunar_dd();
      uint32_t get_lunar_mm();
      uint32_t get_lunar_yy();
      char* get_canchi_ngay();
      char* get_canchi_thang();
      char* get_canchi_nam();
    private:
      uint32_t _dd;
      uint32_t _mm;
      uint32_t _yy;
      uint32_t k;
      uint32_t jdn;
      uint32_t a11;
};

#endif