#include "vn_lunar.h"
#include <Arduino.h>
#include <math.h>

uint32_t lunar_dd;
uint32_t lunar_mm;
uint32_t lunar_yy;
char*	CanChiNgay;
char*	CanChiThang;
char*	CanChiNam;

//=====================================================================================================//
//get_Julius
uint64_t vn_lunar::getJulius(uint32_t _dd, uint32_t _mm, uint32_t _yy) {
  uint32_t a, y, m, jd;
  a = floor((14 - _mm) / 12);
  y = _yy + 4800 - a;
  m = _mm + 12 * a - 3;
  jd = _dd + floor((153 * m + 2) / 5) + 365 * y + floor(y / 4) - floor(y / 100) + floor(y / 400) - 32045;
  if (jd < 2299161) {
    jd = _dd + floor((153 * m + 2) / 5) + 365 * y + floor(y / 4) - 32083;
  }
  return jd;
}
//=====================================================================================================//
//get_New_Moon_Day
uint32_t vn_lunar::getNewMoonDay(uint32_t k) {
  double T, T2, T3, dr, Jd1, M, Mpr, F, C1, deltat, JdNew;
  T = k / 1236.85;
  T2 = T * T;
  T3 = T2 * T;
  dr = PI / 180;
  double timeZone = 7.0;
  Jd1 = 2415020.75933 + 29.53058868 * k + 0.0001178 * T2 - 0.000000155 * T3;  // Mean new moon
  Jd1 = Jd1 + 0.00033 * sin((166.56 + 132.87 * T - 0.009173 * T2) * dr);  // Sun's mean anomaly
  M = 359.2242 + 29.10535608 * k - 0.0000333 * T2 - 0.00000347 * T3;  // Moon's mean anomaly
  Mpr = 306.0253 + 385.81691806 * k + 0.0107306 * T2 + 0.00001236 * T3;  // Moon's argument of latitude
  F = 21.2964 + 390.67050646 * k - 0.0016528 * T2 - 0.00000239 * T3;
  C1 = (0.1734 - 0.000393 * T) * sin(M * dr) + 0.0021 * sin(2 * dr * M);
  C1 = C1 - 0.4068 * sin(Mpr * dr) + 0.0161 * sin(dr * 2 * Mpr);
  C1 = C1 - 0.0004 * sin(dr * 3 * Mpr);
  C1 = C1 + 0.0104 * sin(dr * 2 * F) - 0.0051 * sin(dr * (M + Mpr));
  C1 = C1 - 0.0074 * sin(dr * (M - Mpr)) + 0.0004 * sin(dr * (2 * F + M));
  C1 = C1 - 0.0004 * sin(dr * (2 * F - M)) - 0.0006 * sin(dr * (2 * F + Mpr));
  C1 = C1 + 0.0010 * sin(dr * (2 * F - Mpr)) + 0.0005 * sin(dr * (2 * Mpr + M));
  if (T < -11) {
    deltat = 0.001 + 0.000839 * T + 0.0002261 * T2 - 0.00000845 * T3 - 0.000000081 * T * T3;
  } else {
    deltat = -0.000278 + 0.000265 * T + 0.000262 * T2;
  }
  JdNew = Jd1 + C1 - deltat;
  return floor(JdNew + 0.5 + timeZone / 24);
}
//=====================================================================================================//
//get_Sun_Longitude
uint32_t vn_lunar::getSunLongitude(uint32_t jdn) {
  double timeZone = 7.0;
  double T, T2, dr, M, L0, DL, L;
  // Time in Julian centuries from 2000-01-01 12:00:00 GMT
  T = (jdn - 2451545.5 - timeZone / 24) / 36525;
  T2 = T * T;   // degree to radian
  dr = PI / 180;  // mean anomaly, degree
  M = 357.52910 + 35999.05030 * T - 0.0001559 * T2 - 0.00000048 * T * T2;  // mean longitude, degree
  L0 = 280.46645 + 36000.76983 * T + 0.0003032 * T2;
  DL = (1.914600 - 0.004817 * T - 0.000014 * T2) * sin(dr * M);
  DL = DL + (0.019993 - 0.000101 * T) * sin(dr * 2 * M) + 0.000290 * sin(dr * 3 * M);
  L = L0 + DL; // true longitude, degree
  L = L * dr;  // Normalize to (0, 2*PI)
  L = L - PI * 2 * floor(L / (PI * 2));
  return floor(L / PI * 6);
}
//=====================================================================================================//
//get_Lunar_Month_11
uint32_t vn_lunar::getLunarMonthll(uint32_t yy) {
  double k, off, nm, sunLong;
  off = getJulius(31, 12, yy) - 2415021;
  k = floor(off / 29.530588853);
  nm = getNewMoonDay(floor(k));
  // sun longitude at local midnight
  sunLong = getSunLongitude(floor(nm));
  if (sunLong >= 9) {
    nm = getNewMoonDay(floor(k) - 1);
  }
  return floor(nm);
}
//=====================================================================================================//
//get_Leap_Month_Offset
uint32_t vn_lunar::getLeapMonthOffset(int64_t a11) {
  int64_t last, arc;
  uint32_t k, i;
  k = floor((a11 - 2415021.076998695) / 29.530588853 + 0.5);
  last = 0;
  // We start with the month following lunar month 11
  i = 1;
  arc = getSunLongitude(floor(getNewMoonDay(floor(k + i))));
  do {
    last = arc;
    i++;
    arc = getSunLongitude(floor(getNewMoonDay(floor(k + i))));
  } while (arc != last && i < 14);
  return i - 1;
}
//=====================================================================================================//
//convert_Solar_to_Lunar
uint32_t vn_lunar::convertSolar2Lunar(uint32_t _dd, uint32_t _mm, uint32_t _yy) {
  double dayNumber, monthStart, a11, b11, lunarDay, lunarMonth, lunarYear;
  //double lunarLeap;
  uint32_t k, diff;
  dayNumber = getJulius(_dd, _mm, _yy);
  k = floor((dayNumber - 2415021.076998695) / 29.530588853);
  monthStart = getNewMoonDay(k + 1);
  if (monthStart > dayNumber) {
    monthStart = getNewMoonDay(k);
  }
  a11 = getLunarMonthll(_yy);
  b11 = a11;
  if (a11 >= monthStart) {
    lunarYear = _yy;
    a11 = getLunarMonthll(_yy - 1);
  } else {
    lunarYear = _yy + 1;
    b11 = getLunarMonthll(_yy + 1);
  }
  lunarDay = dayNumber - monthStart + 1;
  diff = floor((monthStart - a11) / 29);
  //lunarLeap = 0;
  lunarMonth = diff + 11;
  if (b11 - a11 > 365) {
    uint32_t leapMonthDiff = getLeapMonthOffset(a11);
    if (diff >= leapMonthDiff) {
      lunarMonth = diff + 10;
      if (diff == leapMonthDiff) {
        //lunarLeap = 1;
      }
    }
  }
  if (lunarMonth > 12) {
    lunarMonth = lunarMonth - 12;
  }
  if (lunarMonth >= 11 && diff < 4) {
    lunarYear -= 1;
  }

  lunar_dd = (uint32_t)(lunarDay);
  lunar_mm = (uint32_t)(lunarMonth);
  lunar_yy = (uint32_t)(lunarYear);

//  String lunar_day = String(floorlunar_dd);
//  String lunar_mon = String(floorlunar_mm);
//  String lunar_year = String(floorlunar_yy);
//  String str_lunar = String(lunar_day + "/" + lunar_mon + "/" + lunar_year);

}

char * Chi(uint32_t y)
{
  char * Chi;
  switch(y)
  {
    case 0:
    {
    	Chi= "Ti";
		  break;
    }
    case 1:
    {
      Chi= "Suu";
	    break;
    }
    case 2:
    {
      Chi= "Dan";
	    break;
    }
    case 3:
    {
      Chi= "Mao";
	    break;
    }
    case 4:
    {
      Chi= "Thin";
	    break;
    }
    case 5:
    {
      Chi= "Ty";
	    break;
    }
    case 6:
    {
      Chi= "Ngo";
	    break;
    }
    case 7:
    {
      Chi= "Mui";
	    break;
    }
    case 8:
    {
      Chi= "Than";
	    break;
    }
    case 9:
    {
      Chi= "Dau";
	    break;
    }
    case 10:
    {
      Chi= "Tuat";
	    break;
    }
    default:
    {
      Chi= "Hoi";
	    break;
    }
  }
  return Chi;
}


char * Can(uint32_t x)
{
  char * Can;
  switch(x)
  {
    case 0:
    {
      Can= "Giap";
	    break;
    }
    case 1:
    {
      Can= "At";
	    break;
    }
    case 2:
    {
      Can= "Binh";
	    break;
    }
    case 3:
    {
      Can= "Dinh";
	    break;
    }
    case 4:
    {
      Can= "Mau";
	    break;
    }
    case 5:
    {
      Can= "Ky";
	    break;
    }
    case 6:
    {
      Can= "Canh";
	    break;
    }
    case 7:
    {
      Can= "Tan";
	    break;
    }
    case 8:
    {
      Can= "Nham";
	    break;
    }
    default:
    {
     Can= "Quy";
	   break;
    }
  }
  return Can;
}

char * CharsConnect(char * a, char * b)
{
  size_t newlen = strlen(a) + strlen(b);
  char *r = malloc(newlen + 2);
  strcpy(r, a);
  strcat(r, " ");
  strcat(r, b);
  return r;
}

uint32_t vn_lunar::convert2CanChi(uint32_t dd, uint32_t mm, uint32_t yy, uint32_t lunar_dd, uint32_t lunar_mm, uint32_t lunar_yy)
{

	// Tìm can chi ngày
	char * CAN_Ngay; 
	char * CHI_Ngay; 
	char * CAN_CHI_Ngay;
  uint64_t JD = this->getJulius(dd, mm, yy);
	uint32_t X = (uint32_t)(JD + 9.5) % 10;					// Tìm Can ngày
	CAN_Ngay = Can(X);

	uint32_t Y = (uint32_t)(JD + 1.5) % 12;				// Tìm Chi ngày
	CHI_Ngay = Chi(Y);

	CAN_CHI_Ngay = CharsConnect(CAN_Ngay, CHI_Ngay);					// Ra được can-chi theo ngày.

	// Tìm can chi tháng
	char * CAN_Thang;
	char * CHI_Thang;
	char * CAN_CHI_Thang;
	uint32_t month = lunar_mm;				
	uint32_t year = lunar_yy;	

	// Tìm can chi năm
	char * CAN_Nam;
	char * CHI_Nam;
	char * CAN_CHI_Nam;

	uint32_t T = year % 10;											// Tìm Can năm
	if (T <= 3)
		CAN_Nam = Can(T + 6);
	else 
		CAN_Nam = Can(T - 4);

	uint32_t M = year % 60;											// Tìm Chi năm
	while (M > 11)
	{
		M = M - 12;
	}
	if (M <= 3)
		CHI_Nam = Chi(M + 8);
	else 
		CHI_Nam = Chi(M - 4);

	uint32_t Z;														// Tìm Can tháng
	if(CAN_Nam == "Giap" ||CAN_Nam == "Ky" ){
		 Z = (month + 1)%10	;
	}
	if(CAN_Nam == "At" ||CAN_Nam == "Canh" ){
		 Z = (month + 3)%10	;
	}
	if(CAN_Nam == "Binh" ||CAN_Nam == "Tan" ){
		 Z = (month + 5)%10	;
	}
	if(CAN_Nam == "Dinh" ||CAN_Nam == "Nham" ){
		 Z = (month + 7)%10	;
	}
	if(CAN_Nam == "Mau" ||CAN_Nam == "Quy" ){
		 Z = (month + 9)%10	;
	}
									
	CAN_Thang = Can(Z);

	if (month >= 11)											// Tìm Chi tháng
		CHI_Thang = Chi(month - 11);
	else
		CHI_Thang = Chi(month+1);

	CAN_CHI_Thang = CharsConnect(CAN_Thang, CHI_Thang);				// Ra được can-chi theo tháng.
	CAN_CHI_Nam = CharsConnect(CAN_Nam, CHI_Nam);


	// Chuyển kết quả
	CanChiNgay = CAN_CHI_Ngay;
	CanChiThang = CAN_CHI_Thang;
	CanChiNam = CAN_CHI_Nam;
}

uint32_t vn_lunar::get_lunar_dd() {
  return lunar_dd;
}
uint32_t vn_lunar::get_lunar_mm() {
  return lunar_mm;
}
uint32_t vn_lunar::get_lunar_yy() {
  return lunar_yy;
}
char* vn_lunar::get_canchi_ngay() {
  return CanChiNgay;
}
char* vn_lunar::get_canchi_thang() {
  return CanChiThang;
}
char* vn_lunar::get_canchi_nam() {
  return CanChiNam;
}
