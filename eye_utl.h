#ifndef _EYE_UTL_
#define _EYE_UTL_

#include "imatrix.h"
#include "matrix.h"
#include "defs.h"

#include "geom.h"
#include "vbmp.h"
#include "lwml.h"

#include "hooke_opt.h"

namespace eye {

// nlsdf
using namespace lwml;

// HOOKE AND JEEWS

class err_func_circle : public i_function {
public:
  err_func_circle( const int_matrix& img, int max_sh, int dr, int r0 )
  : _img(img), _r0(r0), _dr(dr), _max_sh(max_sh) 
  {
    // âû÷èñëÿåì ðàäèóñ r2 îáëàñòè, êîòîðàÿ áóäåò âõîäèòü â íàøó â ëþáîì ñëó÷àå
    int r2 =  _r0 - 2*_max_sh - _dr;
    if( r2 < 0 )
      r2 = 0;
  }

  virtual int func( int x, int y, int r ) const;

private:
  int_matrix _img;
  int _r0, _dr;
  int _max_sh;
};

// Ïîñòðîèòü ìàñêó ôèëüòðàöèè äëÿ ÔÂ×/ÔÍ× â âèäå Ãàóñà
void get_gauss_fhf( matrix& gauss_flt, int sigma );
void get_gauss_flf( matrix& gauss_flt, int sigma );

// Ôóíêöèÿ äëÿ ôèëüòðàöèè, â ñïåêòðàëüíîé îáëàñòè èäåò äîìíîæåíèå íàìàñêó gauss_flt
void freq_flt( matrix& dst, const matrix& gauss_flt );

// ôèëüòðàöèÿ èçîáðàæåíèÿ ÔÂ×
void filt_hight_pass( const matrix& src, matrix& dst, int sigma = 20 );
// ôèëüòðàöèÿ èçîáðàæåíèÿ ÔÍ×
void filt_low_pass( const matrix& src, matrix& dst, int sigma = 20 );

// Áèíàðèçàöèÿ èçîáðàæåíèÿ, ïîðîã äëÿ áèíàðèçàöèè âûáèðàåòñÿ òàê ÷òîáû ïèêñåëåé îñòàëîñü
// îò down_rate äî up_rate îò ðàçìåðà èçîáðàæåíèÿ
void bin_by_object_fraction( matrix& flt, int_matrix& binar, real down_rate, real up_rate );

// áèíàðèçàöèÿ èçîáðàæåíèÿ ïî äâóì ïîðîãàì, ãèñòåðåçèñ
// ïðåäïîëàãàåòñÿ ÷òî çíà÷åíèÿ ìàòðèöû íåîòðèöàòåëüíû
void binar_hysteresis( matrix& image, int_matrix& bin, int val1, int val2 );

// Íàõîæäåíèå ãîðèçîíòàëüíûõ ñäâèãîâ ïî áèíàðíîé ñåòè ñîñóäîâ 
// Îñü óãëîâ ðàçáèâàåì íà cell_num ÷àñòåé è èùåì îïòèìàëüíûé ñäâèã äëÿ êàæäîãî ôðàãìåíòà
void find_shifts( int_matrix& first_eye, int_matrix& second_eye, int cell_num, int min_hro, int max_hro, int_vector& shifts );

// Ïî âåêòîðó ñìåùåíèé íàõîäèò óãîë ïîâîðîòà
// âîçâðàùàåò ñóììó ìîäóëåé ñäâèãîâ.
int find_rotation( const int_vector& shifts, double* alpha );

// Âñïîìîãàòåëüíàÿ ôóíêöèÿ ïîìåòêè êîìïîíåíò ñâÿçíîñòè
int mark_clusters( int_matrix& bin );

// Íàõîäåíèå öåíòðà ëàçåðà (êðàñíîé òî÷êè) ïî öâåòîâîìó ïðîñòðàíñòâó xyz,
// áèíàðíîé ìàñêå "êðàñíîãî" è ìàñêå áëèêîâ glare
int find_laser( const matrix& xyz, const int_matrix& bin, const int_matrix &glare, int_matrix &mask,
                fpoint center, fpoint radius, bool is_big, int *x, int *y );

// Íàõîæäåíèå öåíòðà êðàñíîé òî÷êè ïî áóôåðó ñ êàìåðû (ôîðìàò bayer)
// Åñëè íå ïîëó÷èëîñü íàéòè êðàñíóþ òî÷êó,
// òî âîçâðàùàå 1 è çàíîñèò â x, y öåíòð ìàññ çàñâåòîâ.
int get_red_target( bool is_bayerGR, uchar* ptr, int lx, int ly, int *x, int *y );

// Íàõîæäåíèå áîëåå òî÷íî öåíòðà êðóãà è åãî ðàäèóñà (äëÿ îïðåäåëåíèÿ êà÷åñòâà ëàçåðà)
static void circle_correction( int_matrix& img, int r, int *x0, int *y0, int *r0 );

// Àíàëèçóåò ðàâíîìåðíîñòü öâåòà â êðóãå ñ öåíòðîì â x0, y0 è ðàäèóñà r0
real check_uniformity( bool is_bayerGR, uchar* ptr, int lx, int ly, int x0, int y0, int r0 );

}; // namespace eye

#endif // _EYE_UTL_
