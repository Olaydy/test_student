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

using namespace lwml;

// HOOKE AND othe person

class err_func_circle : public i_function {
public:
  err_func_circle( const int_matrix& img, int max_sh, int dr, int r0 )
  : _img(img), _r0(r0), _dr(dr), _max_sh(max_sh) 
  {
    // вычисляем радиус r2 области, которая будет входить в нашу в любом случае
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

// Построить маску фильтрации для ФВЧ/ФНЧ в виде Гауса
void get_gauss_fhf( matrix& gauss_flt, int sigma );
void get_gauss_flf( matrix& gauss_flt, int sigma );

// Функция для фильтрации, в спектральной области идет домножение намаску gauss_flt
void freq_flt( matrix& dst, const matrix& gauss_flt );

// фильтрация изображения ФВЧ
void filt_hight_pass( const matrix& src, matrix& dst, int sigma = 20 );
// фильтрация изображения ФНЧ
void filt_low_pass( const matrix& src, matrix& dst, int sigma = 20 );

// Бинаризация изображения, порог для бинаризации выбирается так чтобы пикселей осталось
// от down_rate до up_rate от размера изображения
void bin_by_object_fraction( matrix& flt, int_matrix& binar, real down_rate, real up_rate );

// бинаризация изображения по двум порогам, гистерезис
// предполагается что значения матрицы неотрицательны
void binar_hysteresis( matrix& image, int_matrix& bin, int val1, int val2 );

// Нахождение горизонтальных сдвигов по бинарной сети сосудов 
// Ось углов разбиваем на cell_num частей и ищем оптимальный сдвиг для каждого фрагмента
void find_shifts( int_matrix& first_eye, int_matrix& second_eye, int cell_num, int min_hro, int max_hro, int_vector& shifts );

// По вектору смещений находит угол поворота
// возвращает сумму модулей сдвигов.
int find_rotation( const int_vector& shifts, double* alpha );

// Вспомогательная функция пометки компонент связности
int mark_clusters( int_matrix& bin );

// Находение центра лазера (красной точки) по цветовому пространству xyz,
// бинарной маске "красного" и маске бликов glare
int find_laser( const matrix& xyz, const int_matrix& bin, const int_matrix &glare, int_matrix &mask,
                fpoint center, fpoint radius, bool is_big, int *x, int *y );

// Нахождение центра красной точки по буферу с камеры (формат bayer)
// Если не получилось найти красную точку,
// то возвращае 1 и заносит в x, y центр масс засветов.
int get_red_target( bool is_bayerGR, uchar* ptr, int lx, int ly, int *x, int *y );

// Нахождение более точно центра круга и его радиуса (для определения качества лазера)
static void circle_correction( int_matrix& img, int r, int *x0, int *y0, int *r0 );

// Анализует равномерность цвета в круге с центром в x0, y0 и радиуса r0
real check_uniformity( bool is_bayerGR, uchar* ptr, int lx, int ly, int x0, int y0, int r0 );

}; // namespace eye

#endif // _EYE_UTL_
