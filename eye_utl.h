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
    // ��������� ������ r2 �������, ������� ����� ������� � ���� � ����� ������
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

// ��������� ����� ���������� ��� ���/��� � ���� �����
void get_gauss_fhf( matrix& gauss_flt, int sigma );
void get_gauss_flf( matrix& gauss_flt, int sigma );

// ������� ��� ����������, � ������������ ������� ���� ���������� ������� gauss_flt
void freq_flt( matrix& dst, const matrix& gauss_flt );

// ���������� ����������� ���
void filt_hight_pass( const matrix& src, matrix& dst, int sigma = 20 );
// ���������� ����������� ���
void filt_low_pass( const matrix& src, matrix& dst, int sigma = 20 );

// ����������� �����������, ����� ��� ����������� ���������� ��� ����� �������� ��������
// �� down_rate �� up_rate �� ������� �����������
void bin_by_object_fraction( matrix& flt, int_matrix& binar, real down_rate, real up_rate );

// ����������� ����������� �� ���� �������, ����������
// �������������� ��� �������� ������� ��������������
void binar_hysteresis( matrix& image, int_matrix& bin, int val1, int val2 );

// ���������� �������������� ������� �� �������� ���� ������� 
// ��� ����� ��������� �� cell_num ������ � ���� ����������� ����� ��� ������� ���������
void find_shifts( int_matrix& first_eye, int_matrix& second_eye, int cell_num, int min_hro, int max_hro, int_vector& shifts );

// �� ������� �������� ������� ���� ��������
// ���������� ����� ������� �������.
int find_rotation( const int_vector& shifts, double* alpha );

// ��������������� ������� ������� ��������� ���������
int mark_clusters( int_matrix& bin );

// ��������� ������ ������ (������� �����) �� ��������� ������������ xyz,
// �������� ����� "��������" � ����� ������ glare
int find_laser( const matrix& xyz, const int_matrix& bin, const int_matrix &glare, int_matrix &mask,
                fpoint center, fpoint radius, bool is_big, int *x, int *y );

// ���������� ������ ������� ����� �� ������ � ������ (������ bayer)
// ���� �� ���������� ����� ������� �����,
// �� ��������� 1 � ������� � x, y ����� ���� ��������.
int get_red_target( bool is_bayerGR, uchar* ptr, int lx, int ly, int *x, int *y );

// ���������� ����� ����� ������ ����� � ��� ������� (��� ����������� �������� ������)
static void circle_correction( int_matrix& img, int r, int *x0, int *y0, int *r0 );

// ��������� ������������� ����� � ����� � ������� � x0, y0 � ������� r0
real check_uniformity( bool is_bayerGR, uchar* ptr, int lx, int ly, int x0, int y0, int r0 );

}; // namespace eye

#endif // _EYE_UTL_
