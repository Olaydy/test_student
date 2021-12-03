#include "eye_utl.h"

#include "lwml.h"

#include "vbmp.h"
#include "filename.h"
#include "stream.h"
#include "morfolog.h"

#include "geom.h"
#include "igeom.h"
#include "fft2d.h"
#include "trig_gen.h"
#include "medianfilt.h"
#include "fft.h"
#include "filler.h"

#include "gauss_blur.h"
#include "morfolog.h"
#include "hooke_opt.h"

#define SIZE 512

namespace eye {

using namespace lwml;

void get_gauss_fhf( matrix& gauss_flt, int sigma )
{
  for( int y = 0; y < gauss_flt.str(); y++ ){
    for( int x = 0; x < gauss_flt.col(); x++ ){
      gauss_flt(y, x) = gauss(0) * gauss(0) - gauss(x / sigma) * gauss(y / sigma);
    }
  }
}

void get_gauss_flf( matrix& gauss_flt, int sigma )
{
  for( int y = 0; y < gauss_flt.str(); y++ ){
    for( int x = 0; x < gauss_flt.col(); x++ ){
       gauss_flt(y, x) = gauss(x / sigma) * gauss(y / sigma) /( gauss(0) * gauss(0));
    }
  }
}

void binrevers_cols( matrix& dst )
{
  int i, j, k;
  int n = dst.str();
  for( j = i = 0; i < n - 1; i++ ){
    if( i < j ){
      for( int col = 0; col < dst.col(); col++ )
        t_swap(dst(i, col), dst(j, col));
    }
    for( k = n / 2; k <= j; k /= 2 )
      j -= k;
    j += k;
  }
}

void scale_cols( matrix& dst )
{
  int n = dst.str();
  for( int col = 0; col < dst.col(); col++ ){
    for( int i = 0; i < n; i++ ){
      dst.at(i, col) /= n;
    }
  }
}

void resort_cols( matrix& dst, matrix& im )
{
  int n = dst.str();
  for( int col = 0; col < dst.col(); col++ ){
    for( int i = 1; i < n / 2; i++ ){
      t_swap(dst(i, col), dst(n-i, col));
      t_swap(im(i, col), im(n-i, col));
    }
  }
}

void freq_flt( matrix& dst, const matrix& gauss_flt )
{
  // вычисляем спектр, в результате близке к нулю значения, находятся в углах изображения
  
  int ly = dst.str();
  int lx = dst.col();

  int ns = dst.str();
  int nc = dst.col();
  matrix im(ns, nc, 0);
  
  // обработка столбцов
  real ur, ui, wr, wi, tr, ti, ur2;
  int i, j, l, le1, le2, ip;
  int rs = intlog2(ns);

  // fft for columns
  scale_cols(dst);
  binrevers_cols(dst);
  for( le2 = l = 1; l <= rs; l++ ){
    le1 = le2;
    le2 *= 2;
    ur = 1.0;
    ui = 0.0;
    wr = cos(M_PI / le1);
    wi = -sin(M_PI / le1);
    for( j = 0; j < le1; j++ ){
      for( i = j; i < ns; i += le2 ){
        ip = i + le1;
        for( int c = 0; c < nc; c++ ){
          tr = dst.at(ip, c) * ur - im.at(ip, c) * ui;
          ti = dst.at(ip, c) * ui + im.at(ip, c) * ur;
          dst.at(ip, c) = dst.at(i, c);
          im.at(ip, c) = im.at(i, c);
          dst.at(i, c) += tr;
          im.at(i, c) += ti;
          dst.at(ip, c) -= tr;
          im.at(ip, c) -= ti;
        }
      }
      ur2 = ur * wr - ui * wi;
      ui = ur * wi + ui * wr;
      ur = ur2;
    }
  }

  // fft for strings
  vector rsv(nc);
  vector isv(nc);
  for( int s = 0; s < ns; s++ ){
    dst.get_str(rsv, s);
    im.get_str(isv, s);
    fft::cfft(rsv, isv);
    dst.put_str(rsv, s);
    im.put_str(isv, s);
  }

  // Умножаем на маску фильтра
  for( int y = 0; y < ly/2; y++ ){
    for( int x = 0; x < lx/2; x++ ){
      dst(y, x) *= gauss_flt(y, x);
      im(y, x) *= gauss_flt(y, x);
      dst(ly-1-y, x) *= gauss_flt(y, x);
      im(ly-1-y, x) *= gauss_flt(y, x);
      dst(y, lx-1-x) *= gauss_flt(y, x);
      im(y, lx-1-x) *= gauss_flt(y, x);
      dst(ly-1-y, lx-1-x) *= gauss_flt(y, x);
      im(ly-1-y, lx-1-x) *= gauss_flt(y, x);
    }
  }
  
  // обратное БПФ
  //fft2d::cifft(dst, im);
  for( int s = 0; s < ns; s++ ){
    dst.get_str(rsv, s);
    im.get_str(isv, s);
    fft::cifft(rsv, isv);
    dst.put_str(rsv, s);
    im.put_str(isv, s);
  }

  // ifft for columns
  resort_cols(dst, im);
  binrevers_cols(dst);
  binrevers_cols(im);
  for( le2 = l = 1; l <= rs; l++ ){
    le1 = le2;
    le2 *= 2;
    ur = 1.0;
    ui = 0.0;
    wr = cos(M_PI / le1);
    wi = -sin(M_PI / le1);
    for( j = 0; j < le1; j++ ){
      for( i = j; i < ns; i += le2 ){
        ip = i + le1;
        for( int c = 0; c < nc; c++ ){
          tr = dst.at(ip, c) * ur - im.at(ip, c) * ui;
          ti = dst.at(ip, c) * ui + im.at(ip, c) * ur;
          dst.at(ip, c) = dst.at(i, c);
          im.at(ip, c) = im.at(i, c);
          dst.at(i, c) += tr;
          im.at(i, c) += ti;
          dst.at(ip, c) -= tr;
          im.at(ip, c) -= ti;
        }
      }
      ur2 = ur * wr - ui * wi;
      ui = ur * wi + ui * wr;
      ur = ur2;
    }
  }
}

void filt_hight_pass( const matrix& src, matrix& dst, int sigma )
{
  int ly = dst.str();
  int lx = dst.col();
  matrix gauss_flt(ly/2, lx/2, 0);
  get_gauss_fhf(gauss_flt, sigma);
  
  for( int y = 0; y < ly; y++ ){
    for( int x = 0; x < lx; x++ ){
      dst(y, x) = src(y, x);
    }
  }
  freq_flt(dst, gauss_flt);
}

void filt_low_pass( const matrix& src, matrix& dst, int sigma )
{
  int ly = dst.str();
  int lx = dst.col();
  matrix gauss_flt(ly/2, lx/2, 0);
  get_gauss_flf(gauss_flt, sigma);
  
  for( int y = 0; y < ly; y++ ){
    for( int x = 0; x < lx; x++ ){
      dst(y, x) = src(y, x);
    }
  }
  freq_flt(dst, gauss_flt);
}

void bin_by_object_fraction( matrix& flt, int_matrix& binar, real down_rate, real up_rate )
{
  const int hist_len = 100;
  int_vector hist(hist_len, 0);
  
  int lx = flt.col(), ly = flt.str();
  double min = flt.min();
  double max = flt.max();
  double k = (hist_len-1) / (max - min);
  for( int y = 0; y < ly; y++ ){
    for( int x = 0; x < lx; x++ ){
      flt.at(y, x) = (flt.at(y, x) - min)*k;
      // best thing for time improving
      hist[(int_cast)(flt.at(y, x))]++;
    }
  }
//  if( zzz_dump() )
//    vbmp::save(zzz_dump_name("flt.bmp").ascstr(), flt);

  // Уровни бинаризации выбираем так, чтобы было достаточно "сосудов" (1/15 < num < 1/6)
  int th1 = 0;
  int th2 = 0;
  int sum = 0;
  for( int k = 0; k < hist_len; k++ ){
    sum += hist[k];
    if( sum < lx*ly * down_rate )
      th1 = k;
    if( sum > lx*ly * up_rate ){
      th2 = k;
      break;
    }
  }
  //printf("lev1=%d lev2=%d\n", th1, th2);
  // Бинаризуем изображение
  binar_hysteresis(flt, binar, th1, th2);
}

// предполагается что значения матрицы неотрицательны
void binar_hysteresis( matrix& image, int_matrix& bin, int val1, int val2 )
{
  int lx = image.col();
  int ly = image.str();

  for( int x = 0; x < lx; x++ ){
    for( int y = 0; y < ly; y++ ){
      // пометка -2, соответствует объекту, -1: не объект, точка просмотрена
      if( image.at(y, x) < val1 && image.at(y, x) >= 0 ){
        image.at(y, x) = -2;
      }
      else if( image.at(y, x) > val2 )
        image.at(y, x) = -1;
      else{
         bool fl = false; // есть ли соседи > val2
        // Смотрим соседей
        for( int k = -1; !fl && k < 2; k++ ){
          for( int l = -1; !fl && l < 2; l++ ){
            if( k && l ){
              int yy = y + k;
              int xx = x + l;

              if( !( xx >= 0 && yy >= 0 && yy < ly && xx < lx ) )
                continue;

              if( !fl && image.at(yy, xx) > val2 ){
                image.at(yy, xx) = -1;
                fl = true;
              }
            }
          }
        }
        if( !fl )
          image.at(y, x) = -2;
      }
    }
  }

  for( int x = 0; x < image.col(); x++ ){
    for( int y = 0; y < image.str(); y++ ){
      if( image.at(y, x) == -2 )
        bin.at(y, x) = 1;
      else
        bin.at(y, x) = 0;
    }
  }
  //if( zzz_dump() ){
  //  vbmp::save(zzz_dump_name("bin_ves.bmp").ascstr(), bin);
  //}
}

// Нахождение горизонтальных сдвигов по бинарной сети сосудов 
// Ось углов разбиваем на cell_num частей и ищем оптимальный сдвиг для каждого фрагмента
void find_shifts( int_matrix& first_eye, int_matrix& second_eye, int cell_num, int min_hro, int max_hro, int_vector& shifts )
{ 
  int phi_len = first_eye.col();
  int sz = phi_len/cell_num;
  int max_sh = 20;
  int_vector err_func(max_sh*2+1);
  int_vector sh_vec(cell_num);
  
  // Бежим по всем ячейкам. В которых ищем сдвиг
  for( int k = 0; k < cell_num; k++ ){
    err_func.set_zero();
    // Бежим по ячейке
    for( int i = min_hro; i < max_hro; i++ ){
      for( int j = k*sz; j < (k+1) * sz; j++ ){
        // Смотрим все возможные сдвиги
        for( int shift = -max_sh; shift <= max_sh; shift++ ){
            //printf("%d %d    ", j, (phi_len+j + shift)%phi_len);
          if( ( first_eye(i, (phi_len+j + shift)%phi_len) == 0 && second_eye(i, j) == 0 ) ||
              ( first_eye(i, (phi_len+j + shift)%phi_len) == 1 && second_eye(i, j) == 1 ) ){
            err_func[shift+max_sh] += 1;
          }
          else
            err_func[shift+max_sh] -= 5;
        }
      }
    }
   
    int opt = err_func.max_idx();
    shifts[k] = opt - max_sh;
  }
}

// возвращает сумму модулей сдвигов.
int find_rotation( const int_vector& shifts, double* alpha )
{
  int cell_num = shifts.len();
  int min_rot_cells = cell_num* 3/4;
  int flag = 1; // 0 - means it's not rotation
  // Ищем последовательность из отрицательных сдвигов 
  int start = -1;
  int end = -1;
  // Сначала ищем отрицательный сдвиг
  for( int k = 0; k < cell_num; k++ ){
    if( shifts[k] < 0 ){
      start = k;
      break;
    }
  }
  // Не нашли отрицательных сдвигов
  if( start == -1 ){
    // Ищем поворот против часовой, из неотрицательных сдвигов
    int sum = 0;
    int abs_sum = 0;
    int num = 0;
    for( int k = 0; k < cell_num; k++ ){
      // Считаем сумму модулей сдвигов
      abs_sum += abs(shifts[k]);
      if( shifts[k] > 0 && flag != 0 ){
        sum += shifts[k];
        num++;
      }
      // Если встретили отрицательный сдвиг, то не поворот
      else if( shifts[k] < 0){ // !! Возможно есть исключения, когда ротация + закатывание
        flag = 0;
      }
      
    }
    if( num >= min_rot_cells )
      *alpha = 1.0*sum/cell_num;
    else
      *alpha = 0.0;
    return abs_sum;
  }
  
  int count = 1;
  // Теперь ищем последовательность отрицательных сдвигов
  for( int k = start+1; k < cell_num; k++ ){
    if( shifts[k] >= 0){
      end = k-1;
      break;
    }
    else
      count++;
  }
  
  // Если надо, то ищем с конца (т.к замыкание по кругу)
  if( start == 0 && end != -1 ){
    for( int k = cell_num-1; k > 1; k-- ){
      if( shifts[k] < 0){
        count++;
        start = k;
      }
      else
        break;
    }
  }
  
  // проверяем на поворот, все сдвиги должны быть одного знака, кроме может быть 2х нулей
  // поворот по часовой, отрицательные сдвиги
  if( count >= min_rot_cells ){
    int sum = 0;
    int abs_sum = 0;
    for( int k = 0; k < cell_num; k++ ){
      abs_sum += abs(shifts[k]);
      if( shifts[k] <= 0 && flag != 0 )
        sum += shifts[k];
      else{// Если положительный сдвиг, то не поворот
        flag = 0;
      }
    }
    *alpha = 1.0*sum/cell_num;
    return abs_sum;
  }
  // Если отрицательных половина, то это может соответствовать общему сдвигу в плоскости
  // Если отрицательных 4-7 это может быть закатывание, надо проверять
  int abs_sum = 0;
  for( int k = 0; k < cell_num; k++ ){
    abs_sum += abs(shifts[k]);
  }
  *alpha = 0.0;
  // Отдельно обрабатываем случай, когда один отрицательный сдвиг -1,
  // считаем этот сдвиг шумовым, и вычисляем величину поворота
  if( count == 1 && shifts[start] == -1){
    // Делаем вид что на месте -1 был 0, т.е. сумму модулей уменьшаем на 2
    abs_sum -= 2; 
    *alpha = 1.0*abs_sum/cell_num;
  }
  
  return abs_sum;
}

int mark_clusters( int_matrix& bin )
{
  // помечаем компоненты связности индексами начиная с 1
  int cnum = 0;
  for( int i = 0; i < bin.str(); i++ ){
    for( int j = 0; j < bin.col(); j++ ){
      if( bin(i, j) == -1 )
        filler::mark_4n(bin, j, i, ++cnum);
    }
  }
  return cnum;
}

int find_laser( const matrix& xyz, const int_matrix& bin, const int_matrix &glare, int_matrix &mask,
                fpoint center, fpoint radius, bool is_big, int *x, int *y )
{
  // помечаем компоненты связности индексами начиная с 1
  // пометки начинаем из пикселей bin == -1, но распространяем и на точки засветов
  // т.к. красная точка лазера может содержать блик в центре
  if( zzz_dump() ){
    vbmp::save(zzz_dump_name("bin__.bmp").ascstr(), bin);
    vbmp::save(zzz_dump_name("glare__.bmp").ascstr(), mask);
    vbmp::save(zzz_dump_name("glare2__.bmp").ascstr(), glare);
  }

  int cnum = 0;
  for( int i = 0; i < bin.str(); i++ ){
    for( int j = 0; j < bin.col(); j++ ){
      if( bin(i, j) == -1 )
        filler::mark_8n(mask, j, i, ++cnum);
    }
  }
  if( zzz_dump() )
    vbmp::save(zzz_dump_name("filler.bmp").ascstr(), mask);

  // нет ни одной компонетны - лазер не нашли
  if( cnum == 0 )
    return -1;
    
  // строим список найденных компонент и их площади, по засчетам, красному и xyz
  t_list<int> list;
  t_list<int> sz;
  t_list<int> sz_red;
  t_list<int> sz_xyz;
  int min_sz = bin.str()*bin.str();
  for( int i = 0; i < bin.str(); i++ ){
    for( int j = 0; j < bin.col(); j++ ){
      // крайние кластера делаем нереально большими, чтоб не подошли
      if( i == 0 || i == bin.str()-1 || j == 0 || j == bin.col()-1 ){
        int k;
        for( k = 0; k < list.len(); k++ )
          if( list[k] == mask(i, j) )
            break;
        if( k != list.len() )
          sz[k] = min_sz;
      }      
      if( mask(i, j) == 0 ) 
        mask(i, j) = -1;
      else if( mask(i, j) > 0 ){
        int k;
        for( k = 0; k < list.len(); k++ )
          if( list[k] == mask(i, j) )
            break;
        if( k == list.len() ){
          list.put(mask(i, j));
          sz.put(0);
          sz_red.put(0);
          sz_xyz.put(0);
          if( bin(i, j) == -1 )
            sz_red[k]++;
          if( xyz(i, j) != 0 )
            sz_xyz[k]++;          
          if( glare(i, j) == -1 )
            sz[k]++;
        }
        else{
           if( glare(i, j) == -1 )
            sz[k]++;
          if( bin(i, j) == -1 )
            sz_red[k]++;
            if( xyz(i, j) != 0 )
            sz_xyz[k]++;          
        }
        mask(i, j) = k+1;
      }
    }
  }
 
  // отбираем подходящие кластера, у которых есть красное, xyz  и смотрим соотношение красного и засветов 
  t_list<int> good_cl;
  for( int k = 0; k < list.len(); k++ ){    
    if( sz_red[k] == 0 || sz_xyz[k] == 0 )
      continue;
    //printf("! k=%d cl=%d sz=%d sz_red=%d sz_xyz=%d f=%.3f  ", k, list[k], sz[k], sz_red[k], sz_xyz[k], 1.0*sz[k] / sz_red[k]);    
    if( is_big && 1.0*sz[k] / sz_red[k] > 0.7 )
      continue;
    if( is_big && 1.0*sz[k] / sz_red[k] < 0.1 && sz_xyz[k] < 3 )
      continue;
    // если ищем красное, без засвета внутри, то требование более жесткое
    if( !is_big && 1.0*sz[k] / sz_red[k] > 0.9 )
      continue;
    good_cl.put(k);
    //printf("\n put %d   len=%d", list[k], good_cl.len());
  }
  // не нашли подходящий
  if( good_cl.len() == 0 )
    return -1;
      
  int_matrix mask1(bin.str(), bin.col(), 0);
  int_vector center_x(good_cl.len(), 0);
  int_vector center_y(good_cl.len(), 0);
  int_vector counts(good_cl.len(), 0);
  for( int i = 0; i < bin.str(); i++ ){
    for( int j = 0; j < bin.col(); j++ ){
      for( int k = 0; k < good_cl.len(); k++ ){
        if( mask(i, j) == good_cl[k]+1 ){
          mask1(i, j) = good_cl[k]+1;
          center_y[k] += i;
          center_x[k] += j;
          counts[k]++;
        }
      }
    }
  }
  
  //for( int k = 0; k < good_cl.len(); k++ ){
  //  printf("%d: %d %d %d\n  ", good_cl[k]+1, center_x[k], center_y[k], counts[k]);
  //}

  if( zzz_dump() )
    vbmp::save(zzz_dump_name("res1.bmp").ascstr(), mask1);

  // Если нашли несколько подходящих кластеров, то вычисляем среднее
  if( good_cl.len() >= 2 ){
    int in = 0;
    int idx_in = -1;
    for( int k = 0; k < good_cl.len(); k++ ){
      if( abs(center_x[k] / counts[k] - center.x()) < radius.x() && abs(center_y[k] / counts[k] - center.y()) < radius.y() ){  
        in++;
        idx_in = k;
      }
    }    
    // Нашли единственный подходящий кластер в районе центра засветов
    if( in == 1 ){
      *x = center_x[idx_in] / counts[idx_in];
      *y = center_y[idx_in] / counts[idx_in];
      return 0;
    }
    // нашли несколько внутри
    else if( in > 1 ){
        return -1;
    }
    // все снаружи, то считаем среднее арифметическое
    *x = 0;
    *y = 0;
    for( int k =0; k < good_cl.len(); k++ ){
      *x += center_x[k] / counts[k];
      *y += center_y[k] / counts[k];
    }
    *x /= good_cl.len();
    *y /= good_cl.len();   
  }
  
  // Если нашли один подходящий кластер, то проверяем его на близость к центру засветов
  if( good_cl.len() == 1 ){
    *x = center_x[0] / counts[0];
    *y = center_y[0] / counts[0];
  }
  
  return 0;
}

int get_red_target( bool is_bayerGR, uchar* ptr, int lx, int ly, int *x, int *y )
{
  int size = 100;
  int shift_x = lx/2 - size;
  int shift_y = ly/2 - size;
  
  // Нахождение матрицы красного
  int_matrix bin(2*size, 2*size, 0);
  int_matrix glare(2*size, 2*size, 0);
  matrix xyz(2*size, 2*size, 0);
  
  double min = 10000;
  int r, g, b;
  // бежим в координатах вырезанной картинки
  // формируем матрицу засветов, пометок красного и цветового пространства xyz
  // Сразу ищем центр масс засветов
  double str = 0, col = 0;
  int count = 0;  
  for( int yy = 0; yy < 2*size; yy++ ){
    for( int xx = 0; xx < 2*size; xx++ ){
      // переходим к координатам исходной картинки
      int x = xx + shift_x;
      // переворачиваем y чтобы 0 был сверху.
      int y = ly - 1 - (yy + shift_y);
      int xPos = x - (x - 1) % 2;
      int yPos = (y - (y - 1) % 2) * lx;
      // сразу инвертируем значения
      if( is_bayerGR ){
        r = ptr[xPos + yPos - lx];
        g = (ptr[xPos + yPos] + ptr[xPos - 1 + yPos - lx]) /2;
        b = ptr[xPos - 1  + yPos];
      }
      else{
        r = ptr[xPos + yPos];
        g = (ptr[xPos + yPos - lx] + ptr[xPos - 1 + yPos]) /2;
        b = ptr[xPos - 1 + yPos - lx];
      }      
      rgbcolor rgb(1.0/255*r, 1.0/255*g, 1.0/255*b);
      hsbcolor hsb = color_conv::rgb2hsb(rgb);
      
      // формируем матрицу засветов      
      if( hsb.brightness() > 0.3 && hsb.saturation() <= 0.2 )
        glare(yy, xx) = -1;      
      if( hsb.brightness() > 0.5 && hsb.saturation() > 0.2 && hsb.saturation() < 0.98 && (hsb.hue() * 360.0 > 20 && hsb.hue() * 360.0 < 310) )
        glare(yy, xx) = -1;
      
      // вычисление центра масс засветов
      if( glare(yy, xx) == -1 ){        
        count++;
        str += yy;
        col += xx;
      }
        
      // делаем пороговую бинаризацию по каналу красного и насыщенности 
      if( hsb.saturation() > 0.6 && r > 70 && (hsb.hue() * 360.0 < 30 || hsb.hue() * 360.0 > 300) )
        bin(yy, xx) = -1;
      
      // переходим в цветовое пространство xyz, берем Y (http://unick-soft.ru/article.php?id=32)
      xyz(yy, xx) = -0.158 * r + 0.252 * g - 0.003 * b;
      // сразу считаем min по Y, не на засветах
      if( glare(yy, xx) != 0 )
        continue;
      if( xyz(yy, xx) < min )
        min = xyz(yy, xx);
    }
  }
  // Находим центр масс засветов, в дальнейшем искать около него     
  if( count == 0 )
    return -1;
  
  double xx_center = col / count;
  double yy_center = str / count;  
  // также находим дисперсию
  str = 0, col = 0;
  for( int i = 0; i < glare.str(); i++ ){
    for( int j = 0; j < glare.col(); j++ ){
      if( glare(i, j) == -1 ){
        str += (i-yy_center)*(i-yy_center);
        col += (j-xx_center)*(j-xx_center);
      }
    }
  }
  
  // вычисляем максимальную дисперсию 
  double disp_x = sqrt(col / count);
  double disp_y = sqrt(str / count);
  double disp_xx = disp_x;
  double disp_yy = disp_y;
  double disp_frac = disp_y / disp_x;
  if( disp_frac < 0.8 )
    disp_yy = disp_x + (disp_x - disp_y) / 2;
    
  if( disp_y > disp_x ){
    disp_frac = disp_x / disp_y;
    if( disp_frac < 0.8 )
      disp_xx = disp_y + (disp_y - disp_x) / 2;     
  }
  
  //printf("GLARE center %.2f %.2f dispersion %.2f %.2f (new disp %.2f %.2f) disp_frac=%.2f  ", xx_center, yy_center, disp_x, disp_y, disp_xx, disp_yy, disp_frac);
  if( disp_frac < 0.8 || disp_x < 15 || disp_y < 15 )
    return -1;
  for( int i = 0; i < glare.str(); i++ ){
    for( int j = 0; j < glare.col(); j++ ){
      // отсеиваем точки xyz по порогу min*0.5
      if( xyz(i, j) > min*0.5 )
        xyz(i, j) = 0;  
      if( abs(j - xx_center) > disp_xx*0.65 || abs(i - yy_center) > disp_yy*0.65 )
        bin(i, j) = 0;
    }
  }
  
  if( zzz_dump() ){
    vbmp::save(zzz_dump_name("glare.bmp").ascstr(), glare);
    vbmp::save(zzz_dump_name("bin_first.bmp").ascstr(), bin);
    vbmp::save(zzz_dump_name("xyz.bmp").ascstr(), xyz);
  }

  int_matrix bin2(2*size, 2*size);
  int_matrix glare2(2*size, 2*size);
  bin2.copy(bin);  
  
  // гасим красные точки рядом с засветами
  int_matrix mask_dil(2*size, 2*size);
  int_matrix mask_dil2(2*size, 2*size);
  dilatation <primitive_manhattan>(glare, mask_dil2);
  glare2.copy(mask_dil2);
  dilatation <primitive_manhattan>(mask_dil2, mask_dil);
  dilatation <primitive_manhattan>(mask_dil, mask_dil2);
  dilatation <primitive_manhattan>(mask_dil2, mask_dil);
  for( int yy = 0; yy < 2*size; yy++ ){
    for( int xx = 0; xx < 2*size; xx++ ){
      if( mask_dil(yy, xx) == -1 )
        bin2(yy, xx) = 0;
    }
  }
  
  // расширяем красное и поглощаем засветы
  erosion <primitive_infty>(bin2, mask_dil);
  dilatation <primitive_manhattan>(mask_dil, mask_dil2);
  dilatation <primitive_manhattan>(mask_dil2, mask_dil);
  for( int yy = 0; yy < 2*size; yy++ ){
    for( int xx = 0; xx < 2*size; xx++ ){
      if( mask_dil(yy, xx) == -1 )
        glare(yy, xx) = -1;
    }
  }
  
  int_matrix glare_copy(2*size, 2*size);
  glare_copy.copy(glare);
    
  int_matrix mask_dil_copy(2*size, 2*size);
  mask_dil_copy.copy(mask_dil);
  
  // если не нашли красного лазера без засветов, то будем искать красное с малым засветом
    
  // гасим засветы рядом с красными точками
  closing <primitive_infty>(bin, mask_dil2);
  //erosion <primitive_infty>(mask_dil, mask_dil2);
  //dilatation <primitive_manhattan>(mask_dil, mask_dil2);
  dilatation <primitive_manhattan>(mask_dil2, mask_dil);
  if( zzz_dump() ){   
    vbmp::save(zzz_dump_name("bin_preproc.bmp").ascstr(), mask_dil);
    vbmp::save(zzz_dump_name("bin__after.bmp").ascstr(), bin);
  }  
  // расширяем красное и поглощаем засветы
  for( int yy = 0; yy < 2*size; yy++ ){
    for( int xx = 0; xx < 2*size; xx++ ){
      if( mask_dil(yy, xx) == -1 )
        glare_copy(yy, xx) = -1;     
    }
  }      
    
  // Нахождение лазера
  int res = find_laser(xyz, mask_dil, glare2, glare_copy, fpoint(xx_center, yy_center), fpoint(disp_xx*0.65, disp_yy*0.65), true, x, y);
  // Если найденный центр близок к центру засветов, то выдаем ответ, что нашли
  if( res != -1 && abs(*x - xx_center) < disp_xx*0.65 && abs(*y - yy_center) < disp_yy*0.65 ){
    *x += shift_x;
    *y += shift_y;
    return res;  
  }
    
  // Нахождение лазера
  res = find_laser(xyz, mask_dil_copy, glare2, glare, fpoint(xx_center, yy_center), fpoint(disp_xx, disp_yy), false, x, y);
  // Если найденный центр близок к центру засветов, то выдаем ответ, что нашли
  if( res != -1 && abs(*x - xx_center) < disp_xx*0.65 && abs(*y - yy_center) < disp_yy*0.65 ){  
    *x += shift_x;
    *y += shift_y;
    return res;
  }
  
  // если разброс дисперсии для x y отличается значительно, то центр засветов находим плохо 
  if( disp_frac < 0.8 )
    return -1;
  // возвращаем 1 если не нашли красную точку, но нашли примерный центр бликов  
  else{
    *x = shift_x + xx_center;
    *y = shift_y + yy_center;
    return 1;
  }
}

// Нахождение более точно центра круга и его радиуса (для определения качества лазера)
static void circle_correction( int_matrix& img, int r, int *x0, int *y0, int *r0 )
{
  int lx = img.col();
  int ly = img.str();

  int_matrix bin(ly, lx);
  // строим бинарную картинку, на которой более точно будем находить центр и радиус
  for( int y = 0; y < ly; y++ ){
    for( int x = 0; x < lx; x++ ){
      if( img(y, x) < 80 )
        bin(y, x) = -2;
      else
        bin(y, x) = 1;        
    }
  }

  r -= 0.15*r;
  
  int dr = 10;
  int max_sh = 10;
  err_func_circle hook_circle(bin, max_sh, dr, r);
  int value_old;
  int value = hooke_opt::calc(hook_circle, lx/2, ly/2, r, 16, 16, 4, 0.5, 10, x0, y0, r0, &value_old);
  printf("\n!!! %d, %d, %d  -> %d, %d, %d \n", lx/2, ly/2, r, *x0, *y0, *r0);
  if(zzz_dump()){
    vbmp::save(zzz_dump_name("bin.bmp").ascstr(), bin);
  }
}

double local_entropy( int_matrix& img, matrix& dst, int_matrix& glare, int r )
{
  int lx = img.col();
  int ly = img.str();

  int area_sz = 9;
  int area_r = area_sz / 2;
  const int max = 81; //area_sz*area_sz
  double entroy_table[max+1];
  const float log2 = log(2.0f);
  entroy_table[0] = 0.0;
  float frequency = 0;
  for (int i = 1; i < max+1; i++){
    frequency = (float)i / max;
    entroy_table[i] = frequency * (log(frequency) / log2);
  }

  int x0, y0, r0;
  circle_correction(img, r, &x0, &y0, &r0);
  
  // отступаем от края из-за неровностей
  r0 -= 10;
  matrix diff(ly, lx, 0);
    
  // для вычисления энтропии строим каждый раз гистограмму и по ней вычисляем энтропию  
  // также строим матрицу diff, разность переженых и недо жженых точек
  int hist_count[256] = {0};
  int flag; // попал ли в область засвет?
  for( int y = area_sz; y < ly - area_sz; y++ ){
    for( int x = area_sz; x < lx - area_sz; x++ ){
      flag = 0;
      double rr = (x-x0)*(x-x0) + (y-y0)*(y-y0);      
      // смотрим только точки внутри окружности
      if( rr > (r0 - area_r)*(r0 - area_r) )
        continue;     
        
      // для каждой точки смотрим окрестность, строим гистограмму встречаемости значения яркости 
      for( int sh_y = -area_r; sh_y <= area_r; sh_y++ ){
        for( int sh_x = -area_r; sh_x <= area_r; sh_x++ ){
          if( glare(y + sh_y, x + sh_x) == -1 )
            flag = 1;
          hist_count[img(y + sh_y, x + sh_x)]++;         
        }
      }

      // по гистограмме считаем энтропию, и моду
      double e = 0;
      int moda = 0;
      for( int k = 0; k < 256; k++ ){
        if( hist_count[k] != 0 ){
          // ищем самую часто встречающуюся яркость
          if( hist_count[moda] < hist_count[k] )
            moda = k;
          //if( flag == 0 && hist_count[k] <= 10 )
          //  printf("%d %d (back =%d)  (%d %d %d) moda=%d \n", k, hist_count[k], hist_count[40], x, y, img(y, x), moda);          
          e -= entroy_table[hist_count[k]];     
        }
      }
      // Ищем как много пикселей и с какой яркостью "пережженых"
      double down = 0.0, up = 0.0;
      int down_num = 0, up_num = 0;
      for( int k = 0; k < moda; k++ ){
        if( hist_count[k] != 0 ){
          down_num += hist_count[k];
          down += hist_count[k] * (moda - k);
          hist_count[k] = 0;
        }
      }
      
      // не забываем обнулить в гистограмме значение с модой
      hist_count[moda] = 0;
      
      // Ищем как много пикселей и с какой яркостью "недо жженых"
      for( int k = moda+1; k < 256; k++ ){
        if( hist_count[k] != 0 ){
          up_num += hist_count[k];
          up += hist_count[k] * (k - moda);
          hist_count[k] = 0;
        }
      }
      
      // вычисляем долю якрости недо и пере жженых пикселей
      if( down_num != 0 )
        down /= down_num;
      else
        down = 0;
      if( up_num != 0 )
        up /= up_num;
      else 
        up = 0;  
      //if( down != 0 || up != 0 )  
        //printf(" moda=%d: %f %d %f %d \n", moda, down, down_num, up, up_num);
      
      if( flag == 1 )
        dst(y, x)= -1*entroy_table[max];
      else{
        dst(y, x)= e;
        diff(y, x) = up - down;
      }
    }
  }
  
  printf("lx=%d ly=%d   ", lx, ly);
  if(zzz_dump()){
    vbmp::save(zzz_dump_name("res.bmp").ascstr(), dst);    
    vbmp::save(zzz_dump_name("diff.bmp").ascstr(), diff);    
  }
  
  // Вычисляем количество ошибок лазера и их среднюю яркость
  double wrongs = 0.0;
  int wrongs_num = 0;
  for( int y = area_sz; y < ly - area_sz; y++ ){
    for( int x = area_sz; x < lx - area_sz; x++ ){  
      if( diff(y, x) != 0 ){
        wrongs += diff(y, x);
        wrongs_num++;
      }
    }
  }
  // wrongs_num/81 - количество плохих пикселей, в среднем яркость ошибки wrongs/wrongs_num
  printf(" !! %f %f %d %d  %f\n", wrongs, wrongs/81, wrongs_num, wrongs_num/81, wrongs/wrongs_num);
  
  // Оценка равномерности ошибок
  // Сглаживаем сильно изображение, затем уменьшаем размер до 9 на 9
  // Затем смотрим пространственное распределение ошибок
  
  // Вырезаем из изображения объемлющий квадрат, содержащий лишь нашу окружность, отступив от края
  printf("\n\n %d, %d, x=%d y=%d r=%d (%d %d, %d %d)\n",  lx, img.col(), x0, y0, 2*r0, y0-r0, y0+r0, x0-r0, x0+r0);
  matrix cut(2*r0, 2*r0, 0);
  cut.copy_sub(dst, y0-r0, x0-r0);
  if(zzz_dump()){
    vbmp::save(zzz_dump_name("cut.bmp").ascstr(), cut);
  }
  
  gauss_blur bl(2*r0, 2*r0, 0);
  // После уменьшения картинки, на радиус попадут 4 точкки, r0/4
  // Но сглаживаем еще посильнее, чтобы было влияние на соседние точки
  // !! Вероятно не нужно столь сильно учитывать яркость точек при размытии, а сделать их более равнозначными
  int sigma = r0/3;
  bl.proc(cut, sigma, sigma);
  bl.get(cut);
  if(zzz_dump()){
    vbmp::save(zzz_dump_name("res_blur.bmp").ascstr(), cut);
  }
  
  // Умньшаем картинку до размера 9 на 9, запоминаем min max
  double scale = 2*r0/8;
  matrix scaled(2*r0/scale+1, 2*r0/scale + 1, 0);
  double min_val = cut(0, 0), max_val = cut(0, 0);
  for( int y = 0; y < 2*r0; y+=scale ){
    for( int x = 0; x < 2*r0; x+=scale ){
      scaled(y/scale, x/scale) = cut(y, x);
      if( cut(y, x) > max_val )
        max_val = cut(y, x);
      if( cut(y, x) < min_val )
        min_val = cut(y, x);
    }
  }
  printf("scaled :scale=%f %d %d (%f, %f) min/max %f %f\n", scale, lx, ly, ly/scale+1, lx/scale + 1, min_val, max_val);
  if(zzz_dump())
    vbmp::save(zzz_dump_name("scale.bmp").ascstr(), scaled);
  
  
  // Считаем энтропию получившейся картинки 9 на 9
  printf("\n");
  for( int y = 0; y < 9; y++ ){
    for( int x = 0; x < 9; x++ ){
      int val = 5*(scaled(y, x) - min_val)/max_val;
      printf("%3d  ", val);
      hist_count[val]++;
    }
    printf("\n");
  }
  double e = 0;
  for( int k = 0; k < 256; k++ ){
    if( hist_count[k] != 0 )
      e -= entroy_table[hist_count[k]];
  }
  printf("\n      entr = %f\n", e);

  
  // вычисляем среднее значение в области 3 на 3 для 9 точек, чтобы учеть пространственную связь
  double points[9] = {0};
  double sum = 0;
  for( int k = 0; k < 9; k++ ){
    int xx = 3*(k%3) + 1;
    int yy = 3*(k/3) + 1;
    for( int y = -1; y <= 1; y++ ){
      for( int x = -1; x <= 1; x++ ){
        int val = 30*(scaled(yy + y, xx + x) - min_val)/max_val;
        points[k] += val;
        //printf("%2d  ", val);
      }
    }
    //printf(" !!%d(%d) %d(%d) %.2f  ", xx, k%3, yy, k/3, points[k]/9);
    printf(" %.2f  ", points[k]/9);
    sum += points[k]/9;
    if( k%3 == 2 )
      printf("\n");
  }
  
  // считаем среднее дисперсию, и симметричность
  double mid = sum / 9;
  sum = 0;
  double m3 = 0.0;
  for( int k = 0; k < 9; k++ ){
    points[k] /= 9;
    sum += (points[k] - mid) * (points[k] - mid);
    m3 += (points[k] - mid) * (points[k] - mid) * (points[k] - mid);
  }
  
  sum /= 9;
  m3 /= 9;
  printf("mid=%.2f d=%.2f %.2f,   %.2f, %.2f \n", mid, sum, sqrt(sum), m3, m3/sqrt(sum));  
}

real check_uniformity( bool is_bayerGR, uchar* ptr, int lx, int ly, int x0, int y0, int r0 )
{
  int center_x = x0;
  int center_y = y0;
  
  int shift_x = center_x - r0;
  int shift_y = center_y - r0;
  int size_x = r0 + 9;
  int size_y = r0 + 9;
  int_matrix glare(2*size_y, 2*size_x, 0);
  matrix xyz(2*size_y, 2*size_x);
 
  printf("\n%d %d shift %d %d\n", lx, ly, shift_x, shift_y);
  
  double max = INT_MIN;
  double min = INT_MAX;
  double r, g, b;
  // бежим в координатах вырезанной картинки
  for( int yy = 0; yy < 2*size_y; yy++ ){
    for( int xx = 0; xx < 2*size_x; xx++ ){
      // переходим к координатам исходной картинки
      int x = xx + shift_x;
      // переворачиваем y чтобы 0 был сверху.
      int y = (ly - 1 - (yy + shift_y));
          
      int xPos = x - (x - 1) % 2;
      int yPos = (y - (y - 1) % 2) * lx;
      // сразу инвертируем значения
      if( is_bayerGR ){
        r = 255 - ptr[xPos + yPos - lx];
        g = 255 - (ptr[xPos + yPos] + ptr[xPos - 1 + yPos - lx]) /2;
        b = 255 - ptr[xPos - 1  + yPos];
      }
      else{
        r = 255 - ptr[xPos + yPos];
        g = 255 - (ptr[xPos + yPos - lx] + ptr[xPos - 1 + yPos]) /2;
        b = 255 - ptr[xPos - 1 + yPos - lx];
      }
        
      rgbcolor rgb(1.0/255*r, 1.0/255*g, 1.0/255*b);
      hsbcolor hsb = color_conv::rgb2hsb(rgb);
      // формируем матрицу засветов
      if( hsb.brightness() < 0.75 && hsb.saturation() > 0.7 )
        glare(yy, xx) = -1;
        //glare(yy, xx) = hsb.brightness()*255;
      
      // переходим в цветовое пространство xyz, берем Y (http://unick-soft.ru/article.php?id=32)
      xyz(yy, xx) = (-0.158 * r + 0.252 * g - 0.003 * b);      
      
      if( min > xyz(yy, xx) )
        min = xyz(yy, xx);
      if( max < xyz(yy, xx) )
        max = xyz(yy, xx);
    }
  }
    
  // Нормируем значения в диапазон 0-255
  int_matrix gray(2*size_y, 2*size_x);
  for( int yy = 0; yy < 2*size_y; yy++ ){
    for( int xx = 0; xx < 2*size_x; xx++ )
      gray(yy, xx) = (xyz(yy, xx) - min) *255 / (max - min);
  }
  
  int_matrix tmp(2*size_y, 2*size_x, 0);
  
  erosion<primitive_manhattan>(glare, tmp);
  dilatation<primitive_manhattan>(tmp, glare);
  dilatation<primitive_manhattan>(glare, tmp);
  dilatation<primitive_manhattan>(tmp, glare);
  dilatation<primitive_manhattan>(glare, tmp);
  dilatation<primitive_manhattan>(tmp, glare);
  dilatation<primitive_manhattan>(glare, tmp);
  dilatation<primitive_manhattan>(tmp, glare);
  dilatation<primitive_manhattan>(glare, tmp);
  dilatation<primitive_manhattan>(tmp, glare);  
    
  if(zzz_dump()){
    vbmp::save(zzz_dump_name("xyz.bmp").ascstr(), xyz);
    vbmp::save(zzz_dump_name("glare.bmp").ascstr(), glare);
  }
  
  matrix entr(2*size_y, 2*size_x, 0.0);
  local_entropy(gray, entr, glare, r0);

  return 0;
}

}; // namespace eye
