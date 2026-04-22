# -*- coding: utf-8 -*-
# cython: boundscheck=False, wraparound=False, cdivision=True
# distutils: define_macros=NPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION
"""
bedcmmPitch class(Cython)

Author: YASUHARA Wataru
Copyright (c) 2026, Feel a Piece of the World
"""
import numpy as np
cimport numpy as cnp
cimport cython

ctypedef cnp.float64_t DTYPE_d_t
ctypedef cnp.int32_t DTYPE_i_t

# Cレベルの関数を利用するために math をインポート
from libc.math cimport fmin,log,INFINITY,fabs

cdef double EPS = 1e-300

cdef double _parabolic_peak(double[:] y, Py_ssize_t i):
    cdef double ym1, y0, yp1, denom
    if i <= 0 or i >= y.shape[0] - 1:
        return 0.0
    ym1 = y[i - 1]
    y0 = y[i]
    yp1 = y[i + 1]
    denom = ym1 - 2.0 * y0 + yp1
    if fabs(denom) < EPS:
        return 0.0
    return 0.5 * (ym1 - yp1) / denom


cdef double _gaussian_peak(double[:] y, Py_ssize_t i):
    cdef double ym1, y0, yp1, denom
    if i <= 0 or i >= y.shape[0] - 1:
        return 0.0
    ym1 = log(max(y[i - 1], EPS))
    y0 = log(max(y[i], EPS))
    yp1 = log(max(y[i + 1], EPS))
    denom = ym1 - 2.0 * y0 + yp1
    if fabs(denom) < 1e-300:
        return 0.0
    return 0.5 * (ym1 - yp1) / denom

cdef double _centroid_peak(double[:] y, Py_ssize_t i, int half_window=1):
    cdef Py_ssize_t lo, hi, n, j
    cdef double s = 0.0
    cdef double sw = 0.0
    cdef double w, x

    n = y.shape[0]
    lo = i - half_window
    if lo < 0:
        lo = 0
    hi = i + half_window
    if hi >= n:
        hi = n - 1

    for j in range(lo, hi + 1):
        w = y[j]
        if w > 0.0:
            s += w
            sw += (<double>j) * w

    if s <= 0.0:
        return 0.0

    x = sw / s
    return x - (<double>i)


cdef cnp.ndarray[DTYPE_d_t, ndim=1] _periodicity_1d_core_cy(double[:] data, Py_ssize_t[:] periods):
    """
    1次元データの周期性解析コアロジック (Cython版)
    """
    cdef Py_ssize_t n_data = data.shape[0]
    cdef Py_ssize_t n_periods = periods.shape[0]
    cdef cnp.ndarray[DTYPE_d_t, ndim=1] result = np.zeros(n_periods, dtype=np.float64)
    
    cdef Py_ssize_t p_idx, i, lag
    cdef double temp_sum
    cdef Py_ssize_t count
    cdef double mean_val = 0

    for p_idx in range(n_periods):
        lag = periods[p_idx]
        
        temp_sum = 0
        count = n_data - lag
        
        # 内側のループを純粋なCで回す
        # スライシングを使わず、インデックス i と i+lag を直接比較
        for i in range(count):
            # libc.math.fmin を使うことで Python の min() 呼び出しを回避
            temp_sum += fmin(data[i], data[i + lag])
                
        if count > 0:
            result[p_idx] = <double> temp_sum / count
        else:
            result[p_idx] = 0
            
    return result

cdef double _mean(double[:] data):
    cdef double mean,sum

    sum = 0
    for i in range(len(data)):
        sum+= data[i]    
    mean = sum/len(data)

    return mean

cdef double _min(double[:] data):
    cdef double mean,sum

    min = INFINITY
    for i in range(len(data)):
        if data[i] < min:
            min = data[i]

    return min


cdef Py_ssize_t _peak_detect_threshould_cy(double[:] bedcmm_result,
                                           DTYPE_d_t threshould,
                                           Py_ssize_t[:] search_sample,
                                           DTYPE_i_t bedcmm_smooth,
                                           str interpolator_mode):

    cdef int first_peak_idx,i,prev_sign,sign
    cdef double peak_ym1,peak_yp1,delta_x,ip_x

    # 最初のピークの取得
    prev_sign = 0
    sign = 0
    first_peak_idx = 0
    for i in range(1,len(bedcmm_result)):
        if (bedcmm_result[i] - bedcmm_result[i-1]) > 0:
            sign = 1
        else:
            sign = -1
        
        if (prev_sign == 1) and (sign == -1) and (bedcmm_result[i] > threshould):
            first_peak_idx = i - 1
            break

        prev_sign = sign

    return first_peak_idx

cdef Py_ssize_t _peak_detect_maximum_cy(double[:] bedcmm_result,
                                        Py_ssize_t[:] search_sample,
                                        DTYPE_i_t bedcmm_smooth,
                                        str interpolator_mode):
                                       
    cdef int max_peak_idx,i
    cdef double max_temp
    cdef double peak_ym1,peak_yp1,delta_x,ip_x

    # 最大値のインデックスの取得
    max_temp = -INFINITY
    max_peak_idx = -1
    for i in range(len(bedcmm_result)):
        if bedcmm_result[i] > max_temp:
            max_temp = bedcmm_result[i]
            max_peak_idx = i

    return max_peak_idx


cpdef cnp.ndarray[DTYPE_d_t, ndim=1] calc_Pitch_core_cy(double[:] data,
                                                       DTYPE_d_t fs,
                                                       DTYPE_i_t window_size,
                                                       DTYPE_i_t hop_size,
                                                       Py_ssize_t[:] search_sample,
                                                       str pp_mode,
                                                       DTYPE_i_t bedcmm_smooth,
                                                       str pitch_detect_mode,
                                                       DTYPE_d_t pitch_detect_thre,
                                                       str interpolator_mode):
    cdef list Pitch = []
    cdef int i,j
    cdef double smooth_temp
    cdef cnp.ndarray[DTYPE_d_t, ndim=1] smooth_temp_array = np.zeros(len(search_sample)-bedcmm_smooth+1)
    cdef double threshould,delta_x,ip_x,min
    cdef Py_ssize_t max_idx_int
    cdef cnp.ndarray[DTYPE_d_t, ndim=1] bedcmm_result = np.zeros(len(search_sample))
    cdef cnp.ndarray[DTYPE_d_t, ndim=1] calc_data = np.zeros(window_size)

    for i in range(window_size, len(data),hop_size):
        calc_data = np.asarray(data[i-window_size:i])
        bedcmm_result = _periodicity_1d_core_cy(calc_data,search_sample)

        if bedcmm_smooth > 1:
            for j in range(len(search_sample)-bedcmm_smooth+1):
                # smoothing
                smooth_temp = 0
                for k in range(bedcmm_smooth):
                    smooth_temp+=<double> bedcmm_result[j+k]/bedcmm_smooth
                smooth_temp_array[j] = smooth_temp
            bedcmm_result = smooth_temp_array
        elif bedcmm_smooth == 1:
            bedcmm_result = bedcmm_result
        else:
            raise Exception('bedcmm_smooth > 0 and int')
        
        if pitch_detect_mode == 'dynamic':
            threshould = _mean(calc_data)*pitch_detect_thre
            max_idx_int = _peak_detect_threshould_cy(bedcmm_result,threshould,search_sample,bedcmm_smooth,interpolator_mode)
        elif pitch_detect_mode == 'static':
            max_idx_int = _peak_detect_threshould_cy(bedcmm_result,pitch_detect_thre,search_sample,bedcmm_smooth,interpolator_mode)
        elif pitch_detect_mode == 'maximum':
            max_idx_int = _peak_detect_maximum_cy(bedcmm_result,search_sample,bedcmm_smooth,interpolator_mode)

        if max_idx_int != 0:
            if interpolator_mode == 'parabolic':
                delta_x = _parabolic_peak(bedcmm_result,max_idx_int)
            elif interpolator_mode == 'gaussian':

                if pp_mode == 'threshould_diff':
                    min = _min(bedcmm_result)
                    for j in range(len(bedcmm_result)):
                        bedcmm_result[i] = bedcmm_result[i] - min

                delta_x = _gaussian_peak(bedcmm_result,max_idx_int)
            elif interpolator_mode == 'centroid':
            
                if pp_mode == 'threshould_diff':
                    min = _min(bedcmm_result)
                    for j in range(len(bedcmm_result)):
                        bedcmm_result[i] = bedcmm_result[i] - min

                delta_x = _centroid_peak(bedcmm_result,max_idx_int)
            elif interpolator_mode == 'no':
                delta_x = 0.0
            else:
                raise Exception('interpolator_mode is quadratic,centroid,gaussian or no')

            ip_x = <double> search_sample[max_idx_int] + delta_x + ((bedcmm_smooth-1)/2)

        else:
            ip_x = np.nan

        if np.isnan(ip_x):
            Pitch.append(np.nan)
        else:
            Pitch.append(fs/ip_x)

    return np.array(Pitch,dtype=np.float64)

cpdef cnp.ndarray[DTYPE_d_t, ndim=1] calc_Pitch_negaposi_core_cy(double[:] data_posi,
                                                                double[:] data_nega,
                                                                DTYPE_d_t fs,
                                                                DTYPE_i_t window_size,
                                                                DTYPE_i_t hop_size,
                                                                Py_ssize_t[:] search_sample,
                                                                str pp_mode,
                                                                DTYPE_i_t bedcmm_smooth,
                                                                str pitch_detect_mode,
                                                                DTYPE_d_t pitch_detect_thre,
                                                                str interpolator_mode):
    cdef list Pitch = []
    cdef int i,j
    cdef double smooth_temp
    cdef cnp.ndarray[DTYPE_d_t, ndim=1] smooth_temp_array = np.zeros(len(search_sample)-bedcmm_smooth+1)
    cdef double threshould,delta_x,ip_x
    cdef Py_ssize_t max_idx_int
    cdef cnp.ndarray[DTYPE_d_t, ndim=1] bedcmm_result = np.zeros(len(search_sample))
    cdef cnp.ndarray[DTYPE_d_t, ndim=1] bedcmm_result_posi = np.zeros(len(search_sample))
    cdef cnp.ndarray[DTYPE_d_t, ndim=1] bedcmm_result_nega = np.zeros(len(search_sample))
    cdef cnp.ndarray[DTYPE_d_t, ndim=1] calc_data_posi = np.zeros(window_size)
    cdef cnp.ndarray[DTYPE_d_t, ndim=1] calc_data_nega = np.zeros(window_size)

    for i in range(window_size, len(data_posi),hop_size):
        calc_data_posi = np.asarray(data_posi[i-window_size:i])
        calc_data_nega = np.asarray(data_nega[i-window_size:i])
        bedcmm_result = _periodicity_1d_core_cy(calc_data_posi,search_sample) + _periodicity_1d_core_cy(calc_data_nega,search_sample)

        if bedcmm_smooth > 1:
            for j in range(len(search_sample)-bedcmm_smooth+1):
                # smoothing
                smooth_temp = 0
                for k in range(bedcmm_smooth):
                    smooth_temp+=bedcmm_result[j+k]/bedcmm_smooth
                smooth_temp_array[j] = smooth_temp
            bedcmm_result = smooth_temp_array
        elif bedcmm_smooth == 1:
            bedcmm_result = bedcmm_result
        else:
            raise Exception('bedcmm_smooth > 0 and int')
        
        if pitch_detect_mode == 'dynamic':
            threshould = (_mean(calc_data_posi)+_mean(calc_data_nega))*pitch_detect_thre
            max_idx_int = _peak_detect_threshould_cy(bedcmm_result,threshould,search_sample,bedcmm_smooth,interpolator_mode)
        elif pitch_detect_mode == 'static':
            max_idx_int = _peak_detect_threshould_cy(bedcmm_result,pitch_detect_thre,search_sample,bedcmm_smooth,interpolator_mode)
        elif pitch_detect_mode == 'maximum':
            max_idx_int = _peak_detect_maximum_cy(bedcmm_result,search_sample,bedcmm_smooth,interpolator_mode)

        if max_idx_int != 0:
            if interpolator_mode == 'parabolic':
                delta_x = _parabolic_peak(bedcmm_result,max_idx_int)

            elif interpolator_mode == 'gaussian':
            
                if pp_mode == 'threshould_diff':
                    min = _min(bedcmm_result)
                    for j in range(len(bedcmm_result)):
                        bedcmm_result[i] = bedcmm_result[i] - min

                delta_x = _gaussian_peak(bedcmm_result,max_idx_int)
            elif interpolator_mode == 'centroid':
            
                if pp_mode == 'threshould_diff':
                    min = _min(bedcmm_result)
                    for j in range(len(bedcmm_result)):
                        bedcmm_result[i] = bedcmm_result[i] - min

                delta_x = _centroid_peak(bedcmm_result,max_idx_int)
            elif interpolator_mode == 'no':
                delta_x = 0.0
            else:
                raise Exception('interpolator_mode is quadratic,centroid,gaussian or no')

            ip_x = <double> search_sample[max_idx_int] + delta_x + ((bedcmm_smooth-1)/2)

        else:
            ip_x = np.nan

        if np.isnan(ip_x):
            Pitch.append(np.nan)
        else:
            Pitch.append(fs/ip_x)

    return np.array(Pitch,dtype=np.float64)

cpdef cnp.ndarray[DTYPE_d_t, ndim=1] calc_bedcmm_core_cy(double[:] data,
                                                        DTYPE_i_t window_size,
                                                        DTYPE_i_t hop_size,
                                                        Py_ssize_t[:] search_sample):

    cdef cnp.ndarray[DTYPE_d_t, ndim=1] calc_data = np.zeros(window_size)

    bedcmm_result_list = []
    for i in range(window_size, len(data),hop_size):
        calc_data = np.asarray(data[i-window_size:i])
        bedcmm_result_list.append(_periodicity_1d_core_cy(calc_data,search_sample))

    bedcmm_result = np.array(bedcmm_result_list)

    return bedcmm_result

cpdef cnp.ndarray[DTYPE_d_t, ndim=1] calc_bedcmm_negaposi_core_cy(double[:] data_posi,
                                                                 double[:] data_nega,
                                                                 DTYPE_i_t window_size,
                                                                 DTYPE_i_t hop_size,
                                                                 Py_ssize_t[:] search_sample):
    bedcmm_result_list = []
    cdef cnp.ndarray[DTYPE_d_t, ndim=1] calc_data_posi = np.zeros(window_size)
    cdef cnp.ndarray[DTYPE_d_t, ndim=1] calc_data_nega = np.zeros(window_size)

    for i in range(window_size, len(data_posi),hop_size):
        calc_data_posi = np.asarray(data_posi[i-window_size:i])
        calc_data_nega = np.asarray(data_nega[i-window_size:i])
        bedcmm_result_list.append(_periodicity_1d_core_cy(calc_data_posi,search_sample) + _periodicity_1d_core_cy(calc_data_nega,search_sample))

    bedcmm_result = np.array(bedcmm_result_list)

    return bedcmm_result

