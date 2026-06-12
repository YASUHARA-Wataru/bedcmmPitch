# -*- coding: utf-8 -*-
# cython: boundscheck=False, wraparound=False, cdivision=True
# distutils: define_macros=NPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION
# cython: language_level=3
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
ctypedef cnp.intp_t DTYPE_ip_t

# Cレベルの関数を利用するために math をインポート
from libc.math cimport fmin,log,INFINITY,fabs,exp,pow,log2

cdef double EPS = 1e-300

cdef (double,double)_parabolic_peak(double[:] y, Py_ssize_t i):
    cdef double ym1, y0, yp1, denom,delta,y_peak
    if i <= 0 or i >= y.shape[0] - 1:
        return 0.0,0.0
    ym1 = y[i - 1]
    y0 = y[i]
    yp1 = y[i + 1]
    denom = ym1 - 2.0 * y0 + yp1

    if fabs(denom) < EPS:
        return 0.0,y0

    delta = 0.5 * (ym1 - yp1) / denom
    y_peak = y0 - 0.25 * (ym1 - yp1) * delta    
 
    return delta,y_peak


cdef (double,double) _gaussian_peak(double[:] y, Py_ssize_t i):
    cdef double ym1, y0, yp1, denom, delta

    if i <= 0 or i >= y.shape[0] - 1:
        return 0.0,0.0

    ym1 = log(max(y[i - 1], EPS))
    y0 = log(max(y[i], EPS))
    yp1 = log(max(y[i + 1], EPS))
    denom = ym1 - 2.0 * y0 + yp1
    if fabs(denom) < 1e-300:
        return 0.0,exp(y0)
    delta = 0.5 * (ym1 - yp1) / denom 

    return delta ,exp(y0 - 0.25 * (ym1 - yp1) * delta)

cdef (double,double) _centroid_peak(double[:] y, Py_ssize_t i, int half_window=1):
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
        return 0.0,0.0

    x = sw / s
    return x - (<double>i),y[i]


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
    cdef Py_ssize_t i

    sum = 0
    for i in range(len(data)):
        sum+= data[i]    
    mean = sum/len(data)

    return mean

cdef double _min(double[:] data):
    cdef double mean,sum
    cdef Py_ssize_t i

    min = INFINITY
    for i in range(len(data)):
        if data[i] < min:
            min = data[i]

    return min

cdef double _sum(double[:] data):
    cdef double sum
    cdef Py_ssize_t i

    sum = 0
    for i in range(len(data)):
        sum+= data[i]    

    return sum


cdef Py_ssize_t _argmax(double[:] x):

    cdef:
        Py_ssize_t i
        Py_ssize_t idx=0

    cdef double v=x[0]

    for i in range(1,x.shape[0]):
        if x[i] > v:
            v = x[i]
            idx = i

    return idx

cdef Py_ssize_t _peak_detect_threshold_cy(double[:] bedcmm_result,
                                           DTYPE_d_t threshold):

    cdef Py_ssize_t first_peak_idx,i
    cdef int prev_sign,sign
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
        
        if (prev_sign == 1) and (sign == -1) and (bedcmm_result[i] > threshold):
            first_peak_idx = i - 1
            break

        prev_sign = sign

    return first_peak_idx

cdef Py_ssize_t _peak_detect_maximum_cy(double[:] bedcmm_result):

    cdef list peaks
    cdef double max_value
    cdef Py_ssize_t max_index
    # ピークの取得
    cdef int prev_sign = 0
    cdef int sign = 0
    cdef Py_ssize_t i

    peaks = []
    for i in range(1,len(bedcmm_result)):
        if (bedcmm_result[i] - bedcmm_result[i-1]) > 0:
            sign = 1
        else:
            sign = -1
        
        if (prev_sign == 1) and (sign == -1) and (bedcmm_result[i]):
            peaks.append(i - 1)

        prev_sign = sign

    if len(peaks) == 0:
        max_index = np.nan
    else:
        # 最大値の取得
        max_value = -INFINITY
        for i in peaks:
            if max_value < bedcmm_result[i]:
                max_value = bedcmm_result[i]
                max_index = i

    return max_index

cdef _calc_peak_max_value(double[:] bedcmm_result):

    cdef list peaks
    cdef double max_value
    # ピークの取得
    cdef int prev_sign = 0
    cdef int sign = 0
    cdef Py_ssize_t i

    peaks = []
    for i in range(1,len(bedcmm_result)):
        if (bedcmm_result[i] - bedcmm_result[i-1]) > 0:
            sign = 1
        else:
            sign = -1
        
        if (prev_sign == 1) and (sign == -1):
            peaks.append(i - 1)

        prev_sign = sign

    if len(peaks) == 0:
        max_value = np.nan
    else:
        # 最大値の取得
        max_value = -INFINITY
        for i in peaks:
            if max_value < bedcmm_result[i]:
                max_value = bedcmm_result[i]

    return max_value

cpdef cnp.ndarray[DTYPE_d_t, ndim=2] calc_Pitch_core_cy(double[:] data,
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
    cdef double threshold,delta_x,ip_x,min,peak_value,mean_data,score
    cdef Py_ssize_t max_idx_int
    cdef cnp.ndarray[DTYPE_d_t, ndim=1] bedcmm_result = np.zeros(len(search_sample))
    cdef cnp.ndarray[DTYPE_d_t, ndim=1] calc_data = np.zeros(window_size)

    for i in range(window_size, len(data),hop_size):
        calc_data = np.asarray(data[i-window_size:i])
        bedcmm_result = _periodicity_1d_core_cy(calc_data,search_sample)
        mean_data = _mean(calc_data)

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
        
        if pitch_detect_mode == 'score':
            threshold = mean_data*pitch_detect_thre
            max_idx_int = _peak_detect_threshold_cy(bedcmm_result,threshold)
        elif pitch_detect_mode == 'static':
            max_idx_int = _peak_detect_threshold_cy(bedcmm_result,pitch_detect_thre)
        elif pitch_detect_mode == 'maximum':
            max_idx_int = _peak_detect_maximum_cy(bedcmm_result)
        elif pitch_detect_mode == 'peak':
            peak_value = _calc_peak_max_value(bedcmm_result)
            threshold = peak_value*pitch_detect_thre
            max_idx_int = _peak_detect_threshold_cy(bedcmm_result,threshold)
        else:
            raise Exception('pitch_detect_mode is score,static,maximum,peak.')

        if max_idx_int != 0:
            if interpolator_mode == 'parabolic':
                delta_x,peak_value = _parabolic_peak(bedcmm_result,max_idx_int)
            elif interpolator_mode == 'gaussian':

                if pp_mode == 'threshold_diff':
                    min = _min(bedcmm_result)
                    for j in range(len(bedcmm_result)):
                        bedcmm_result[i] = bedcmm_result[i] - min

                delta_x,peak_value = _gaussian_peak(bedcmm_result,max_idx_int)
            elif interpolator_mode == 'centroid':
            
                if pp_mode == 'threshold_diff':
                    min = _min(bedcmm_result)
                    for j in range(len(bedcmm_result)):
                        bedcmm_result[j] = bedcmm_result[j] - min

                delta_x,peak_value = _centroid_peak(bedcmm_result,max_idx_int)
            elif interpolator_mode == 'no':
                peak_value = bedcmm_result[max_idx_int]
                delta_x = 0.0
            else:
                raise Exception('interpolator_mode is parabolic,centroid,gaussian or no')

            ip_x = <double> search_sample[max_idx_int] + delta_x + ((bedcmm_smooth-1)/2)
            score = peak_value/mean_data
        else:
            ip_x = np.nan
            score = np.nan

        if np.isnan(ip_x):
            Pitch.append([np.nan,np.nan])
        else:
            Pitch.append([fs/ip_x,score])

    return np.array(Pitch,dtype=np.float64)

cpdef cnp.ndarray[DTYPE_d_t, ndim=2] calc_Pitch_negaposi_core_cy(double[:] data_posi,
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
    cdef double threshold,delta_x,ip_x,peak_value,mean_data,score
    cdef Py_ssize_t max_idx_int
    cdef cnp.ndarray[DTYPE_d_t, ndim=1] bedcmm_result = np.zeros(len(search_sample))
    cdef cnp.ndarray[DTYPE_d_t, ndim=1] calc_data_posi = np.zeros(window_size)
    cdef cnp.ndarray[DTYPE_d_t, ndim=1] calc_data_nega = np.zeros(window_size)

    for i in range(window_size, len(data_posi),hop_size):
        calc_data_posi = np.asarray(data_posi[i-window_size:i])
        calc_data_nega = np.asarray(data_nega[i-window_size:i])
        mean_data = _mean(calc_data_posi)+_mean(calc_data_nega)
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
        
        if pitch_detect_mode == 'score':
            threshold = mean_data*pitch_detect_thre
            max_idx_int = _peak_detect_threshold_cy(bedcmm_result,threshold)
        elif pitch_detect_mode == 'static':
            max_idx_int = _peak_detect_threshold_cy(bedcmm_result,pitch_detect_thre)
        elif pitch_detect_mode == 'maximum':
            max_idx_int = _peak_detect_maximum_cy(bedcmm_result)
        elif pitch_detect_mode == 'peak':
            peak_value = _calc_peak_max_value(bedcmm_result)
            threshold = peak_value*pitch_detect_thre
            max_idx_int = _peak_detect_threshold_cy(bedcmm_result,threshold)
        else:
            raise Exception('pitch_detect_mode is score,static,maximum,peak.')

        if max_idx_int != 0:
            if interpolator_mode == 'parabolic':
                delta_x,peak_value = _parabolic_peak(bedcmm_result,max_idx_int)

            elif interpolator_mode == 'gaussian':
            
                if pp_mode == 'threshold_diff':
                    min = _min(bedcmm_result)
                    for j in range(len(bedcmm_result)):
                        bedcmm_result[j] = bedcmm_result[j] - min

                delta_x,peak_value = _gaussian_peak(bedcmm_result,max_idx_int)
            elif interpolator_mode == 'centroid':
            
                if pp_mode == 'threshold_diff':
                    min = _min(bedcmm_result)
                    for j in range(len(bedcmm_result)):
                        bedcmm_result[i] = bedcmm_result[i] - min

                delta_x,peak_value = _centroid_peak(bedcmm_result,max_idx_int)
            elif interpolator_mode == 'no':
                peak_value = bedcmm_result[max_idx_int]
                delta_x = 0.0
            else:
                raise Exception('interpolator_mode is quadratic,centroid,gaussian or no')

            ip_x = <double> search_sample[max_idx_int] + delta_x + ((bedcmm_smooth-1)/2)
            score = peak_value/mean_data
        else:
            ip_x = np.nan
            score = np.nan

        if np.isnan(ip_x):
            Pitch.append([np.nan,np.nan])
        else:
            Pitch.append([fs/ip_x,score])

    return np.array(Pitch,dtype=np.float64)

cpdef calc_bedcmm_core_cy(double[:] data,
                          DTYPE_i_t window_size,
                          DTYPE_i_t hop_size,
                          Py_ssize_t[:] search_sample):
 
    cdef cnp.ndarray[DTYPE_d_t, ndim=1] calc_data = np.zeros(window_size)
    cdef list bedcmm_result_list,mean_data_list

    bedcmm_result_list = []
    mean_data_list =[]

    for i in range(window_size, len(data),hop_size):
        calc_data = np.asarray(data[i-window_size:i])
        bedcmm_result_list.append(_periodicity_1d_core_cy(calc_data,search_sample))
        mean_data_list.append(_mean(calc_data))

    bedcmm_result = np.array(bedcmm_result_list)
    mean_data = np.array(mean_data_list)

    return bedcmm_result,mean_data

cpdef calc_bedcmm_negaposi_core_cy(double[:] data_posi,
                                   double[:] data_nega,
                                   DTYPE_i_t window_size,
                                   DTYPE_i_t hop_size,
                                   Py_ssize_t[:] search_sample):

    cdef cnp.ndarray[DTYPE_d_t, ndim=1] calc_data_posi = np.zeros(window_size)
    cdef cnp.ndarray[DTYPE_d_t, ndim=1] calc_data_nega = np.zeros(window_size)
    cdef list bedcmm_result_list
    cdef list mean_data_list

    bedcmm_result_list = []
    mean_data_list = []

    for i in range(window_size, len(data_posi),hop_size):
        calc_data_posi = np.asarray(data_posi[i-window_size:i])
        calc_data_nega = np.asarray(data_nega[i-window_size:i])
        bedcmm_result_list.append(_periodicity_1d_core_cy(calc_data_posi,search_sample) + _periodicity_1d_core_cy(calc_data_nega,search_sample))
        mean_data_list.append(_mean(calc_data_posi)+_mean(calc_data_nega))

    bedcmm_result = np.array(bedcmm_result_list)
    mean_data = np.array(mean_data_list)

    return bedcmm_result,mean_data

cpdef cnp.ndarray[DTYPE_d_t, ndim=1] _prev_first_peak_zero(double[:] data):

    cdef cnp.ndarray[DTYPE_d_t, ndim=1] result = np.copy(data)
    cdef Py_ssize_t first_peak_idx = 0
    # ピークの取得
    cdef int prev_sign = 0
    cdef int sign = 0
    cdef Py_ssize_t i

    for i in range(1,len(data)):
        if (data[i] - data[i-1]) > 0:
            sign = 1
        else:
            sign = -1
        
        if prev_sign*sign == -1:
            first_peak_idx = i
            break

        prev_sign = sign

    if (first_peak_idx != 0):
        for i in range(first_peak_idx):
            result[i] = 0
    else:
        for i in range(len(data)):
            result[i] = 1.0 / data.shape[0]

    return result

cpdef cnp.ndarray[DTYPE_d_t, ndim=2] create_transition_matrix_cy(
    int n_periods,
    double fs,
    double sigma
):
    cdef:
        cnp.ndarray[DTYPE_d_t, ndim=2] A
        cnp.ndarray[DTYPE_d_t, ndim=1] log_freqs
        int i,j
        double diff
        double denom

    log_freqs = np.empty(n_periods,dtype=np.float64)

    for i in range(n_periods):
        log_freqs[i] = log2(fs/(i+1))

    A = np.empty((n_periods,n_periods),dtype=np.float64)

    denom = 2.0*sigma*sigma

    for i in range(n_periods):
        for j in range(n_periods):
            diff = log_freqs[i]-log_freqs[j]
            A[i,j] = exp(-(diff*diff)/denom)

    cdef double s

    for j in range(n_periods):
        s = 0.0

        for i in range(n_periods):
            s += A[i,j]

        for i in range(n_periods):
            A[i,j] /= s

    return A

cdef void matvec(
    double[:, :] A,
    double[:] x,
    double[:] y
):
    cdef:
        Py_ssize_t i,j
        Py_ssize_t n=A.shape[0]
        double s

    for i in range(n):

        s = 0.0

        for j in range(n):
            s += A[i,j]*x[j]

        y[i] = s

cdef cnp.ndarray[DTYPE_d_t, ndim=1] _probability(double[:] data):
    cdef double sum
    cdef cnp.ndarray[DTYPE_d_t, ndim=1] result = np.zeros_like(data)

    cdef Py_ssize_t i

    sum = _sum(data)

    for i in range(len(data)):
        result[i] = data[i]/sum

    return result

cpdef cnp.ndarray[DTYPE_d_t, ndim=2] calc_Pitch_bayes_negaposi_core_cy(double[:] data_posi,
                                                                       double[:] data_nega,
                                                                       DTYPE_d_t fs,
                                                                       DTYPE_i_t window_size,
                                                                       DTYPE_i_t hop_size,
                                                                       Py_ssize_t[:] search_sample,
                                                                       str pp_mode,
                                                                       DTYPE_d_t alpha,
                                                                       DTYPE_d_t sigma,
                                                                       str interpolator_mode):
    cdef list Pitch = []
    cdef int i,j
    cdef double delta_x,ip_x,peak_value,mean_data,score,probability
    cdef Py_ssize_t max_idx_int
    cdef cnp.ndarray[DTYPE_d_t, ndim=1] bedcmm_result = np.zeros(len(search_sample))
    cdef cnp.ndarray[DTYPE_d_t, ndim=1] calc_data_posi = np.zeros(window_size)
    cdef cnp.ndarray[DTYPE_d_t, ndim=1] calc_data_nega = np.zeros(window_size)
    cdef cnp.ndarray[DTYPE_d_t, ndim=1] likelihoods = np.zeros(len(search_sample))
    cdef cnp.ndarray[DTYPE_d_t, ndim=1] posterior = np.zeros(len(search_sample))
    cdef cnp.ndarray[DTYPE_d_t, ndim=1] prediction = np.zeros(len(search_sample))
    cdef cnp.ndarray[DTYPE_d_t, ndim=2] A_org = np.zeros((window_size//2,window_size//2))
    cdef cnp.ndarray[DTYPE_d_t, ndim=2] A = np.zeros((len(search_sample),len(search_sample)))
    
    A_org = create_transition_matrix_cy(window_size//2,fs=fs,sigma=sigma)
    A = A_org[search_sample[0]:search_sample[len(search_sample)-1]+1,search_sample[0]:search_sample[len(search_sample)-1]+1]

    for i in range(window_size, len(data_posi),hop_size):
        calc_data_posi = np.asarray(data_posi[i-window_size:i])
        calc_data_nega = np.asarray(data_nega[i-window_size:i])
        mean_data = _mean(calc_data_posi)+_mean(calc_data_nega)
        bedcmm_result = _periodicity_1d_core_cy(calc_data_posi,search_sample) + _periodicity_1d_core_cy(calc_data_nega,search_sample)
        likelihoods = _prev_first_peak_zero(bedcmm_result)
        likelihoods = _probability(likelihoods)

        if _sum(posterior) < 0.1:
            posterior = likelihoods
        else:
            matvec(A,posterior,prediction)
            for j in range(len(prediction)):
                posterior[j] = pow(prediction[j],alpha) * likelihoods[j] 

        posterior = _probability(posterior)

        max_idx_int = _argmax(posterior)

        if max_idx_int != 0:
            if interpolator_mode == 'parabolic':
                delta_x,peak_value = _parabolic_peak(bedcmm_result,max_idx_int)

            elif interpolator_mode == 'gaussian':
            
                if pp_mode == 'threshold_diff':
                    min = _min(bedcmm_result)
                    for j in range(len(bedcmm_result)):
                        bedcmm_result[j] = bedcmm_result[j] - min

                delta_x,peak_value = _gaussian_peak(bedcmm_result,max_idx_int)
            elif interpolator_mode == 'centroid':
            
                if pp_mode == 'threshold_diff':
                    min = _min(bedcmm_result)
                    for j in range(len(bedcmm_result)):
                        bedcmm_result[i] = bedcmm_result[i] - min

                delta_x,peak_value = _centroid_peak(bedcmm_result,max_idx_int)
            elif interpolator_mode == 'no':
                peak_value = bedcmm_result[max_idx_int]
                delta_x = 0.0
            else:
                raise Exception('interpolator_mode is quadratic,centroid,gaussian or no')

            ip_x = <double> search_sample[max_idx_int] + delta_x
            score = peak_value/mean_data
            probability = posterior[max_idx_int]
        else:
            ip_x = np.nan
            score = np.nan
            probability = np.nan

        if np.isnan(ip_x):
            Pitch.append([np.nan,np.nan,np.nan])
        else:
            Pitch.append([fs/ip_x,score,probability])

    return np.array(Pitch,dtype=np.float64)

cpdef cnp.ndarray[DTYPE_d_t, ndim=2] calc_Pitch_bayes_core_cy(double[:] data,
                                                              DTYPE_d_t fs,
                                                              DTYPE_i_t window_size,
                                                              DTYPE_i_t hop_size,
                                                              Py_ssize_t[:] search_sample,
                                                              str pp_mode,
                                                              DTYPE_d_t alpha,
                                                              DTYPE_d_t sigma,
                                                              str interpolator_mode):
    cdef list Pitch = []
    cdef int i,j
    cdef double delta_x,ip_x,peak_value,mean_data,score,probability
    cdef Py_ssize_t max_idx_int
    cdef cnp.ndarray[DTYPE_d_t, ndim=1] bedcmm_result = np.zeros(len(search_sample))
    cdef cnp.ndarray[DTYPE_d_t, ndim=1] bedcmm_result_posi = np.zeros(len(search_sample))
    cdef cnp.ndarray[DTYPE_d_t, ndim=1] bedcmm_result_nega = np.zeros(len(search_sample))
    cdef cnp.ndarray[DTYPE_d_t, ndim=1] calc_data = np.zeros(window_size)
    cdef cnp.ndarray[DTYPE_d_t, ndim=1] likelihoods = np.zeros(len(search_sample))
    cdef cnp.ndarray[DTYPE_d_t, ndim=1] posterior = np.zeros(len(search_sample))
    cdef cnp.ndarray[DTYPE_d_t, ndim=1] prediction = np.zeros(len(search_sample))
    cdef cnp.ndarray[DTYPE_d_t, ndim=2] A_org = np.zeros((window_size//2,window_size//2))
    cdef cnp.ndarray[DTYPE_d_t, ndim=2] A = np.zeros((len(search_sample),len(search_sample)))
    
    A_org = create_transition_matrix_cy(window_size//2,fs=fs,sigma=sigma)
    A = A_org[search_sample[0]:search_sample[len(search_sample)-1]+1,search_sample[0]:search_sample[len(search_sample)-1]+1]

    for i in range(window_size, len(data),hop_size):
        calc_data = np.asarray(data[i-window_size:i])
        mean_data = _mean(calc_data)
        bedcmm_result = _periodicity_1d_core_cy(calc_data,search_sample)
        likelihoods = _prev_first_peak_zero(bedcmm_result)
        likelihoods = _probability(likelihoods)

        if _sum(posterior) < 0.5:
            posterior = likelihoods
        else:
            matvec(A,posterior,prediction)
            for j in range(len(prediction)):
                posterior[j] = pow(prediction[j],alpha) * likelihoods[j] 

        posterior = _probability(posterior)

        max_idx_int = _argmax(posterior)

        if max_idx_int != 0:
            if interpolator_mode == 'parabolic':
                delta_x,peak_value = _parabolic_peak(bedcmm_result,max_idx_int)

            elif interpolator_mode == 'gaussian':
            
                if pp_mode == 'threshold_diff':
                    min = _min(bedcmm_result)
                    for j in range(len(bedcmm_result)):
                        bedcmm_result[i] = bedcmm_result[i] - min

                delta_x,peak_value = _gaussian_peak(bedcmm_result,max_idx_int)
            elif interpolator_mode == 'centroid':
            
                if pp_mode == 'threshold_diff':
                    min = _min(bedcmm_result)
                    for j in range(len(bedcmm_result)):
                        bedcmm_result[i] = bedcmm_result[i] - min

                delta_x,peak_value = _centroid_peak(bedcmm_result,max_idx_int)
            elif interpolator_mode == 'no':
                peak_value = bedcmm_result[max_idx_int]
                delta_x = 0.0
            else:
                raise Exception('interpolator_mode is quadratic,centroid,gaussian or no')

            ip_x = <double> search_sample[max_idx_int] + delta_x
            score = peak_value/mean_data
            probability = posterior[max_idx_int]
        else:
            ip_x = np.nan
            score = np.nan
            probability = np.nan

        if np.isnan(ip_x):
            Pitch.append([np.nan,np.nan,np.nan])
        else:
            Pitch.append([fs/ip_x,score,probability])

    return np.array(Pitch,dtype=np.float64)


cpdef cnp.ndarray[DTYPE_i_t, ndim=1] viterbi_cy(
    double[:, :] log_likelihood,
    double[:, :] log_A,
    double beta
):

    cdef:
        Py_ssize_t n_frames = log_likelihood.shape[0]
        Py_ssize_t n_states = log_likelihood.shape[1]

        cnp.ndarray[DTYPE_d_t, ndim=2] score_np
        cnp.ndarray[DTYPE_ip_t, ndim=2] back_ptr_np
        cnp.ndarray[DTYPE_ip_t, ndim=1] path_np

        double[:, :] score
        Py_ssize_t[:, :] back_ptr
        Py_ssize_t[:] path

        Py_ssize_t t,j,i
        Py_ssize_t best_prev

        double best_score
        double candidate

    score_np = np.zeros(
        (n_frames,n_states),
        dtype=np.float64
    )

    back_ptr_np = np.zeros(
        (n_frames,n_states),
        dtype=np.intp
    )

    path_np = np.zeros(
        n_frames,
        dtype=np.intp
    )

    score = score_np
    back_ptr = back_ptr_np
    path = path_np

    for j in range(n_states):
        score[0,j] = log_likelihood[0,j]

    for t in range(1,n_frames):

        for j in range(n_states):

            best_prev = 0
            best_score = (
                score[t-1,0]
                + log_A[0,j]
            )

            for i in range(1,n_states):

                candidate = (
                    score[t-1,i]
                    + log_A[i,j]
                )

                if candidate > best_score:
                    best_score = candidate
                    best_prev = i

            score[t,j] = (
                best_score
                + beta*log_likelihood[t,j]
            )

            back_ptr[t,j] = best_prev

    best_prev = 0
    best_score = score[n_frames-1,0]

    for i in range(1,n_states):

        if score[n_frames-1,i] > best_score:
            best_score = score[n_frames-1,i]
            best_prev = i

    path[n_frames-1] = best_prev

    for t in range(n_frames-2,-1,-1):

        path[t] = back_ptr[
            t+1,
            path[t+1]
        ]

    return path_np


cdef cnp.ndarray[DTYPE_d_t, ndim=2] _log_dim2(double[:,:] A):

    cdef:
        Py_ssize_t max_i = A.shape[0]
        Py_ssize_t max_j = A.shape[1]
        Py_ssize_t i,j
        cnp.ndarray[DTYPE_d_t, ndim=2] result = np.zeros((max_i,max_j))


    for i in range(max_i):
        for j in range(max_j):
            result[i,j] = log(A[i,j] + EPS)

    return result

cdef cnp.ndarray[DTYPE_d_t, ndim=1] _log_dim1(double[:] A):

    cdef:
        Py_ssize_t max_i = A.shape[0]
        Py_ssize_t i
        cnp.ndarray[DTYPE_d_t, ndim=1] result = np.zeros(max_i)


    for i in range(max_i):
        result[i] = log(A[i] + EPS)

    return result


cpdef cnp.ndarray[DTYPE_d_t, ndim=2] calc_Pitch_viterbi_negaposi_core_cy(double[:] data_posi,
                                                                         double[:] data_nega,
                                                                         DTYPE_d_t fs,
                                                                         DTYPE_i_t window_size,
                                                                         DTYPE_i_t hop_size,
                                                                         Py_ssize_t[:] search_sample,
                                                                         str pp_mode,
                                                                         DTYPE_d_t beta,
                                                                         DTYPE_d_t sigma,
                                                                         str interpolator_mode):
    cdef list Pitch = []
    cdef list likelihood_list = []
    cdef list likelihood_list_log = []
    cdef list bedcmm_result_list = []
    cdef list mean_data_list = []
    cdef int i,frame_num
    cdef double delta_x,ip_x,peak_value,mean_data,score,probability
    cdef Py_ssize_t max_idx_int
    cdef cnp.ndarray[DTYPE_d_t, ndim=1] bedcmm_result = np.zeros(len(search_sample))
    cdef cnp.ndarray[DTYPE_d_t, ndim=1] calc_data_posi = np.zeros(window_size)
    cdef cnp.ndarray[DTYPE_d_t, ndim=1] calc_data_nega = np.zeros(window_size)
    cdef cnp.ndarray[DTYPE_d_t, ndim=1] likelihoods = np.zeros(len(search_sample))
    cdef cnp.ndarray[DTYPE_d_t, ndim=1] posterior = np.zeros(len(search_sample))
    cdef cnp.ndarray[DTYPE_d_t, ndim=1] prediction = np.zeros(len(search_sample))
    cdef cnp.ndarray[DTYPE_d_t, ndim=2] A_org = np.zeros((window_size//2,window_size//2))
    cdef cnp.ndarray[DTYPE_d_t, ndim=2] A = np.zeros((len(search_sample),len(search_sample)))
    cdef cnp.ndarray[DTYPE_ip_t, ndim=1] path_list = np.zeros((len(range(window_size, len(data_posi),hop_size))),dtype=np.intp)
    cdef cnp.ndarray[DTYPE_d_t, ndim=2] likelihood_array = np.zeros((len(range(window_size, len(data_posi),hop_size)),len(search_sample)),dtype=np.float64)

    A_org = create_transition_matrix_cy(window_size//2,fs=fs,sigma=sigma)
    A = A_org[search_sample[0]:search_sample[len(search_sample)-1]+1,search_sample[0]:search_sample[len(search_sample)-1]+1]
    A = _log_dim2(A)

    for i in range(window_size, len(data_posi),hop_size):
        calc_data_posi = np.asarray(data_posi[i-window_size:i])
        calc_data_nega = np.asarray(data_nega[i-window_size:i])
        mean_data = _mean(calc_data_posi)+_mean(calc_data_nega)
        mean_data_list.append(mean_data)

        bedcmm_result = _periodicity_1d_core_cy(calc_data_posi,search_sample) + _periodicity_1d_core_cy(calc_data_nega,search_sample)
        bedcmm_result_list.append(bedcmm_result)

        likelihoods = _prev_first_peak_zero(bedcmm_result)
        likelihoods = _probability(likelihoods)
        likelihood_list.append(likelihoods)

    for likelihoods in likelihood_list:
        likelihood_list_log.append(_log_dim1(likelihoods))

    if len(likelihood_list) == 0:
        return np.array([])
    
    likelihood_array = np.array(likelihood_list_log)

    path_list = viterbi_cy(likelihood_array,A,beta=beta)

    for frame_num,max_idx_int in enumerate(path_list):
        bedcmm_result = bedcmm_result_list[frame_num]
        if max_idx_int != 0:
            if interpolator_mode == 'parabolic':
                delta_x,peak_value = _parabolic_peak(bedcmm_result,max_idx_int)

            elif interpolator_mode == 'gaussian':
            
                if pp_mode == 'threshold_diff':
                    min = _min(bedcmm_result)
                    for j in range(len(bedcmm_result)):
                        bedcmm_result[i] = bedcmm_result[i] - min

                delta_x,peak_value = _gaussian_peak(bedcmm_result,max_idx_int)
            elif interpolator_mode == 'centroid':
            
                if pp_mode == 'threshold_diff':
                    min = _min(bedcmm_result)
                    for j in range(len(bedcmm_result)):
                        bedcmm_result[i] = bedcmm_result[i] - min

                delta_x,peak_value = _centroid_peak(bedcmm_result,max_idx_int)
            elif interpolator_mode == 'no':
                peak_value = bedcmm_result[max_idx_int]
                delta_x = 0.0
            else:
                raise Exception('interpolator_mode is quadratic,centroid,gaussian or no')

            ip_x = <double> search_sample[max_idx_int] + delta_x
            score = peak_value/mean_data_list[frame_num]
            probability = likelihood_list[frame_num][max_idx_int]
        else:
            ip_x = np.nan
            score = np.nan
            probability = np.nan

        if np.isnan(ip_x):
            Pitch.append([np.nan,np.nan,np.nan])
        else:
            Pitch.append([fs/ip_x,score,probability])

    return np.array(Pitch,dtype=np.float64)


cpdef cnp.ndarray[DTYPE_d_t, ndim=2] calc_Pitch_viterbi_core_cy(double[:] data,
                                                                DTYPE_d_t fs,
                                                                DTYPE_i_t window_size,
                                                                DTYPE_i_t hop_size,
                                                                Py_ssize_t[:] search_sample,
                                                                str pp_mode,
                                                                DTYPE_d_t beta,
                                                                DTYPE_d_t sigma,
                                                                str interpolator_mode):
    cdef list Pitch = []
    cdef list likelihood_list = []
    cdef list likelihood_list_log = []
    cdef list bedcmm_result_list = []
    cdef list mean_data_list = []
    cdef int i,frame_num
    cdef double delta_x,ip_x,peak_value,mean_data,score,probability
    cdef Py_ssize_t max_idx_int
    cdef cnp.ndarray[DTYPE_d_t, ndim=1] bedcmm_result = np.zeros(len(search_sample))
    cdef cnp.ndarray[DTYPE_d_t, ndim=1] calc_data = np.zeros(window_size)
    cdef cnp.ndarray[DTYPE_d_t, ndim=1] likelihoods = np.zeros(len(search_sample))
    cdef cnp.ndarray[DTYPE_d_t, ndim=1] posterior = np.zeros(len(search_sample))
    cdef cnp.ndarray[DTYPE_d_t, ndim=1] prediction = np.zeros(len(search_sample))
    cdef cnp.ndarray[DTYPE_d_t, ndim=2] A_org = np.zeros((window_size//2,window_size//2))
    cdef cnp.ndarray[DTYPE_d_t, ndim=2] A = np.zeros((len(search_sample),len(search_sample)))
    cdef cnp.ndarray[DTYPE_ip_t, ndim=1] path_list = np.zeros((len(range(window_size, len(data),hop_size))),dtype=np.intp)
    cdef cnp.ndarray[DTYPE_d_t, ndim=2] likelihood_array = np.zeros((len(range(window_size, len(data),hop_size)),len(search_sample)),dtype=np.float64)

    A_org = create_transition_matrix_cy(window_size//2,fs=fs,sigma=sigma)
    A = A_org[search_sample[0]:search_sample[len(search_sample)-1]+1,search_sample[0]:search_sample[len(search_sample)-1]+1]
    A = _log_dim2(A)

    for i in range(window_size, len(data),hop_size):
        calc_data = np.asarray(data[i-window_size:i])
        mean_data = _mean(calc_data)
        mean_data_list.append(mean_data)

        bedcmm_result = _periodicity_1d_core_cy(calc_data,search_sample)
        bedcmm_result_list.append(bedcmm_result)

        likelihoods = _prev_first_peak_zero(bedcmm_result)
        likelihoods = _probability(likelihoods)
        likelihood_list.append(likelihoods)

    for likelihoods in likelihood_list:
        likelihood_list_log.append(_log_dim1(likelihoods))

    likelihood_array = np.array(likelihood_list_log)

    path_list = viterbi_cy(likelihood_array,A,beta=beta)

    for frame_num,max_idx_int in enumerate(path_list):
        bedcmm_result = bedcmm_result_list[frame_num]
        if max_idx_int != 0:
            if interpolator_mode == 'parabolic':
                delta_x,peak_value = _parabolic_peak(bedcmm_result,max_idx_int)

            elif interpolator_mode == 'gaussian':
            
                if pp_mode == 'threshold_diff':
                    min = _min(bedcmm_result)
                    for j in range(len(bedcmm_result)):
                        bedcmm_result[i] = bedcmm_result[i] - min

                delta_x,peak_value = _gaussian_peak(bedcmm_result,max_idx_int)
            elif interpolator_mode == 'centroid':
            
                if pp_mode == 'threshold_diff':
                    min = _min(bedcmm_result)
                    for j in range(len(bedcmm_result)):
                        bedcmm_result[i] = bedcmm_result[i] - min

                delta_x,peak_value = _centroid_peak(bedcmm_result,max_idx_int)
            elif interpolator_mode == 'no':
                peak_value = bedcmm_result[max_idx_int]
                delta_x = 0.0
            else:
                raise Exception('interpolator_mode is quadratic,centroid,gaussian or no')

            ip_x = <double> search_sample[max_idx_int] + delta_x
            score = peak_value/mean_data_list[frame_num]
            probability = likelihood_list[frame_num][max_idx_int]
        else:
            ip_x = np.nan
            score = np.nan
            probability = np.nan

        if np.isnan(ip_x):
            Pitch.append([np.nan,np.nan,np.nan])
        else:
            Pitch.append([fs/ip_x,score,probability])

    return np.array(Pitch,dtype=np.float64)

