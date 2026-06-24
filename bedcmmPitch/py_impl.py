# -*- coding: utf-8 -*-
"""
bedcmmPitch method class(Python)

Author: YASUHARA Wataru
Copyright (c) 2026, Feel a Piece of the World
"""
import numpy as np
import math
import warnings
from ._config import implementation
if implementation == 'Cython':
    from .cy_impl import calc_Pitch_core_cy,calc_Pitch_negaposi_core_cy,calc_bedcmm_core_cy,calc_bedcmm_negaposi_core_cy
    from .cy_impl import calc_Pitch_bayes_negaposi_core_cy,calc_Pitch_bayes_core_cy
    from .cy_impl import calc_Pitch_viterbi_negaposi_core_cy,calc_Pitch_viterbi_core_cy

EPS = 1e-300

def _parabolic_peak(y, i: int) -> float:
    """
    3点パラボラ補完: 返り値はピークの相対オフセット(delta)
    y[i-1], y[i], y[i+1] を使う。
    """
    y = np.asarray(y, dtype=np.float64)
    if i <= 0 or i >= len(y) - 1:
        return 0.0

    ym1 = y[i - 1]
    y0 = y[i]
    yp1 = y[i + 1]
    denom = ym1 - 2.0 * y0 + yp1

    if abs(denom) < EPS:
        return 0.0,y0

    delta = 0.5 * (ym1 - yp1) / denom
    y_peak = y0 - 0.25 * (ym1 - yp1) * delta    
    
    return  delta, y_peak 


def _gaussian_peak(y, i: int, eps: float = 1e-12) -> float:
    """
    3点 Gaussian 補完（log振幅をパラボラ補完）。
    y は非負振幅を想定。
    """
    y = np.asarray(y, dtype=np.float64)
    if i <= 0 or i >= len(y) - 1:
        return 0.0

    ym1 = math.log(max(float(y[i - 1]), eps))
    y0 = math.log(max(float(y[i]), eps))
    yp1 = math.log(max(float(y[i + 1]), eps))

    denom = (ym1 - 2*y0 + yp1)
    if abs(denom) < 1e-12:
        return 0.0, y0

    delta = 0.5 * (ym1 - yp1) / denom
    l_peak = y0 - 0.25 * (ym1 - yp1) * delta

    y_peak = math.exp(l_peak)

    return delta, y_peak

def _centroid_peak(y , i: int, half_window: int = 1) -> float:
    """
    重心法。i を中心に [i-half_window, i+half_window] の重心を返す。
    y は非負値を想定。
    最大値は、iの値をそのまま採用する事とする。
    """
    y = np.asarray(y, dtype=np.float64)
    n = len(y)
    lo = max(0, i - half_window)
    hi = min(n - 1, i + half_window)

    idx = np.arange(lo, hi + 1, dtype=np.float64)
    w = np.clip(y[lo : hi + 1], 0.0, None)

    s = float(np.sum(w))
    if s <= 0.0:
        return 0.0

    x_bar = float(np.sum(idx * w) / s)

    return x_bar - float(i),y[i]

def _periodicity(data,period):

    result = np.zeros(len(period))
    for p_idx,a_preiod in enumerate(period):

        temp_data = np.zeros(data.shape[0] - a_preiod)
        for index in range(data.shape[0] - a_preiod):
            temp_data[index] = min([data[index],data[index+a_preiod]])
                
        result[p_idx]=np.mean((temp_data))

    return result    

def _peak_detect_threshold(bedcmm_result,threshold):

    priod_thre_ind = (bedcmm_result > threshold)[:-2]
    priod_diff = np.diff(bedcmm_result)
    plus_peaks = np.where((priod_diff[:-1] > 0) & (priod_diff[1:] * priod_diff[:-1] < 0 ) & priod_thre_ind)[0]

    if len(plus_peaks) > 0:
        plus_peak_idx = plus_peaks[0]+1
    else:
        plus_peak_idx = np.nan

    return plus_peak_idx


def _peak_detect_maximum(bedcmm_result):

    plus_peak_idx = np.argmax(bedcmm_result)
        
    return plus_peak_idx

def _calc_peak_max_value(bedcmm_result):

    priod_diff = np.diff(bedcmm_result)
    plus_peaks = np.where((priod_diff[:-1] > 0) & (priod_diff[1:] * priod_diff[:-1] < 0 ))[0]
    max_value = np.max(bedcmm_result[plus_peaks])
    return max_value

def calc_Pitch_core(data,
                    fs,
                    window_size,
                    hop_size,
                    search_sample,
                    pp_mode,
                    bedcmm_smooth,
                    pitch_detect_mode,
                    pitch_detect_thre,
                    interpolator_mode):
    
    Pitch = []
    for i in range(window_size, len(data),hop_size):
        calc_data = data[i-window_size:i]
        bedcmm_result = _periodicity(calc_data,search_sample)
        mean_data = np.mean(calc_data)

        if bedcmm_smooth > 1:
            filt = np.ones(bedcmm_smooth)/bedcmm_smooth
            bedcmm_result = np.convolve(bedcmm_result,filt,mode='valid')
        elif bedcmm_smooth == 1:
            bedcmm_result = bedcmm_result
        else:
            raise Exception('bedcmm_smooth > 0 and int')
        
        if pitch_detect_mode == 'score':
            threshold = mean_data*pitch_detect_thre
            max_idx_int = _peak_detect_threshold(bedcmm_result,threshold)
        elif pitch_detect_mode == 'static':
            max_idx_int = _peak_detect_threshold(bedcmm_result,pitch_detect_thre)
        elif pitch_detect_mode == 'maximum':
            max_idx_int = _peak_detect_maximum(bedcmm_result)
        elif pitch_detect_mode == 'peak':
            peak_value = _calc_peak_max_value(bedcmm_result)
            threshold = peak_value*pitch_detect_thre
            max_idx_int = _peak_detect_threshold(bedcmm_result,threshold)
        else:
            raise Exception('pitch_detect_mode is score,static,maximum,peak.')

        if ~np.isnan(max_idx_int):
            if max_idx_int != 0:
                if interpolator_mode == 'parabolic':
                    delta,peak_value = _parabolic_peak(bedcmm_result,max_idx_int)
                elif interpolator_mode == 'gaussian':
                    if pp_mode == 'threshold_diff':
                        bedcmm_result = bedcmm_result - min(bedcmm_result)
                    delta,peak_value = _gaussian_peak(bedcmm_result,max_idx_int)
                elif interpolator_mode == 'centroid':
                    if pp_mode == 'threshold_diff':
                        bedcmm_result = bedcmm_result - min(bedcmm_result)
                    delta,peak_value = _centroid_peak(bedcmm_result,max_idx_int)
                elif interpolator_mode == 'no':
                    peak_value = bedcmm_result[max_idx_int]
                    delta = 0
                else:
                    raise Exception('interpolator_mode is parabolic,centroid,gaussian or no')

                if delta < -0.5:
                    delta = -0.5
                if delta > 0.5:
                    delta = 0.5

                peak_idx = search_sample[max_idx_int]+delta

                peak_idx = peak_idx + ((bedcmm_smooth-1)/2)
                peak_score = peak_value/mean_data
            else:
                peak_idx = np.nan
                peak_score = np.nan
        else:
            peak_idx = np.nan
            peak_score =np.nan

        if np.isnan(peak_idx):
            Pitch.append([np.nan,np.nan])
        else:
            Pitch.append([fs/peak_idx,peak_score])

    Pitch = np.array(Pitch)

    return Pitch

def calc_Pitch_negaposi_core(data_posi,data_nega,
                             fs,
                             window_size,
                             hop_size,
                             search_sample,
                             pp_mode,
                             bedcmm_smooth,
                             pitch_detect_mode,
                             pitch_detect_thre,
                             interpolator_mode):

    Pitch = []
    for i in range(window_size, len(data_posi),hop_size):
        calc_data_posi = data_posi[i-window_size:i]
        calc_data_nega = data_nega[i-window_size:i]
        bedcmm_result = _periodicity(calc_data_posi,search_sample) + _periodicity(calc_data_nega,search_sample)
        mean_data = np.mean(calc_data_posi)+np.mean(calc_data_nega)

        if bedcmm_smooth > 1:
            filt = np.ones(bedcmm_smooth)/bedcmm_smooth
            bedcmm_result = np.convolve(bedcmm_result,filt,mode='valid')
        elif bedcmm_smooth == 1:
            bedcmm_result = bedcmm_result
        else:
            raise Exception('bedcmm_smooth > 0 and int')
        
        if pitch_detect_mode == 'score':
            threshold = mean_data*pitch_detect_thre
            max_idx_int = _peak_detect_threshold(bedcmm_result,threshold)
        elif pitch_detect_mode == 'static':
            max_idx_int = _peak_detect_threshold(bedcmm_result,pitch_detect_thre)
        elif pitch_detect_mode == 'maximum':
            max_idx_int = _peak_detect_maximum(bedcmm_result)
        elif pitch_detect_mode == 'peak':
            peak_value = _calc_peak_max_value(bedcmm_result)
            threshold = peak_value*pitch_detect_thre
            max_idx_int = _peak_detect_threshold(bedcmm_result,threshold)
        else:
            raise Exception('pitch_detect_mode is score,static,maximum,peak.')

        if ~np.isnan(max_idx_int):
            if max_idx_int != 0:
                if interpolator_mode == 'parabolic':
                    delta,peak_value = _parabolic_peak(bedcmm_result,max_idx_int)
                elif interpolator_mode == 'gaussian':
                    if pp_mode == 'threshold_diff':
                        bedcmm_result = bedcmm_result - min(bedcmm_result)
                    delta,peak_value = _gaussian_peak(bedcmm_result,max_idx_int)
                elif interpolator_mode == 'centroid':
                    if pp_mode == 'threshold_diff':
                        bedcmm_result = bedcmm_result - min(bedcmm_result)
                    delta,peak_value = _centroid_peak(bedcmm_result,max_idx_int)
                elif interpolator_mode == 'no':
                    peak_value = bedcmm_result[max_idx_int]
                    delta = 0
                else:
                    raise Exception('interpolator_mode is quadratic,centroid,gaussian or no')

                if delta < -0.5:
                    delta = -0.5
                if delta > 0.5:
                    delta = 0.5

                peak_idx = search_sample[max_idx_int]+delta

                peak_idx = peak_idx + ((bedcmm_smooth-1)/2)
                peak_score = peak_value/mean_data
            else:
                peak_idx = np.nan
                peak_score =np.nan
        else:
            peak_idx = np.nan
            peak_score = np.nan

        if np.isnan(peak_idx):
            Pitch.append([np.nan,np.nan])
        else:
            Pitch.append([fs/peak_idx,peak_score])

    Pitch = np.array(Pitch)

    return Pitch


def calc_Pitch(data,
               fs=44100,
               window_size=2048,
               hop_size=256,
               fmin=65,
               fmax=2000,
               pp_mode='positive+negative',
               pp_threshold=0,
               bedcmm_smooth=3,
               pitch_detect_mode='peak',
               pitch_detect_thre=0.85,
               interpolator_mode='parabolic'):
    
    data = data.copy()
    data = np.ascontiguousarray(data, dtype=np.float64)

    if data.ndim != 1:
        raise Exception('data must be 1D array.')

    # データ前処理
    if pp_mode == 'positive':
        data[data < 0] = 0
    elif pp_mode == 'negative':
        data[data > 0] = 0
        data[data < 0] = -data[data < 0]
    elif pp_mode == 'positive+negative':
        data_pos = np.zeros_like(data)
        data_neg = np.zeros_like(data)
        data_pos[data > 0] = data[data > 0]
        data_neg[data < 0] = -data[data < 0]
    elif pp_mode == 'threshold_diff':
        data = data - pp_threshold
    else:
        raise Exception('pp_mode is only positive,negative,positive+negative,threshold_diff.')

    if fmin is None:
        if fmax is None:
            search_sample = np.arange(int(window_size/2), dtype=np.intp)
        else:
            start_range = int(np.floor(1/fmax*fs))
            search_sample = np.arange(start_range,int(window_size/2), dtype=np.intp)
    else:
        if fmax is None:
            end_range = int(np.ceil(1/fmin*fs))
            search_sample = np.arange(end_range+1, dtype=np.intp)
        else:
            start_range = int(np.floor(1/fmax*fs))
            end_range = int(np.ceil(1/fmin*fs))
            search_sample = np.arange(start_range,end_range+1, dtype=np.intp)
    
        if end_range > window_size:
            raise Exception(f'fmin must be lager than {fs/window_size} Hz')
        elif end_range > (window_size//2):
            warnings.warn(f'fmin might be lager than {fs/(window_size//2)} Hz')

    # 処理実行
    if pp_mode == 'positive+negative':
        if implementation == 'Cython':
            Pitch = calc_Pitch_negaposi_core_cy(data_pos,data_neg,
                                                fs,
                                                window_size,
                                                hop_size,
                                                search_sample,
                                                pp_mode,
                                                bedcmm_smooth,
                                                pitch_detect_mode,
                                                pitch_detect_thre,
                                                interpolator_mode)
        else:
            Pitch = calc_Pitch_negaposi_core(data_pos,data_neg,
                                            fs,
                                            window_size,
                                            hop_size,
                                            search_sample,
                                            pp_mode,
                                            bedcmm_smooth,
                                            pitch_detect_mode,
                                            pitch_detect_thre,
                                            interpolator_mode)
    else:
        if implementation == 'Cython':
            Pitch = calc_Pitch_core_cy(data,
                                       fs,
                                       window_size,
                                       hop_size,
                                       search_sample,
                                       pp_mode,
                                       bedcmm_smooth,
                                       pitch_detect_mode,
                                       pitch_detect_thre,
                                       interpolator_mode)
        else:
            Pitch = calc_Pitch_core(data,
                                    fs,
                                    window_size,
                                    hop_size,
                                    search_sample,
                                    pp_mode,
                                    bedcmm_smooth,
                                    pitch_detect_mode,
                                    pitch_detect_thre,
                                    interpolator_mode)

    Pitch = np.array(Pitch)
    if Pitch.size == 0:
        return np.array([]), np.array([])

    Pitch_data = Pitch[:,0]
    Pitch_score = Pitch[:,1]

    return Pitch_data,Pitch_score

def calc_bedcmm(data,
                fs=44100,
                window_size=2048,
                hop_size=256,
                fmax=None,
                fmin=None,
                pp_mode='positive+negative',
                pp_threshold=0):

    data = data.copy()
    data = np.ascontiguousarray(data, dtype=np.float64)

    # データ前処理
    if pp_mode == 'positive':
        data[data < 0] = 0
    elif pp_mode == 'negative':
        data[data > 0] = 0
        data[data < 0] = -data[data < 0]
    elif pp_mode == 'positive+negative':
        data_pos = np.zeros_like(data)
        data_neg = np.zeros_like(data)
        data_pos[data > 0] = data[data > 0]
        data_neg[data < 0] = -data[data < 0]
    elif pp_mode == 'threshold_diff':
        data = data - pp_threshold
    else:
        raise Exception('pp_mode is only positive,negative,positive+negative,threshold_diff.')

    if fmin is None:
        if fmax is None:
            search_sample = np.arange(int(window_size/2), dtype=np.intp)
        else:
            start_range = int(np.floor(1/fmax*fs))
            search_sample = np.arange(start_range,int(window_size/2), dtype=np.intp)
    else:
        if fmax is None:
            end_range = int(np.ceil(1/fmin*fs))
            search_sample = np.arange(end_range+1, dtype=np.intp)
        else:
            start_range = int(np.floor(1/fmax*fs))
            end_range = int(np.ceil(1/fmin*fs))
            search_sample = np.arange(start_range,end_range+1, dtype=np.intp)

        if end_range > window_size:
            raise Exception(f'fmin must be lager than {fs/window_size} Hz')
        elif end_range > (window_size//2):
            warnings.warn(f'fmin might be lager than {fs/(window_size//2)} Hz')

    if pp_mode == 'positive+negative':
        if implementation == 'Cython':
            result = calc_bedcmm_negaposi_core_cy(data_pos,
                                                         data_neg,
                                                         window_size,
                                                         hop_size,
                                                         search_sample)

            bedcmm_result = result[0]
            mean_data = result[1]
        else:
            bedcmm_result,mean_data = calc_bedcmm_negaposi_core(data_pos,
                                                      data_neg,
                                                      window_size,
                                                      hop_size,
                                                      search_sample)
    else:
        if implementation == 'Cython':
            result = calc_bedcmm_core_cy(data,
                                                window_size,
                                                hop_size,
                                                search_sample)
            bedcmm_result = result[0]
            mean_data = result[1]
        else:
            bedcmm_result,mean_data = calc_bedcmm_core(data,
                                             window_size,
                                             hop_size,
                                             search_sample)
    return bedcmm_result,mean_data


def calc_bedcmm_core(data,
                     window_size,
                     hop_size,
                     search_sample):
 
    bedcmm_result_list = []
    mean_data_list = []
    for i in range(window_size, len(data),hop_size):
        calc_data = data[i-window_size:i]
        bedcmm_result_list.append(_periodicity(calc_data,search_sample))
        mean_data_list.append(np.mean(calc_data))

    bedcmm_result = np.array(bedcmm_result_list)
    mean_data = np.array(mean_data_list)

    return bedcmm_result,mean_data

def calc_bedcmm_negaposi_core(data_pos,
                              data_neg,
                              window_size,
                              hop_size,
                              search_sample):
 
    bedcmm_result_list = []
    mean_data_list = []
    for i in range(window_size, len(data_pos),hop_size):
        calc_data_posi = data_pos[i-window_size:i]
        calc_data_nega = data_neg[i-window_size:i]
        bedcmm_result_list.append(_periodicity(calc_data_posi,search_sample) + _periodicity(calc_data_nega,search_sample))
        mean_data_list.append(np.mean(calc_data_posi)+np.mean(calc_data_nega))

    bedcmm_result = np.array(bedcmm_result_list)
    mean_data = np.array(mean_data_list)

    return bedcmm_result,mean_data

def create_transition_matrix(
    n_periods,
    fs,
    sigma=0.1
):
    periods = np.arange(
        1,
        n_periods + 1
    )

    freqs = fs / periods

    log_freqs = np.log2(freqs)

    diff = (
        log_freqs[:, None]
        - log_freqs[None, :]
    )

    A = np.exp(
        -(diff**2)
        /(2*sigma**2)
    )

    A /= np.sum(
        A,
        axis=0,
        keepdims=True
    )

    return A

def calc_Pitch_bayes_negaposi_core(data_posi,data_nega,
                                   fs,
                                   window_size,
                                   hop_size,
                                   search_sample,
                                   pp_mode,
                                   alpha,
                                   sigma,
                                   interpolator_mode):

    Pitch = []
    # 遷移行列の取得
    A = create_transition_matrix(
        window_size//2,
        fs=fs,
        sigma=sigma
    )        
    A = A[search_sample[0]:search_sample[-1]+1,search_sample[0]:search_sample[-1]+1]
    # ベイズ推定に必要な変数の初期化
    posterior = None

    for i in range(window_size, len(data_posi),hop_size):
        calc_data_posi = data_posi[i-window_size:i]
        calc_data_nega = data_nega[i-window_size:i]
        bedcmm_result = _periodicity(calc_data_posi,search_sample) + _periodicity(calc_data_nega,search_sample)
        mean_data = np.mean(calc_data_posi)+np.mean(calc_data_nega)

        priod_diff = np.diff(bedcmm_result)
        up_inds = np.where(priod_diff > 0 )[0]

        likelihoods = np.zeros_like(bedcmm_result)
        if len(up_inds) > 0:
            likelihoods[up_inds[0]:] = bedcmm_result[up_inds[0]:]
        else:
            likelihoods = np.ones_like(bedcmm_result)/len(bedcmm_result)

        likelihoods /= np.sum(likelihoods)
        if posterior is None:
            posterior = likelihoods
        else:
            prediction = A @ posterior
            posterior = prediction**alpha * likelihoods
        
        posterior /= np.sum(posterior)

        # MAP推定
        max_idx_int = np.argmax(posterior)

        if ~np.isnan(max_idx_int):
            if max_idx_int != 0:
                if interpolator_mode == 'parabolic':
                    delta,peak_value = _parabolic_peak(bedcmm_result,max_idx_int)
                elif interpolator_mode == 'gaussian':
                    if pp_mode == 'threshold_diff':
                        bedcmm_result = bedcmm_result - min(bedcmm_result)
                    delta,peak_value = _gaussian_peak(bedcmm_result,max_idx_int)
                elif interpolator_mode == 'centroid':
                    if pp_mode == 'threshold_diff':
                        bedcmm_result = bedcmm_result - min(bedcmm_result)
                    delta,peak_value = _centroid_peak(bedcmm_result,max_idx_int)
                elif interpolator_mode == 'no':
                    peak_value = bedcmm_result[max_idx_int]
                    delta = 0
                else:
                    raise Exception('interpolator_mode is quadratic,centroid,gaussian or no')

                if delta < -0.5:
                    delta = -0.5
                if delta > 0.5:
                    delta = 0.5

                peak_idx = search_sample[max_idx_int]+delta

                peak_score = peak_value/mean_data
                peak_posister = posterior[max_idx_int]
            else:
                peak_idx = np.nan
                peak_score =np.nan
                peak_posister = np.nan
        else:
            peak_idx = np.nan
            peak_score = np.nan
            peak_posister = np.nan

        if np.isnan(peak_idx):
            Pitch.append([np.nan,np.nan,np.nan])
        else:
            Pitch.append([fs/peak_idx,peak_score,peak_posister])

    Pitch = np.array(Pitch)

    return Pitch


def calc_Pitch_bayes_core(data,
                          fs,
                          window_size,
                          hop_size,
                          search_sample,
                          pp_mode,
                          alpha,
                          sigma,
                          interpolator_mode):

    Pitch = []
    # 遷移行列の取得
    A = create_transition_matrix(
        window_size,
        fs=fs,
        sigma=sigma
    )        
    A = A[search_sample[0]:search_sample[-1]+1,search_sample[0]:search_sample[-1]+1]
    # ベイズ推定に必要な変数の初期化
    posterior = None

    for i in range(window_size, len(data),hop_size):
        calc_data = data[i-window_size:i]
        bedcmm_result = _periodicity(calc_data,search_sample)
        mean_data = np.mean(calc_data)

        priod_diff = np.diff(bedcmm_result)
        up_inds = np.where(priod_diff > 0 )[0]

        likelihoods = np.zeros_like(bedcmm_result)
        if len(up_inds) > 0:
            likelihoods[up_inds[0]:] = bedcmm_result[up_inds[0]:]
        else:
            likelihoods = np.ones_like(bedcmm_result)/len(bedcmm_result)

        likelihoods /= np.sum(likelihoods)
        if posterior is None:
            posterior = likelihoods
        else:
            prediction = A @ posterior
            posterior = prediction**alpha * likelihoods
        
        posterior /= np.sum(posterior)

        # MAP推定
        max_idx_int = np.argmax(posterior)

        if ~np.isnan(max_idx_int):
            if max_idx_int != 0:
                if interpolator_mode == 'parabolic':
                    delta,peak_value = _parabolic_peak(bedcmm_result,max_idx_int)
                elif interpolator_mode == 'gaussian':
                    if pp_mode == 'threshold_diff':
                        bedcmm_result = bedcmm_result - min(bedcmm_result)
                    delta,peak_value = _gaussian_peak(bedcmm_result,max_idx_int)
                elif interpolator_mode == 'centroid':
                    if pp_mode == 'threshold_diff':
                        bedcmm_result = bedcmm_result - min(bedcmm_result)
                    delta,peak_value = _centroid_peak(bedcmm_result,max_idx_int)
                elif interpolator_mode == 'no':
                    peak_value = bedcmm_result[max_idx_int]
                    delta = 0
                else:
                    raise Exception('interpolator_mode is quadratic,centroid,gaussian or no')

                if delta < -0.5:
                    delta = -0.5
                if delta > 0.5:
                    delta = 0.5

                peak_idx = search_sample[max_idx_int]+delta

                peak_score = peak_value/mean_data
                peak_posister = posterior[max_idx_int]
            else:
                peak_idx = np.nan
                peak_score =np.nan
                peak_posister = np.nan
        else:
            peak_idx = np.nan
            peak_score = np.nan
            peak_posister = np.nan

        if np.isnan(peak_idx):
            Pitch.append([np.nan,np.nan,np.nan])
        else:
            Pitch.append([fs/peak_idx,peak_score,peak_posister])

    Pitch = np.array(Pitch)

    return Pitch


def calc_Pitch_with_bayes(data,
                          fs=44100,
                          window_size=2048,
                          hop_size=256,
                          fmin=65,
                          fmax=2000,
                          pp_mode='positive+negative',
                          pp_threshold=0,
                          alpha = 0.7,
                          sigma = 0.1,
                          interpolator_mode='parabolic'):
    
    data = data.copy()
    data = np.ascontiguousarray(data, dtype=np.float64)

    if data.ndim != 1:
        raise Exception('data must be 1D array.')

    # データ前処理
    if pp_mode == 'positive':
        data[data < 0] = 0
    elif pp_mode == 'negative':
        data[data > 0] = 0
        data[data < 0] = -data[data < 0]
    elif pp_mode == 'positive+negative':
        data_pos = np.zeros_like(data)
        data_neg = np.zeros_like(data)
        data_pos[data > 0] = data[data > 0]
        data_neg[data < 0] = -data[data < 0]
    elif pp_mode == 'threshold_diff':
        data = data - pp_threshold
    else:
        raise Exception('pp_mode is only positive,negative,positive+negative,threshold_diff.')

    if fmin is None:
        if fmax is None:
            search_sample = np.arange(int(window_size/2), dtype=np.intp)
        else:
            start_range = int(np.floor(1/fmax*fs))
            search_sample = np.arange(start_range,int(window_size/2), dtype=np.intp)
    else:
        if fmax is None:
            end_range = int(np.ceil(1/fmin*fs))
            search_sample = np.arange(end_range+1, dtype=np.intp)
        else:
            start_range = int(np.floor(1/fmax*fs))
            end_range = int(np.ceil(1/fmin*fs))
            search_sample = np.arange(start_range,end_range+1, dtype=np.intp)
    
        if end_range > (window_size//2):
            raise Exception(f'fmin must be lager than {fs/(window_size//2)} Hz')


    # 処理実行
    if pp_mode == 'positive+negative':
        if implementation == 'Cython':
            Pitch = calc_Pitch_bayes_negaposi_core_cy(data_pos,data_neg,
                                                      fs,
                                                      window_size,
                                                      hop_size,
                                                      search_sample,
                                                      pp_mode,
                                                      alpha,
                                                      sigma,
                                                      interpolator_mode)
        else:
            Pitch = calc_Pitch_bayes_negaposi_core(data_pos,data_neg,
                                                   fs,
                                                   window_size,
                                                   hop_size,
                                                   search_sample,
                                                   pp_mode,
                                                   alpha,
                                                   sigma,
                                                   interpolator_mode)
    else:
        if implementation == 'Cython':
            Pitch = calc_Pitch_bayes_core_cy(data,
                                             fs,
                                             window_size,
                                             hop_size,
                                             search_sample,
                                             pp_mode,
                                             alpha,
                                             sigma,
                                             interpolator_mode)
        else:
            Pitch = calc_Pitch_bayes_core(data,
                                          fs,
                                          window_size,
                                          hop_size,
                                          search_sample,
                                          pp_mode,
                                          alpha,
                                          sigma,
                                          interpolator_mode)

    Pitch = np.array(Pitch)
    if Pitch.size == 0:
        return np.array([]), np.array([]),np.array([])

    Pitch_data = Pitch[:,0]
    Pitch_score = Pitch[:,1]
    Pitch_prob = Pitch[:,2]

    return Pitch_data,Pitch_score,Pitch_prob


def viterbi(log_likelihood, log_A, beta):

    n_frames, n_states = log_likelihood.shape

    score = np.zeros(
        (n_frames, n_states)
    )

    back_ptr = np.zeros(
        (n_frames, n_states),
        dtype=np.int32
    )

    # 初期化
    score[0] = log_likelihood[0]

    # Forward
    for t in range(1, n_frames):

        for j in range(n_states):

            candidates = (
                score[t-1]
                + log_A[:, j]
            )

            best_prev = np.argmax(
                candidates
            )

            score[t, j] = (
                candidates[best_prev]
                + beta * log_likelihood[t, j]
            )

            back_ptr[t, j] = best_prev

    # 終端
    path = np.zeros(
        n_frames,
        dtype=np.int32
    )

    path[-1] = np.argmax(
        score[-1]
    )

    # Backtracking
    for t in range(
        n_frames-2,
        -1,
        -1
    ):
        path[t] = back_ptr[
            t+1,
            path[t+1]
        ]

    return path

def calc_Pitch_viterbi_negaposi_core(data_posi,data_nega,
                                     fs,
                                     window_size,
                                     hop_size,
                                     search_sample,
                                     pp_mode,
                                     beta,
                                     sigma,
                                     interpolator_mode):

    # 遷移行列の取得
    A = create_transition_matrix(
        window_size//2,
        fs=fs,
        sigma=sigma
    )        
    A = A[search_sample[0]:search_sample[-1]+1,search_sample[0]:search_sample[-1]+1]
    log_A = np.log(A + 1e-300)
    likelihood_list = []
    bedcmm_result_list = []
    mean_data_list = []

    for i in range(window_size, len(data_posi),hop_size):
        calc_data_posi = data_posi[i-window_size:i]
        calc_data_nega = data_nega[i-window_size:i]
        bedcmm_result = _periodicity(calc_data_posi,search_sample) + _periodicity(calc_data_nega,search_sample)
        bedcmm_result_list.append(bedcmm_result)
        mean_data = np.mean(calc_data_posi)+np.mean(calc_data_nega)
        mean_data_list.append(mean_data)

        priod_diff = np.diff(bedcmm_result)
        up_inds = np.where(priod_diff > 0 )[0]

        likelihoods = np.zeros_like(bedcmm_result)
        if len(up_inds) > 0:
            likelihoods[up_inds[0]:] = bedcmm_result[up_inds[0]:]
        else:
            likelihoods = np.ones_like(bedcmm_result)/len(bedcmm_result)

        likelihoods /= np.sum(likelihoods)
        likelihood_list.append(likelihoods)

    log_likelihood_list = np.log(np.array(likelihood_list)+ 1e-300)
    path_list = viterbi(log_likelihood_list,log_A,beta=beta)

    Pitch = []
    for frame_num,max_idx_int in enumerate(path_list):
        if ~np.isnan(max_idx_int):
            bedcmm_result = bedcmm_result_list[frame_num]
            if max_idx_int != 0:
                if interpolator_mode == 'parabolic':
                    delta,peak_value = _parabolic_peak(bedcmm_result,max_idx_int)
                elif interpolator_mode == 'gaussian':
                    if pp_mode == 'threshold_diff':
                        bedcmm_result = bedcmm_result - min(bedcmm_result)
                    delta,peak_value = _gaussian_peak(bedcmm_result,max_idx_int)
                elif interpolator_mode == 'centroid':
                    if pp_mode == 'threshold_diff':
                        bedcmm_result = bedcmm_result - min(bedcmm_result)
                    delta,peak_value = _centroid_peak(bedcmm_result,max_idx_int)
                elif interpolator_mode == 'no':
                    peak_value = bedcmm_result[max_idx_int]
                    delta = 0
                else:
                    raise Exception('interpolator_mode is quadratic,centroid,gaussian or no')

                if delta < -0.5:
                    delta = -0.5
                if delta > 0.5:
                    delta = 0.5

                peak_idx = search_sample[max_idx_int]+delta
                peak_score = peak_value/mean_data_list[frame_num]
                peak_likelihood= likelihood_list[frame_num][max_idx_int]
            else:
                peak_idx = np.nan
                peak_score =np.nan
                peak_likelihood = np.nan
        else:
            peak_idx = np.nan
            peak_score = np.nan
            peak_likelihood = np.nan

        if np.isnan(peak_idx):
            Pitch.append([np.nan,np.nan,np.nan])
        else:
            Pitch.append([fs/peak_idx,peak_score,peak_likelihood])

    Pitch = np.array(Pitch)

    return Pitch


def calc_Pitch_viterbi_core(data,
                            fs,
                            window_size,
                            hop_size,
                            search_sample,
                            pp_mode,
                            beta,
                            sigma,
                            interpolator_mode):

    # 遷移行列の取得
    A = create_transition_matrix(
        window_size//2,
        fs=fs,
        sigma=sigma
    )        
    A = A[search_sample[0]:search_sample[-1]+1,search_sample[0]:search_sample[-1]+1]
    log_A = np.log(A + 1e-300)
    likelihood_list = []
    bedcmm_result_list = []
    mean_data_list = []

    for i in range(window_size, len(data),hop_size):
        calc_data = data[i-window_size:i]
        bedcmm_result = _periodicity(calc_data,search_sample)
        bedcmm_result_list.append(bedcmm_result)
        mean_data = np.mean(calc_data)
        mean_data_list.append(mean_data)

        priod_diff = np.diff(bedcmm_result)
        up_inds = np.where(priod_diff > 0 )[0]

        likelihoods = np.zeros_like(bedcmm_result)
        if len(up_inds) > 0:
            likelihoods[up_inds[0]:] = bedcmm_result[up_inds[0]:]
        else:
            likelihoods = np.ones_like(bedcmm_result)/len(bedcmm_result)

        likelihoods /= np.sum(likelihoods)
        likelihood_list.append(likelihoods)

    log_likelihood_list = np.log(np.array(likelihood_list)+ 1e-300)
    path_list = viterbi(log_likelihood_list,log_A,beta=beta)

    Pitch = []
    for frame_num,max_idx_int in enumerate(path_list):
        if ~np.isnan(max_idx_int):
            bedcmm_result = bedcmm_result_list[frame_num]
            if max_idx_int != 0:
                if interpolator_mode == 'parabolic':
                    delta,peak_value = _parabolic_peak(bedcmm_result,max_idx_int)
                elif interpolator_mode == 'gaussian':
                    if pp_mode == 'threshold_diff':
                        bedcmm_result = bedcmm_result - min(bedcmm_result)
                    delta,peak_value = _gaussian_peak(bedcmm_result,max_idx_int)
                elif interpolator_mode == 'centroid':
                    if pp_mode == 'threshold_diff':
                        bedcmm_result = bedcmm_result - min(bedcmm_result)
                    delta,peak_value = _centroid_peak(bedcmm_result,max_idx_int)
                elif interpolator_mode == 'no':
                    peak_value = bedcmm_result[max_idx_int]
                    delta = 0
                else:
                    raise Exception('interpolator_mode is quadratic,centroid,gaussian or no')

                if delta < -0.5:
                    delta = -0.5
                if delta > 0.5:
                    delta = 0.5

                peak_idx = search_sample[max_idx_int]+delta

                peak_score = peak_value/mean_data_list[frame_num]
                peak_likelihood= likelihood_list[frame_num][max_idx_int]
            else:
                peak_idx = np.nan
                peak_score =np.nan
                peak_likelihood = np.nan
        else:
            peak_idx = np.nan
            peak_score = np.nan
            peak_likelihood = np.nan

        if np.isnan(peak_idx):
            Pitch.append([np.nan,np.nan,np.nan])
        else:
            Pitch.append([fs/peak_idx,peak_score,peak_likelihood])

    Pitch = np.array(Pitch)

    return Pitch



def calc_Pitch_with_viterbi(data,
                            fs=44100,
                            window_size=2048,
                            hop_size=256,
                            fmin=65,
                            fmax=2000,
                            pp_mode='positive+negative',
                            pp_threshold=0,
                            beta = 10,
                            sigma = 0.1,
                            interpolator_mode='parabolic'):
    
    data = data.copy()
    data = np.ascontiguousarray(data, dtype=np.float64)

    if data.ndim != 1:
        raise Exception('data must be 1D array.')

    # データ前処理
    if pp_mode == 'positive':
        data[data < 0] = 0
    elif pp_mode == 'negative':
        data[data > 0] = 0
        data[data < 0] = -data[data < 0]
    elif pp_mode == 'positive+negative':
        data_pos = np.zeros_like(data)
        data_neg = np.zeros_like(data)
        data_pos[data > 0] = data[data > 0]
        data_neg[data < 0] = -data[data < 0]
    elif pp_mode == 'threshold_diff':
        data = data - pp_threshold
    else:
        raise Exception('pp_mode is only positive,negative,positive+negative,threshold_diff.')

    if fmin is None:
        if fmax is None:
            search_sample = np.arange(int(window_size/2), dtype=np.intp)
        else:
            start_range = int(np.floor(1/fmax*fs))
            search_sample = np.arange(start_range,int(window_size/2), dtype=np.intp)
    else:
        if fmax is None:
            end_range = int(np.ceil(1/fmin*fs))
            search_sample = np.arange(end_range+1, dtype=np.intp)
        else:
            start_range = int(np.floor(1/fmax*fs))
            end_range = int(np.ceil(1/fmin*fs))
            search_sample = np.arange(start_range,end_range+1, dtype=np.intp)
    
        if end_range > (window_size//2):
            raise Exception(f'fmin must be lager than {fs/(window_size//2)} Hz')

    # 処理実行
    if pp_mode == 'positive+negative':
        if implementation == 'Cython':
            Pitch = calc_Pitch_viterbi_negaposi_core_cy(data_pos,data_neg,
                                                      fs,
                                                      window_size,
                                                      hop_size,
                                                      search_sample,
                                                      pp_mode,
                                                      beta,
                                                      sigma,
                                                      interpolator_mode)
        else:
            Pitch = calc_Pitch_viterbi_negaposi_core(data_pos,data_neg,
                                                   fs,
                                                   window_size,
                                                   hop_size,
                                                   search_sample,
                                                   pp_mode,
                                                   beta,
                                                   sigma,
                                                   interpolator_mode)
    else:
        if implementation == 'Cython':
            Pitch = calc_Pitch_viterbi_core_cy(data,
                                             fs,
                                             window_size,
                                             hop_size,
                                             search_sample,
                                             pp_mode,
                                             beta,
                                             sigma,
                                             interpolator_mode)
        else:
            Pitch = calc_Pitch_viterbi_core(data,
                                          fs,
                                          window_size,
                                          hop_size,
                                          search_sample,
                                          pp_mode,
                                          beta,
                                          sigma,
                                          interpolator_mode)

    Pitch = np.array(Pitch)
    
    if Pitch.size == 0:
        return np.array([]), np.array([]),np.array([])

    Pitch_data = Pitch[:,0]
    Pitch_score = Pitch[:,1]
    Pitch_prob = Pitch[:,2]

    return Pitch_data,Pitch_score,Pitch_prob


def main():
    pass

if __name__ == "__main__":
    main()
