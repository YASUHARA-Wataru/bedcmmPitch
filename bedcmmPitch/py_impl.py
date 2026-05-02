# -*- coding: utf-8 -*-
"""
bedcmmPitch method class(Python)

Author: YASUHARA Wataru
Copyright (c) 2026, Feel a Piece of the World
"""
import numpy as np
import math
from ._config import implementation
if implementation == 'Cython':
    from .cy_impl import calc_Pitch_core_cy,calc_Pitch_negaposi_core_cy,calc_bedcmm_core_cy,calc_bedcmm_negaposi_core_cy

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
            if max_idx_int != search_sample[0]:
                if interpolator_mode == 'parabolic':
                    delta,peak_value = _parabolic_peak(bedcmm_result,max_idx_int)
                    peak_idx = search_sample[max_idx_int]+delta
                elif interpolator_mode == 'gaussian':
                    if pp_mode == 'threshold_diff':
                        bedcmm_result = bedcmm_result - min(bedcmm_result)
                    delta,peak_value = _gaussian_peak(bedcmm_result,max_idx_int)
                    peak_idx = search_sample[max_idx_int]+delta
                elif interpolator_mode == 'centroid':
                    if pp_mode == 'threshold_diff':
                        bedcmm_result = bedcmm_result - min(bedcmm_result)
                    delta,peak_value = _centroid_peak(bedcmm_result,max_idx_int)
                    peak_idx = search_sample[max_idx_int]+delta
                elif interpolator_mode == 'no':
                    peak_value = bedcmm_result[max_idx_int]
                    peak_idx = float(search_sample[max_idx_int])
                else:
                    raise Exception('interpolator_mode is parabolic,centroid,gaussian or no')
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
            if max_idx_int != search_sample[0]:
                if interpolator_mode == 'parabolic':
                    delta,peak_value = _parabolic_peak(bedcmm_result,max_idx_int)
                    peak_idx = search_sample[max_idx_int]+delta
                elif interpolator_mode == 'gaussian':
                    if pp_mode == 'threshold_diff':
                        bedcmm_result = bedcmm_result - min(bedcmm_result)
                    delta,peak_value = _gaussian_peak(bedcmm_result,max_idx_int)
                    peak_idx = search_sample[max_idx_int]+delta
                elif interpolator_mode == 'centroid':
                    if pp_mode == 'threshold_diff':
                        bedcmm_result = bedcmm_result - min(bedcmm_result)
                    delta,peak_value = _centroid_peak(bedcmm_result,max_idx_int)
                    peak_idx = search_sample[max_idx_int]+delta
                elif interpolator_mode == 'no':
                    peak_value = bedcmm_result[max_idx_int]
                    peak_idx = float(search_sample[max_idx_int])
                else:
                    raise Exception('interpolator_mode is quadratic,centroid,gaussian or no')
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
               pitch_range=[65,2000],
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

    if pitch_range is None:
        search_sample = np.arange(int(window_size/2), dtype=np.intp)
    else:
        if len(pitch_range) != 2:
            raise Exception('pitch_rage is [min freq(low), max freq(high)].')
        start_range = int(np.floor(1/pitch_range[1]*fs))
        end_range = int(np.ceil(1/pitch_range[0]*fs))
        search_sample = np.arange(start_range,end_range+1, dtype=np.intp)

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

    Pitch_data = Pitch[:,0]
    Pitch_score = Pitch[:,1]

    return Pitch_data,Pitch_score

def calc_bedcmm(data,
                fs=44100,
                window_size=2048,
                hop_size=256,
                pitch_range=None,
                pp_mode='positive+negative',
                pp_threshold=0):

    data = data.copy()
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

    if pitch_range is None:
        search_sample = np.arange(int(window_size/2), dtype=np.intp)
    else:
        if len(pitch_range) != 2:
            raise Exception('pitch_rage is [start freq(low), end freq(high)].')
        start_range = int(np.floor(1/pitch_range[1]*fs))
        end_range = int(np.ceil(1/pitch_range[0]*fs))
        search_sample = np.arange(start_range,end_range+1, dtype=np.intp)

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


def main():
    pass

if __name__ == "__main__":
    main()
