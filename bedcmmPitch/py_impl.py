# -*- coding: utf-8 -*-
"""
bedcmmPitch method class(Python)

Author: YASUHARA Wataru
Copyright (c) 2026, Feel a Piece of the World
"""
import numpy as np
from ._config import implementation
if implementation == 'Cython':
    from .cy_impl import calc_Pitch_core_cy,calc_Pitch_negaposi_core_cy,calc_bedcmm_core_cy,calc_bedcmm_negaposi_core_cy

def _periodicity(data,period):

    result = np.zeros(len(period))
    for p_idx,a_preiod in enumerate(period):

        temp_data = np.zeros(data.shape[0] - a_preiod)
        for index in range(data.shape[0] - a_preiod):
            temp_data[index] = min([data[index],data[index+a_preiod]])
                
        result[p_idx]=np.mean((temp_data))

    return result    

def _peak_detect_threshould(bedcmm_result,threshould,search_sample,bedcmm_smooth,interpolator_mode):
    priod_thre_ind = (bedcmm_result > threshould)[:-2]
    priod_diff = np.diff(bedcmm_result)
    plus_peaks = np.where((priod_diff[:-1] > 0) & (priod_diff[1:] * priod_diff[:-1] < 0 ) & priod_thre_ind)[0]
    if len(plus_peaks) > 0:
        plus_peak_idx = plus_peaks[0]+1
        plus_peak = int(search_sample[plus_peak_idx])

        if interpolator_mode == 'parabolic':
            c = np.polyfit(np.arange(plus_peak-1,plus_peak+2),bedcmm_result[plus_peak_idx-1:plus_peak_idx+2],2)
            max_idx = -c[1]/(2*c[0]) + ((bedcmm_smooth-1)/2)
        elif interpolator_mode == 'no':
            max_idx = plus_peak + ((bedcmm_smooth-1)/2)
        else:
            raise Exception('interpolator_mode is parabolic or no')
    else:
        max_idx = np.nan
        
    return max_idx


def _peak_detect_maximum(bedcmm_result,search_sample,bedcmm_smooth,interpolator_mode):

    plus_peak_idx = np.argmax(bedcmm_result)
    plus_peak = search_sample[plus_peak_idx]
    if plus_peak != search_sample[0]:
        if interpolator_mode == 'parabolic':
            c = np.polyfit(np.arange(plus_peak-1,plus_peak+2),bedcmm_result[plus_peak_idx-1:plus_peak_idx+2],2)
            max_idx = -c[1]/(2*c[0]) + ((bedcmm_smooth-1)/2)
        elif interpolator_mode == 'no':
            max_idx = plus_peak + ((bedcmm_smooth-1)/2)
        else:
            raise Exception('interpolator_mode is quadratic or no')
    else:
        max_idx = np.nan
        
    return max_idx


def calc_Pitch_core(data,
                    fs,
                    window_size,
                    hop_size,
                    search_sample,
                    bedcmm_smooth,
                    pitch_detect_mode,
                    pitch_detect_thre,
                    interpolator_mode):
    
    Pitch = []
    for i in range(window_size, len(data),hop_size):
        calc_data = data[i-window_size:i]
        bedcmm_result = _periodicity(calc_data,search_sample)

        if bedcmm_smooth > 1:
            filt = np.ones(bedcmm_smooth)/bedcmm_smooth
            bedcmm_result = np.convolve(bedcmm_result,filt,mode='valid')
        elif bedcmm_smooth == 1:
            bedcmm_result = bedcmm_result
        else:
            raise Exception('bedcmm_smooth > 0 and int')
        
        if pitch_detect_mode == 'dynamic':
            threshould = np.mean(calc_data)*pitch_detect_thre
            max_idx = _peak_detect_threshould(bedcmm_result,threshould,search_sample,bedcmm_smooth,interpolator_mode)
        elif pitch_detect_mode == 'static':
            max_idx = _peak_detect_threshould(bedcmm_result,pitch_detect_thre,search_sample,bedcmm_smooth,interpolator_mode)
        elif pitch_detect_mode == 'maximum':
            max_idx = _peak_detect_maximum(bedcmm_result,search_sample,bedcmm_smooth,interpolator_mode)

        if np.isnan(max_idx):
            Pitch.append(np.nan)
        else:
            Pitch.append(fs/max_idx)

    Pitch = np.array(Pitch)

    return Pitch

def calc_Pitch_negaposi_core(data_posi,data_nega,
                             fs,
                             window_size,
                             hop_size,
                             search_sample,
                             bedcmm_smooth,
                             pitch_detect_mode,
                             pitch_detect_thre,
                             interpolator_mode):

    Pitch = []
    for i in range(window_size, len(data_posi),hop_size):
        calc_data_posi = data_posi[i-window_size:i]
        calc_data_nega = data_nega[i-window_size:i]
        bedcmm_result = _periodicity(calc_data_posi,search_sample) + _periodicity(calc_data_nega,search_sample)

        if bedcmm_smooth != 0:
            filt = np.ones(bedcmm_smooth)/bedcmm_smooth
            bedcmm_result = np.convolve(bedcmm_result,filt,mode='valid')
        
        if pitch_detect_mode == 'dynamic':
            threshould = (np.mean(calc_data_posi) + np.mean(calc_data_nega))*pitch_detect_thre
            max_idx = _peak_detect_threshould(bedcmm_result,threshould,search_sample,bedcmm_smooth,interpolator_mode)
        elif pitch_detect_mode == 'static':
            max_idx = _peak_detect_threshould(bedcmm_result,pitch_detect_thre,search_sample,bedcmm_smooth,interpolator_mode)
        elif pitch_detect_mode == 'maximum':
            max_idx = _peak_detect_maximum(bedcmm_result,bedcmm_smooth,interpolator_mode)

        if np.isnan(max_idx):
            Pitch.append(np.nan)
        else:
            Pitch.append(fs/max_idx)

    Pitch = np.array(Pitch)

    return Pitch


def calc_Pitch(data,
               fs=44100,
               window_size=2048,
               hop_size=256,
               pitch_range=None,
               pp_mode='positive',
               pp_threshould=0,
               bedcmm_smooth=3,
               pitch_detect_mode='dynamic',
               pitch_detect_thre=0.33333,
               interpolator_mode='parabolic'):
    
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
    elif pp_mode == 'threshould_diff':
        data = data - pp_threshould
    else:
        raise Exception('pp_mode is only positive,negative,positive+negative,threshould_diff.')

    if pitch_range is None:
        search_sample = np.arange(int(window_size/2))
    else:
        if len(pitch_range) != 2:
            raise Exception('pitch_rage is [start freq(low), end freq(high)].')
        start_range = int(np.floor(1/pitch_range[1]*fs))
        end_range = int(np.ceil(1/pitch_range[0]*fs))
        search_sample = np.arange(start_range,end_range+1)

    # 処理実行
    if pp_mode == 'positive+negative':
        if implementation == 'Cython':
            Pitch = calc_Pitch_negaposi_core_cy(data_pos,data_neg,
                                                fs,
                                                window_size,
                                                hop_size,
                                                search_sample,
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
                                    bedcmm_smooth,
                                    pitch_detect_mode,
                                    pitch_detect_thre,
                                    interpolator_mode)

    return Pitch

def calc_bedcmm(data,
                fs=44100,
                window_size=2048,
                hop_size=256,
                pitch_range=None,
                pp_mode='positive',
                pp_threshould=0):

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
    elif pp_mode == 'threshould_diff':
        data = data - pp_threshould
    else:
        raise Exception('pp_mode is only positive,negative,positive+negative,threshould_diff.')

    if pitch_range is None:
        search_sample = np.arange(int(window_size/2))
    else:
        if len(pitch_range) != 2:
            raise Exception('pitch_rage is [start freq(low), end freq(high)].')
        start_range = int(np.floor(1/pitch_range[1]*fs))
        end_range = int(np.ceil(1/pitch_range[0]*fs))
        search_sample = np.arange(start_range,end_range+1)

    if pp_mode == 'positive+negative':
        if implementation == 'Cython':
            bedcmm_result = calc_bedcmm_negaposi_core_cy(data,
                                                         window_size,
                                                         hop_size,
                                                         search_sample)
        else:
            bedcmm_result = calc_bedcmm_negaposi_core(data,
                                                      window_size,
                                                      hop_size,
                                                      search_sample)
    else:
        if implementation == 'Cython':
            bedcmm_result = calc_bedcmm_core_cy(data,
                                                window_size,
                                                hop_size,
                                                search_sample)
        else:
            bedcmm_result = calc_bedcmm_core(data,
                                             window_size,
                                             hop_size,
                                             search_sample)
    return bedcmm_result


def calc_bedcmm_core(data,
                     window_size,
                     hop_size,
                     search_sample):
 
    bedcmm_result_list = []
    for i in range(window_size, len(data),hop_size):
        calc_data = data[i-window_size:i]
        bedcmm_result_list.append(_periodicity(calc_data,search_sample))

    bedcmm_result = np.array(bedcmm_result_list)

    return bedcmm_result

def calc_bedcmm_negaposi_core(data_pos,
                              data_neg,
                              window_size,
                              hop_size,
                              search_sample):
 
    bedcmm_result_list = []
    for i in range(window_size, len(data_pos),hop_size):
        calc_data_posi = data_pos[i-window_size:i]
        calc_data_nega = data_neg[i-window_size:i]
        bedcmm_result_list.append(_periodicity(calc_data_posi,search_sample) + _periodicity(calc_data_nega,search_sample))

    bedcmm_result = np.array(bedcmm_result_list)

    return bedcmm_result


def main():
    pass

if __name__ == "__main__":
    main()
