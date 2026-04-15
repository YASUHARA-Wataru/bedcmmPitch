# bedcmmPitch
## Algorism of Pitch detection
This is Pitch detection algrorism using periodicity analysis.
They have two functions, calc_Pitch(calculate Pitch) and calc_bedcmm(periodicity analysis).

## 関数
- calc_Pitch
    - input(default) : datade description
        - data               : 1D array data(signal data)
        - fs(44100)          : float(sampling rate)
        - window_size(2048)  : int(window size)
        - hop_size(256)      : int(hop size)
        - pitch_range(None)  : [start_freq, end_freq](range of search)
        - pp_mode('positive'): str(perprocessing mode('positive','negative','positive+negative','threshould_diff'))
        - pp_threshould(0)   : float(perprocessing threshould(using in 'threshould_diff mode'))
        - bedcmm_smooth(3)   : int(bedcmm result smoothing size(if size 1 means do not smoothing))
        - pitch_detect_mode('dynamic'): str(bedcmm result peak detect mode(using in 'dynamic','static','maximum'))
        - pitch_detect_thre(0.33333): float(bedcmm result peak detect threshould(using in 'dynamic','static'))
        - interpolator_mode('parabolic'): str(peak index interpolator mode('parabolic','no'))
    - output
        - Pitch data:1D array data(pitch data)
- calc_bedcmm
    - input
        - data               : 1D array data(signal data)
        - fs(44100)          : float(sampling rate)
        - window_size(2048)  : int(window size)
        - hop_size(256)      : int(hop size)
        - pitch_range(None)  : [start_freq, end_freq](range of search)
        - pp_mode('positive'): str(perprocessing mode('positive','negative','positive+negative','threshould_diff'))
        - pp_threshould(0)   : float(perprocessing threshould(using in 'threshould_diff mode'))
    - output
        - bedcmm data:2D array data(time,bedcmm)

## Caution

## easy Install
1. 
2. 
3. 
    1. 
4. 
5. 

※4.
