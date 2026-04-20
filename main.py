import numpy as np
import matplotlib.pyplot as plt
import bedcmmPitch

def plot_pitch(t,times,Pitch,f_start,f_end):
    plt.figure()
    plt.plot(t[times], Pitch, label='bedcmmPitch',c='c',)
    plt.axhline(f_start, color='red', alpha=0.5, label='Actual Start Freqency')
    plt.axhline(f_end, color='red', alpha=0.5, label='Actual End Freqency')
    plt.xlabel('Time[s]')
    plt.ylabel('Freqency[Hz]')
    plt.show()
    pass

def main():
    #TODO:bedcmm scoreの導入
    print("Hello from bedcmmPitch!")

    fs = 44100  # サンプリングレート
    duration = 0.2  # 秒
    change_time = 0.1
    amp = 10
    t = np.linspace(0, duration, int(fs * duration), endpoint=False)
    #f_start, f_end = 440.0, 466.16
    f_start, f_end = 440.0, 880.0
    pitch_range = [380,1000]
    signal = np.where(t < change_time, amp*np.sin(2 * np.pi * f_start * t), amp*np.sin(2 * np.pi * f_end * t))

    window_size = 2048
    hop_size = 256
    times = np.arange(window_size, len(signal), hop_size)

    # default check
    print('default check')
    Pitch = bedcmmPitch.calc_Pitch(signal,
                                   window_size=window_size,
                                   hop_size=hop_size)
    plot_pitch(t,times,Pitch,f_start,f_end)
    # ppmode check
    print("ppmode check posi")
    Pitch = bedcmmPitch.calc_Pitch(signal,
                                   window_size=window_size,
                                   hop_size=hop_size,
                                   pitch_range=pitch_range,
                                   pp_mode='positive',
                                   pitch_detect_mode='dynamic',
                                   interpolator_mode='parabolic')

    plot_pitch(t,times,Pitch,f_start,f_end)

    print("ppmode check nega")
    Pitch = bedcmmPitch.calc_Pitch(signal,
                                   window_size=window_size,
                                   hop_size=hop_size,
                                   pitch_range=pitch_range,
                                   pp_mode='negative',
                                   pitch_detect_mode='dynamic',
                                   interpolator_mode='parabolic')
    plot_pitch(t,times,Pitch,f_start,f_end)


    print("ppmode check thrediff")
    Pitch = bedcmmPitch.calc_Pitch(signal,
                                   window_size=window_size,
                                   hop_size=hop_size,
                                   pitch_range=pitch_range,
                                   pp_mode='threshould_diff',
                                   pitch_detect_mode='maximum',
                                   interpolator_mode='parabolic')

    plot_pitch(t,times,Pitch,f_start,f_end)

    # pitch_range check
    print("pitch range check")
    Pitch = bedcmmPitch.calc_Pitch(signal,
                                   window_size=window_size,
                                   hop_size=hop_size,
                                   pitch_range=pitch_range,
                                   pp_mode='positive',
                                   pitch_detect_mode='maximum',
                                   interpolator_mode='parabolic')

    plot_pitch(t,times,Pitch,f_start,f_end)
    

    # detect mode check
    print("detect mode check")
    Pitch = bedcmmPitch.calc_Pitch(signal,
                                   window_size=window_size,
                                   hop_size=hop_size,
                                   pitch_range=pitch_range,
                                   pp_mode='positive',
                                   pitch_detect_mode='static',
                                   pitch_detect_thre=0.05,
                                   interpolator_mode='parabolic')
    plot_pitch(t,times,Pitch,f_start,f_end)
    Pitch = bedcmmPitch.calc_Pitch(signal,
                                   window_size=window_size,
                                   hop_size=hop_size,
                                   pitch_range=pitch_range,
                                   pp_mode='positive',
                                   pitch_detect_mode='maximum',
                                   interpolator_mode='parabolic')
    plot_pitch(t,times,Pitch,f_start,f_end)

    # interpolator mode check
    print("interpolator check")
    print("interpolator gaussian")
    Pitch = bedcmmPitch.calc_Pitch(signal,
                                   window_size=window_size,
                                   hop_size=hop_size,
                                   pitch_range=pitch_range,
                                   pp_mode='positive',
                                   pitch_detect_mode='dynamic',
                                   interpolator_mode='gaussian')
    plot_pitch(t,times,Pitch,f_start,f_end)
    print("interpolator threshould gaussian")
    Pitch = bedcmmPitch.calc_Pitch(signal,
                                   window_size=window_size,
                                   hop_size=hop_size,
                                   pitch_range=pitch_range,
                                   pp_mode='threshould_diff',
                                   pitch_detect_mode='dynamic',
                                   interpolator_mode='gaussian')
    print("interpolator centroid")
    Pitch = bedcmmPitch.calc_Pitch(signal,
                                   window_size=window_size,
                                   hop_size=hop_size,
                                   pitch_range=pitch_range,
                                   pp_mode='positive',
                                   pitch_detect_mode='dynamic',
                                   interpolator_mode='centroid')
    plot_pitch(t,times,Pitch,f_start,f_end)
    print("interpolator threshould gaussian")
    Pitch = bedcmmPitch.calc_Pitch(signal,
                                   window_size=window_size,
                                   hop_size=hop_size,
                                   pitch_range=pitch_range,
                                   pp_mode='threshould_diff',
                                   pitch_detect_mode='dynamic',
                                   interpolator_mode='centroid')

    plot_pitch(t,times,Pitch,f_start,f_end)
    print("interpolator no")
    Pitch = bedcmmPitch.calc_Pitch(signal,
                                   window_size=window_size,
                                   hop_size=hop_size,
                                   pitch_range=pitch_range,
                                   pp_mode='positive',
                                   pitch_detect_mode='dynamic',
                                   interpolator_mode='no')
    plot_pitch(t,times,Pitch,f_start,f_end)

    print("bedcmm result check")
    bedcmm_result = bedcmmPitch.calc_bedcmm(signal,
                                            window_size=window_size,
                                            hop_size=hop_size)
    plt.figure()
    plt.pcolormesh(t[times],np.arange(0,int(window_size/2)/fs,1/fs),bedcmm_result.T, shading="auto")
    plt.show()

    bedcmm_result = bedcmmPitch.calc_bedcmm(signal,
                                            window_size=window_size,
                                            hop_size=hop_size,
                                            pitch_range=pitch_range)
    start_range = int(np.floor(1/pitch_range[1]*fs))
    end_range = int(np.ceil(1/pitch_range[0]*fs))
    search_time = np.arange(start_range,end_range+1)/fs
    plt.figure()
    plt.pcolormesh(t[times],search_time,bedcmm_result.T, shading="auto")
    plt.show()
    
    pass
if __name__ == "__main__":
    main()
