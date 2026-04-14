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
    print("Hello from bedcmmpitch!")

    fs = 44100  # サンプリングレート
    duration = 0.2  # 秒
    change_time = 0.1
    amp = 1
    t = np.linspace(0, duration, int(fs * duration), endpoint=False)
    f_start, f_end = 440.0, 466.16
    signal = np.where(t < change_time, amp*np.sin(2 * np.pi * f_start * t), amp*np.sin(2 * np.pi * f_end * t))

    window_size = 2048
    hop_size = 256
    times = np.arange(window_size, len(signal), hop_size)

    # default check
    print('default check')
    Pitch = bedcmmPitch.calc_Pitch(signal)
    plot_pitch(t,times,Pitch,f_start,f_end)
    # ppmode check
    print("ppmode check posinega")
    Pitch = bedcmmPitch.calc_Pitch(signal,
                                   pitch_range=[400,500],
                                   pp_mode='positive+negative',
                                   pitch_detect_mode='dynamic',
                                   interpolator_mode='parabolic')

    plot_pitch(t,times,Pitch,f_start,f_end)

    print("ppmode check thrediff")
    Pitch = bedcmmPitch.calc_Pitch(signal,
                                   pitch_range=[380,500],
                                   pp_mode='threshould_diff',
                                   pitch_detect_mode='maximum',
                                   interpolator_mode='parabolic')

    plot_pitch(t,times,Pitch,f_start,f_end)

    # pitch_range check
    print("pitch range check")
    Pitch = bedcmmPitch.calc_Pitch(signal,
                                   pitch_range=[400,500],
                                   pp_mode='positive',
                                   pitch_detect_mode='maximum',
                                   interpolator_mode='parabolic')

    plot_pitch(t,times,Pitch,f_start,f_end)
    

    # detect mode check
    print("detect mode check")
    Pitch = bedcmmPitch.calc_Pitch(signal,
                                   pitch_range=[400,500],
                                   pp_mode='positive',
                                   pitch_detect_mode='static',
                                   pitch_detect_thre=0.05,
                                   interpolator_mode='parabolic')
    plot_pitch(t,times,Pitch,f_start,f_end)
    Pitch = bedcmmPitch.calc_Pitch(signal,
                                   pitch_range=[400,500],
                                   pp_mode='positive',
                                   pitch_detect_mode='maximum',
                                   interpolator_mode='parabolic')
    plot_pitch(t,times,Pitch,f_start,f_end)
    # interpolator mode check
    print("interpolator check")
    Pitch = bedcmmPitch.calc_Pitch(signal,
                                   pitch_range=[400,500],
                                   pp_mode='positive',
                                   pitch_detect_mode='dynamic',
                                   interpolator_mode='parabolic')
    plot_pitch(t,times,Pitch,f_start,f_end)
    Pitch = bedcmmPitch.calc_Pitch(signal,
                                   pitch_range=[400,500],
                                   pp_mode='positive',
                                   pitch_detect_mode='dynamic',
                                   interpolator_mode='no')
    plot_pitch(t,times,Pitch,f_start,f_end)

    pass
if __name__ == "__main__":
    main()
