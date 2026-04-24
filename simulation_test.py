import numpy as np
import matplotlib.pyplot as plt
import bedcmmPitch

def breathy_attack(
    fs=44100,
    duration=0.5,
    f0=200,
    attack_time=0.1,       # アタック時間
    noise_level=0.8,       # 初期ノイズ強度
    harmonics=5
):
    t = np.arange(0, duration, 1/fs)

    signal = np.zeros_like(t)

    # エンベロープ（ノイズ→周期）
    env_voiced = 1 / (1 + np.exp(-10 * (t - attack_time/2)))
    env_noise = 1 - env_voiced

    # 周波数の揺れ（アタック時のみ不安定）
    jitter = 1 + 0.05 * np.random.randn(len(t))
    f_inst = f0 * jitter

    phase = 2 * np.pi * np.cumsum(f_inst) / fs

    # 周期成分
    voiced = np.zeros_like(t)
    for k in range(1, harmonics+1):
        voiced += (1.0 / k) * np.sin(k * phase)

    # ノイズ成分（息）
    noise = np.random.randn(len(t))

    # 合成
    signal = env_voiced * voiced + env_noise * noise_level * noise

    return signal


def generate_signal(fs=44100, duration=2.0, f_base=220.0):
    t = np.arange(0, duration, 1/fs)

    # --- ピッチ揺らぎ ---
    f0 = f_base * (1 + 0.01 * np.sin(2*np.pi*5*t))
    phase = 2*np.pi * np.cumsum(f0) / fs

    # --- 倍音構成 ---
    signal = np.zeros_like(t)
    for k in range(1, 6):
        signal += (1.0/k) * np.sin(k * phase)

    # --- 振幅ゆらぎ ---
    envelope = 0.8 + 0.2 * np.sin(2*np.pi*0.5*t)
    signal *= envelope

    return signal


def add_spike_noise(signal, fs, spike_rate=50, spike_amp=5.0):
    noisy = signal.copy()
    n_samples = len(signal)

    n_spikes = int(spike_rate * (n_samples / fs))

    spike_positions = np.random.randint(0, n_samples, n_spikes)

    for pos in spike_positions:
        width = np.random.randint(1, 4)  # 1〜3サンプル
        amp = spike_amp * (0.5 + np.random.rand())

        noisy[pos:pos+width] += amp * (2*np.random.rand() - 1)

    return noisy


def add_white_noise(signal, snr_db=20):
    power = np.mean(signal**2)
    noise_power = power / (10**(snr_db/10))
    noise = np.sqrt(noise_power) * np.random.randn(len(signal))
    return signal + noise



def main():
    np.random.seed(101)
    #pattern = 'breathy'
    pattern = 'spike_noise'
    # --- 生成 ---
    fs = 44100
    if pattern == 'spike_noise':
        # スパイクノイズパターン
        sig = generate_signal(fs=fs)
        # ノイズ付加
        sig_spike = add_spike_noise(sig, fs, spike_rate=500)
        signal = add_white_noise(sig_spike, snr_db=10)
    elif pattern == 'breathy':
        # 吐息が混ざったパターン
        signal = breathy_attack(fs=fs)

    t = np.arange(0,len(signal)/fs,1/fs)

    # パラメータ設定
    window_size = 1024*2
    hop_size = 256
    times = np.arange(window_size, len(signal), hop_size)

    fo_bedcmm, bedcmm_score = bedcmmPitch.calc_Pitch(signal,pp_mode='positive+negative',pitch_range=[65,1220],window_size=window_size)

    fo_bedcmm[bedcmm_score < 0.5] = np.nan

    # 可視化
    fig,ax =plt.subplots(2,1,figsize=(8, 6),sharex=True)
    plt.suptitle("Response Time Comparison: (YIN vs pYIN) vs bedcmm")
    ax[0].plot(t[times], fo_bedcmm, label='bedcmm', linewidth=2,c='g', alpha=0.7)
    f_start, f_end = 100.0,300.0
    f_min = min(f_start-(f_end-f_start)*0.2,f_end+(f_end-f_start)*0.2)
    f_max = max(f_start-(f_end-f_start)*0.2,f_end+(f_end-f_start)*0.2)
    ax[0].set_ylim(f_min,f_max)
    ax[0].legend()
    ax[0].grid(True)
    ax[0].set_xlabel('Time[s]')
    ax[0].set_ylabel('Freqency[Hz]')
    ax[1].plot(t,signal)
    ax[1].set_xlabel('Time[s]')
    ax[1].set_ylabel('Amplitude')
    plt.show()
    pass

if __name__ == "__main__":
    main()