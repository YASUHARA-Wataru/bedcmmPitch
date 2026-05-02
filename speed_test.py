import numpy as np
import time
import pandas as pd
import bedcmmPitch

# ====== あなたの関数に置き換える ======
def pitch_detect(signal, sr):
    """
    ここを自分のピッチ検出関数に置き換えてください
    """
    #pitch,score = bedcmmPitch.calc_Pitch(signal,fs=sr,pp_mode='positive+negative')
    pitch,score = bedcmmPitch.calc_Pitch(signal,fs=sr,pp_mode='positive')

    return pitch


# ====== テスト用音声生成 ======
def generate_audio(duration_sec, sr=44100, freq=220, noise_std=0.01, seed=0):
    rng = np.random.default_rng(seed)
    t = np.arange(int(sr * duration_sec)) / sr

    # 正弦波（0〜1にシフト）
    signal = 0.5 * (np.sin(2 * np.pi * freq * t) + 1.0)

    # ノイズ追加
    noise = rng.normal(0, noise_std, len(signal))
    signal = signal + noise

    # クリップ（非負にする）
    signal = np.clip(signal, 0, None)

    return signal


# ====== ベンチマーク ======
def benchmark(durations=[10, 20, 30], sr=44100, repeat=3):
    results = []

    for duration in durations:
        signal = generate_audio(duration, sr)

        times = []
        for _ in range(repeat):
            start = time.perf_counter()
            pitch_detect(signal, sr)
            end = time.perf_counter()
            times.append(end - start)

        avg_time = np.mean(times)
        rt_factor = avg_time / duration

        results.append({
            "audio_sec": duration,
            "proc_sec": avg_time,
            "real_time_factor": rt_factor
        })

    return pd.DataFrame(results)


# ====== 実行 ======
if __name__ == "__main__":
    df = benchmark(repeat=5)

    print("\n=== Benchmark Result ===")
    print(df.to_string(index=False, float_format="%.4f"))