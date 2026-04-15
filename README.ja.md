# bedcmmPitch
## ピッチ検出アルゴリズム
周期性解析手法を用いたピッチ検出のアルゴリズム
Pitchを計算するcalc_Pitch関数と途中の周期性解析の情報を出力するcalc_bedcmmの二つの関数をもつ

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

## 注意点
基本的にデフォルトパラメータで動かしてください。
差異はあると思うので、興味があればパラメータを変更してみてください。
pitch_detect_threは、pitch_detect_modeがdynamicの時は、入力信号の平均値と閾値の掛け算で、bedcmmの結果の閾値を作成します。また、pitch_detect_modeがstaticの時は、そのまま閾値として用います。
pitch_detect_modeがmaxの時は、pitch_rangeを絞らないとうまく機能しません。

## 簡易的なインストール方法
1. フォルダ「bedcmmPitch」を作業フォルダにコピー
2. 「setup.py」を作業フォルダににコピー
3. numpy、cython、setuptools(、matplotlib(main表示用のため))がライブラリに入っている事を確認
    1. 上記ライブラリが無ければ、uvやpipで導入
4. 作業フォルダにいて、環境が動く状態でプロンプトから「python setup.py build_ext --inplace」の実行
5. main.pyの実行(動作確認)

※4.は飛ばしてもpure python実装が動きます。


