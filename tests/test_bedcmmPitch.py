# tests/test_bedcmmpitch.py

import numpy as np
import pytest

from bedcmmPitch import calc_Pitch, calc_bedcmm


def generate_sine(
    freq=440.0,
    fs=44100,
    duration=1.0,
    amplitude=1.0
):
    t = np.arange(int(fs * duration)) / fs
    return amplitude * np.sin(2 * np.pi * freq * t)


def test_import():
    """
    import確認
    """
    from bedcmmPitch import calc_Pitch
    assert calc_Pitch is not None


def test_calc_pitch_runs():
    """
    calc_Pitch が正常終了する
    """
    x = generate_sine()

    pitch, score = calc_Pitch(
        x,
        fs=44100,
        window_size=2048,
        hop_size=256,
    )

    assert pitch is not None
    assert score is not None


def test_calc_pitch_shape():
    """
    pitch と score の shape が一致
    """
    x = generate_sine()

    pitch, score = calc_Pitch(x)

    assert pitch.shape == score.shape
    assert pitch.ndim == 1


def test_calc_pitch_detect_440hz():
    """
    440Hz を大まかに検出できる
    """
    x = generate_sine(freq=440.0)

    pitch, score = calc_Pitch(
        x,
        fs=44100,
        window_size=4096,
        hop_size=512,
        fmin=100,
        fmax=1000,
    )

    valid_pitch = pitch[~np.isnan(pitch)]

    assert len(valid_pitch) > 0

    mean_pitch = np.mean(valid_pitch)

    # ±10Hzくらいの緩い判定
    assert abs(mean_pitch - 440.0) < 10.0


def test_calc_bedcmm_runs():
    """
    calc_bedcmm が正常終了する
    """
    x = generate_sine()

    result, mean_data = calc_bedcmm(x)

    assert result is not None
    assert mean_data is not None


def test_calc_bedcmm_shape():
    """
    bedcmm出力shape確認
    """
    x = generate_sine()

    result, mean_data = calc_bedcmm(x)

    assert result.ndim == 2
    assert mean_data.ndim == 1
    assert result.shape[0] == mean_data.shape[0]


def test_invalid_dimension():
    """
    2次元入力で例外
    """
    x = np.zeros((10, 10))

    with pytest.raises(Exception):
        calc_Pitch(x)


@pytest.mark.parametrize(
    "pp_mode",
    [
        "positive",
        "negative",
        "positive+negative",
        "threshold_diff",
    ]
)
def test_pp_modes(pp_mode):
    """
    pp_mode が全て動作する
    """
    x = generate_sine()

    pitch, score = calc_Pitch(
        x,
        pp_mode=pp_mode
    )

    assert pitch is not None
    assert score is not None


def test_invalid_pp_mode():
    """
    不正 pp_mode で例外
    """
    x = generate_sine()

    with pytest.raises(Exception):
        calc_Pitch(
            x,
            pp_mode="invalid_mode"
        )

@pytest.mark.parametrize(
    "mode",
    [
        "parabolic",
        "gaussian",
        "centroid",
        "no",
    ]
)

def test_interpolator_modes(mode):

    x = generate_sine()

    pitch, score = calc_Pitch(
        x,
        interpolator_mode=mode
    )

    assert pitch is not None

@pytest.mark.parametrize(
    "mode",
    [
        "score",
        "static",
        "maximum",
        "peak",
    ]
)
def test_pitch_detect_modes(mode):

    x = generate_sine()

    pitch, score = calc_Pitch(
        x,
        pitch_detect_mode=mode
    )

    assert pitch is not None

def test_fmin_only():

    x = generate_sine()

    pitch, score = calc_Pitch(
        x,
        fmin=80,
        fmax=None
    )

    assert pitch is not None


def test_fmax_only():

    x = generate_sine()

    pitch, score = calc_Pitch(
        x,
        fmin=None,
        fmax=1000
    )

    assert pitch is not None


def test_no_fmin_fmax():

    x = generate_sine()

    pitch, score = calc_Pitch(
        x,
        fmin=None,
        fmax=None
    )

    assert pitch is not None

def test_bedcmm_smooth():

    x = generate_sine()

    pitch, score = calc_Pitch(
        x,
        bedcmm_smooth=1
    )

    assert pitch is not None


def test_invalid_bedcmm_smooth():

    x = generate_sine()

    with pytest.raises(Exception):
        calc_Pitch(
            x,
            bedcmm_smooth=0
        )

def test_nan_input():

    x = generate_sine()

    x[100:200] = np.nan

    pitch, score = calc_Pitch(x)

    assert pitch is not None

def test_short_input():

    x = np.zeros(100)

    pitch, score = calc_Pitch(
        x,
        window_size=2048
    )

    assert len(pitch) == 0

def test_silence():

    x = np.zeros(44100)

    pitch, score = calc_Pitch(x)

    assert pitch is not None