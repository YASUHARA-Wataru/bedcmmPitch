# bedcmmPitch
## Pitch Detection Algorithm

This repository provides a pitch detection algorithm based on a periodicity analysis method.  
It is intended for research and Proof of Concept (PoC) use.
Robust to spike-like noise and impulsive artifacts.

The implementation includes two main functions:
- `calc_Pitch`: computes pitch values
- `calc_bedcmm`: outputs intermediate periodicity analysis results

## About bedcmm

`bedcmm` is a periodicity analysis method used as the core of this pitch detection algorithm.

For more details, please refer to the original repository:
https://github.com/YASUHARA-Wataru/bedcmm
---

## Requirements
- numpy
- cython
- matplotlib

Install with:

pip install numpy cython matplotlib

## Functions

### calc_Pitch

**Input (default):**
- `data` : 1D array (signal data)
- `fs (44100)` : float (sampling rate)
- `window_size (2048)` : int (window size)
- `hop_size (256)` : int (hop size)
- `pitch_range (None)` : [start_freq, end_freq] (search range)
- `pp_mode ('positive+negative')` : str (preprocessing mode: `'positive'`, `'negative'`, `'positive+negative'`, `'threshold_diff'`)
- `pp_threshold (0)` : float (preprocessing threshold, used in `'threshold_diff'` mode)
- `bedcmm_smooth (3)` : int (smoothing size for bedcmm result; 1 means no smoothing)
- `pitch_detect_mode ('peak')` : str (peak detection mode: `'score'`, `'static'`, `'maximum'`, `'peak'`)
- `pitch_detect_threshold (0.85)` : float (threshold for peak detection)
- `interpolator_mode ('parabolic')` : str (peak interpolation mode: `'parabolic'`, `'centroid'`, `'gaussian'`, `'no'`)

**Output:**
- Pitch data: 1D array
- Score data: 1D array (bedcmm score)

---

### calc_bedcmm

**Input:**
- `data` : 1D array (signal data)
- `fs (44100)` : float (sampling rate)
- `window_size (2048)` : int (window size)
- `hop_size (256)` : int (hop size)
- `pitch_range (None)` : [start_freq, end_freq] (search range)
- `pp_mode ('positive')` : str (preprocessing mode)
- `pp_threshold (0)` : float (used in `'threshold_diff'` mode)

**Output:**
- bedcmm data: 2D array (time, bedcmm)
- mean data: 1D array (time)

---

## Notes

It is recommended to use the default parameters.

You may observe differences depending on parameter settings, so feel free to experiment if needed.

- When `pitch_detect_mode = 'score'`, the threshold is calculated as:  
  (mean of input signal) × `pitch_detect_threshold`

- When `pitch_detect_mode = 'static'`, `pitch_detect_threshold` is used directly as a threshold.

- When `pitch_detect_mode = 'peak'`, the threshold is:  
  (maximum peak value of bedcmm) × `pitch_detect_threshold`

- When `pitch_detect_mode = 'maximum'`, it may not function properly unless `pitch_range` is restricted.

- The default `pp_mode` prioritizes accuracy but is computationally expensive.  
  For faster computation, you can use `'positive'` or `'negative'`.

- When using `pp_mode = 'threshold_diff'`, all values must be non-negative; otherwise, the scoring mechanism may not function correctly.

---

## Quick Installation

1. Copy the `bedcmmPitch` folder into your working directory  
2. Copy `setup.py` into your working directory  
3. Ensure the following libraries are installed:
   - numpy
   - cython
   - setuptools  
   - (optional) matplotlib (for visualization)

   If not installed, use `pip` or `uv` to install them.

4. In your working directory, run:```python setup.py build_ext --inplace```


5. Run `main.py` to verify functionality

*Note:* Step 4 can be skipped; the pure Python implementation will still work.

---

## License

This repository is available **for research, evaluation, and Proof of Concept (PoC) purposes only**.

### Permitted Use

- Research use
- Algorithm evaluation and validation
- Proof of Concept (PoC)

PoC usage is allowed **without time limitation**.

### Restrictions

The following uses are prohibited:

- Commercial use (including integration into products or services, or any monetized usage)
- Use in production environments (including integration into operational systems)

### Commercial Use

If you wish to use this software for commercial purposes or in a production environment,  
a **separate commercial license agreement is required**.

Please contact via Issues or direct email for inquiries.

### Disclaimer

This software is provided "AS IS", without warranty of any kind.  
The author shall not be held liable for any damages arising from the use of this software.

---

*For detailed terms, please refer to the `LICENSE` file.*

---

## Contact
For commercial use, please contact us with a brief description of your use case.
fapow.contact[at]gmail.com