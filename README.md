# fft

Public-domain single-header library implementing radix-2 decimation-in-time FFT
(i.e. FFT for powers of 2). Pascal port of
[wareya/fft](https://github.com/wareya/fft).

## api

```pascal
procedure fft_core(
  input_real: PDouble; // pointer to real-valued spatial samples (for audio,
                       // this is where your entire audio signal goes)
  input_imag: PDouble; // pointer to imaginary-valued ones (not useful for
                       // audio) in_imag is allowed to be nil. If so, it will
                       // be treated as if it were all zeroes.
  size: UInt64; // number of complex samples per domain. for audio, this is the
                // number of real samples you have. must be a power of 2.
                // Algorithm will definitely fail and possibly crash otherwise,
                // not tested.
  gap: UInt64;  // must be 1 for outside callers. used for recursion.
  output_real: PDouble; // pointer to space for real-valued output. does not
                        // need to be initialized. must be allocated.
  output_imag: PDouble; // same as above, for imaginary. not optional.
                        // out_real and out_imag work together to store a
                        // complex number (2d vector) representing the phase
                        // and amplitude of the given frequency band, even for
                        // wholly real inputs.
  forwards: Boolean // if true, transform is forwards (fft). if false,
                    // transform is backwards (ifft).
  );

{$IF not Defined(FFT_CORE_ONLY)}
procedure normalize_fft(input_real, input_imag: PDouble; size: UInt64);
  // divide the amplitude of each bin by the number of bins. obligatory
  // after fft() for audio. modifies the input.

procedure half_normalize_fft(input_real, input_imag: PDouble; size: UInt64);

procedure fft(input_real, input_imag: PDouble; size: UInt64; output_real, output_imag: PDouble);
  // compute forwards fft, args are same as fft_core.

procedure ifft(input_real, input_imag: PDouble; size: UInt64; output_real, output_imag: PDouble);
  // compute backwards fft (inverse fft, ifft), args are same as fft_core

procedure sanitize_fft(input_real, input_imag: PDouble; size: UInt64);
  // moves all data to positive-frequency bins.  yes, FFTs have negative
  // frequencies for some reason. they're used to retain correlation data for
  // complex inputs. for real inputs, the negative frequencies just mirror the
  // positive ones and sap half their amplitude, therefore this function. for
  // an explanation of what negative frequencies mean, see:
  // http://dsp.stackexchange.com/questions/431/what-is-the-physical-significance-of-negative-frequencies

procedure unsanitize_fft(input_real, input_imag: PDouble; size: UInt64);
  // undo the above. note again that these two fuctions are not sensical
  // for complex inputs.
{$ENDIF}
```
