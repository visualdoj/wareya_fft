unit wareya_fft;

//===== fft.hpp =========================
//      Public-domain single-header library
//      implementing radix-2 decimation-in-time FFT (i.e. FFT for powers of 2)
//
//  This software is dual-licensed to the public domain and under the following
//  license: you are granted a perpetual, irrevocable license to copy, modify,
//  publish, and distribute this file as you see fit.
//
//  THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
//  WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
//  MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY
//  SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
//  WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
//  ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR
//  IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
//

{$MODE FPC}
{$MODESWITCH DEFAULTPARAMETERS}
{$MODESWITCH OUT}
{$MODESWITCH RESULT}

interface

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

implementation

// address of cell if base adderss not nil, nil otherwise
function fft_private_safe_addrof(ptr: PDouble; i: PtrInt): Pointer; inline;
begin
  if ptr <> nil then begin
    Exit(@(ptr[i]));
  end else
    Exit(nil);
end;

// For a 8-sample input, the FFT's last three bins contain "negative"
// frequencies. (So, the last (size/2)-1 bins.) They are only meaningful for
// complex inputs.
procedure fft_core(input_real, input_imag: PDouble; size, gap: UInt64; output_real, output_imag: PDouble; forwards: Boolean);
var
  i: UInt64;
  a_real, a_imag, b_real, b_imag: Double;
  twiddle_real, twiddle_imag, bias_real, bias_imag: Double;
begin
  if size = 1 then begin
    output_real[0] := input_real[0];
    if input_imag <> nil then begin
      output_imag[0] := input_imag[0];
    end else
      output_imag[0] := 0;
  end else begin
    size := size div 2;
    // This algorithm works by extending the concept of how two-bin DFTs
    // (discrete fourier transform) work, in order to correlate decimated DFTs,
    // recursively.  No, I'm not your guy if you want a proof of why it works,
    // but it does.
    fft_core(input_real, input_imag, size, gap * 2, output_real, output_imag, forwards);
    fft_core(@input_real[gap], fft_private_safe_addrof(input_imag,gap), size, gap*2, @output_real[size], @output_imag[size], forwards);
    // non-combed decimated output to non-combed correlated output
    i := 0;
    while i < size do begin
      a_real := output_real[i];
      a_imag := output_imag[i];
      b_real := output_real[i + size];
      b_imag := output_imag[i + size];

      twiddle_real := cos(pi * i / size);
      twiddle_imag := sin(pi * i / size) * (1 - 2 * Ord(forwards));
      // complex multiplication (vector angle summing and length multiplication)
      bias_real := b_real * twiddle_real - b_imag * twiddle_imag;
      bias_imag := b_imag * twiddle_real + b_real * twiddle_imag;
      // real output (sum of real parts)
      output_real[i       ] := a_real + bias_real;
      output_real[i + size] := a_real - bias_real;
      // imag output (sum of imaginary parts)
      output_imag[i       ] := a_imag + bias_imag;
      output_imag[i + size] := a_imag - bias_imag;
      Inc(i);
    end;
  end;
end;

{$IF not Defined(FFT_CORE_ONLY)}
procedure normalize_fft(input_real, input_imag: PDouble; size: UInt64);
var
  i: UInt64;
begin
  i := 0;
  while i < size do begin
    input_real[i] := input_real[i] / size;
    input_imag[i] := input_imag[i] / size;
    Inc(i);
  end;
end;
procedure half_normalize_fft(input_real, input_imag: PDouble; size: UInt64);
var
  i: UInt64;
begin
  i := 0;
  while i < size do begin
    input_real[i] := input_real[i] / Sqrt(size);
    input_imag[i] := input_imag[i] / Sqrt(size);
    Inc(i);
  end;
end;
procedure fft(input_real, input_imag: PDouble; size: UInt64; output_real, output_imag: PDouble);
begin
  fft_core(input_real, input_imag, size, 1, output_real, output_imag, True);
  half_normalize_fft(output_real, output_imag, size); // allows calling fft() four times to result in the original signal with no amplitude change
end;
procedure ifft(input_real, input_imag: PDouble; size: UInt64; output_real, output_imag: PDouble);
begin
  fft_core(input_real, input_imag, size, 1, output_real, output_imag, False);
  half_normalize_fft(output_real, output_imag, size); // see above, also causes ifft(fft(x)) to result in the original signal with no amplitude change
end;

// boost bins that are split into positive (A-handed spin) and negative
// (B-handed spin) parts only useful if former input signal was not complex,
// for only needing to look at one bin to get the magnitude
// FIXME or HELPME: How come the nyquist frequency is quiet in saw waves, but
// loud in pure signal?
procedure sanitize_fft(input_real, input_imag: PDouble; size: UInt64);
var
  i: UInt64;
begin
  i := 1;
  while i < size div 2 do begin
    input_real[i] := input_real[i] * 2;
    input_imag[i] := input_imag[i] * 2;
    input_real[size - i] := input_real[size - i] * 2;
    input_imag[size - i] := input_imag[size - i] * 2;
    Inc(i);
  end;
end;
// opposite of above
procedure unsanitize_fft(input_real, input_imag: PDouble; size: UInt64);
var
  i: UInt64;
begin
  i := 1;
  while i < size div 2 do begin
    input_real[i] := input_real[i] / 2;
    input_imag[i] := input_imag[i] / 2;
    input_real[size - i] := input_real[size - i] / 2;
    input_imag[size - i] := input_imag[size - i] / 2;
    Inc(i);
  end;
end;
{$ENDIF}

end.
