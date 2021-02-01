{$MODE FPC}
{$MODESWITCH DEFAULTPARAMETERS}
{$MODESWITCH OUT}
{$MODESWITCH RESULT}

uses
  wareya_fft;

function lerp(a, b, i: Double): Double; inline;
begin
  Exit(a * (1 - i) + b * i);
end;

type
TTestEnum = (
    TEST_IMPULSE
   ,TEST_DC
   ,TEST_NYQUIST
   ,TEST_HALF
   ,TEST_SQUARE
   ,TEST_SQUAREDOUBLE
   ,TEST_SAW
   ,TEST_NOISE
   ,TEST_COMPLEXIMPULSE
   ,TEST_COMPLEXDC
   ,TEST_COMPLEXNOISE
   ,TEST_BLANK
);

var
  TEST: TTestEnum = TEST_NYQUIST;

procedure init(in_real, in_imag: PDouble; size: UInt64);
var
  i: Int32;
begin
  for i := 0 to size - 1 do begin
    case TEST of
    TEST_IMPULSE: begin
        in_real[i] := 0;
        in_imag[i] := 0;
        in_real[0] := 1;
      end;
    TEST_DC: begin
        in_real[i] := 1;
        in_imag[i] := 0;
      end;
    TEST_NYQUIST: begin
        in_real[i] := cos(pi * i);
        in_imag[i] := 0;
      end;
    TEST_HALF: begin
        in_real[i] := cos(pi * i / 2);
        in_imag[i] := 0;
      end;
    TEST_SQUARE: begin
        if i < size / 2 then begin
          in_real[i] :=   1;
        end else
          in_real[i] := - 1;
        in_imag[i] := 0;
      end;
    TEST_SQUAREDOUBLE: begin
        if (i mod (size div 2)) < (size div 4) then begin
          in_real[i] :=   1;
        end else
          in_real[i] := - 1;
        in_imag[i] := 0;
      end;
    TEST_SAW: begin
        in_real[i] := lerp(-1, 1, i / double(size - 1)); // subtract 1 from size to avoid DC from incorrect starting phase
        in_imag[i] := 0;
      end;
    TEST_NOISE: begin
        in_real[i] := Random(1024) / 512.0 - 1;
        in_imag[i] := 0;
      end;
    TEST_COMPLEXIMPULSE: begin
        in_real[i] := 0;
        in_imag[i] := 0;
        in_real[0] := 1;
        in_imag[0] := 1;
      end;
    TEST_COMPLEXDC: begin
        in_real[i] := 1;
        in_imag[i] := 1;
      end;
    TEST_COMPLEXNOISE: begin
        in_real[i] := Random(1024) / 512.0 - 1;
        in_imag[i] := Random(1024) / 512.0 - 1;
      end;
    else // TEST_BLANK
      in_real[i] := 0;
      in_imag[i] := 0;
    end;
  end;
end;

const
  SIZE = 16;
var
  in_real: PDouble;
  in_imag: PDouble;
  out_real: PDouble;
  out_imag: PDouble;
  i: Int32;
  signal_was_real: Boolean;
  factor: Double;
  S: ShortString;
begin
  if ParamCount > 0 then begin
    S := LowerCase(ParamStr(1));
    if      S = 'impulse'        then TEST := TEST_IMPULSE
    else if S = 'dc'             then TEST := TEST_DC
    else if S = 'nyquist'        then TEST := TEST_NYQUIST
    else if S = 'half'           then TEST := TEST_HALF
    else if S = 'square'         then TEST := TEST_SQUARE
    else if S = 'squaredouble'   then TEST := TEST_SQUAREDOUBLE
    else if S = 'saw'            then TEST := TEST_SAW
    else if S = 'noise'          then TEST := TEST_NOISE
    else if S = 'compleximpulse' then TEST := TEST_COMPLEXIMPULSE
    else if S = 'complexdc'      then TEST := TEST_COMPLEXDC
    else if S = 'complexnoise'   then TEST := TEST_COMPLEXNOISE
    else if S = 'blank'          then TEST := TEST_BLANK;
  end;

  Randomize;
  in_real  := GetMem(SizeOf(Double) * SIZE);
  in_imag  := GetMem(SizeOf(Double) * SIZE);
  out_real := GetMem(SizeOf(Double) * SIZE);
  out_imag := GetMem(SizeOf(Double) * SIZE);
  init(in_real, in_imag, SIZE);

  Writeln('data');
  Writeln('sample'#9'real'#9'imag'#9'mag');
  for i := 0 to SIZE - 1 do
    Writeln(i, #9, in_real[i]:0:2, #9, in_imag[i]:0:2, #9, Sqrt(Sqr(in_real[i]) + Sqr(in_imag[i])), ' ');
  Writeln('');

  fft(in_real, in_imag, SIZE, out_real, out_imag);

  // end FFT display for real inputs early
  signal_was_real := (TEST < TEST_COMPLEXIMPULSE) or (TEST > TEST_COMPLEXNOISE);

  Writeln('transform');
  Writeln('bin'#9'real'#9'imag'#9'mag');
  for i := 0 to SIZE - 1 do begin
    // boost nonunique frequencies for real inputs on display
    if (i > 0) and (i < SIZE div 2) and signal_was_real then begin
      factor := 2;
    end else
      factor := 1;
    Writeln(i, #9, out_real[i]*factor:0:2, #9, out_imag[i]*factor:0:2, #9, Sqrt(Sqr(out_real[i]) + Sqr(out_imag[i]))*factor);
    // end FFT display for real inputs early
    if (i >= SIZE div 2) and signal_was_real then
      break;
  end;
  Writeln('');

  ifft(out_real, out_imag, SIZE, in_real, in_imag);

  Writeln('inverse');
  Writeln('sample'#9'real'#9'imag'#9'mag');
  for i := 0 to SIZE - 1 do
    Writeln(i, #9, in_real[i]:0:2, #9, in_imag[i]:0:2, #9, Sqrt(Sqr(in_real[i]) + Sqr(in_imag[i])));
  Writeln('');
end.
