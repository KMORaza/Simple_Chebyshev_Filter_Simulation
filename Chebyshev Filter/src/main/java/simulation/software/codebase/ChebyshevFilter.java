package simulation.software.codebase;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Implements a Chebyshev filter, which provides a sharper transition between passband and stopband
 * compared to Butterworth filters, at the cost of ripple in the passband (Type I) or stopband (Type II).
 * The filter is designed using the pole-zero placement method in the s-plane, with transformations
 * for different filter types (low-pass, high-pass, band-pass, band-stop).
 */
public class ChebyshevFilter {
    /**
     * Enumerates Chebyshev filter types:
     * - TYPE_I: Passband ripple, flat stopband.
     * - TYPE_II: Flat passband, stopband ripple (inverse Chebyshev).
     */
    public enum FilterType { TYPE_I, TYPE_II }

    /**
     * Enumerates filter design types:
     * - LOW_PASS: Passes low frequencies, attenuates high frequencies.
     * - HIGH_PASS: Passes high frequencies, attenuates low frequencies.
     * - BAND_PASS: Passes a frequency band, attenuates others.
     * - BAND_STOP: Attenuates a frequency band, passes others.
     */
    public enum DesignType { LOW_PASS, HIGH_PASS, BAND_PASS, BAND_STOP }

    /**
     * Enumerates design representation methods:
     * - TRANSFER_FUNCTION: Polynomial ratio H(s) = N(s)/D(s).
     * - ZPK: Zero-pole-gain form.
     * - SOS: Second-order sections form for numerical stability.
     */
    public enum DesignMethod { TRANSFER_FUNCTION, ZPK, SOS }

    private FilterType filterType;
    private DesignType designType;
    private int order;
    private double passbandRipple; // in dB for Type I
    private double stopbandAttenuation; // in dB for Type II
    private double cutoffFrequency1; // rad/s, lower cutoff for band filters
    private double cutoffFrequency2; // rad/s, upper cutoff for band filters
    private double passbandEdge; // rad/s, edge of passband
    private double stopbandEdge; // rad/s, edge of stopband
    private double epsilon; // Ripple factor for Chebyshev response
    private List<Complex> poles; // Filter poles in s-plane
    private List<Complex> zeros; // Filter zeros in s-plane
    private double gain; // Filter gain

    /**
     * Constructs a Chebyshev filter with specified parameters.
     *
     * @param filterType Filter type (Type I or Type II).
     * @param designType Filter design (low-pass, high-pass, band-pass, band-stop).
     * @param order Filter order (number of poles).
     * @param passbandRipple Passband ripple in dB (Type I).
     * @param stopbandAttenuation Stopband attenuation in dB (Type II).
     * @param cutoffFrequency1 First cutoff frequency in rad/s.
     * @param cutoffFrequency2 Second cutoff frequency in rad/s (for band filters).
     * @param passbandEdge Passband edge frequency in rad/s.
     * @param stopbandEdge Stopband edge frequency in rad/s.
     */
    public ChebyshevFilter(FilterType filterType, DesignType designType, int order, double passbandRipple,
                           double stopbandAttenuation, double cutoffFrequency1, double cutoffFrequency2,
                           double passbandEdge, double stopbandEdge) {
        this.filterType = filterType;
        this.designType = designType;
        this.order = order;
        this.passbandRipple = passbandRipple;
        this.stopbandAttenuation = stopbandAttenuation;
        this.cutoffFrequency1 = cutoffFrequency1;
        this.cutoffFrequency2 = cutoffFrequency2;
        this.passbandEdge = passbandEdge;
        this.stopbandEdge = stopbandEdge;
        this.epsilon = calculateEpsilon();
        this.poles = new ArrayList<>();
        this.zeros = new ArrayList<>();
        // Gain for Type I: normalizes passband; Type II: unity gain
        this.gain = filterType == FilterType.TYPE_I ? 1.0 / (Math.pow(2, order - 1) * epsilon) : 1.0;
        calculatePolesAndZeros();
    }

    /**
     * Calculates the ripple factor epsilon based on filter type.
     * For Type I: ε = sqrt(10^(Rp/10) - 1), where Rp is passband ripple in dB.
     * For Type II: ε = 1 / sqrt(10^(As/10) - 1), where As is stopband attenuation in dB.
     *
     * @return Epsilon value for Chebyshev response.
     */
    private double calculateEpsilon() {
        if (filterType == FilterType.TYPE_I) {
            return Math.sqrt(Math.pow(10, passbandRipple / 10.0) - 1);
        } else {
            return 1.0 / Math.sqrt(Math.pow(10, stopbandAttenuation / 10.0) - 1);
        }
    }

    /**
     * Computes poles and zeros for the Chebyshev filter in the s-plane.
     * For Type I: Poles lie on an ellipse defined by sinh and cosh functions.
     * For Type II: Poles are scaled reciprocals, zeros lie on the imaginary axis.
     * Poles are calculated as: s_m = sinh(σ) * sin(θ_m) + j * cosh(σ) * cos(θ_m),
     * where σ = asinh(1/ε)/n, θ_m = π/2 * (2m-1)/n.
     */
    private void calculatePolesAndZeros() {
        // asinh term defines the shape of the pole ellipse
        double asinhTerm = Math.log(1.0 / epsilon + Math.sqrt(1.0 / (epsilon * epsilon) + 1)) / order;
        // Center frequency for band filters
        double omega_0 = designType == DesignType.BAND_PASS || designType == DesignType.BAND_STOP ?
                Math.sqrt(cutoffFrequency1 * cutoffFrequency2) : cutoffFrequency1;
        // Bandwidth for band filters
        double bw = designType == DesignType.BAND_PASS || designType == DesignType.BAND_STOP ?
                cutoffFrequency2 - cutoffFrequency1 : 0;

        for (int m = 1; m <= order; m++) {
            // Pole angle: distributes poles symmetrically
            double theta_m = Math.PI / 2.0 * (2 * m - 1) / order;
            double realPart = Math.sinh(asinhTerm) * Math.sin(theta_m);
            double imagPart = Math.cosh(asinhTerm) * Math.cos(theta_m);

            if (filterType == FilterType.TYPE_I) {
                if (realPart < 0) { // Only left-half plane poles for stability
                    Complex pole = new Complex(-realPart, imagPart);
                    transformPole(pole);
                }
            } else {
                // Type II: Scale poles by 1/|s|^2, add zeros at 1/cos(θ_m)
                double scale = 1.0 / (realPart * realPart + imagPart * imagPart);
                if (realPart < 0) {
                    Complex pole = new Complex(-realPart * scale, imagPart * scale);
                    transformPole(pole);
                }
                Complex zero = new Complex(0, 1.0 / Math.cos(theta_m));
                transformZero(zero);
            }
        }

        if (filterType == FilterType.TYPE_I && order % 2 == 0) {
            adjustEvenOrderPoles();
        }
    }

    /**
     * Transforms a pole based on the filter design type.
     * - Low-pass: Direct mapping.
     * - High-pass: s → 1/s.
     * - Band-pass/stop: s → (s^2 + ω_0^2)/(s * BW), where ω_0 is center frequency, BW is bandwidth.
     *
     * @param pole Pole to transform.
     */
    private void transformPole(Complex pole) {
        if (designType == DesignType.LOW_PASS) {
            poles.add(pole);
        } else if (designType == DesignType.HIGH_PASS) {
            double denominator = pole.re * pole.re + pole.im * pole.im;
            poles.add(new Complex(-pole.re / denominator, pole.im / denominator));
        } else {
            double omega_0 = Math.sqrt(cutoffFrequency1 * cutoffFrequency2);
            double bw = cutoffFrequency2 - cutoffFrequency1;
            Complex s = pole;
            Complex term = new Complex(0, omega_0).multiply(s).add(new Complex(omega_0 * omega_0, 0));
            Complex sqrtTerm = term.sqrt();
            Complex s1 = s.multiply(bw / 2.0).add(sqrtTerm);
            Complex s2 = s.multiply(bw / 2.0).subtract(sqrtTerm);
            if (s1.re < 0) poles.add(s1);
            if (s2.re < 0) poles.add(s2);
        }
    }

    /**
     * Transforms a zero based on the filter design type.
     * - Low-pass: Direct mapping.
     * - High-pass: s → 1/s.
     * - Band-pass: Zeros at cutoff frequencies.
     * - Band-stop: Zero at center frequency.
     *
     * @param zero Zero to transform.
     */
    private void transformZero(Complex zero) {
        if (designType == DesignType.LOW_PASS) {
            zeros.add(zero);
        } else if (designType == DesignType.HIGH_PASS) {
            double denominator = zero.re * zero.re + zero.im * zero.im;
            zeros.add(new Complex(-zero.re / denominator, zero.im / denominator));
        } else if (designType == DesignType.BAND_PASS) {
            zeros.add(new Complex(0, cutoffFrequency1));
            zeros.add(new Complex(0, cutoffFrequency2));
        } else if (designType == DesignType.BAND_STOP) {
            double omega_0 = Math.sqrt(cutoffFrequency1 * cutoffFrequency2);
            zeros.add(new Complex(0, omega_0));
        }
    }

    /**
     * Adjusts poles for even-order Type I filters to ensure correct passband ripple.
     * Scales poles to match the ripple specification using cos(π*(n-1)/(2n)).
     */
    private void adjustEvenOrderPoles() {
        List<Complex> newPoles = new ArrayList<>();
        double angle = Math.cos(Math.PI * (order - 1) / (2.0 * order));
        for (Complex pole : poles) {
            double scale = Math.sqrt((pole.re * pole.re + pole.im * pole.im + angle * angle) / (1 - angle * angle));
            newPoles.add(new Complex(pole.re * scale, pole.im * scale));
        }
        poles = newPoles;
    }

    /**
     * Checks filter stability by ensuring all poles have negative real parts.
     *
     * @return True if all poles are in the left-half s-plane.
     */
    public boolean isStable() {
        for (Complex pole : poles) {
            if (pole.getReal() >= 0) {
                return false;
            }
        }
        return true;
    }

    /**
     * Calculates minimum filter order to meet passband and stopband requirements.
     * Uses: n ≥ log(√(As/Ap) + √(As/Ap - 1)) / log(ωs/ωp + √((ωs/ωp)^2 - 1)),
     * where Ap is passband ripple factor, As is stopband attenuation factor,
     * ωp is passband edge, ωs is stopband edge.
     *
     * @param passbandRipple Passband ripple in dB.
     * @param stopbandAttenuation Stopband attenuation in dB.
     * @param passbandEdge Passband edge frequency in rad/s.
     * @param stopbandEdge Stopband edge frequency in rad/s.
     * @return Minimum required order.
     */
    public static int calculateMinOrder(double passbandRipple, double stopbandAttenuation, double passbandEdge, double stopbandEdge) {
        double Ap = Math.pow(10, passbandRipple / 10.0) - 1;
        double As = Math.pow(10, stopbandAttenuation / 10.0) - 1;
        double ratio = stopbandEdge / passbandEdge;
        double numerator = Math.log(Math.sqrt(As / Ap) + Math.sqrt(As / Ap - 1));
        double denominator = Math.log(ratio + Math.sqrt(ratio * ratio - 1));
        return (int) Math.ceil(numerator / denominator);
    }

    /**
     * Computes the frequency response (magnitude in dB) from startFreq to endFreq.
     * H(jω) = gain * Π(s - z_i) / Π(s - p_i), where s = jω/ω_0.
     * Magnitude: 20 * log10(|H(jω)|).
     *
     * @param startFreq Starting frequency in rad/s.
     * @param endFreq Ending frequency in rad/s.
     * @param points Number of frequency points.
     * @return Array of gain values in dB.
     */
    public double[] getFrequencyResponse(double startFreq, double endFreq, int points) {
        double[] frequencies = new double[points];
        double[] gain = new double[points];
        double step = (endFreq - startFreq) / (points - 1);
        double omega_0 = designType == DesignType.BAND_PASS || designType == DesignType.BAND_STOP ?
                Math.sqrt(cutoffFrequency1 * cutoffFrequency2) : cutoffFrequency1;

        for (int i = 0; i < points; i++) {
            double omega = startFreq + i * step;
            Complex s = new Complex(0, omega / omega_0);
            Complex H = computeTransferFunction(s);
            gain[i] = 20 * Math.log10(H.abs());
            frequencies[i] = omega;
        }
        return gain;
    }

    /**
     * Computes the linear magnitude response from startFreq to endFreq.
     * |H(jω)| = |gain * Π(jω/ω_0 - z_i) / Π(jω/ω_0 - p_i)|.
     *
     * @param startFreq Starting frequency in rad/s.
     * @param endFreq Ending frequency in rad/s.
     * @param points Number of frequency points.
     * @return Array of linear magnitude values.
     */
    public double[] getMagnitudeResponseLinear(double startFreq, double endFreq, int points) {
        double[] frequencies = new double[points];
        double[] magnitude = new double[points];
        double step = (endFreq - startFreq) / (points - 1);
        double omega_0 = designType == DesignType.BAND_PASS || designType == DesignType.BAND_STOP ?
                Math.sqrt(cutoffFrequency1 * cutoffFrequency2) : cutoffFrequency1;

        for (int i = 0; i < points; i++) {
            double omega = startFreq + i * step;
            Complex s = new Complex(0, omega / omega_0);
            Complex H = computeTransferFunction(s);
            magnitude[i] = H.abs();
            frequencies[i] = omega;
        }
        return magnitude;
    }

    /**
     * Computes the phase response from startFreq to endFreq.
     * Phase: arg(H(jω)) = atan2(Im(H), Re(H)).
     *
     * @param startFreq Starting frequency in rad/s.
     * @param endFreq Ending frequency in rad/s.
     * @param points Number of frequency points.
     * @return Array of phase values in degrees.
     */
    public double[] getPhaseResponse(double startFreq, double endFreq, int points) {
        double[] frequencies = new double[points];
        double[] phase = new double[points];
        double step = (endFreq - startFreq) / (points - 1);
        double omega_0 = designType == DesignType.BAND_PASS || designType == DesignType.BAND_STOP ?
                Math.sqrt(cutoffFrequency1 * cutoffFrequency2) : cutoffFrequency1;

        for (int i = 0; i < points; i++) {
            double omega = startFreq + i * step;
            Complex s = new Complex(0, omega / omega_0);
            Complex H = computeTransferFunction(s);
            phase[i] = Math.toDegrees(Math.atan2(H.im, H.re));
            frequencies[i] = omega;
        }
        return phase;
    }

    /**
     * Computes the group delay from startFreq to endFreq.
     * Group delay: -d[arg(H(jω))]/dω ≈ -Σ[Re((s - p_i)^{-1})], where s = jω/ω_0.
     *
     * @param startFreq Starting frequency in rad/s.
     * @param endFreq Ending frequency in rad/s.
     * @param points Number of frequency points.
     * @return Array of group delay values in seconds.
     */
    public double[] getGroupDelay(double startFreq, double endFreq, int points) {
        double[] frequencies = new double[points];
        double[] groupDelay = new double[points];
        double step = (endFreq - startFreq) / (points - 1);
        double omega_0 = designType == DesignType.BAND_PASS || designType == DesignType.BAND_STOP ?
                Math.sqrt(cutoffFrequency1 * cutoffFrequency2) : cutoffFrequency1;

        for (int i = 0; i < points; i++) {
            double omega = startFreq + i * step;
            Complex s = new Complex(0, omega / omega_0);
            groupDelay[i] = computeGroupDelay(s) * omega_0;
            frequencies[i] = omega;
        }
        return groupDelay;
    }

    /**
     * Computes group delay for a given complex frequency s.
     * Formula: -Σ[Re((s - p_i)^{-1})].
     *
     * @param s Complex frequency.
     * @return Group delay contribution.
     */
    private double computeGroupDelay(Complex s) {
        double sum = 0;
        for (Complex pole : poles) {
            Complex diff = s.subtract(pole);
            sum += diff.re / (diff.re * diff.re + diff.im * diff.im);
        }
        return -sum;
    }

    /**
     * Computes the transfer function H(s) = gain * Π(s - z_i) / Π(s - p_i).
     *
     * @param s Complex frequency.
     * @return Transfer function value.
     */
    private Complex computeTransferFunction(Complex s) {
        Complex numerator = new Complex(gain, 0);
        Complex denominator = new Complex(1, 0);
        for (Complex pole : poles) {
            denominator = denominator.multiply(s.subtract(pole));
        }
        for (Complex zero : zeros) {
            numerator = numerator.multiply(s.subtract(zero));
        }
        return numerator.divide(denominator);
    }

    /**
     * Computes the impulse response using the inverse Laplace transform.
     * h(t) = Σ[gain * p_i * e^(p_i * t)] for real poles,
     * or 2 * gain * e^(Re(p_i) * t) * [Re(p_i) * cos(Im(p_i) * t) - Im(p_i) * sin(Im(p_i) * t)] for complex poles.
     *
     * @param duration Time duration in seconds.
     * @param dt Time step in seconds.
     * @return Impulse response array.
     */
    public double[] getImpulseResponse(double duration, double dt) {
        int points = (int) (duration / dt) + 1;
        double[] impulse = new double[points];
        double[] time = new double[points];
        for (int i = 0; i < points; i++) {
            double t = i * dt;
            time[i] = t;
            double sum = 0;
            for (Complex pole : poles) {
                if (pole.im == 0) {
                    sum += pole.re * Math.exp(pole.re * t);
                } else {
                    double real = pole.re;
                    double imag = pole.im;
                    sum += 2 * Math.exp(real * t) * (real * Math.cos(imag * t) - imag * Math.sin(imag * t));
                }
            }
            impulse[i] = gain * sum;
        }
        return impulse;
    }

    /**
     * Computes the step response using the inverse Laplace transform.
     * Step response: ∫h(τ)dτ = Σ[gain * (e^(p_i * t) - 1)/p_i] for real poles,
     * or Σ[gain * (Re(p_i) * (e^(Re(p_i) * t) * cos(Im(p_i) * t) - 1) + Im(p_i) * e^(Re(p_i) * t) * sin(Im(p_i) * t)) / |p_i|^2] for complex poles.
     *
     * @param duration Time duration in seconds.
     * @param dt Time step in seconds.
     * @return Step response array.
     */
    public double[] getStepResponse(double duration, double dt) {
        int points = (int) (duration / dt) + 1;
        double[] step = new double[points];
        double[] time = new double[points];
        for (int i = 0; i < points; i++) {
            double t = i * dt;
            time[i] = t;
            double sum = 0;
            for (Complex pole : poles) {
                if (pole.im == 0) {
                    sum += (Math.exp(pole.re * t) - 1) / pole.re;
                } else {
                    double real = pole.re;
                    double imag = pole.im;
                    double denom = real * real + imag * imag;
                    sum += (real * (Math.exp(real * t) * Math.cos(imag * t) - 1) + imag * Math.exp(real * t) * Math.sin(imag * t)) / denom;
                }
            }
            step[i] = gain * sum;
        }
        return step;
    }

    /**
     * Applies the filter to an input signal using convolution with the impulse response.
     * y[n] = Σ[h[k] * x[n-k]], where h is the impulse response, x is the input signal.
     *
     * @param input Input signal array.
     * @param dt Time step in seconds.
     * @return Filtered output signal.
     */
    public double[] applyFilter(double[] input, double dt) {
        double[] impulse = getImpulseResponse(1.0, dt);
        double[] output = new double[input.length];
        for (int n = 0; n < input.length; n++) {
            double sum = 0;
            for (int k = 0; k < impulse.length && n - k >= 0; k++) {
                sum += impulse[k] * input[n - k];
            }
            output[n] = sum;
        }
        return output;
    }

    /**
     * Generates a frequency vector from startFreq to endFreq.
     *
     * @param startFreq Starting frequency in rad/s.
     * @param endFreq Ending frequency in rad/s.
     * @param points Number of points.
     * @return Frequency array.
     */
    public double[] getFrequencies(double startFreq, double endFreq, int points) {
        double[] frequencies = new double[points];
        double step = (endFreq - startFreq) / (points - 1);
        for (int i = 0; i < points; i++) {
            frequencies[i] = startFreq + i * step;
        }
        return frequencies;
    }

    /**
     * Returns the list of filter poles.
     *
     * @return List of complex poles.
     */
    public List<Complex> getPoles() {
        return poles;
    }

    /**
     * Returns the list of filter zeros.
     *
     * @return List of complex zeros.
     */
    public List<Complex> getZeros() {
        return zeros;
    }

    /**
     * Returns the filter gain.
     *
     * @return Gain value.
     */
    public double getGain() {
        return gain;
    }

    /**
     * Returns the filter design in the specified representation.
     *
     * @param method Design representation method.
     * @return String representation of the filter design.
     */
    public String getDesignRepresentation(DesignMethod method) {
        switch (method) {
            case TRANSFER_FUNCTION:
                return getTransferFunction();
            case ZPK:
                return getZPK();
            case SOS:
                return getSOS();
            default:
                return "Unknown method";
        }
    }

    /**
     * Returns the transfer function as numerator and denominator polynomials.
     * H(s) = gain * Π(s - z_i) / Π(s - p_i).
     *
     * @return String of polynomial coefficients.
     */
    private String getTransferFunction() {
        double[] num = polyFromRoots(zeros);
        double[] den = polyFromRoots(poles);
        for (int i = 0; i < num.length; i++) {
            num[i] *= gain;
        }
        return "Numerator: " + Arrays.toString(num) + "\nDenominator: " + Arrays.toString(den);
    }

    /**
     * Returns the zero-pole-gain representation.
     *
     * @return String of zeros, poles, and gain.
     */
    private String getZPK() {
        StringBuilder sb = new StringBuilder();
        sb.append("Zeros:\n");
        for (Complex z : zeros) {
            sb.append(String.format("(%.4f, %.4f)\n", z.re, z.im));
        }
        sb.append("Poles:\n");
        for (Complex p : poles) {
            sb.append(String.format("(%.4f, %.4f)\n", p.re, p.im));
        }
        sb.append(String.format("Gain: %.4f\n", gain));
        return sb.toString();
    }

    /**
     * Returns the second-order sections (SOS) representation.
     * Each section: [b0, b1, b2, a0, a1, a2], where H(s) = Π(b0 + b1s + b2s^2)/(a0 + a1s + a2s^2).
     *
     * @return String of SOS matrix.
     */
    private String getSOS() {
        List<double[]> sos = new ArrayList<>();
        List<Complex> remainingPoles = new ArrayList<>(poles);
        List<Complex> remainingZeros = new ArrayList<>(zeros);

        while (!remainingPoles.isEmpty()) {
            double[] section = new double[6];
            section[3] = 1.0; // a0 = 1

            Complex pole1 = remainingPoles.remove(0);
            Complex pole2 = remainingPoles.stream()
                    .filter(p -> Math.abs(p.re - pole1.re) < 1e-10 && Math.abs(p.im + pole1.im) < 1e-10)
                    .findFirst().orElse(null);
            if (pole2 != null) {
                remainingPoles.remove(pole2);
                section[4] = -2 * pole1.re; // a1 = -2 * Re(p)
                section[5] = pole1.re * pole1.re + pole1.im * pole1.im; // a2 = |p|^2
            } else {
                section[4] = -pole1.re; // a1 = -Re(p)
                section[5] = 0; // a2 = 0 for real pole
            }

            if (!remainingZeros.isEmpty()) {
                Complex zero1 = remainingZeros.remove(0);
                Complex zero2 = remainingZeros.stream()
                        .filter(z -> Math.abs(z.re - zero1.re) < 1e-10 && Math.abs(z.im + zero1.im) < 1e-10)
                        .findFirst().orElse(null);
                if (zero2 != null) {
                    remainingZeros.remove(zero2);
                    section[1] = -2 * zero1.re; // b1 = -2 * Re(z)
                    section[2] = zero1.re * zero1.re + zero1.im * zero1.im; // b2 = |z|^2
                } else {
                    section[1] = -zero1.re; // b1 = -Re(z)
                    section[2] = 0; // b2 = 0 for real zero
                }
            }
            section[0] = 1.0; // b0 = 1
            sos.add(section);
        }

        double sectionGain = Math.pow(gain, 1.0 / sos.size());
        for (double[] section : sos) {
            section[0] *= sectionGain;
            section[1] *= sectionGain;
            section[2] *= sectionGain;
        }

        StringBuilder sb = new StringBuilder("SOS Matrix:\n");
        for (int i = 0; i < sos.size(); i++) {
            sb.append(String.format("Section %d: [%.4f, %.4f, %.4f, %.4f, %.4f, %.4f]\n",
                    i + 1, sos.get(i)[0], sos.get(i)[1], sos.get(i)[2], sos.get(i)[3], sos.get(i)[4], sos.get(i)[5]));
        }
        return sb.toString();
    }

    /**
     * Converts roots to polynomial coefficients.
     * For roots r_i, computes Π(s - r_i) using iterative polynomial multiplication.
     *
     * @param roots List of complex roots.
     * @return Polynomial coefficients.
     */
    private double[] polyFromRoots(List<Complex> roots) {
        double[] poly = new double[roots.size() + 1];
        poly[0] = 1.0;
        for (Complex root : roots) {
            double[] newPoly = new double[poly.length + 1];
            for (int i = 0; i < poly.length; i++) {
                newPoly[i] += poly[i];
                newPoly[i + 1] -= root.re * poly[i];
                if (Math.abs(root.im) > 1e-10) {
                    newPoly[i + 1] -= root.im * poly[i];
                }
            }
            poly = newPoly;
        }
        return poly;
    }

    // Getters
    public FilterType getFilterType() { return filterType; }
    public DesignType getDesignType() { return designType; }
    public int getOrder() { return order; }
    public double getPassbandRipple() { return passbandRipple; }
    public double getStopbandAttenuation() { return stopbandAttenuation; }
    public double getCutoffFrequency1() { return cutoffFrequency1; }
    public double getCutoffFrequency2() { return cutoffFrequency2; }
    public double getPassbandEdge() { return passbandEdge; }
    public double getStopbandEdge() { return stopbandEdge; }

    /**
     * Complex number class for s-plane calculations.
     */
    public static class Complex {
        private final double re;
        private final double im;

        public Complex(double real, double imag) {
            this.re = real;
            this.im = imag;
        }

        public Complex add(Complex other) {
            return new Complex(re + other.re, im + other.im);
        }

        public Complex subtract(Complex other) {
            return new Complex(re - other.re, im - other.im);
        }

        public Complex multiply(Complex other) {
            double real = re * other.re - im * other.im;
            double imag = re * other.im + im * other.re;
            return new Complex(real, imag);
        }

        public Complex multiply(double scalar) {
            return new Complex(re * scalar, im * scalar);
        }

        public Complex divide(Complex other) {
            double denominator = other.re * other.re + other.im * other.im;
            double real = (re * other.re + im * other.im) / denominator;
            double imag = (im * other.re - re * other.im) / denominator;
            return new Complex(real, imag);
        }

        public Complex sqrt() {
            double r = Math.sqrt(re * re + im * im);
            double theta = Math.atan2(im, re) / 2.0;
            return new Complex(Math.sqrt(r) * Math.cos(theta), Math.sqrt(r) * Math.sin(theta));
        }

        public double abs() {
            return Math.sqrt(re * re + im * im);
        }

        public double getReal() { return re; }
        public double getImag() { return im; }
    }
}