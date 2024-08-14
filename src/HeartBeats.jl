module HeartBeats

export detect_heartbeats, example_ecg

using DelimitedFiles, DSP, Statistics

"""
    example_ecg()

Return example ECG data as obtained from `scipy.misc.electrocardiogram()`.
"""
function example_ecg()
    vec(readdlm(joinpath(pkgdir(HeartBeats, "data"), "ecg.txt")))
end

"""
    detect_heartbeats(ecg, fs)

Detect heartbeats in ECG signal `ecg` (recorded with sampling frequency `fs`).
"""
function detect_heartbeats(ecg::AbstractVector, fs::Real)
    # Pan-Tompkins detector based on https://github.com/cbrnr/sleepecg
    # bandpass filter between 5-30Hz
    filtered = filtfilt(digitalfilter(Bandpass(5, 30, fs=fs), Butterworth(2)), ecg)

    # set samples before first zero crossing (within first two seconds) to zero
    start = findfirst(!iszero, diff(signbit.(filtered[1:2*fs])))
    !isnothing(start) && filtered[1:start] .= 0

    # 5-point derivative
    derivative = xcorr(filtered, [-1, -2, 0, 2, 1])[3:end-2]

    integrated = conv(derivative.^2, ones(trunc(Int, 0.15 * fs)))
    middle, half = length(integrated) ÷ 2, length(derivative) ÷ 2
    integrated = integrated[middle - half + 1 : middle + half]  # mode="same"

    signal_len = length(filtered)
    beat_mask = zeros(Bool, signal_len)
    refractory_samples = trunc(Int, 0.2 * fs)  # 200ms
    t_wave_samples = trunc(Int, 0.36 * fs)  # 360ms

    n = trunc(Int, 2 * fs)
    s_f = maximum(filtered[1:n])  # filtered signal
    n_f = mean(filtered[1:n])  # filtered noise
    s_i = maximum(integrated[1:n])  # integrated signal
    n_i = mean(integrated[1:n])  # integrated noise
    threshold_i1 = n_i + 0.25 * (s_i - n_i)
    threshold_f1 = n_f + 0.25 * (s_f - n_f)

    rr_intervals = zeros(signal_len ÷ refractory_samples)
    num_peaks_found = 0

    peak_index = -refractory_samples + 2
    previous_peak_index = -refractory_samples + 2
    RR_missed_limit = nothing
    RR_low_limit = nothing
    RR_high_limit = nothing

    # searchback
    do_searchback = true
    index = 2
    while index < signal_len
        PEAKF = nothing
        PEAKI = nothing
        signal_peak_found = false
        noise_peak_found = false
        if (num_peaks_found > 1 && index - previous_peak_index > RR_missed_limit && do_searchback) ||
           (num_peaks_found == 0 && index > fs) ||
           (num_peaks_found == 1 && index - previous_peak_index > 1.5 * fs)
            for i = 1:16
                found_a_candidate = false
                searchback_divisor = 2^i
                best_searchback_index = previous_peak_index + refractory_samples
                best_candidate_amplitude = -1
                searchback_index = best_searchback_index
                while searchback_index < index
                    PEAKF = filtered[searchback_index]
                    if PEAKF > filtered[searchback_index + 1]
                        if PEAKF > filtered[searchback_index - 1]
                            PEAKI = integrated[searchback_index]
                            if (threshold_f1 / searchback_divisor < PEAKF && PEAKF < threshold_f1 && threshold_i1 / searchback_divisor < PEAKI) ||
                               (threshold_i1 / searchback_divisor < PEAKI && PEAKI < threshold_i1 && threshold_f1 / searchback_divisor < PEAKF)
                                if PEAKF > best_candidate_amplitude
                                    best_searchback_index = searchback_index
                                    best_candidate_amplitude = filtered[searchback_index]
                                    found_a_candidate = true
                                end
                            end
                        end
                        searchback_index += 1
                    end
                    searchback_index += 1
                end
                if found_a_candidate
                    s_i = 0.25 * PEAKI + 0.75 * s_i
                    s_f = 0.25 * PEAKF + 0.75 * s_f
                    signal_peak_found = true
                    peak_index = best_searchback_index
                    do_searchback = false
                    break
                end
            end
        elseif filtered[index] > filtered[index + 1]
            if filtered[index] > filtered[index - 1]
                PEAKF = filtered[index]
                PEAKI = integrated[index]
                if PEAKF > threshold_f1 && PEAKI > threshold_i1
                    s_f = 0.125 * PEAKF + 0.875 * s_f
                    s_i = 0.125 * PEAKI + 0.875 * s_i

                    signal_peak_found = true
                    peak_index = index
                else
                    noise_peak_found = true
                end
            end
            index += 1
        end
        if signal_peak_found && num_peaks_found > 0
            RR = peak_index - previous_peak_index

            # t-wave identification
            if RR < t_wave_samples
                reverse_index = peak_index
                max_slope_in_this_peak = -1
                while reverse_index > 0
                    amplitude_here = filtered[reverse_index]
                    amplitude_before = filtered[reverse_index - 1]
                    if amplitude_before > amplitude_here
                        break
                    end
                    slope = amplitude_here - amplitude_before
                    if slope > max_slope_in_this_peak
                        max_slope_in_this_peak = slope
                    end
                    reverse_index -= 1
                end

                reverse_index = previous_peak_index
                max_slope_in_previous_peak = -1
                while reverse_index > 0
                    amplitude_here = filtered[reverse_index]
                    amplitude_before = filtered[reverse_index - 1]
                    if amplitude_before > amplitude_here
                        break
                    end
                    slope = amplitude_here - amplitude_before
                    if (slope > max_slope_in_previous_peak)
                        max_slope_in_previous_peak = slope
                    end
                    reverse_index -= 1
                end

                if max_slope_in_this_peak < max_slope_in_previous_peak / 2.0
                    signal_peak_found = false
                    noise_peak_found = true
                end
            end
        end

        if signal_peak_found
            num_peaks_found += 1
            beat_mask[peak_index] = true

            threshold_i1 = n_i + 0.25 * (s_i - n_i)
            threshold_f1 = n_f + 0.25 * (s_f - n_f)

            if num_peaks_found > 1
                rr_intervals[num_peaks_found] = peak_index - previous_peak_index

                # learning phase 2
                if num_peaks_found == 2
                    RR_low_limit = 0.92 * rr_intervals[num_peaks_found]
                    RR_high_limit = 1.16 * rr_intervals[num_peaks_found]
                end

                # RR average 1 / RR average 2
                RR_sum = 0
                RR_count = 0
                irregular = false
                for i = num_peaks_found:-1:1
                    RR_n = rr_intervals[i]
                    if RR_low_limit < RR_n && RR_n < RR_high_limit
                        RR_sum += RR_n
                        RR_count += 1
                        if RR_count >= 8
                            break
                        end
                    else
                        irregular = true
                    end
                end

                RR_average = RR_sum / RR_count
                RR_low_limit = 0.92 * RR_average
                RR_high_limit = 1.16 * RR_average
                RR_missed_limit = 1.66 * RR_average

                if irregular
                    threshold_f1 /= 2
                    threshold_i1 /= 2
                end
            end

            do_searchback = true
            previous_peak_index = peak_index
            index = peak_index + refractory_samples

        elseif noise_peak_found
            n_i = 0.125 * PEAKI + 0.875 * n_i
            n_f = 0.125 * PEAKF + 0.875 * n_f
            threshold_i1 = n_i + 0.25 * (s_i - n_i)
            threshold_f1 = n_f + 0.25 * (s_f - n_f)
        end
        index += 1
    end
    findall(beat_mask)
end

end
