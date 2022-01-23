![License](https://img.shields.io/github/license/cbrnr/HeartBeats.jl)

HeartBeats.jl
=============
HeartBeats.jl provides a heartbeat detector based on the approach described by [Pan & Tompkins (1985)](https://ieeexplore.ieee.org/document/4122029). It is based on the implementation available in the Python package [SleepECG](https://github.com/cbrnr/sleepecg).

## Installation
Use the package manager to add HeartBeats.jl by typing `] add HeartBeats` in the Julia REPL.


## Example
HeartBeats.jl contains a short example ECG dataset taken from [`scipy.misc.electrocardiogram()`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.misc.electrocardiogram.html). The function `example_ecg()` returns this data, which was sampled with a sampling frequency of 360&nbsp;Hz, as a `Vector{Float64}`. We can use this dataset to showcase the `detect_heartbeats()` function:

```julia
using HeartBeats

ecg = example_ecg()
fs = 360  # sampling frequency

beats = detect_heartbeats(ecg, fs)
```

The `beats` array will contain all detected R peak locations.

## Benchmark
The detector is based on the Python implementation available in [SleepECG](https://github.com/cbrnr/sleepecg). It is about 18× faster than the Python implementation and only 2× slower than the C implementation. Follow these steps to reproduce the benchmark:

1. Export all data records used in the `'runtime'` benchmark by including `export_records = True` in `config.yml` (refer to the [SleepECG documentation](https://sleepecg.readthedocs.io/en/stable/heartbeat_detection.html) for details on how to set up and run the benchmarks). This will generate 15 text files.
2. Move those text files to a folder that you can access from Julia (i.e. set the Julia working directory accordingly).
3. Run the following code snippet to benchmark the runtime for 60 minute data segments (note that you need to add the `CSV` package):

```julia
using CSV, HeartBeats

function run_benchmark()
    total = 0
    fs = 128
    files = filter(x -> endswith(x, ".txt"), readdir(".", join=true))
    for file in files
        println(file)
        data = CSV.File(file, comment="#")
        ecg = view(data[:ecg], 1:60*60*fs)  # 60 minutes
        stats = @timed detect_heartbeats(ecg, fs)
        total += stats.time
    end
    println("Total: $total, Average: $(total / length(files))")
    total
end

run_benchmark()
```