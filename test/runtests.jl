using HeartBeats, Test

ecg = example_ecg()
fs = 360

@test length(detect_heartbeats(ecg, fs)) == 478
