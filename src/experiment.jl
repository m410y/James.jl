struct Experiment{M,N,K<:Spectrum,D<:Detector,S<:Sample}
    spectrum::K
    detector::Motorized{M,D}
    sample::Motorized{N,S}
end

abstract type FrameSetting end

struct StillFrameSetting{M,N} <: FrameSetting
    angles_d::Vec{M}
    angles_s::Vec{N}
end

struct ScanFrameSetting{M,N} <: FrameSetting
    angles_d::Vec{M}
    angles_s::Vec{N}
    n_axis::Integer
    increment::Number
end

function model(experiment::Experiment, setting::StillFrameSetting)
    StillModel(
        experiment.spectrum,
        experiment.detector(setting.angles_d...),
        experiment.sample(setting.angles_s...),
    )
end

function model(experiment::Experiment, setting::ScanFrameSetting)
    ScanModel(
        experiment.spectrum,
        experiment.detector(setting.angles_d...),
        experiment.sample(setting.angles_s...),
        experiment.sample.goniometer.axes[setting.n_axis],
        setting.increment,
    )
end

function Base.show(io::IO, ::MIME"text/plain", experiment::Experiment)
    show(io, "text/plain", experiment.spectrum)
    println(io)
    show(io, "text/plain", experiment.detector.goniometer)
    println(io)
    show(io, "text/plain", experiment.detector.object)
    println(io)
    show(io, "text/plain", experiment.sample.goniometer)
    println(io)
    show(io, "text/plain", experiment.sample.object)
end

function Base.show(io::IO, ::MIME"text/plain", setting::StillFrameSetting)
    print(io, "StillFrameSetting:\n")
    print(io, "  detector angles: $(Tuple(setting.angles_d))\n")
    print(io, "  sample angles: $(Tuple(setting.angles_s))")
end

function Base.show(io::IO, ::MIME"text/plain", setting::ScanFrameSetting)
    print(io, "ScanFrameSetting:\n")
    print(io, "  detector angles: $(Tuple(setting.angles_d))\n")
    print(io, "  sample angles: $(Tuple(setting.angles_s))\n")
    print(io, "  sample axis $(setting.n_axis) increment: $(setting.increment)")
end
