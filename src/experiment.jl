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
