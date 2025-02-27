struct Experiment{M,N,K<:Spectrum,D<:Detector,S<:Sample}
    spectrum::K
    detector::Motorized{M,D}
    sample::Motorized{N,S}
end

abstract type FrameSettings end

struct StillFrameSettings{M,N} <: FrameSettings
    angles_d::Vec{M}
    angles_s::Vec{N}
end

struct ScanFrameSettings{M,N} <: FrameSettings
    angles_d::Vec{M}
    angles_s::Vec{N}
    n_axis::Integer
    increment::Number
end

function model(model::Experiment, settings::StillFrameSettings)
    StillModel(
        model.spectrum,
        model.detector(settings.angles_d...),
        model.sample(settings.angles_s...),
    )
end

function model(model::Experiment, settings::ScanFrameSettings)
    ScanModel(
        model.spectrum,
        model.detector(settings.angles_d...),
        model.sample(settings.angles_s...),
        model.sample.goniometer.axes[settings.n_axis],
        settings.increment
    )
end
